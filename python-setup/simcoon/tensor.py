"""
Unified Tensor2 and Tensor4 classes — single and batch in one class (scipy.Rotation style).

Single tensor: wraps ``_CppTensor2`` / ``_CppTensor4`` exactly as before.
Batch of N tensors: stores ``(N, 6)`` / ``(N, 6, 6)`` numpy arrays, batch ops loop in C++.

``t.single`` is True for a single tensor, False for a batch.

Examples
--------
>>> import simcoon as smc
>>> import numpy as np

>>> # Single tensor (unchanged API)
>>> L = smc.Tensor4.stiffness(smc.L_iso([70000, 0.3], 'Enu'))
>>> eps = smc.Tensor2.strain(np.array([0.01, -0.003, -0.003, 0.005, 0.002, 0.001]))
>>> sigma = L @ eps

>>> # Batch from (N, 6) array
>>> eps_batch = smc.Tensor2.strain(np.random.randn(100, 6) * 0.01)
>>> eps_batch.single  # False
>>> len(eps_batch)     # 100
>>> eps_batch[0]       # single Tensor2

>>> # Batch operations loop in C++
>>> sigma_batch = L @ eps_batch
"""

import collections.abc
import numpy as np

from simcoon._core import (
    _CppTensor2,
    _CppTensor4,
    _dyadic,
    _auto_dyadic,
    VoigtType,
    Tensor4Type,
)

# C++ batch functions
from simcoon._core import (
    _batch_t2_rotate,
    _batch_t2_push_forward,
    _batch_t2_pull_back,
    _batch_t2_mises,
    _batch_t2_trace,
    _batch_t4_contract,
    _batch_t4_rotate,
    _batch_t4_push_forward,
    _batch_t4_pull_back,
    _batch_t4_inverse,
)


# ======================================================================
# Voigt conversion helpers (for batch numpy ↔ mat)
# ======================================================================

def _mat_to_voigt(m, vtype):
    """Convert (..., 3, 3) matrices to (..., 6) Voigt vectors."""
    v = np.empty((*m.shape[:-2], 6), dtype=np.float64)
    v[..., 0] = m[..., 0, 0]
    v[..., 1] = m[..., 1, 1]
    v[..., 2] = m[..., 2, 2]
    if vtype == VoigtType.strain:
        v[..., 3] = m[..., 0, 1] + m[..., 1, 0]
        v[..., 4] = m[..., 0, 2] + m[..., 2, 0]
        v[..., 5] = m[..., 1, 2] + m[..., 2, 1]
    else:
        v[..., 3] = 0.5 * (m[..., 0, 1] + m[..., 1, 0])
        v[..., 4] = 0.5 * (m[..., 0, 2] + m[..., 2, 0])
        v[..., 5] = 0.5 * (m[..., 1, 2] + m[..., 2, 1])
    return v


def _voigt_to_mat(v, vtype):
    """Convert (..., 6) Voigt vectors to (..., 3, 3) matrices."""
    m = np.empty((*v.shape[:-1], 3, 3), dtype=np.float64)
    m[..., 0, 0] = v[..., 0]
    m[..., 1, 1] = v[..., 1]
    m[..., 2, 2] = v[..., 2]
    if vtype == VoigtType.strain:
        m[..., 0, 1] = m[..., 1, 0] = 0.5 * v[..., 3]
        m[..., 0, 2] = m[..., 2, 0] = 0.5 * v[..., 4]
        m[..., 1, 2] = m[..., 2, 1] = 0.5 * v[..., 5]
    else:
        m[..., 0, 1] = m[..., 1, 0] = v[..., 3]
        m[..., 0, 2] = m[..., 2, 0] = v[..., 4]
        m[..., 1, 2] = m[..., 2, 1] = v[..., 5]
    return m


def _get_rotation_matrices(R, N):
    """Extract (N, 3, 3) rotation matrices from a scipy Rotation."""
    from scipy.spatial.transform import Rotation as ScipyRotation
    if not isinstance(R, ScipyRotation):
        raise TypeError(f"Expected scipy Rotation, got {type(R)}")
    mats = np.atleast_3d(R.as_matrix())
    if mats.ndim == 2:
        mats = mats[np.newaxis]
    if mats.shape == (3, 3, 1):
        mats = mats.transpose(2, 0, 1)
    if mats.shape[0] == 1:
        mats = np.broadcast_to(mats, (N, 3, 3)).copy()
    if mats.shape[0] != N:
        raise ValueError(
            f"Rotation batch size {mats.shape[0]} != tensor batch size {N}"
        )
    return mats


# ======================================================================
# Tensor2 — unified single / batch
# ======================================================================

class Tensor2:
    """A 2nd-order tensor with type tag for Voigt convention and rotation dispatch.

    Handles both single tensors and batches of N tensors.
    Use ``.single`` to check which mode.

    Single: wraps one ``_CppTensor2`` internally.
    Batch: stores ``(N, 6)`` Voigt array, dispatches to C++ batch loops.

    Do not construct directly — use ``stress()``, ``strain()``, etc.
    """

    __slots__ = ("_cpp", "_voigt_data", "_vtype", "_single")

    def __init__(self, data):
        """Create a Tensor2.

        Parameters
        ----------
        data : _CppTensor2, list of Tensor2, or internal
        """
        if isinstance(data, _CppTensor2):
            # Single from C++ object
            self._cpp = data
            self._voigt_data = None
            self._vtype = data.vtype
            self._single = True
        elif isinstance(data, list) and data and isinstance(data[0], Tensor2):
            # Batch from list of single Tensor2
            if not all(t._single for t in data):
                raise ValueError("Cannot nest batches")
            ref = data[0].vtype
            v = np.empty((len(data), 6), dtype=np.float64)
            for i, t in enumerate(data):
                if t.vtype != ref:
                    raise ValueError(f"Mixed VoigtType: {ref} vs {t.vtype} at index {i}")
                v[i] = t.voigt
            self._cpp = None
            self._voigt_data = v
            self._vtype = ref
            self._single = False
        else:
            raise TypeError("Use Tensor2.stress(), Tensor2.strain(), etc.")

    @classmethod
    def _from_batch_voigt(cls, voigt_N6, vtype):
        """Internal: create batch from (N, 6) array + vtype, no copy."""
        obj = object.__new__(cls)
        obj._cpp = None
        obj._voigt_data = np.ascontiguousarray(voigt_N6, dtype=np.float64)
        obj._vtype = vtype
        obj._single = False
        return obj

    @classmethod
    def _from_single_cpp(cls, cpp_obj):
        """Internal: create single from _CppTensor2."""
        obj = object.__new__(cls)
        obj._cpp = cpp_obj
        obj._voigt_data = None
        obj._vtype = cpp_obj.vtype
        obj._single = True
        return obj

    # ------------------------------------------------------------------
    # Factory methods — shape detection for single vs batch
    # ------------------------------------------------------------------

    @classmethod
    def stress(cls, data):
        """Create stress tensor(s) from array.

        Parameters
        ----------
        data : array_like
            (6,) or (3,3) for single, (N,6) or (N,3,3) for batch.
        """
        return cls._from_data(data, VoigtType.stress)

    @classmethod
    def strain(cls, data):
        """Create strain tensor(s) from array.

        Parameters
        ----------
        data : array_like
            (6,) or (3,3) for single, (N,6) or (N,3,3) for batch.
        """
        return cls._from_data(data, VoigtType.strain)

    @classmethod
    def from_mat(cls, m, vtype):
        """Create from 3x3 matrix or (N,3,3) batch with explicit VoigtType."""
        m = np.asarray(m, dtype=np.float64)
        if m.shape == (3, 3):
            return cls._from_single_cpp(_CppTensor2.from_mat(m, vtype))
        elif m.ndim == 3 and m.shape[1:] == (3, 3):
            return cls._from_batch_voigt(_mat_to_voigt(m, vtype), vtype)
        raise ValueError(f"Expected (3,3) or (N,3,3), got {m.shape}")

    @classmethod
    def from_voigt(cls, v, vtype):
        """Create from Voigt vector (6,) or batch (N,6) with explicit VoigtType."""
        v = np.asarray(v, dtype=np.float64)
        if v.ndim == 1 and v.size == 6:
            return cls._from_single_cpp(_CppTensor2.from_voigt(v.ravel(), vtype))
        elif v.ndim == 2 and v.shape[1] == 6:
            return cls._from_batch_voigt(v.copy(), vtype)
        raise ValueError(f"Expected (6,) or (N,6), got {v.shape}")

    @classmethod
    def zeros(cls, vtype=VoigtType.stress):
        """Create a single zero tensor."""
        return cls._from_single_cpp(_CppTensor2.zeros(vtype))

    @classmethod
    def identity(cls, vtype=VoigtType.stress):
        """Create a single identity tensor (3x3 eye)."""
        return cls._from_single_cpp(_CppTensor2.identity(vtype))

    @classmethod
    def from_tensor(cls, t, n):
        """Broadcast a single Tensor2 to a batch of size n."""
        if not t._single:
            raise ValueError("from_tensor requires a single tensor")
        v = np.broadcast_to(t.voigt[np.newaxis, :], (n, 6)).copy()
        return cls._from_batch_voigt(v, t.vtype)

    @classmethod
    def concatenate(cls, batches):
        """Join multiple batches (or singles) into one batch."""
        parts = list(batches)
        if not parts:
            raise ValueError("Nothing to concatenate")
        vtype = parts[0].vtype
        arrays = []
        for b in parts:
            if b.vtype != vtype:
                raise ValueError(f"Mixed VoigtType: {vtype} vs {b.vtype}")
            arrays.append(b.voigt if b._single else b._voigt_data)
        combined = np.concatenate(
            [a[np.newaxis] if a.ndim == 1 else a for a in arrays], axis=0
        )
        return cls._from_batch_voigt(combined, vtype)

    # ------------------------------------------------------------------
    # Properties
    # ------------------------------------------------------------------

    @property
    def single(self):
        """True if this is a single tensor, False if batch."""
        return self._single

    @property
    def voigt(self):
        """Voigt vector: (6,) for single, (N,6) for batch."""
        if self._single:
            return self._cpp.voigt.ravel()
        return self._voigt_data

    @property
    def mat(self):
        """Matrix: (3,3) for single, (N,3,3) for batch."""
        if self._single:
            return self._cpp.mat
        return _voigt_to_mat(self._voigt_data, self._vtype)

    @property
    def vtype(self):
        """VoigtType tag."""
        return self._vtype

    @property
    def voigt_T(self):
        """(6, N) transposed Voigt array (batch only, for C++ interop)."""
        if self._single:
            raise AttributeError("voigt_T only available on batch")
        return self._voigt_data.T

    # ------------------------------------------------------------------
    # Sequence protocol (batch only)
    # ------------------------------------------------------------------

    def __len__(self):
        if self._single:
            raise TypeError("single Tensor2 has no len()")
        return self._voigt_data.shape[0]

    def __getitem__(self, key):
        if self._single:
            raise TypeError("single Tensor2 is not subscriptable")
        if isinstance(key, (int, np.integer)):
            if key < 0:
                key += len(self)
            if key < 0 or key >= len(self):
                raise IndexError(f"index {key} out of range for batch of size {len(self)}")
            return Tensor2.from_voigt(self._voigt_data[key].copy(), self._vtype)
        v = self._voigt_data[key]
        if v.ndim == 1:
            return Tensor2.from_voigt(v.copy(), self._vtype)
        return Tensor2._from_batch_voigt(v.copy(), self._vtype)

    def __iter__(self):
        if self._single:
            raise TypeError("single Tensor2 is not iterable")
        for i in range(len(self)):
            yield self[i]

    def __reversed__(self):
        if self._single:
            raise TypeError("single Tensor2 is not reversible")
        for i in range(len(self) - 1, -1, -1):
            yield self[i]

    def __contains__(self, item):
        if self._single:
            raise TypeError("single Tensor2 does not support 'in'")
        if isinstance(item, Tensor2) and item._single:
            return np.any(np.all(self._voigt_data == item.voigt[np.newaxis, :], axis=1))
        return False

    def count(self, item):
        """Count occurrences of a single Tensor2 in batch."""
        if self._single:
            raise TypeError("count only on batch")
        if isinstance(item, Tensor2) and item._single:
            return int(np.sum(np.all(self._voigt_data == item.voigt[np.newaxis, :], axis=1)))
        return 0

    # ------------------------------------------------------------------
    # Numpy array protocol
    # ------------------------------------------------------------------

    def __array__(self, dtype=None):
        v = self.voigt
        return v.astype(dtype) if dtype else v

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        if method != "__call__":
            return NotImplemented
        raw_inputs = []
        batch_vtype = None
        for inp in inputs:
            if isinstance(inp, Tensor2):
                if inp._single:
                    raw_inputs.append(inp.voigt[np.newaxis, :])
                else:
                    raw_inputs.append(inp._voigt_data)
                batch_vtype = batch_vtype or inp.vtype
            else:
                raw_inputs.append(inp)
        if batch_vtype is None:
            return NotImplemented
        _ALLOWED = {np.add, np.subtract, np.multiply, np.negative,
                    np.true_divide, np.positive}
        if ufunc not in _ALLOWED:
            return NotImplemented
        result = ufunc(*raw_inputs, **kwargs)
        if isinstance(result, np.ndarray) and result.ndim == 2 and result.shape[1] == 6:
            return Tensor2._from_batch_voigt(result, batch_vtype)
        return result

    # ------------------------------------------------------------------
    # Methods — single dispatches to C++, batch to C++ batch functions
    # ------------------------------------------------------------------

    def is_symmetric(self, tol=1e-12):
        """Check symmetry (single only)."""
        if not self._single:
            raise NotImplementedError("is_symmetric not supported on batch")
        return self._cpp.is_symmetric(tol)

    def rotate(self, R, active=True):
        """Rotate tensor(s).

        Parameters
        ----------
        R : simcoon.Rotation or scipy.spatial.transform.Rotation
            Single or batch rotation.
        active : bool
            Active (True) or passive rotation.
        """
        from simcoon.rotation import Rotation as SmcRotation
        from scipy.spatial.transform import Rotation as ScipyRotation

        if self._single:
            if isinstance(R, SmcRotation):
                cpp_rot = R._to_cpp()
            else:
                cpp_rot = R
            return Tensor2._from_single_cpp(self._cpp.rotate(cpp_rot, active))

        # Batch: get rotation matrices
        N = len(self)
        if isinstance(R, ScipyRotation):
            mats = _get_rotation_matrices(R, N)
        elif isinstance(R, SmcRotation):
            mats = _get_rotation_matrices(R, N)
        else:
            raise TypeError(f"Expected Rotation, got {type(R)}")

        result_voigt = _batch_t2_rotate(
            self._voigt_data, self._vtype, mats, active
        )
        return Tensor2._from_batch_voigt(result_voigt, self._vtype)

    def push_forward(self, F):
        """Push-forward via deformation gradient F."""
        F = np.asarray(F, dtype=np.float64)
        if self._single:
            return Tensor2._from_single_cpp(self._cpp.push_forward(F))
        F_batch = F[np.newaxis] if F.ndim == 2 else F
        result = _batch_t2_push_forward(self._voigt_data, self._vtype, F_batch)
        return Tensor2._from_batch_voigt(result, self._vtype)

    def pull_back(self, F):
        """Pull-back via deformation gradient F."""
        F = np.asarray(F, dtype=np.float64)
        if self._single:
            return Tensor2._from_single_cpp(self._cpp.pull_back(F))
        F_batch = F[np.newaxis] if F.ndim == 2 else F
        result = _batch_t2_pull_back(self._voigt_data, self._vtype, F_batch)
        return Tensor2._from_batch_voigt(result, self._vtype)

    def mises(self):
        """Von Mises equivalent: float for single, (N,) for batch."""
        if self._single:
            return _mises_impl(self)
        return _batch_t2_mises(self._voigt_data, self._vtype)

    def trace(self):
        """Trace: float for single, (N,) for batch."""
        if self._single:
            m = self.mat
            return m[0, 0] + m[1, 1] + m[2, 2]
        return _batch_t2_trace(self._voigt_data, self._vtype)

    def dev(self):
        """Deviatoric part → Tensor2 (single or batch)."""
        if self._single:
            v = self.voigt.copy()
            tr = v[0] + v[1] + v[2]
            v[0] -= tr / 3.0
            v[1] -= tr / 3.0
            v[2] -= tr / 3.0
            return Tensor2.from_voigt(v, self._vtype)
        v = self._voigt_data.copy()
        tr = v[:, 0] + v[:, 1] + v[:, 2]
        v[:, 0] -= tr / 3.0
        v[:, 1] -= tr / 3.0
        v[:, 2] -= tr / 3.0
        return Tensor2._from_batch_voigt(v, self._vtype)

    def norm(self):
        """Frobenius norm: float for single, (N,) for batch."""
        m = self.mat
        if self._single:
            return np.sqrt(np.sum(m * m))
        return np.sqrt(np.sum(m * m, axis=(-2, -1)))

    # ------------------------------------------------------------------
    # Arithmetic
    # ------------------------------------------------------------------

    def __add__(self, other):
        if not isinstance(other, Tensor2):
            return NotImplemented
        if self._single and other._single:
            return Tensor2._from_single_cpp(self._cpp + other._cpp)
        sv = self.voigt if self._single else self._voigt_data
        ov = other.voigt if other._single else other._voigt_data
        if sv.ndim == 1:
            sv = sv[np.newaxis, :]
        if ov.ndim == 1:
            ov = ov[np.newaxis, :]
        return Tensor2._from_batch_voigt(sv + ov, self._vtype)

    def __radd__(self, other):
        if isinstance(other, Tensor2):
            return other.__add__(self)
        return NotImplemented

    def __sub__(self, other):
        if not isinstance(other, Tensor2):
            return NotImplemented
        if self._single and other._single:
            return Tensor2._from_single_cpp(self._cpp - other._cpp)
        sv = self.voigt if self._single else self._voigt_data
        ov = other.voigt if other._single else other._voigt_data
        if sv.ndim == 1:
            sv = sv[np.newaxis, :]
        if ov.ndim == 1:
            ov = ov[np.newaxis, :]
        return Tensor2._from_batch_voigt(sv - ov, self._vtype)

    def __rsub__(self, other):
        if isinstance(other, Tensor2):
            return other.__sub__(self)
        return NotImplemented

    def __neg__(self):
        if self._single:
            return Tensor2._from_single_cpp(-self._cpp)
        return Tensor2._from_batch_voigt(-self._voigt_data, self._vtype)

    def __mul__(self, other):
        if isinstance(other, (int, float)):
            if self._single:
                return Tensor2._from_single_cpp(self._cpp * float(other))
            return Tensor2._from_batch_voigt(self._voigt_data * other, self._vtype)
        if not self._single:
            other = np.asarray(other)
            if other.ndim <= 1:
                return Tensor2._from_batch_voigt(
                    self._voigt_data * other[..., np.newaxis], self._vtype
                )
        return NotImplemented

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        if isinstance(other, (int, float)):
            if self._single:
                return Tensor2._from_single_cpp(self._cpp / float(other))
            return Tensor2._from_batch_voigt(self._voigt_data / other, self._vtype)
        if not self._single:
            other = np.asarray(other)
            if other.ndim <= 1:
                return Tensor2._from_batch_voigt(
                    self._voigt_data / other[..., np.newaxis], self._vtype
                )
        return NotImplemented

    def __mod__(self, other):
        """Double contraction A_ij B_ij → float or (N,)."""
        if not isinstance(other, Tensor2):
            return NotImplemented
        ma = self.mat
        mb = other.mat
        if self._single and other._single:
            return np.sum(ma * mb)
        if ma.ndim == 2:
            ma = ma[np.newaxis]
        if mb.ndim == 2:
            mb = mb[np.newaxis]
        return np.sum(ma * mb, axis=(-2, -1))

    def __eq__(self, other):
        if isinstance(other, Tensor2):
            if self._single and other._single:
                return self._cpp == other._cpp
            sv = self.voigt if self._single else self._voigt_data
            ov = other.voigt if other._single else other._voigt_data
            if sv.ndim == 1 and ov.ndim == 1:
                return np.array_equal(sv, ov)
            if sv.ndim == 1:
                return np.all(ov == sv[np.newaxis, :], axis=1)
            if ov.ndim == 1:
                return np.all(sv == ov[np.newaxis, :], axis=1)
            if sv.shape != ov.shape:
                return False
            return np.array_equal(sv, ov)
        return NotImplemented

    def __hash__(self):
        return id(self)

    def __repr__(self):
        if self._single:
            return f"Tensor2(vtype={self.vtype.name})"
        return f"Tensor2(N={len(self)}, vtype={self._vtype.name})"

    # ------------------------------------------------------------------
    # Internal
    # ------------------------------------------------------------------

    @classmethod
    def _from_data(cls, data, vtype):
        data = np.asarray(data, dtype=np.float64)
        # Single
        if data.shape == (3, 3):
            return cls._from_single_cpp(_CppTensor2.from_mat(data, vtype))
        if data.shape == (6,) or (data.ndim == 1 and data.size == 6):
            return cls._from_single_cpp(_CppTensor2.from_voigt(data.ravel(), vtype))
        # Batch
        if data.ndim == 2 and data.shape[1] == 6:
            return cls._from_batch_voigt(data.copy(), vtype)
        if data.ndim == 3 and data.shape[1:] == (3, 3):
            return cls._from_batch_voigt(_mat_to_voigt(data, vtype), vtype)
        raise ValueError(
            f"Expected (6,), (3,3), (N,6), or (N,3,3), got {data.shape}"
        )

    @classmethod
    def from_list(cls, tensors):
        """Stack a list of single Tensor2 into a batch."""
        return cls(list(tensors))

    @classmethod
    def from_columns(cls, arr, vtype):
        """Create batch from (6, N) column-major array (C++ interop)."""
        arr = np.asarray(arr, dtype=np.float64)
        if arr.ndim != 2 or arr.shape[0] != 6:
            raise ValueError(f"Expected (6, N), got {arr.shape}")
        return cls._from_batch_voigt(arr.T.copy(), vtype)


# ======================================================================
# Tensor4 — unified single / batch
# ======================================================================

class Tensor4:
    """A 4th-order tensor with type tag for rotation dispatch and Voigt storage.

    Handles both single tensors and batches of N tensors.
    Use ``.single`` to check which mode.
    """

    __slots__ = ("_cpp", "_voigt_data", "_type", "_single")

    def __init__(self, data):
        if isinstance(data, _CppTensor4):
            self._cpp = data
            self._voigt_data = None
            self._type = data.type
            self._single = True
        elif isinstance(data, list) and data and isinstance(data[0], Tensor4):
            if not all(t._single for t in data):
                raise ValueError("Cannot nest batches")
            ref = data[0].type
            v = np.empty((len(data), 6, 6), dtype=np.float64)
            for i, t in enumerate(data):
                if t.type != ref:
                    raise ValueError(f"Mixed Tensor4Type: {ref} vs {t.type} at index {i}")
                v[i] = t.mat
            self._cpp = None
            self._voigt_data = v
            self._type = ref
            self._single = False
        else:
            raise TypeError("Use Tensor4.stiffness(), Tensor4.compliance(), etc.")

    @classmethod
    def _from_batch_voigt(cls, voigt_N66, t4type):
        obj = object.__new__(cls)
        obj._cpp = None
        obj._voigt_data = np.ascontiguousarray(voigt_N66, dtype=np.float64)
        obj._type = t4type
        obj._single = False
        return obj

    @classmethod
    def _from_single_cpp(cls, cpp_obj):
        obj = object.__new__(cls)
        obj._cpp = cpp_obj
        obj._voigt_data = None
        obj._type = cpp_obj.type
        obj._single = True
        return obj

    # ------------------------------------------------------------------
    # Factory methods
    # ------------------------------------------------------------------

    @classmethod
    def stiffness(cls, data):
        """Create stiffness tensor(s): (6,6) single or (N,6,6) batch."""
        data = np.asarray(data, dtype=np.float64)
        if data.shape == (6, 6):
            return cls._from_single_cpp(_CppTensor4.from_mat(data, Tensor4Type.stiffness))
        if data.ndim == 3 and data.shape[1:] == (6, 6):
            return cls._from_batch_voigt(data.copy(), Tensor4Type.stiffness)
        raise ValueError(f"Expected (6,6) or (N,6,6), got {data.shape}")

    @classmethod
    def compliance(cls, data):
        """Create compliance tensor(s): (6,6) single or (N,6,6) batch."""
        data = np.asarray(data, dtype=np.float64)
        if data.shape == (6, 6):
            return cls._from_single_cpp(_CppTensor4.from_mat(data, Tensor4Type.compliance))
        if data.ndim == 3 and data.shape[1:] == (6, 6):
            return cls._from_batch_voigt(data.copy(), Tensor4Type.compliance)
        raise ValueError(f"Expected (6,6) or (N,6,6), got {data.shape}")

    @classmethod
    def strain_concentration(cls, data):
        """Create strain concentration tensor(s)."""
        data = np.asarray(data, dtype=np.float64)
        if data.shape == (6, 6):
            return cls._from_single_cpp(_CppTensor4.from_mat(data, Tensor4Type.strain_concentration))
        if data.ndim == 3 and data.shape[1:] == (6, 6):
            return cls._from_batch_voigt(data.copy(), Tensor4Type.strain_concentration)
        raise ValueError(f"Expected (6,6) or (N,6,6), got {data.shape}")

    @classmethod
    def stress_concentration(cls, data):
        """Create stress concentration tensor(s)."""
        data = np.asarray(data, dtype=np.float64)
        if data.shape == (6, 6):
            return cls._from_single_cpp(_CppTensor4.from_mat(data, Tensor4Type.stress_concentration))
        if data.ndim == 3 and data.shape[1:] == (6, 6):
            return cls._from_batch_voigt(data.copy(), Tensor4Type.stress_concentration)
        raise ValueError(f"Expected (6,6) or (N,6,6), got {data.shape}")

    @classmethod
    def from_mat(cls, m, tensor_type):
        """Create from (6,6) or (N,6,6) with explicit Tensor4Type."""
        m = np.asarray(m, dtype=np.float64)
        if m.shape == (6, 6):
            return cls._from_single_cpp(_CppTensor4.from_mat(m, tensor_type))
        if m.ndim == 3 and m.shape[1:] == (6, 6):
            return cls._from_batch_voigt(m.copy(), tensor_type)
        raise ValueError(f"Expected (6,6) or (N,6,6), got {m.shape}")

    @classmethod
    def from_voigt(cls, v, tensor_type):
        """Alias for from_mat (Tensor4 stores Voigt matrix)."""
        return cls.from_mat(v, tensor_type)

    @classmethod
    def identity(cls, tensor_type=Tensor4Type.stiffness):
        return cls._from_single_cpp(_CppTensor4.identity(tensor_type))

    @classmethod
    def identity2(cls, tensor_type=Tensor4Type.stiffness):
        return cls._from_single_cpp(_CppTensor4.identity2(tensor_type))

    @classmethod
    def volumetric(cls, tensor_type=Tensor4Type.stiffness):
        return cls._from_single_cpp(_CppTensor4.volumetric(tensor_type))

    @classmethod
    def deviatoric(cls, tensor_type=Tensor4Type.stiffness):
        return cls._from_single_cpp(_CppTensor4.deviatoric(tensor_type))

    @classmethod
    def deviatoric2(cls, tensor_type=Tensor4Type.stiffness):
        return cls._from_single_cpp(_CppTensor4.deviatoric2(tensor_type))

    @classmethod
    def zeros(cls, tensor_type=Tensor4Type.stiffness):
        return cls._from_single_cpp(_CppTensor4.zeros(tensor_type))

    @classmethod
    def from_tensor(cls, t, n):
        """Broadcast a single Tensor4 to a batch of size n."""
        if not t._single:
            raise ValueError("from_tensor requires a single tensor")
        v = np.broadcast_to(t.mat[np.newaxis, :, :], (n, 6, 6)).copy()
        return cls._from_batch_voigt(v, t.type)

    @classmethod
    def from_list(cls, tensors):
        """Stack a list of single Tensor4 into a batch."""
        return cls(list(tensors))

    @classmethod
    def concatenate(cls, batches):
        """Join multiple batches into one."""
        parts = list(batches)
        if not parts:
            raise ValueError("Nothing to concatenate")
        t4type = parts[0].type
        arrays = []
        for b in parts:
            if b.type != t4type:
                raise ValueError(f"Mixed Tensor4Type: {t4type} vs {b.type}")
            if b._single:
                arrays.append(b.mat[np.newaxis])
            else:
                arrays.append(b._voigt_data)
        return cls._from_batch_voigt(np.concatenate(arrays, axis=0), t4type)

    # ------------------------------------------------------------------
    # Properties
    # ------------------------------------------------------------------

    @property
    def single(self):
        return self._single

    @property
    def mat(self):
        """6x6 matrix: (6,6) single, (N,6,6) batch."""
        if self._single:
            return self._cpp.mat
        return self._voigt_data

    @property
    def voigt(self):
        """Alias for mat (Tensor4 stores Voigt matrix)."""
        return self.mat

    @property
    def type(self):
        return self._type

    # ------------------------------------------------------------------
    # Sequence protocol (batch only)
    # ------------------------------------------------------------------

    def __len__(self):
        if self._single:
            raise TypeError("single Tensor4 has no len()")
        return self._voigt_data.shape[0]

    def __getitem__(self, key):
        if self._single:
            raise TypeError("single Tensor4 is not subscriptable")
        if isinstance(key, (int, np.integer)):
            if key < 0:
                key += len(self)
            if key < 0 or key >= len(self):
                raise IndexError(f"index {key} out of range for batch of size {len(self)}")
            return Tensor4.from_mat(self._voigt_data[key].copy(), self._type)
        v = self._voigt_data[key]
        if v.ndim == 2:
            return Tensor4.from_mat(v.copy(), self._type)
        return Tensor4._from_batch_voigt(v.copy(), self._type)

    def __iter__(self):
        if self._single:
            raise TypeError("single Tensor4 is not iterable")
        for i in range(len(self)):
            yield self[i]

    def __reversed__(self):
        if self._single:
            raise TypeError("single Tensor4 is not reversible")
        for i in range(len(self) - 1, -1, -1):
            yield self[i]

    # ------------------------------------------------------------------
    # Numpy array protocol
    # ------------------------------------------------------------------

    def __array__(self, dtype=None):
        v = self.mat
        return v.astype(dtype) if dtype else v

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        if method != "__call__":
            return NotImplemented
        raw_inputs = []
        batch_type = None
        for inp in inputs:
            if isinstance(inp, Tensor4):
                if inp._single:
                    raw_inputs.append(inp.mat[np.newaxis])
                else:
                    raw_inputs.append(inp._voigt_data)
                batch_type = batch_type or inp.type
            else:
                raw_inputs.append(inp)
        if batch_type is None:
            return NotImplemented
        _ALLOWED = {np.add, np.subtract, np.multiply, np.negative,
                    np.true_divide, np.positive}
        if ufunc not in _ALLOWED:
            return NotImplemented
        result = ufunc(*raw_inputs, **kwargs)
        if isinstance(result, np.ndarray) and result.ndim == 3 and result.shape[1:] == (6, 6):
            return Tensor4._from_batch_voigt(result, batch_type)
        return result

    # ------------------------------------------------------------------
    # Methods
    # ------------------------------------------------------------------

    def contract(self, t):
        """Contract with Tensor2: result = mat @ t.voigt."""
        if self._single and t._single:
            return Tensor2._from_single_cpp(self._cpp.contract(t._cpp))
        # At least one is batch
        t4_data = self.mat[np.newaxis] if self._single else self._voigt_data
        t2_data = t.voigt[np.newaxis] if t._single else t._voigt_data
        t2_vtype = t._vtype
        result_voigt, out_vtype = _batch_t4_contract(
            t4_data, self._type, t2_data, t2_vtype
        )
        return Tensor2._from_batch_voigt(result_voigt, out_vtype)

    def rotate(self, R, active=True):
        """Rotate tensor(s)."""
        from simcoon.rotation import Rotation as SmcRotation
        from scipy.spatial.transform import Rotation as ScipyRotation

        if self._single:
            if isinstance(R, SmcRotation):
                cpp_rot = R._to_cpp()
            else:
                cpp_rot = R
            return Tensor4._from_single_cpp(self._cpp.rotate(cpp_rot, active))

        N = len(self)
        if isinstance(R, (ScipyRotation, SmcRotation)):
            mats = _get_rotation_matrices(R, N)
        else:
            raise TypeError(f"Expected Rotation, got {type(R)}")

        result = _batch_t4_rotate(self._voigt_data, self._type, mats, active)
        return Tensor4._from_batch_voigt(result, self._type)

    def push_forward(self, F):
        """Push-forward via deformation gradient F."""
        F = np.asarray(F, dtype=np.float64)
        if self._single:
            return Tensor4._from_single_cpp(self._cpp.push_forward(F))
        F_batch = F[np.newaxis] if F.ndim == 2 else F
        result = _batch_t4_push_forward(self._voigt_data, self._type, F_batch)
        return Tensor4._from_batch_voigt(result, self._type)

    def pull_back(self, F):
        """Pull-back via deformation gradient F."""
        F = np.asarray(F, dtype=np.float64)
        if self._single:
            return Tensor4._from_single_cpp(self._cpp.pull_back(F))
        F_batch = F[np.newaxis] if F.ndim == 2 else F
        result = _batch_t4_pull_back(self._voigt_data, self._type, F_batch)
        return Tensor4._from_batch_voigt(result, self._type)

    def inverse(self):
        """Invert the 6x6 Voigt matrix. stiffness <-> compliance."""
        if self._single:
            return Tensor4._from_single_cpp(self._cpp.inverse())
        result, inv_type = _batch_t4_inverse(self._voigt_data, self._type)
        return Tensor4._from_batch_voigt(result, inv_type)

    # ------------------------------------------------------------------
    # Arithmetic
    # ------------------------------------------------------------------

    def __add__(self, other):
        if not isinstance(other, Tensor4):
            return NotImplemented
        if self._single and other._single:
            return Tensor4._from_single_cpp(self._cpp + other._cpp)
        sv = self.mat if self._single else self._voigt_data
        ov = other.mat if other._single else other._voigt_data
        if sv.ndim == 2:
            sv = sv[np.newaxis]
        if ov.ndim == 2:
            ov = ov[np.newaxis]
        return Tensor4._from_batch_voigt(sv + ov, self._type)

    def __radd__(self, other):
        if isinstance(other, Tensor4):
            return other.__add__(self)
        return NotImplemented

    def __sub__(self, other):
        if not isinstance(other, Tensor4):
            return NotImplemented
        if self._single and other._single:
            return Tensor4._from_single_cpp(self._cpp - other._cpp)
        sv = self.mat if self._single else self._voigt_data
        ov = other.mat if other._single else other._voigt_data
        if sv.ndim == 2:
            sv = sv[np.newaxis]
        if ov.ndim == 2:
            ov = ov[np.newaxis]
        return Tensor4._from_batch_voigt(sv - ov, self._type)

    def __rsub__(self, other):
        if isinstance(other, Tensor4):
            return other.__sub__(self)
        return NotImplemented

    def __neg__(self):
        if self._single:
            return Tensor4._from_single_cpp(-self._cpp)
        return Tensor4._from_batch_voigt(-self._voigt_data, self._type)

    def __mul__(self, other):
        if isinstance(other, (int, float)):
            if self._single:
                return Tensor4._from_single_cpp(self._cpp * float(other))
            return Tensor4._from_batch_voigt(self._voigt_data * other, self._type)
        if not self._single:
            other = np.asarray(other)
            if other.ndim <= 1:
                return Tensor4._from_batch_voigt(
                    self._voigt_data * other[..., np.newaxis, np.newaxis], self._type
                )
        return NotImplemented

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        if isinstance(other, (int, float)):
            if self._single:
                return Tensor4._from_single_cpp(self._cpp / float(other))
            return Tensor4._from_batch_voigt(self._voigt_data / other, self._type)
        if not self._single:
            other = np.asarray(other)
            if other.ndim <= 1:
                return Tensor4._from_batch_voigt(
                    self._voigt_data / other[..., np.newaxis, np.newaxis], self._type
                )
        return NotImplemented

    def __matmul__(self, other):
        if isinstance(other, Tensor2):
            return self.contract(other)
        return NotImplemented

    def __eq__(self, other):
        if isinstance(other, Tensor4):
            if self._single and other._single:
                return self._cpp == other._cpp
            sv = self.mat if self._single else self._voigt_data
            ov = other.mat if other._single else other._voigt_data
            if sv.ndim == 2 and ov.ndim == 2:
                return np.array_equal(sv, ov)
            if sv.ndim == 2:
                return np.all(ov == sv[np.newaxis], axis=(1, 2))
            if ov.ndim == 2:
                return np.all(sv == ov[np.newaxis], axis=(1, 2))
            if sv.shape != ov.shape:
                return False
            return np.array_equal(sv, ov)
        return NotImplemented

    def __hash__(self):
        return id(self)

    def __repr__(self):
        if self._single:
            return f"Tensor4(type={self.type.name})"
        return f"Tensor4(N={len(self)}, type={self._type.name})"


# ======================================================================
# Module-level free functions
# ======================================================================

def dyadic(a, b):
    """Dyadic product of two single Tensor2 → Tensor4(stiffness)."""
    return Tensor4._from_single_cpp(_dyadic(a._cpp, b._cpp))


def auto_dyadic(a):
    """Dyadic product of a single Tensor2 with itself → Tensor4(stiffness)."""
    return Tensor4._from_single_cpp(_auto_dyadic(a._cpp))


def double_contract(a, b):
    """Batch double contraction A_ij B_ij → (N,)."""
    ma = a.mat
    mb = b.mat
    if ma.ndim == 2:
        ma = ma[np.newaxis]
    if mb.ndim == 2:
        mb = mb[np.newaxis]
    return np.sum(ma * mb, axis=(-2, -1))


def _mises_impl(t):
    """Compute von Mises equivalent for a single Tensor2."""
    v = t.voigt
    tr = v[0] + v[1] + v[2]
    d = v.copy()
    d[0] -= tr / 3.0
    d[1] -= tr / 3.0
    d[2] -= tr / 3.0

    if t.vtype == VoigtType.stress or t.vtype == VoigtType.generic:
        return np.sqrt(1.5 * (d[0]**2 + d[1]**2 + d[2]**2 +
                              2.0 * (d[3]**2 + d[4]**2 + d[5]**2)))
    elif t.vtype == VoigtType.strain:
        return np.sqrt(2.0/3.0 * (d[0]**2 + d[1]**2 + d[2]**2 +
                                   0.5 * (d[3]**2 + d[4]**2 + d[5]**2)))
    raise ValueError("Mises not defined for VoigtType.none")
