"""
Unified Tensor2 and Tensor4 classes — scipy.Rotation-style, numpy-only storage.

Single tensor: stores ``(6,)`` / ``(6, 6)`` numpy array.
Batch of N tensors: stores ``(N, 6)`` / ``(N, 6, 6)`` numpy array.
Type tags are Python strings: ``"stress"``, ``"strain"``, ``"stiffness"``, etc.

``t.single`` is True for a single tensor, False for a batch.

Examples
--------
>>> import simcoon as smc
>>> import numpy as np

>>> # Single tensor
>>> L = smc.Tensor4.stiffness(smc.L_iso([70000, 0.3], 'Enu'))
>>> eps = smc.Tensor2.strain(np.array([0.01, -0.003, -0.003, 0.005, 0.002, 0.001]))
>>> sigma = L @ eps

>>> # Batch from (N, 6) array
>>> eps_batch = smc.Tensor2.strain(np.random.randn(100, 6) * 0.01)
>>> eps_batch.single  # False
>>> len(eps_batch)     # 100
>>> eps_batch[0]       # single Tensor2

>>> # Batch operations
>>> sigma_batch = L @ eps_batch
"""

import numpy as np

from simcoon._core import (
    _CppTensor2,
    _CppTensor4,
    _dyadic,
    _auto_dyadic,
    VoigtType as _VoigtType,
    Tensor4Type as _T4Type,
)

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
# Type mappings (string <-> C++ enum)
# ======================================================================

_VTYPE_MAP = {
    "stress": _VoigtType.stress,
    "strain": _VoigtType.strain,
    "generic": _VoigtType.generic,
    "none": _VoigtType.none,
}
_VTYPE_RMAP = {v: k for k, v in _VTYPE_MAP.items()}

_T4TYPE_MAP = {
    "stiffness": _T4Type.stiffness,
    "compliance": _T4Type.compliance,
    "strain_concentration": _T4Type.strain_concentration,
    "stress_concentration": _T4Type.stress_concentration,
    "generic": _T4Type.generic,
}
_T4TYPE_RMAP = {v: k for k, v in _T4TYPE_MAP.items()}

_T2_TYPES = frozenset(_VTYPE_MAP)
_T4_TYPES = frozenset(_T4TYPE_MAP)


def _check_t2_type(ts):
    if ts not in _T2_TYPES:
        raise ValueError(f"Invalid Tensor2 type '{ts}', expected one of {sorted(_T2_TYPES)}")


def _check_t4_type(ts):
    if ts not in _T4_TYPES:
        raise ValueError(f"Invalid Tensor4 type '{ts}', expected one of {sorted(_T4_TYPES)}")


# ======================================================================
# Voigt conversion helpers
# ======================================================================

def _mat_to_voigt(m, type_str):
    """Convert (..., 3, 3) matrices to (..., 6) Voigt vectors."""
    v = np.empty((*m.shape[:-2], 6), dtype=np.float64)
    v[..., 0] = m[..., 0, 0]
    v[..., 1] = m[..., 1, 1]
    v[..., 2] = m[..., 2, 2]
    if type_str == "strain":
        v[..., 3] = m[..., 0, 1] + m[..., 1, 0]
        v[..., 4] = m[..., 0, 2] + m[..., 2, 0]
        v[..., 5] = m[..., 1, 2] + m[..., 2, 1]
    else:
        v[..., 3] = 0.5 * (m[..., 0, 1] + m[..., 1, 0])
        v[..., 4] = 0.5 * (m[..., 0, 2] + m[..., 2, 0])
        v[..., 5] = 0.5 * (m[..., 1, 2] + m[..., 2, 1])
    return v


def _voigt_to_mat(v, type_str):
    """Convert (..., 6) Voigt vectors to (..., 3, 3) matrices."""
    m = np.empty((*v.shape[:-1], 3, 3), dtype=np.float64)
    m[..., 0, 0] = v[..., 0]
    m[..., 1, 1] = v[..., 1]
    m[..., 2, 2] = v[..., 2]
    if type_str == "strain":
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
    mats = np.asarray(R.as_matrix(), dtype=np.float64)
    if mats.ndim == 2:
        mats = mats[np.newaxis]
    if mats.shape[0] == 1:
        mats = np.broadcast_to(mats, (N, 3, 3)).copy()
    if mats.shape[0] != N:
        raise ValueError(
            f"Rotation batch size {mats.shape[0]} != tensor batch size {N}"
        )
    return mats


# ======================================================================
# _TensorBase — shared logic for Tensor2 and Tensor4
# ======================================================================

class _TensorBase:
    """Private base class for Tensor2 and Tensor4.

    Parameterized by ``_single_ndim``: 1 for Tensor2, 2 for Tensor4.
    Subclasses set this as a class variable.
    """

    __slots__ = ("_data", "_type_str")

    _single_ndim = None  # set by subclass

    # ------------------------------------------------------------------
    # Internal constructors
    # ------------------------------------------------------------------

    @classmethod
    def _create(cls, data, type_str):
        """Internal: create from numpy array + type string."""
        obj = object.__new__(cls)
        obj._data = np.ascontiguousarray(data, dtype=np.float64)
        obj._type_str = type_str
        return obj

    def _ensure_batch(self):
        """Return _data with a leading batch axis if single."""
        if self._data.ndim == self._single_ndim:
            return self._data[np.newaxis]
        return self._data

    def _rewrap(self, result, type_str=None):
        """Wrap batch result, squeezing if original was single."""
        ts = type_str or self._type_str
        if self._data.ndim == self._single_ndim:
            return type(self)._create(result[0], ts)
        return type(self)._create(result, ts)

    # ------------------------------------------------------------------
    # Properties
    # ------------------------------------------------------------------

    @property
    def single(self):
        """True if this is a single tensor, False if batch."""
        return self._data.ndim == self._single_ndim

    @property
    def type(self):
        """Type string."""
        return self._type_str

    # ------------------------------------------------------------------
    # Sequence protocol
    # ------------------------------------------------------------------

    def __len__(self):
        if self.single:
            return 1
        return self._data.shape[0]

    def __getitem__(self, key):
        if self.single:
            raise TypeError(f"single {type(self).__name__} is not subscriptable")
        if isinstance(key, (int, np.integer)):
            if key < 0:
                key += len(self)
            if key < 0 or key >= len(self):
                raise IndexError(f"index {key} out of range for batch of size {len(self)}")
        return type(self)._create(self._data[key].copy(), self._type_str)

    def __iter__(self):
        if self.single:
            raise TypeError(f"single {type(self).__name__} is not iterable")
        for i in range(len(self)):
            yield self[i]

    def __reversed__(self):
        if self.single:
            raise TypeError(f"single {type(self).__name__} is not reversible")
        for i in range(len(self) - 1, -1, -1):
            yield self[i]

    def __contains__(self, item):
        if self.single:
            raise TypeError(f"single {type(self).__name__} does not support 'in'")
        if isinstance(item, type(self)) and item.single:
            axes = tuple(range(1, self._data.ndim))
            return np.any(np.all(self._data == item._data, axis=axes))
        return False

    def count(self, item):
        """Count occurrences of a single tensor in batch."""
        if self.single:
            raise TypeError("count only on batch")
        if isinstance(item, type(self)) and item.single:
            axes = tuple(range(1, self._data.ndim))
            return int(np.sum(np.all(self._data == item._data, axis=axes)))
        return 0

    # ------------------------------------------------------------------
    # Numpy array protocol
    # ------------------------------------------------------------------

    def __array__(self, dtype=None):
        v = self._data.copy()
        return v.astype(dtype) if dtype else v

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        if method != "__call__":
            return NotImplemented
        cls = type(self)
        raw_inputs = []
        result_ts = None
        all_single = True
        for inp in inputs:
            if isinstance(inp, cls):
                d = inp._data
                if d.ndim == self._single_ndim:
                    d = d[np.newaxis]
                else:
                    all_single = False
                raw_inputs.append(d)
                result_ts = result_ts or inp._type_str
            else:
                raw_inputs.append(inp)
        if result_ts is None:
            return NotImplemented
        _ALLOWED = {np.add, np.subtract, np.multiply, np.negative,
                    np.true_divide, np.positive}
        if ufunc not in _ALLOWED:
            return NotImplemented
        result = ufunc(*raw_inputs, **kwargs)
        batch_ndim = self._single_ndim + 1
        if isinstance(result, np.ndarray) and result.ndim == batch_ndim:
            if all_single and result.shape[0] == 1:
                return cls._create(result[0], result_ts)
            return cls._create(result, result_ts)
        return result

    # ------------------------------------------------------------------
    # Arithmetic
    # ------------------------------------------------------------------

    def __add__(self, other):
        if not isinstance(other, type(self)):
            return NotImplemented
        return type(self)._create(self._data + other._data, self._type_str)

    def __radd__(self, other):
        if isinstance(other, type(self)):
            return other.__add__(self)
        return NotImplemented

    def __sub__(self, other):
        if not isinstance(other, type(self)):
            return NotImplemented
        return type(self)._create(self._data - other._data, self._type_str)

    def __rsub__(self, other):
        if isinstance(other, type(self)):
            return other.__sub__(self)
        return NotImplemented

    def __neg__(self):
        return type(self)._create(-self._data, self._type_str)

    def __mul__(self, other):
        if isinstance(other, (int, float)):
            return type(self)._create(self._data * other, self._type_str)
        other = np.asarray(other)
        if other.ndim <= 1:
            expand = (Ellipsis,) + (np.newaxis,) * self._single_ndim
            return type(self)._create(self._data * other[expand], self._type_str)
        return NotImplemented

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        if isinstance(other, (int, float)):
            return type(self)._create(self._data / other, self._type_str)
        other = np.asarray(other)
        if other.ndim <= 1:
            expand = (Ellipsis,) + (np.newaxis,) * self._single_ndim
            return type(self)._create(self._data / other[expand], self._type_str)
        return NotImplemented

    def __eq__(self, other):
        if isinstance(other, type(self)):
            sd, od = self._data, other._data
            if sd.ndim == od.ndim:
                return np.array_equal(sd, od)
            # mixed single/batch
            axes = tuple(range(1, max(sd.ndim, od.ndim)))
            if sd.ndim < od.ndim:
                return np.all(od == sd[np.newaxis], axis=axes)
            return np.all(sd == od[np.newaxis], axis=axes)
        return NotImplemented

    def __hash__(self):
        return id(self)

    def __repr__(self):
        name = type(self).__name__
        if self.single:
            return f"{name}(type='{self._type_str}')"
        return f"{name}(N={len(self)}, type='{self._type_str}')"

    # ------------------------------------------------------------------
    # Shared factory helpers
    # ------------------------------------------------------------------

    @classmethod
    def from_tensor(cls, t, n):
        """Broadcast a single tensor to a batch of size n."""
        if not t.single:
            raise ValueError("from_tensor requires a single tensor")
        expand = (np.newaxis,) + (slice(None),) * t._data.ndim
        shape = (n,) + t._data.shape
        v = np.broadcast_to(t._data[expand], shape).copy()
        return cls._create(v, t._type_str)

    @classmethod
    def from_list(cls, tensors):
        """Stack a list of single tensors into a batch."""
        return cls(list(tensors))

    @classmethod
    def concatenate(cls, batches):
        """Join multiple batches (or singles) into one batch."""
        parts = list(batches)
        if not parts:
            raise ValueError("Nothing to concatenate")
        type_str = parts[0]._type_str
        sn = parts[0]._single_ndim
        arrays = []
        for b in parts:
            if b._type_str != type_str:
                raise ValueError(f"Mixed type: {type_str} vs {b._type_str}")
            d = b._data
            arrays.append(d[np.newaxis] if d.ndim == sn else d)
        return cls._create(np.concatenate(arrays, axis=0), type_str)


# ======================================================================
# Tensor2
# ======================================================================

class Tensor2(_TensorBase):
    """A 2nd-order tensor with type tag for Voigt convention and rotation dispatch.

    Always stores a numpy array: ``(6,)`` for single, ``(N, 6)`` for batch.
    Type is a string: ``"stress"``, ``"strain"``, ``"generic"``, or ``"none"``.
    """

    _single_ndim = 1

    def __init__(self, data):
        if isinstance(data, list) and data and isinstance(data[0], Tensor2):
            if not all(t.single for t in data):
                raise ValueError("Cannot nest batches")
            ref = data[0]._type_str
            v = np.empty((len(data), 6), dtype=np.float64)
            for i, t in enumerate(data):
                if t._type_str != ref:
                    raise ValueError(f"Mixed type: {ref} vs {t._type_str} at index {i}")
                v[i] = t._data
            self._data = v
            self._type_str = ref
        else:
            raise TypeError("Use Tensor2.stress(), Tensor2.strain(), etc.")

    # ------------------------------------------------------------------
    # Factory methods
    # ------------------------------------------------------------------

    @classmethod
    def stress(cls, data):
        """Create stress tensor(s): (6,), (3,3), (N,6), or (N,3,3)."""
        return cls._from_data(data, "stress")

    @classmethod
    def strain(cls, data):
        """Create strain tensor(s): (6,), (3,3), (N,6), or (N,3,3)."""
        return cls._from_data(data, "strain")

    @classmethod
    def from_mat(cls, m, type_str):
        """Create from (3,3) or (N,3,3) matrix with explicit type string."""
        _check_t2_type(type_str)
        m = np.asarray(m, dtype=np.float64)
        if m.shape == (3, 3):
            return cls._create(_mat_to_voigt(m, type_str).ravel(), type_str)
        if m.ndim == 3 and m.shape[1:] == (3, 3):
            return cls._create(_mat_to_voigt(m, type_str), type_str)
        raise ValueError(f"Expected (3,3) or (N,3,3), got {m.shape}")

    @classmethod
    def from_voigt(cls, v, type_str):
        """Create from Voigt vector (6,) or batch (N,6) with explicit type string."""
        _check_t2_type(type_str)
        v = np.asarray(v, dtype=np.float64)
        if v.ndim == 1 and v.size == 6:
            return cls._create(v.copy(), type_str)
        if v.ndim == 2 and v.shape[1] == 6:
            return cls._create(v.copy(), type_str)
        raise ValueError(f"Expected (6,) or (N,6), got {v.shape}")

    @classmethod
    def zeros(cls, type_str="stress"):
        _check_t2_type(type_str)
        return cls._create(np.zeros(6, dtype=np.float64), type_str)

    @classmethod
    def identity(cls, type_str="stress"):
        _check_t2_type(type_str)
        return cls._create(np.array([1, 1, 1, 0, 0, 0], dtype=np.float64), type_str)

    @classmethod
    def from_columns(cls, arr, type_str):
        """Create batch from (6, N) column-major array (C++ interop)."""
        _check_t2_type(type_str)
        arr = np.asarray(arr, dtype=np.float64)
        if arr.ndim != 2 or arr.shape[0] != 6:
            raise ValueError(f"Expected (6, N), got {arr.shape}")
        return cls._create(arr.T.copy(), type_str)

    # ------------------------------------------------------------------
    # Properties
    # ------------------------------------------------------------------

    @property
    def voigt(self):
        """Voigt vector: (6,) for single, (N,6) for batch. Returns a copy."""
        return self._data.copy()

    @property
    def mat(self):
        """Matrix: (3,3) for single, (N,3,3) for batch."""
        return _voigt_to_mat(self._data, self._type_str)

    @property
    def vtype(self):
        """Alias for type (backward compatibility)."""
        return self._type_str

    @property
    def voigt_T(self):
        """(6, N) transposed Voigt array (batch only, for C++ interop)."""
        if self.single:
            raise AttributeError("voigt_T only available on batch")
        return self._data.T.copy()

    # ------------------------------------------------------------------
    # Domain methods
    # ------------------------------------------------------------------

    def is_symmetric(self, tol=1e-12):
        """Check symmetry of the 3x3 matrix (single only)."""
        if not self.single:
            raise NotImplementedError("is_symmetric not supported on batch")
        m = self.mat
        return (abs(m[0, 1] - m[1, 0]) < tol and
                abs(m[0, 2] - m[2, 0]) < tol and
                abs(m[1, 2] - m[2, 1]) < tol)

    def rotate(self, R, active=True):
        """Rotate tensor(s) via Rotation."""
        from simcoon.rotation import Rotation as SmcRotation
        from scipy.spatial.transform import Rotation as ScipyRotation

        data_2d = self._ensure_batch()
        N = data_2d.shape[0]

        if isinstance(R, (ScipyRotation, SmcRotation)):
            mats = _get_rotation_matrices(R, N)
        else:
            raise TypeError(f"Expected Rotation, got {type(R)}")

        result = _batch_t2_rotate(data_2d, _VTYPE_MAP[self._type_str], mats, active)
        return self._rewrap(result)

    def push_forward(self, F, metric=True):
        """Push-forward via deformation gradient F."""
        F = np.asarray(F, dtype=np.float64)
        data_2d = self._ensure_batch()
        F_batch = F[np.newaxis] if F.ndim == 2 else F
        result = _batch_t2_push_forward(
            data_2d, _VTYPE_MAP[self._type_str], F_batch, metric)
        return self._rewrap(result)

    def pull_back(self, F, metric=True):
        """Pull-back via deformation gradient F."""
        F = np.asarray(F, dtype=np.float64)
        data_2d = self._ensure_batch()
        F_batch = F[np.newaxis] if F.ndim == 2 else F
        result = _batch_t2_pull_back(
            data_2d, _VTYPE_MAP[self._type_str], F_batch, metric)
        return self._rewrap(result)

    def mises(self):
        """Von Mises equivalent: float for single, (N,) for batch."""
        d = self._data.copy()
        tr = d[..., 0] + d[..., 1] + d[..., 2]
        d[..., 0] -= tr / 3.0
        d[..., 1] -= tr / 3.0
        d[..., 2] -= tr / 3.0
        diag2 = d[..., 0]**2 + d[..., 1]**2 + d[..., 2]**2
        shear2 = d[..., 3]**2 + d[..., 4]**2 + d[..., 5]**2
        if self._type_str == "strain":
            return np.sqrt((2.0 / 3.0) * (diag2 + 0.5 * shear2))
        if self._type_str in ("stress", "generic"):
            return np.sqrt(1.5 * (diag2 + 2.0 * shear2))
        raise ValueError("Mises not defined for type 'none'")

    def trace(self):
        """Trace: float for single, (N,) for batch."""
        return self._data[..., 0] + self._data[..., 1] + self._data[..., 2]

    def dev(self):
        """Deviatoric part -> Tensor2 (single or batch)."""
        d = self._data.copy()
        tr = d[..., 0] + d[..., 1] + d[..., 2]
        d[..., 0] -= tr / 3.0
        d[..., 1] -= tr / 3.0
        d[..., 2] -= tr / 3.0
        return Tensor2._create(d, self._type_str)

    def norm(self):
        """Frobenius norm: float for single, (N,) for batch."""
        d = self._data
        diag2 = d[..., 0]**2 + d[..., 1]**2 + d[..., 2]**2
        shear2 = d[..., 3]**2 + d[..., 4]**2 + d[..., 5]**2
        if self._type_str == "strain":
            return np.sqrt(diag2 + 0.5 * shear2)
        return np.sqrt(diag2 + 2.0 * shear2)

    def __mod__(self, other):
        """Double contraction A_ij B_ij -> float or (N,)."""
        if not isinstance(other, Tensor2):
            return NotImplemented
        ma = self.mat
        mb = other.mat
        if ma.ndim == 2 and mb.ndim == 2:
            return np.sum(ma * mb)
        if ma.ndim == 2:
            ma = ma[np.newaxis]
        if mb.ndim == 2:
            mb = mb[np.newaxis]
        return np.sum(ma * mb, axis=(-2, -1))

    # ------------------------------------------------------------------
    # Internal
    # ------------------------------------------------------------------

    @classmethod
    def _from_data(cls, data, type_str):
        _check_t2_type(type_str)
        data = np.asarray(data, dtype=np.float64)
        if data.shape == (3, 3):
            return cls._create(_mat_to_voigt(data, type_str).ravel(), type_str)
        if data.shape == (6,) or (data.ndim == 1 and data.size == 6):
            return cls._create(data.copy(), type_str)
        if data.ndim == 2 and data.shape[1] == 6:
            return cls._create(data.copy(), type_str)
        if data.ndim == 3 and data.shape[1:] == (3, 3):
            return cls._create(_mat_to_voigt(data, type_str), type_str)
        raise ValueError(
            f"Expected (6,), (3,3), (N,6), or (N,3,3), got {data.shape}"
        )


# ======================================================================
# Tensor4
# ======================================================================

class Tensor4(_TensorBase):
    """A 4th-order tensor with type tag for rotation dispatch and Voigt storage.

    Always stores numpy array: ``(6,6)`` for single, ``(N,6,6)`` for batch.
    Type is a string: ``"stiffness"``, ``"compliance"``, etc.
    """

    _single_ndim = 2

    def __init__(self, data):
        if isinstance(data, list) and data and isinstance(data[0], Tensor4):
            if not all(t.single for t in data):
                raise ValueError("Cannot nest batches")
            ref = data[0]._type_str
            v = np.empty((len(data), 6, 6), dtype=np.float64)
            for i, t in enumerate(data):
                if t._type_str != ref:
                    raise ValueError(f"Mixed type: {ref} vs {t._type_str} at index {i}")
                v[i] = t._data
            self._data = v
            self._type_str = ref
        else:
            raise TypeError("Use Tensor4.stiffness(), Tensor4.compliance(), etc.")

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _to_cpp(self):
        """Create a temporary _CppTensor4 for single-point C++ operations."""
        return _CppTensor4.from_mat(self._data, _T4TYPE_MAP[self._type_str])

    # ------------------------------------------------------------------
    # Factory methods
    # ------------------------------------------------------------------

    @classmethod
    def _typed_factory(cls, data, type_str):
        _check_t4_type(type_str)
        data = np.asarray(data, dtype=np.float64)
        if data.shape == (6, 6):
            return cls._create(data.copy(), type_str)
        if data.ndim == 3 and data.shape[1:] == (6, 6):
            return cls._create(data.copy(), type_str)
        raise ValueError(f"Expected (6,6) or (N,6,6), got {data.shape}")

    @classmethod
    def stiffness(cls, data):
        """Create stiffness tensor(s): (6,6) single or (N,6,6) batch."""
        return cls._typed_factory(data, "stiffness")

    @classmethod
    def compliance(cls, data):
        """Create compliance tensor(s): (6,6) single or (N,6,6) batch."""
        return cls._typed_factory(data, "compliance")

    @classmethod
    def strain_concentration(cls, data):
        """Create strain concentration tensor(s)."""
        return cls._typed_factory(data, "strain_concentration")

    @classmethod
    def stress_concentration(cls, data):
        """Create stress concentration tensor(s)."""
        return cls._typed_factory(data, "stress_concentration")

    @classmethod
    def from_mat(cls, m, type_str):
        """Create from (6,6) or (N,6,6) with explicit type string."""
        return cls._typed_factory(m, type_str)

    @classmethod
    def from_voigt(cls, v, type_str):
        """Alias for from_mat."""
        return cls.from_mat(v, type_str)

    @classmethod
    def identity(cls, type_str="stiffness"):
        _check_t4_type(type_str)
        cpp = _CppTensor4.identity(_T4TYPE_MAP[type_str])
        return cls._create(np.array(cpp.mat), type_str)

    @classmethod
    def identity2(cls, type_str="stiffness"):
        _check_t4_type(type_str)
        cpp = _CppTensor4.identity2(_T4TYPE_MAP[type_str])
        return cls._create(np.array(cpp.mat), type_str)

    @classmethod
    def volumetric(cls, type_str="stiffness"):
        _check_t4_type(type_str)
        cpp = _CppTensor4.volumetric(_T4TYPE_MAP[type_str])
        return cls._create(np.array(cpp.mat), type_str)

    @classmethod
    def deviatoric(cls, type_str="stiffness"):
        _check_t4_type(type_str)
        cpp = _CppTensor4.deviatoric(_T4TYPE_MAP[type_str])
        return cls._create(np.array(cpp.mat), type_str)

    @classmethod
    def deviatoric2(cls, type_str="stiffness"):
        _check_t4_type(type_str)
        cpp = _CppTensor4.deviatoric2(_T4TYPE_MAP[type_str])
        return cls._create(np.array(cpp.mat), type_str)

    @classmethod
    def zeros(cls, type_str="stiffness"):
        _check_t4_type(type_str)
        return cls._create(np.zeros((6, 6), dtype=np.float64), type_str)

    # ------------------------------------------------------------------
    # Properties
    # ------------------------------------------------------------------

    @property
    def mat(self):
        """6x6 matrix: (6,6) single, (N,6,6) batch. Returns a copy."""
        return self._data.copy()

    @property
    def voigt(self):
        """Alias for mat."""
        return self.mat

    # ------------------------------------------------------------------
    # Domain methods
    # ------------------------------------------------------------------

    @staticmethod
    def _infer_contraction_vtype(type_str):
        if type_str in ("stiffness", "stress_concentration", "generic"):
            return "stress"
        return "strain"

    def contract(self, t):
        """Contract with Tensor2: result = mat @ t.voigt."""
        out_ts = self._infer_contraction_vtype(self._type_str)

        if self.single:
            result = (self._data @ t._data.T).T
            if t.single:
                return Tensor2._create(result.ravel(), out_ts)
            return Tensor2._create(result, out_ts)

        t2 = t._data[np.newaxis] if t.single else t._data
        t2 = np.broadcast_to(t2, (self._data.shape[0], 6))
        result = np.matmul(self._data, t2[..., np.newaxis]).squeeze(-1)
        return Tensor2._create(result, out_ts)

    def rotate(self, R, active=True):
        """Rotate tensor(s)."""
        from simcoon.rotation import Rotation as SmcRotation
        from scipy.spatial.transform import Rotation as ScipyRotation

        if self.single:
            cpp_t4 = self._to_cpp()
            if isinstance(R, SmcRotation):
                cpp_rot = R._to_cpp()
            else:
                cpp_rot = R
            cpp_result = cpp_t4.rotate(cpp_rot, active)
            return Tensor4._create(np.array(cpp_result.mat),
                                   _T4TYPE_RMAP.get(cpp_result.type, self._type_str))

        N = len(self)
        if isinstance(R, (ScipyRotation, SmcRotation)):
            mats = _get_rotation_matrices(R, N)
        else:
            raise TypeError(f"Expected Rotation, got {type(R)}")
        result = _batch_t4_rotate(self._data, _T4TYPE_MAP[self._type_str], mats, active)
        return self._rewrap(result)

    def push_forward(self, F, metric=True):
        """Push-forward via deformation gradient F."""
        F = np.asarray(F, dtype=np.float64)
        if self.single:
            cpp_result = self._to_cpp().push_forward(F, metric)
            return Tensor4._create(np.array(cpp_result.mat), self._type_str)
        F_batch = F[np.newaxis] if F.ndim == 2 else F
        result = _batch_t4_push_forward(
            self._data, _T4TYPE_MAP[self._type_str], F_batch, metric)
        return self._rewrap(result)

    def pull_back(self, F, metric=True):
        """Pull-back via deformation gradient F."""
        F = np.asarray(F, dtype=np.float64)
        if self.single:
            cpp_result = self._to_cpp().pull_back(F, metric)
            return Tensor4._create(np.array(cpp_result.mat), self._type_str)
        F_batch = F[np.newaxis] if F.ndim == 2 else F
        result = _batch_t4_pull_back(
            self._data, _T4TYPE_MAP[self._type_str], F_batch, metric)
        return self._rewrap(result)

    def inverse(self):
        """Invert the 6x6 Voigt matrix. stiffness <-> compliance."""
        if self.single:
            cpp_result = self._to_cpp().inverse()
            inv_ts = _T4TYPE_RMAP.get(cpp_result.type, self._type_str)
            return Tensor4._create(np.array(cpp_result.mat), inv_ts)
        result, inv_type_enum = _batch_t4_inverse(
            self._data, _T4TYPE_MAP[self._type_str])
        inv_ts = _T4TYPE_RMAP.get(inv_type_enum, self._type_str)
        return Tensor4._create(result, inv_ts)

    # ------------------------------------------------------------------
    # Arithmetic overrides (Tensor4 * Tensor2 = contraction)
    # ------------------------------------------------------------------

    def __mul__(self, other):
        if isinstance(other, Tensor2):
            return self.contract(other)
        return super().__mul__(other)

    def __matmul__(self, other):
        if isinstance(other, Tensor2):
            return self.contract(other)
        return NotImplemented


# ======================================================================
# Module-level free functions
# ======================================================================

def dyadic(a, b):
    """Dyadic product of two single Tensor2 -> Tensor4(stiffness)."""
    if not (isinstance(a, Tensor2) and a.single and isinstance(b, Tensor2) and b.single):
        raise ValueError("dyadic requires two single Tensor2 arguments")
    cpp_a = _CppTensor2.from_voigt(a._data, _VTYPE_MAP[a._type_str])
    cpp_b = _CppTensor2.from_voigt(b._data, _VTYPE_MAP[b._type_str])
    result = _dyadic(cpp_a, cpp_b)
    return Tensor4._create(np.array(result.mat),
                           _T4TYPE_RMAP.get(result.type, "stiffness"))


def auto_dyadic(a):
    """Dyadic product of a single Tensor2 with itself -> Tensor4(stiffness)."""
    if not (isinstance(a, Tensor2) and a.single):
        raise ValueError("auto_dyadic requires a single Tensor2")
    cpp_a = _CppTensor2.from_voigt(a._data, _VTYPE_MAP[a._type_str])
    result = _auto_dyadic(cpp_a)
    return Tensor4._create(np.array(result.mat),
                           _T4TYPE_RMAP.get(result.type, "stiffness"))


def double_contract(a, b):
    """Batch double contraction A_ij B_ij -> float or (N,)."""
    ma = a.mat
    mb = b.mat
    if ma.ndim == 2:
        ma = ma[np.newaxis]
    if mb.ndim == 2:
        mb = mb[np.newaxis]
    return np.sum(ma * mb, axis=(-2, -1))
