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

>>> # Concentration tensors (strain/stress) use the engineering Voigt convention,
>>> # exactly like stiffness/compliance: wrap a 6x6 producer and contract it.
>>> # The identity is eye(6) (it returns the field unchanged), and inverse keeps
>>> # the concentration type. push_forward/pull_back are undefined (mixed indices).
>>> F = np.array([[1.2, 0.15, 0.0], [0.0, 0.95, 0.0], [0.0, 0.0, 1.0 / (1.2 * 0.95)]])
>>> D = np.array([[0.03, 0.01, 0.0], [0.01, -0.02, 0.0], [0.0, 0.0, -0.01]])
>>> A = smc.Tensor4.strain_concentration(smc.A_R(F))   # De = A^R : D
>>> De = A.contract(smc.Tensor2.strain(smc.t2v_strain(D)))
>>> D0 = A.inverse().contract(De)                      # round-trips: D0 == D
"""

import numpy as np

from simcoon._core import (
    _CppTensor2,
    _CppTensor4,
    _dyadic,
    _auto_dyadic,
    _sym_dyadic,
    _auto_sym_dyadic,
    Tensor2Type as _Tensor2Type,
    Tensor4Type as _T4Type,
)

from simcoon._core import (
    _batch_rotate,
    _batch_push_forward,
    _batch_pull_back,
    _batch_mises,
    _batch_trace,
    _batch_contract,
    _batch_inverse,
)


# ======================================================================
# Type mappings (string <-> C++ enum)
# ======================================================================

_VTYPE_MAP = {
    "stress": _Tensor2Type.stress,
    "strain": _Tensor2Type.strain,
    "symmetric": _Tensor2Type.symmetric,
    "none": _Tensor2Type.none,      # no Voigt type: any 3x3, 9 components
}
_VTYPE_RMAP = {v: k for k, v in _VTYPE_MAP.items()}
_T2_ALIASES = {"generic": "symmetric"}   # pre-2.2 name

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
    """Validate a Tensor2 type string and return its canonical form."""
    ts = _T2_ALIASES.get(ts, ts)
    if ts not in _T2_TYPES:
        raise ValueError(f"Invalid Tensor2 type '{ts}', expected one of {sorted(_T2_TYPES)}")
    return ts


def _t2_ncomp(ts):
    """Stored components: 6 (Voigt, symmetric) or 9 ("none": any 3x3, row-major)."""
    return 9 if ts == "none" else 6


def _require_voigt_type(ts):
    if ts == "none":
        raise ValueError("A Tensor2 of type 'none' (any 3x3, no Voigt convention) has no Voigt "
                         "or Mandel vector: use mat")


def _check_t4_type(ts):
    if ts not in _T4_TYPES:
        raise ValueError(f"Invalid Tensor4 type '{ts}', expected one of {sorted(_T4_TYPES)}")


# ======================================================================
# Kelvin-Mandel congruence factors (mirror C++ mandel_factors in tensor.hpp)
# ======================================================================

_SQ2 = float(np.sqrt(2.0))


def _t2_mandel_factor(type_str):
    """eng->Mandel shear factor for a Voigt 6-vector: mandel = sqrt2 * t_ij."""
    return (1.0 / _SQ2) if type_str == "strain" else _SQ2


def _t4_mandel_factors(type_str):
    """eng->Mandel (row, col) shear factors for a 6x6 (Mandel->eng divides)."""
    if type_str in ("stiffness", "generic"):
        return _SQ2, _SQ2
    if type_str == "compliance":
        return 1.0 / _SQ2, 1.0 / _SQ2
    if type_str == "strain_concentration":
        return 1.0 / _SQ2, _SQ2
    return _SQ2, 1.0 / _SQ2  # stress_concentration


def _to_f_cube(arr):
    """Convert (N,R,C) C-order array to (R,C,N) F-order for zero-copy arma cube."""
    return np.asfortranarray(arr.transpose(1, 2, 0))


def _from_f_cube(arr):
    """Convert (R,C,N) F-order array from arma cube to (N,R,C) C-order."""
    return np.ascontiguousarray(arr.transpose(2, 0, 1))


# ======================================================================
# Voigt conversion helpers
# ======================================================================

def _require_symmetric(m, rtol=1e-10):
    """Tensor2 stores 6 components: a non-symmetric 3x3 would be symmetrised silently."""
    skew = np.abs(m - np.swapaxes(m, -1, -2)).max()
    if skew > rtol * max(np.abs(m).max(), 1.0):
        raise ValueError(
            f"Tensor2 stores symmetric tensors only (6 components); the input is not symmetric "
            f"(max |m - m^T| = {skew:.3g}). Non-symmetric or two-point tensors such as F, R, L "
            "or PK1 take the 9-component type 'none'")


def _mat_to_voigt(m, type_str):
    """Convert (..., 3, 3) matrices to the stored components: (..., 6) Voigt vectors for
    the symmetric types, the row-major (..., 9) for type "none"."""
    if type_str == "none":
        return np.ascontiguousarray(m, dtype=np.float64).reshape(*m.shape[:-2], 9)
    _require_symmetric(m)
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
    """Convert the stored components to (..., 3, 3) matrices."""
    if type_str == "none":
        return v.reshape(*v.shape[:-1], 3, 3).copy()
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


def _to_cpp_rotation(R):
    """Single-tensor rotation argument -> C++ Rotation. Accepts a simcoon or scipy Rotation,
    raises TypeError otherwise. (Batch paths use _get_rotation_matrices instead.)"""
    from simcoon.rotation import Rotation as SmcRotation
    from scipy.spatial.transform import Rotation as ScipyRotation
    if isinstance(R, SmcRotation):
        return R._to_cpp()
    if isinstance(R, ScipyRotation):
        return SmcRotation.from_scipy(R)._to_cpp()
    raise TypeError(f"Expected Rotation, got {type(R)}")


def _get_rotation_matrices(R, N):
    """Extract rotation matrices as (3,3,N) F-order for zero-copy to arma cube."""
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
    # Transpose to (3,3,N) F-order: arr_to_cube then borrows it without a copy
    return np.asfortranarray(mats.transpose(1, 2, 0))


# ======================================================================
# Basis — the reference system the components are written in
# ======================================================================

def _as_smc_rotation(R):
    """simcoon or scipy Rotation -> simcoon Rotation (TypeError otherwise)."""
    from simcoon.rotation import Rotation as SmcRotation
    from scipy.spatial.transform import Rotation as ScipyRotation
    if isinstance(R, SmcRotation):
        return R
    if isinstance(R, ScipyRotation):
        return SmcRotation.from_scipy(R)
    raise TypeError(f"Expected Rotation, got {type(R)}")


class Basis:
    r"""The basis the components of a `Tensor2` / `Tensor4` are written in.

    A tensor whose ``basis`` is ``None`` is expressed in the fixed orthonormal
    lab basis :math:`\mathbf{e}_i` (the default, at no cost). A ``Basis`` holds
    the basis vectors :math:`\mathbf{g}_i` through the matrix
    :math:`\mathbf{A}` whose columns are their lab components. Two kinds:

    * **orthonormal** -- built from a rotation, :math:`\mathbf{g}_i =
      \mathbf{R}\,\mathbf{e}_i`. The metric is the identity: every formula of
      the tensor classes holds as in the lab.
    * **natural** -- any three independent vectors, typically the convected
      basis :math:`\mathbf{g}_i = \mathbf{F}\,\mathbf{G}_i` (`from_F`). The
      metric :math:`g_{ij} = \mathbf{g}_i \cdot \mathbf{g}_j`, i.e.
      :math:`\mathbf{g} = \mathbf{A}^T \mathbf{A}`, enters the trace, the
      deviator, the norms and the invariants.

    The variance of the components follows the tensor type: stress is
    contravariant (:math:`\boldsymbol{\sigma} = \sigma^{ij}\,\mathbf{g}_i
    \otimes \mathbf{g}_j`), strain covariant (:math:`\boldsymbol{\varepsilon} =
    \varepsilon_{ij}\,\mathbf{g}^i \otimes \mathbf{g}^j`), stiffness and
    compliance likewise on their four indices, concentration tensors mixed.
    Lab components are :math:`\mathbf{A}\,\hat{\boldsymbol{\sigma}}\,
    \mathbf{A}^T` and :math:`\mathbf{A}^{-T}\,\hat{\boldsymbol{\varepsilon}}\,
    \mathbf{A}^{-1}`.

    One basis is either single (shared by every tensor of a batch, stored
    once) or a batch of N bases, one per tensor. A basis is immutable; tensors
    derived from one another share it by reference.

    The variance is a convention of the components, not a property of the
    tensor: it is a tag (``t.variance``) defaulted from the type, and
    ``t.to_variance(v)`` gives the components of the other variance (with the
    metric in a natural basis). Two-point tensors such as :math:`\mathbf{F}`
    have one index in each configuration and no single basis: ``F`` is always a
    plain array of lab components, consumed by `from_F` and by the transports.

    Parameters
    ----------
    rotation : simcoon.Rotation or scipy.spatial.transform.Rotation, optional
        Orthonormal basis :math:`\mathbf{g}_i = \mathbf{R}\,\mathbf{e}_i`
        (single or batch).
    vectors : array_like, optional
        Natural basis: ``(3,3)`` matrix (or ``(N,3,3)`` batch) whose columns
        are the lab components of the basis vectors. The array is kept by
        reference, not copied: do not modify it afterwards.
    name : str, optional
        A label shown in ``repr`` and in error messages. It is never used to
        compare two bases.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        material = sim.Basis(rotation=sim.Rotation.from_euler('z', 30, degrees=True),
                             name="material")
        sigma = sim.Tensor2.stress(np.array([100., 0., 0., 0., 0., 0.]))
        sigma_m = sigma.to_basis(material)      # same tensor, material components
        sigma_m.to_basis(None)                  # back to lab components

        F = np.array([[1., 0.5, 0.], [0., 1., 0.], [0., 0., 1.]])
        S = sim.Tensor2.stress(np.array([0., 80., 0., 0., 0., 0.]))   # PK2
        tau = S.with_basis(sim.Basis.from_F(F))  # Kirchhoff: same components, convected basis
        tau.trace()                              # uses the metric
    """

    __slots__ = ("_rotation", "_A", "_name", "_cache")

    def __init__(self, rotation=None, vectors=None, name=None):
        if (rotation is None) == (vectors is None):
            raise ValueError("Basis needs exactly one of rotation= or vectors=")
        if rotation is not None:
            self._rotation = _as_smc_rotation(rotation)
            self._A = np.ascontiguousarray(self._rotation.as_matrix(), dtype=np.float64)
        else:
            A = np.asarray(vectors, dtype=np.float64)
            if A.shape != (3, 3) and not (A.ndim == 3 and A.shape[1:] == (3, 3)):
                raise ValueError(f"Expected (3,3) or (N,3,3) basis vectors, got {A.shape}")
            self._rotation = None
            self._A = A
        self._name = name
        self._cache = {}

    @classmethod
    def from_F(cls, F, name=None):
        r"""Convected basis :math:`\mathbf{g}_i = \mathbf{F}\,\mathbf{e}_i` of a deformation gradient.

        Parameters
        ----------
        F : array_like
            ``(3,3)`` deformation gradient or ``(N,3,3)`` batch, kept by
            reference. Its metric is the right Cauchy-Green tensor
            :math:`\mathbf{C} = \mathbf{F}^T \mathbf{F}`.
        name : str, optional
            Label of the basis.

        Returns
        -------
        Basis
            Natural basis.
        """
        return cls(vectors=F, name=name)

    # -- description ------------------------------------------------------

    @property
    def name(self):
        """Label of the basis (``None`` if unnamed)."""
        return self._name

    @property
    def orthonormal(self):
        """True for a basis built from a rotation (identity metric)."""
        return self._rotation is not None

    @property
    def single(self):
        """True for one basis, False for a batch of N bases."""
        return self._A.ndim == 2

    def __len__(self):
        if self.single:
            raise TypeError("single Basis has no len()")
        return self._A.shape[0]

    @property
    def rotation(self):
        """The rotation of an orthonormal basis, ``None`` for a natural basis."""
        return self._rotation

    @property
    def matrix(self):
        r"""Basis matrix :math:`\mathbf{A}` (columns = lab components of the basis vectors):
        ``(3,3)`` or ``(N,3,3)``. Returns a copy."""
        return self._A.copy()

    @property
    def metric(self):
        r"""Covariant metric :math:`g_{ij} = \mathbf{g}_i \cdot \mathbf{g}_j`: ``(3,3)`` or ``(N,3,3)``."""
        return self._metric().copy()

    @property
    def inverse_metric(self):
        r"""Contravariant metric :math:`g^{ij}`: ``(3,3)`` or ``(N,3,3)``."""
        return self._inverse_metric().copy()

    @property
    def det(self):
        r"""Volume of the basis, :math:`\det \mathbf{A} = \sqrt{\det \mathbf{g}}`: float or ``(N,)``."""
        return self._lazy("det", lambda: np.linalg.det(self._A))

    @property
    def reciprocal(self):
        r"""Reciprocal (dual) basis vectors :math:`\mathbf{g}^i`, :math:`\mathbf{g}^i \cdot
        \mathbf{g}_j = \delta^i_j`, as the columns of :math:`\mathbf{A}^{-T}`: ``(3,3)`` or
        ``(N,3,3)``. Covariant components live on them. Returns a copy."""
        return self._reciprocal().copy()

    @property
    def stretch(self):
        r"""Right stretch of the basis, :math:`\mathbf{U} = \sqrt{\mathbf{g}}` (the identity for an
        orthonormal basis): ``(3,3)`` or ``(N,3,3)``. For ``from_F`` it is the right stretch of
        :math:`\mathbf{F}`. Returns a copy."""
        return self._stretch().copy()

    @property
    def polar(self):
        r"""The orthonormal basis closest to this one: the rotation :math:`\mathbf{R} =
        \mathbf{A}\,\mathbf{U}^{-1}` of the polar decomposition :math:`\mathbf{A} = \mathbf{R}\,
        \mathbf{U}`, as a `Basis` (``self`` when already orthonormal). For ``from_F`` it is the
        frame turned by the rotation of :math:`\mathbf{F}`."""
        if self.orthonormal:
            return self
        return self._lazy("polar", self._polar)

    def equals(self, other, tol=1e-12):
        """True if ``other`` is the same basis within ``tol`` (the name is ignored)."""
        if other is self:
            return True
        if not isinstance(other, Basis) or self._A.shape != other._A.shape:
            return False
        return bool(np.allclose(self._A, other._A, rtol=0.0, atol=tol))

    def __repr__(self):
        parts = ["orthonormal" if self.orthonormal else "natural"]
        if not self.single:
            parts.append(f"N={len(self)}")
        if self._name is not None:
            parts.append(f"name={self._name!r}")
        return f"Basis({', '.join(parts)})"

    # -- internal ---------------------------------------------------------

    def _lazy(self, key, compute):
        value = self._cache.get(key)
        if value is None:
            value = self._cache[key] = compute()
        return value

    def _metric(self):
        return self._lazy("g", lambda: np.swapaxes(self._A, -1, -2) @ self._A)

    def _inverse_metric(self):
        return self._lazy("ginv", lambda: np.linalg.inv(self._metric()))

    def _reciprocal(self):
        return self._lazy("recip", lambda: np.swapaxes(np.linalg.inv(self._A), -1, -2))

    def _metric_power(self, p):
        """g^p through the eigen-decomposition of the (symmetric positive) metric."""
        w, Q = np.linalg.eigh(self._metric())
        return (Q * (w ** p)[..., np.newaxis, :]) @ np.swapaxes(Q, -1, -2)

    def _stretch(self):
        if self.orthonormal:
            return self._lazy("U", lambda: np.broadcast_to(np.eye(3), self._A.shape).copy())
        return self._lazy("U", lambda: self._metric_power(0.5))

    def _polar(self):
        from simcoon.rotation import Rotation as SmcRotation
        R = self._A @ self._metric_power(-0.5)          # A U^-1
        return Basis(rotation=SmcRotation.from_matrix(R), name=self._name)

    def _take(self, key):
        """Basis of ``tensor[key]``: shared when single, indexed when per tensor."""
        if self.single:
            return self
        if self._rotation is not None:
            return Basis(rotation=self._rotation[key], name=self._name)
        return Basis(vectors=self._A[key], name=self._name)

    def _post_rotated(self, R):
        """Basis A.Q: the frame turned by R, R being read in this basis (passive rotation)."""
        R = _as_smc_rotation(R)
        if self._rotation is not None:
            return Basis(rotation=self._rotation * R)
        return Basis(vectors=self._A @ R.as_matrix())

    def _pre_rotated(self, R):
        """Basis Q.A: the basis vectors turned by the lab rotation R (active rotation)."""
        R = _as_smc_rotation(R)
        if self._rotation is not None:
            return Basis(rotation=R * self._rotation)
        return Basis(vectors=R.as_matrix() @ self._A)

    def _convected(self, F, forward):
        """Natural basis F.A (push-forward) or F^-1.A (pull-back)."""
        if forward:
            return Basis(vectors=F @ self._A)
        return Basis(vectors=np.linalg.solve(F, self._A))


def _same_basis(a, b):
    """True if two tensor bases (None = lab) are the same reference system."""
    if a is b:
        return True
    if a is None or b is None:
        return False
    return a.equals(b)


def _common_basis(a, b):
    """The basis shared by two operands, or ValueError("Mixed basis")."""
    if a is b or _same_basis(a, b):
        return a
    raise ValueError(f"Mixed basis: {a if a is not None else 'lab'} vs "
                     f"{b if b is not None else 'lab'} (use to_basis() first)")


def _is_orthonormal(basis):
    """Lab (None) or an orthonormal Basis: identity metric."""
    return basis is None or basis._rotation is not None


def _merge_bases(bases, counts):
    """Basis of a batch stacked from parts carrying ``bases`` (``counts`` tensors each).

    The same basis everywhere is kept as one shared basis; different bases are
    stacked into a per-tensor basis; lab and non-lab parts cannot be mixed.
    """
    first = bases[0]
    if all(b is first for b in bases):
        if first is None or first.single:
            return first
    if any(b is None for b in bases):
        if all(b is None for b in bases):
            return None
        raise ValueError("Mixed basis: lab and non-lab tensors (use to_basis() first)")
    if first.single and all(b.single and first.equals(b) for b in bases):
        return first
    mats = [np.broadcast_to(b._A, (n, 3, 3)) for b, n in zip(bases, counts)]
    A = np.ascontiguousarray(np.concatenate(mats, axis=0))
    if all(b.orthonormal for b in bases):
        from simcoon.rotation import Rotation as SmcRotation
        return Basis(rotation=SmcRotation.from_matrix(A))
    return Basis(vectors=A)


# Voigt index pairs: I -> (i, j)
_VI = np.array([0, 1, 2, 0, 0, 1])
_VJ = np.array([0, 1, 2, 1, 2, 2])


def _voigt_operators(M):
    r"""6x6 engineering-Voigt operators of a change of basis.

    ``M`` (``(3,3)`` or ``(N,3,3)``) holds the old basis vectors in the new
    basis. Returns ``(P_sharp, P_flat)``: contravariant components transform as
    ``M X M^T`` (``P_sharp`` on a stress Voigt vector), covariant ones as
    ``M^-T X M^-1`` (``P_flat`` on an engineering strain Voigt vector). They
    are dual: ``P_sharp^T P_flat = I``.
    """
    def pair(B):
        direct = B[..., _VI[:, None], _VI[None, :]] * B[..., _VJ[:, None], _VJ[None, :]]
        cross = B[..., _VI[:, None], _VJ[None, :]] * B[..., _VJ[:, None], _VI[None, :]]
        return direct, cross

    direct, cross = pair(M)
    P_sharp = direct.copy()
    P_sharp[..., :, 3:] += cross[..., :, 3:]

    direct, cross = pair(np.swapaxes(np.linalg.inv(M), -1, -2))
    P_flat = direct.copy()
    P_flat[..., :, 3:] = 0.5 * (direct[..., :, 3:] + cross[..., :, 3:])
    P_flat[..., 3:, :] *= 2.0
    return P_sharp, P_flat


CONTRAVARIANT, COVARIANT = "contravariant", "covariant"
_DUAL = {CONTRAVARIANT: COVARIANT, COVARIANT: CONTRAVARIANT}
_DEFAULT = object()     # sentinel: "the default variance of the type"


def _check_variance(v):
    if v not in _DUAL:
        raise ValueError(f"variance must be 'contravariant' or 'covariant', got {v!r}")
    return v


def _congruence_operator(B, tensorial):
    """6x6 operator of X -> B X B^T on Voigt components with tensorial (stress-like) or
    engineering (strain-like, doubled shear) shear factors."""
    if tensorial:
        return _voigt_operators(B)[0]
    return _voigt_operators(np.swapaxes(np.linalg.inv(B), -1, -2))[1]


def _pair_operator(variance, tensorial, M):
    """Operator of the change of basis M (old vectors in the new basis) on one pair of indices:
    contravariant components transform as M X M^T, covariant ones as M^-T X M^-1."""
    B = M if variance == CONTRAVARIANT else np.swapaxes(np.linalg.inv(M), -1, -2)
    return _congruence_operator(B, tensorial)


def _common_variance(a, b, basis):
    """Variance of a binary result: must agree where it matters (a natural basis)."""
    if a._variance == b._variance:
        return a._variance
    if basis is not None and basis._rotation is None:
        raise ValueError(f"Mixed variance: {a._variance} vs {b._variance} (use to_variance() first)")
    return a._variance


def _relative_matrix(old, new):
    """Old basis vectors written in the new basis: ``A_new^-1 A_old`` (None = lab)."""
    if new is None:
        return old._A
    if old is None:
        return np.linalg.inv(new._A)
    return np.linalg.solve(new._A, old._A)


def _relative_rotation(old, new):
    """Rotation taking components in ``old`` to components in ``new`` (both orthonormal or lab)."""
    if new is None:
        return old._rotation
    if old is None:
        return new._rotation.inv()
    return new._rotation.inv() * old._rotation


def _isotropic_in_basis(t, basis):
    """An isotropic lab tensor (identity, projector) written in ``basis``."""
    if basis is None:
        return t
    if isinstance(basis, Basis) and not basis.single:
        t = type(t).from_tensor(t, len(basis))
    t._check_basis(basis)
    if basis.orthonormal:       # isotropic: same components in every orthonormal basis
        return t.with_basis(basis)
    return t.to_basis(basis)


def _contraction_mats(a, b):
    """3x3 component arrays whose elementwise sum is the double contraction a:b.

    Dual variances (stress:strain) contract as they are; in a natural basis two
    tensors of the same variance need the metric on both indices of one of them.
    """
    basis = _common_basis(a._basis, b._basis)
    ma, mb = a.mat, b.mat
    if basis is not None and basis._rotation is None:
        G = a._metrics()[0]
        b._require_variance()
        if a._variance == b._variance:
            ma = G @ ma @ G
    return ma, mb


# ======================================================================
# _TensorBase — shared logic for Tensor2 and Tensor4
# ======================================================================

class _TensorBase:
    """Private base class for Tensor2 and Tensor4.

    Parameterized by ``_single_ndim``: 1 for Tensor2, 2 for Tensor4.
    Subclasses set this as a class variable.
    """

    __slots__ = ("_data", "_type_str", "_basis", "_variance")

    _single_ndim = None  # set by subclass

    # ------------------------------------------------------------------
    # Internal constructors
    # ------------------------------------------------------------------

    @classmethod
    def _create(cls, data, type_str, basis=None, variance=_DEFAULT):
        """Internal: create from numpy array + type string (+ basis, None = lab; + variance,
        by default the one of the type)."""
        obj = object.__new__(cls)
        obj._data = np.ascontiguousarray(data, dtype=np.float64)
        obj._type_str = type_str
        obj._basis = basis
        obj._variance = cls._DEFAULT_VARIANCE[type_str] if variance is _DEFAULT else variance
        return obj

    def _like(self, data, type_str=None):
        """Same type (unless given), basis and variance as self, other components."""
        return type(self)._create(data, type_str or self._type_str, self._basis, self._variance)

    def _ensure_batch(self):
        """Return _data with a leading batch axis if single."""
        if self._data.ndim == self._single_ndim:
            return self._data[np.newaxis]
        return self._data

    def _rewrap(self, result, type_str=None):
        """Wrap batch result, squeezing if original was single. Keeps the basis."""
        ts = type_str or self._type_str
        if self._data.ndim == self._single_ndim:
            return self._like(result[0], ts)
        return self._like(result, ts)

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

    @property
    def basis(self):
        """The `Basis` the components are written in; ``None`` for the lab basis."""
        return self._basis

    @property
    def variance(self):
        """Variance of the components: ``"contravariant"`` or ``"covariant"`` for a `Tensor2`
        (``None`` when the type has none: ``"symmetric"`` until declared, ``"none"``), a pair
        ``(output, input)`` for a `Tensor4`. Defaults to the one of the type."""
        return self._variance

    def to_variance(self, variance):
        r"""The same tensor with components of another variance (index raising / lowering).

        In the lab or an orthonormal basis the numbers do not change: the tag does.
        In a natural basis one pair of indices is lowered with the metric,
        :math:`T_{ij} = g_{ik}\,T^{kl}\,g_{lj}`, or raised with its inverse. The
        Voigt shear factors stay those of the type. On a ``"symmetric"`` `Tensor2`
        whose variance is not set yet, this declares it.

        Parameters
        ----------
        variance : str or tuple of str
            ``"contravariant"`` / ``"covariant"``; for a `Tensor4` the pair
            ``(output, input)``.

        Returns
        -------
        Tensor2 or Tensor4
            Same tensor, type and basis; components and tag of the new variance.
        """
        variance = self._check_variance_arg(variance)
        if self._type_str == "none":
            raise ValueError("Type 'none' (any 3x3, no Voigt convention) has no variance")
        if variance == self._variance:
            return self
        if self._variance is None or _is_orthonormal(self._basis):
            return type(self)._create(self._data, self._type_str, self._basis, variance)
        return type(self)._create(self._regauged(variance), self._type_str, self._basis, variance)

    # ------------------------------------------------------------------
    # Basis
    # ------------------------------------------------------------------

    def _check_basis(self, basis):
        """Validate a basis for this tensor (kind, batch size, variance)."""
        if basis is None:
            return
        if not isinstance(basis, Basis):
            raise TypeError(f"Expected a Basis or None, got {type(basis)}")
        if not basis.single:
            if self.single:
                raise ValueError("A single tensor needs a single Basis")
            if len(basis) != self._data.shape[0]:
                raise ValueError(
                    f"Basis batch size {len(basis)} != tensor batch size {self._data.shape[0]}")
        if not basis.orthonormal:
            self._require_variance()

    def with_basis(self, basis):
        r"""The tensor with the same components and another basis (a transport).

        The numbers are kept and the basis vectors replaced: the result is a
        different tensor. Giving the convected basis of ``F`` to second
        Piola-Kirchhoff components yields the Kirchhoff stress,
        :math:`\tau^{ij} = S^{IJ}`; likewise :math:`e_{ij} = E_{IJ}` for the
        Green-Lagrange / Almansi pair. No arithmetic, no copy.

        Parameters
        ----------
        basis : Basis or None
            The new basis (``None`` = lab).

        Returns
        -------
        Tensor2 or Tensor4
            Same components and type, new basis.
        """
        self._check_basis(basis)
        return type(self)._create(self._data, self._type_str, basis, self._variance)

    def to_basis(self, basis=None):
        r"""The same tensor, its components re-expressed in another basis.

        With :math:`\mathbf{M} = \mathbf{A}_{new}^{-1} \mathbf{A}_{old}`,
        contravariant components become :math:`\mathbf{M}\,\hat{\mathbf{T}}\,
        \mathbf{M}^T` and covariant ones :math:`\mathbf{M}^{-T}\,
        \hat{\mathbf{T}}\,\mathbf{M}^{-1}`, on each pair of indices.

        Parameters
        ----------
        basis : Basis or None, optional
            Target basis; ``None`` (default) gives the lab components.

        Returns
        -------
        Tensor2 or Tensor4
            Same tensor and type, components in ``basis``.
        """
        old = self._basis
        if _same_basis(old, basis):
            return self
        self._check_basis(basis)
        if _is_orthonormal(old) and _is_orthonormal(basis):
            data = self._rotated_components(_relative_rotation(old, basis), True)
        else:
            self._require_variance()
            data = self._changed_components(_relative_matrix(old, basis))
        return type(self)._create(data, self._type_str, basis, self._variance)

    def _transported(self, F, metric, forward):
        """Push/pull of a tensor that has its own basis: the basis is convected
        (F.A or F^-1.A) and the components are kept, up to the Piola weight."""
        self._require_variance()
        basis = self._basis._convected(F, forward)
        self._check_basis(basis)
        data = self._data
        exponent = self._PIOLA_EXPONENT[self._type_str]
        if metric and exponent:
            J = np.linalg.det(F)
            scale = J ** (exponent if forward else -exponent)
            if np.ndim(scale):
                scale = scale[(Ellipsis,) + (np.newaxis,) * self._single_ndim]
            data = data * scale
        return type(self)._create(data, self._type_str, basis, self._variance)

    def _lab_transport(self, F, metric, forward):
        """Push/pull of a lab tensor whose variance is not the one the C++ kernels assume:
        a change of basis of matrix F (or F^-1) through the Voigt operators, plus the
        Piola weight of the type."""
        self._require_variance()
        M = F if forward else np.linalg.inv(F)
        data = self._changed_components(M)
        exponent = self._PIOLA_EXPONENT[self._type_str]
        if metric and exponent:
            J = np.linalg.det(F)
            scale = J ** (exponent if forward else -exponent)
            if np.ndim(scale):
                scale = scale[(Ellipsis,) + (np.newaxis,) * self._single_ndim]
            data = data * scale
        return self._like(data)

    def _kernel_variance(self):
        """True when the C++ transport kernels apply: default variance of a transportable type."""
        return (self._variance == self._DEFAULT_VARIANCE[self._type_str]
                and self._type_str in self._KERNEL_TYPES)

    def _rotate(self, R, active):
        """rotate() for Tensor2 and Tensor4 (see `Tensor2.rotate`)."""
        basis = self._basis
        if basis is None:
            data = self._rotated_components(R, active)
            return type(self)._create(data, self._type_str,
                                      None if active else Basis(rotation=R), self._variance)
        if active:      # transport: the basis vectors turn, the components stay
            new_basis = basis._pre_rotated(R)
            self._check_basis(new_basis)
            return type(self)._create(self._data, self._type_str, new_basis, self._variance)
        data = self._rotated_components(R, False)
        return type(self)._create(data, self._type_str, basis._post_rotated(R), self._variance)

    # ------------------------------------------------------------------
    # Sequence protocol
    # ------------------------------------------------------------------

    def __len__(self):
        if self.single:
            raise TypeError(f"single {type(self).__name__} has no len()")
        return self._data.shape[0]

    def __getitem__(self, key):
        if self.single:
            raise TypeError(f"single {type(self).__name__} is not subscriptable")
        if isinstance(key, (int, np.integer)):
            if key < 0:
                key += len(self)
            if key < 0 or key >= len(self):
                raise IndexError(f"index {key} out of range for batch of size {len(self)}")
        basis = self._basis
        if basis is not None:
            basis = basis._take(key)
        return type(self)._create(self._data[key].copy(), self._type_str, basis, self._variance)

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
        result_basis = self._basis
        result_variance = self._variance
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
                result_basis = _common_basis(result_basis, inp._basis)
                result_variance = _common_variance(self, inp, result_basis)
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
                return cls._create(result[0], result_ts, result_basis, result_variance)
            return cls._create(result, result_ts, result_basis, result_variance)
        return result

    # ------------------------------------------------------------------
    # Arithmetic
    # ------------------------------------------------------------------

    def __add__(self, other):
        if not isinstance(other, type(self)):
            return NotImplemented
        basis = _common_basis(self._basis, other._basis)
        return type(self)._create(self._data + other._data, self._type_str, basis,
                                  _common_variance(self, other, basis))

    def __radd__(self, other):
        if isinstance(other, type(self)):
            return other.__add__(self)
        return NotImplemented

    def __sub__(self, other):
        if not isinstance(other, type(self)):
            return NotImplemented
        basis = _common_basis(self._basis, other._basis)
        return type(self)._create(self._data - other._data, self._type_str, basis,
                                  _common_variance(self, other, basis))

    def __rsub__(self, other):
        if isinstance(other, type(self)):
            return other.__sub__(self)
        return NotImplemented

    def __neg__(self):
        return self._like(-self._data)

    def __mul__(self, other):
        if isinstance(other, (int, float)):
            return self._like(self._data * other)
        other = np.asarray(other)
        if other.ndim <= 1:
            expand = (Ellipsis,) + (np.newaxis,) * self._single_ndim
            return self._like(self._data * other[expand])
        return NotImplemented

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        if isinstance(other, (int, float)):
            return self._like(self._data / other)
        other = np.asarray(other)
        if other.ndim <= 1:
            expand = (Ellipsis,) + (np.newaxis,) * self._single_ndim
            return self._like(self._data / other[expand])
        return NotImplemented

    def __eq__(self, other):
        if isinstance(other, type(self)):
            if not _same_basis(self._basis, other._basis) or self._variance != other._variance:
                return False  # components in different bases or variances are not comparable
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
        basis = "" if self._basis is None else f", basis={self._basis!r}"
        if self._variance != self._DEFAULT_VARIANCE[self._type_str]:
            basis += f", variance={self._variance!r}"
        if self.single:
            return f"{name}(type='{self._type_str}'{basis})"
        return f"{name}(N={len(self)}, type='{self._type_str}'{basis})"

    # ------------------------------------------------------------------
    # Shared factory helpers
    # ------------------------------------------------------------------

    @classmethod
    def from_tensor(cls, t, n):
        """Broadcast a single tensor to a batch of size ``n`` (copies the data)."""
        if not t.single:
            raise ValueError("from_tensor requires a single tensor")
        expand = (np.newaxis,) + (slice(None),) * t._data.ndim
        shape = (n,) + t._data.shape
        v = np.broadcast_to(t._data[expand], shape).copy()
        return cls._create(v, t._type_str, t._basis, t._variance)

    @classmethod
    def from_list(cls, tensors):
        """Stack a list of single tensors (all of the same type) into a batch."""
        return cls(list(tensors))

    @classmethod
    def concatenate(cls, batches):
        """Join multiple batches and/or singles (all of the same type) into one batch.

        Parts sharing one basis keep it; parts in different bases give a
        per-tensor basis; lab and non-lab parts cannot be mixed."""
        parts = list(batches)
        if not parts:
            raise ValueError("Nothing to concatenate")
        type_str = parts[0]._type_str
        sn = parts[0]._single_ndim
        arrays = []
        variance = parts[0]._variance
        for b in parts:
            if b._type_str != type_str:
                raise ValueError(f"Mixed type: {type_str} vs {b._type_str}")
            if b._variance != variance:
                raise ValueError(f"Mixed variance: {variance} vs {b._variance}")
            d = b._data
            arrays.append(d[np.newaxis] if d.ndim == sn else d)
        basis = _merge_bases([b._basis for b in parts], [a.shape[0] for a in arrays])
        return cls._create(np.concatenate(arrays, axis=0), type_str, basis, variance)


# ======================================================================
# Tensor2
# ======================================================================

class Tensor2(_TensorBase):
    """A 2nd-order tensor with a type tag driving the Voigt convention and rotation dispatch.

    A single object transparently represents either one tensor or a batch
    (scipy ``Rotation`` style): the stored numpy array is ``(6,)`` for a single
    tensor and ``(N, 6)`` for a batch. The type is a string: ``"stress"``,
    ``"strain"`` or ``"symmetric"`` for a symmetric tensor stored as a Voigt
    vector (shear factors ``2*e_ij`` for strain, ``s_ij`` otherwise), or
    ``"none"`` -- no Voigt convention -- for any 3x3 (``F``, ``R``, ``L``,
    ``PK1``) stored as its 9 row-major components, with no Voigt vector, no
    variance and no transport; it accepts an orthonormal `basis` (both legs
    read in that frame), not a natural one. The type also selects the
    rotation/transport rules.

    Construct through the typed factories (`stress`, `strain`, `from_mat`,
    `from_voigt`, `from_mandel`), never through ``Tensor2(array)``.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        sigma = sim.Tensor2.stress(np.array([100., 50., 75., 0., 0., 0.]))
        sigma.mises()          # von Mises equivalent stress
        eps = sim.Tensor2.strain(np.random.randn(1000, 6) * 0.01)  # batch
        len(eps), eps[0]       # 1000, single Tensor2
    """

    _single_ndim = 1
    # power of J = det(F) applied by a push-forward with metric=True (a density weight,
    # tied to the physical kind, not to the variance)
    _PIOLA_EXPONENT = {"stress": -1, "strain": 0, "symmetric": 0}
    # variance of the components unless declared otherwise (None: no variance)
    _DEFAULT_VARIANCE = {"stress": CONTRAVARIANT, "strain": COVARIANT, "symmetric": None, "none": None}
    _KERNEL_TYPES = ("stress", "strain")     # types the C++ transport kernels know

    def _tensorial(self):
        """Voigt shear factors of the type: tensorial (stress-like) or engineering (strain)."""
        return self._type_str != "strain"

    @staticmethod
    def _check_variance_arg(variance):
        return _check_variance(variance)

    def _regauged(self, variance):
        """Components after lowering (to covariant) or raising (to contravariant) both indices
        with the metric of a natural basis."""
        basis = self._basis
        G = basis._metric() if variance == COVARIANT else basis._inverse_metric()
        return _mat_to_voigt(G @ self.mat @ G, self._type_str)

    def __init__(self, data):
        if isinstance(data, list) and data and isinstance(data[0], Tensor2):
            if not all(t.single for t in data):
                raise ValueError("Cannot nest batches")
            ref = data[0]._type_str
            v = np.empty((len(data), _t2_ncomp(ref)), dtype=np.float64)
            for i, t in enumerate(data):
                if t._type_str != ref:
                    raise ValueError(f"Mixed type: {ref} vs {t._type_str} at index {i}")
                if t._variance != data[0]._variance:
                    raise ValueError(f"Mixed variance: {data[0]._variance} vs {t._variance} at index {i}")
                v[i] = t._data
            self._data = v
            self._type_str = ref
            self._basis = _merge_bases([t._basis for t in data], [1] * len(data))
            self._variance = data[0]._variance
        else:
            raise TypeError("Use Tensor2.stress(), Tensor2.strain(), etc.")

    # ------------------------------------------------------------------
    # Factory methods
    # ------------------------------------------------------------------

    @classmethod
    def stress(cls, data):
        """Create stress tensor(s).

        Parameters
        ----------
        data : array_like
            ``(6,)`` Voigt vector, ``(3,3)`` matrix, ``(N,6)`` Voigt batch,
            or ``(N,3,3)`` matrix batch.

        Returns
        -------
        Tensor2
            Stress-typed tensor (single or batch, matching the input shape).
        """
        return cls._from_data(data, "stress")

    @classmethod
    def strain(cls, data):
        """Create strain tensor(s).

        Parameters
        ----------
        data : array_like
            ``(6,)`` Voigt vector (with ``2*e_ij`` shear terms), ``(3,3)``
            matrix, ``(N,6)`` Voigt batch, or ``(N,3,3)`` matrix batch.

        Returns
        -------
        Tensor2
            Strain-typed tensor (single or batch, matching the input shape).
        """
        return cls._from_data(data, "strain")

    @classmethod
    def from_mat(cls, m, type_str):
        """Create from a 3x3 matrix representation with an explicit type.

        Parameters
        ----------
        m : array_like
            ``(3,3)`` matrix or ``(N,3,3)`` batch of matrices.
        type_str : str
            One of ``"stress"``, ``"strain"``, ``"symmetric"`` (symmetric, 6
            components) or ``"none"`` (any 3x3, 9 components).

        Returns
        -------
        Tensor2
            Tensor(s) of the requested type.
        """
        type_str = _check_t2_type(type_str)
        m = np.asarray(m, dtype=np.float64)
        if m.shape == (3, 3):
            return cls._create(_mat_to_voigt(m, type_str).ravel(), type_str)
        if m.ndim == 3 and m.shape[1:] == (3, 3):
            return cls._create(_mat_to_voigt(m, type_str), type_str)
        raise ValueError(f"Expected (3,3) or (N,3,3), got {m.shape}")

    @classmethod
    def from_voigt(cls, v, type_str):
        """Create from a Voigt 6-vector with an explicit type.

        Parameters
        ----------
        v : array_like
            ``(6,)`` Voigt vector or ``(N,6)`` batch. Shear components follow
            the type convention (``2*e_ij`` for strain, ``s_ij`` for stress).
        type_str : str
            One of ``"stress"``, ``"strain"``, ``"symmetric"`` (type ``"none"``
            has no Voigt vector).

        Returns
        -------
        Tensor2
            Tensor(s) of the requested type.
        """
        type_str = _check_t2_type(type_str)
        _require_voigt_type(type_str)
        v = np.asarray(v, dtype=np.float64)
        if v.ndim == 1 and v.size == 6:
            return cls._create(v.copy(), type_str)
        if v.ndim == 2 and v.shape[1] == 6:
            return cls._create(v.copy(), type_str)
        raise ValueError(f"Expected (6,) or (N,6), got {v.shape}")

    @classmethod
    def from_mandel(cls, v, type_str):
        """Create from a Kelvin-Mandel 6-vector (inverse of the `mandel` property).

        Parameters
        ----------
        v : array_like
            ``(6,)`` Kelvin-Mandel vector (``sqrt(2)`` factor on shear terms,
            identical for stress and strain) or ``(N,6)`` batch.
        type_str : str
            One of ``"stress"``, ``"strain"``, ``"symmetric"`` (type ``"none"``
            has no Mandel vector).

        Returns
        -------
        Tensor2
            Tensor(s) of the requested type.
        """
        type_str = _check_t2_type(type_str)
        _require_voigt_type(type_str)
        v = np.asarray(v, dtype=np.float64)
        if not ((v.ndim == 1 and v.size == 6) or (v.ndim == 2 and v.shape[1] == 6)):
            raise ValueError(f"Expected (6,) or (N,6), got {v.shape}")
        v = v.copy()
        v[..., 3:] /= _t2_mandel_factor(type_str)
        return cls._create(v, type_str)

    @classmethod
    def zeros(cls, type_str="stress"):
        """Create a single zero tensor of the given type (default ``"stress"``)."""
        type_str = _check_t2_type(type_str)
        return cls._create(np.zeros(_t2_ncomp(type_str), dtype=np.float64), type_str)

    @classmethod
    def identity(cls, type_str="stress", basis=None):
        r"""Create the identity tensor of the given type.

        In the lab or an orthonormal basis its Voigt components are
        ``[1,1,1,0,0,0]``. In a natural basis the identity is the metric:
        :math:`g^{ij}` for a stress-typed tensor, :math:`g_{ij}` for a
        strain-typed one.

        Parameters
        ----------
        type_str : str, optional
            Tensor type (default ``"stress"``).
        basis : Basis or None, optional
            Basis of the result (default: lab). A batch basis gives a batch.
        """
        type_str = _check_t2_type(type_str)
        if type_str == "none":
            eye = cls._create(np.eye(3).ravel(), type_str)
        else:
            eye = cls._create(np.array([1, 1, 1, 0, 0, 0], dtype=np.float64), type_str)
        return _isotropic_in_basis(eye, basis)

    @classmethod
    def from_columns(cls, arr, type_str):
        """Create a batch from a column-major array (C++ interop).

        Parameters
        ----------
        arr : array_like
            ``(6, N)`` array, one Voigt vector per column (the simcoon C++
            batch convention); ``(9, N)`` for type ``"none"``.
        type_str : str
            One of ``"stress"``, ``"strain"``, ``"symmetric"``, ``"none"``.

        Returns
        -------
        Tensor2
            Batch of N tensors.
        """
        type_str = _check_t2_type(type_str)
        n = _t2_ncomp(type_str)
        arr = np.asarray(arr, dtype=np.float64)
        if arr.ndim != 2 or arr.shape[0] != n:
            raise ValueError(f"Expected ({n}, N), got {arr.shape}")
        return cls._create(arr.T.copy(), type_str)

    # ------------------------------------------------------------------
    # Properties
    # ------------------------------------------------------------------

    @property
    def voigt(self):
        """Voigt vector: (6,) for single, (N,6) for batch. Returns a copy.
        Type ``"none"`` (any 3x3) has no Voigt vector: use `mat`."""
        _require_voigt_type(self._type_str)
        return self._data.copy()

    @property
    def mandel(self):
        """Kelvin-Mandel vector (sqrt2 on shear, identical for stress/strain):
        (6,) for single, (N,6) for batch. Returns a copy."""
        _require_voigt_type(self._type_str)
        v = self._data.copy()
        v[..., 3:] *= _t2_mandel_factor(self._type_str)
        return v

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
        _require_voigt_type(self._type_str)
        if self.single:
            raise AttributeError("voigt_T only available on batch")
        return self._data.T.copy()

    # ------------------------------------------------------------------
    # Domain methods
    # ------------------------------------------------------------------

    def is_symmetric(self, tol=1e-12):
        """Check symmetry of the 3x3 matrix representation (single only).

        Parameters
        ----------
        tol : float, optional
            Absolute tolerance on the off-diagonal differences (default 1e-12).

        Returns
        -------
        bool
        """
        if not self.single:
            raise NotImplementedError("is_symmetric not supported on batch")
        m = self.mat
        return (abs(m[0, 1] - m[1, 0]) < tol and
                abs(m[0, 2] - m[2, 0]) < tol and
                abs(m[1, 2] - m[2, 1]) < tol)

    def _to_cpp(self):
        """Create a temporary _CppTensor2 for single-point C++ operations."""
        return _CppTensor2.from_voigt(self._data, _VTYPE_MAP[self._type_str])

    def rotate(self, R, active=True):
        r"""Rotate the tensor(s) (active) or the frame they are written in (passive).

        Parameters
        ----------
        R : simcoon.Rotation or scipy.spatial.transform.Rotation
            Rotation(s) to apply. On a batch, a single rotation is broadcast
            to all N tensors; a batch of N rotations is applied slice-wise.
        active : bool, optional
            ``True`` (default) rotates the tensor: a new tensor
            :math:`\mathbf{Q}\,\mathbf{X}\,\mathbf{Q}^T`, ``R`` being a lab
            rotation. A lab tensor stays written in the lab; a tensor with its
            own `basis` keeps its components and has its basis vectors turned
            (:math:`\mathbf{A} \to \mathbf{Q}\,\mathbf{A}`).
            ``False`` rotates the frame: the same tensor, with components
            :math:`\mathbf{Q}^T \hat{\mathbf{X}}\,\mathbf{Q}` in the basis
            turned by ``R`` (:math:`\mathbf{A} \to \mathbf{A}\,\mathbf{Q}`,
            ``R`` read in the current basis). The result carries that basis.

        Returns
        -------
        Tensor2
            Same type tag; `basis` as described above.

        Examples
        --------
        .. code-block:: python

            import numpy as np
            import simcoon as sim
            from scipy.spatial.transform import Rotation

            sigma = sim.Tensor2.stress(np.array([100., 0., 0., 0., 0., 0.]))
            R = Rotation.from_euler('z', 45, degrees=True)
            sigma_rot = sigma.rotate(R)                  # another tensor, lab basis
            sigma_45 = sigma.rotate(R, active=False)     # same tensor, basis turned by R
            sigma_45.to_basis(None) == sigma             # True up to round-off
        """
        return self._rotate(R, active)

    def _rotated_components(self, R, active):
        """Components after the orthogonal congruence Q X Q^T (active) or Q^T X Q."""
        if self._type_str == "none":
            # any 3x3: plain congruence, single or batch, any rotation batch
            Q = np.asarray(_as_smc_rotation(R).as_matrix(), dtype=np.float64)
            if not active:
                Q = np.swapaxes(Q, -1, -2)
            m = Q @ self.mat @ np.swapaxes(Q, -1, -2)
            if self.single and m.ndim == 3:
                raise ValueError("A single tensor takes a single rotation")
            return m.reshape(*m.shape[:-2], 9)
        if self.single:
            cpp_result = self._to_cpp().rotate(_to_cpp_rotation(R), active)
            return np.array(cpp_result.voigt).ravel()
        mats = _get_rotation_matrices(R, self._data.shape[0])
        return _batch_rotate(self._data, _VTYPE_MAP[self._type_str], mats, active)

    def _changed_components(self, M):
        """Voigt components after the change of basis of matrix M (old vectors in the new basis).

        With the variance of the type, M X M^T (stress) and M^-T X M^-1 (strain) are
        what the transport kernels compute with F = M and no Piola weight, so those
        are reused; otherwise the 6x6 Voigt operator of the pair is applied."""
        if self._kernel_variance():
            if self.single:
                return np.array(self._to_cpp().push_forward(M, False).voigt).ravel()
            M_batch = M[np.newaxis] if M.ndim == 2 else M
            return _batch_push_forward(
                self._data, _VTYPE_MAP[self._type_str], _to_f_cube(M_batch), False)
        P = _pair_operator(self._variance, self._tensorial(), M)
        return np.einsum("...ij,...j->...i", P, self._data)

    def _require_variance(self):
        if self._variance is None:
            raise ValueError(
                f"Tensor2 type '{self._type_str}' has no variance: a natural basis and "
                "push_forward/pull_back need contravariant (stress-like) or covariant "
                "(strain-like) components. Declare it with to_variance() on a 'symmetric' "
                "tensor; type 'none' (F, R, DR, PK1: two-point or non-symmetric) has no variance "
                "and takes an orthonormal basis only (both legs read in that frame)")

    def _metrics(self):
        """(metric contracting two indices of this tensor, identity of the same
        variance) in a natural basis; None when the metric is the identity."""
        basis = self._basis
        if basis is None or basis._rotation is not None:
            return None
        self._require_variance()
        if self._variance == CONTRAVARIANT:
            return basis._metric(), basis._inverse_metric()
        return basis._inverse_metric(), basis._metric()

    def push_forward(self, F, metric=True):
        r"""Push-forward (reference to current configuration) via the deformation gradient.

        Type-dependent transport: stress is fully contravariant
        (``F s F^T``, Piola with ``1/J``), strain fully covariant
        (``F^-T e F^-1``).

        A push-forward is a transport that carries the basis along: the
        transported tensor has, in the convected basis
        :math:`\mathbf{g}_i = \mathbf{F}\,\mathbf{G}_i`, the components the
        original had in :math:`\mathbf{G}_i`. A lab tensor is returned in lab
        components (the transport is carried out); a tensor with its own
        `basis` keeps its components and gets the convected basis
        :math:`\mathbf{F}\,\mathbf{A}` (only the ``1/J`` weight of
        ``metric=True`` touches the numbers).

        The variance is read from the type tag: the covariant transport is the
        one of the Green-Lagrange / Almansi pair (and the contravariant one of
        the second Piola-Kirchhoff / Kirchhoff pair). A logarithmic strain is
        not the transport of anything: on a strain-typed :math:`\ln \mathbf{V}`
        the result has no meaning. ``F`` is a two-point tensor and is always
        given as lab components.

        Parameters
        ----------
        F : array_like
            ``(3,3)`` deformation gradient, or ``(N,3,3)`` batch (a single F
            is broadcast over a tensor batch).
        metric : bool, optional
            ``True`` (default) includes the ``J = det(F)`` factor (proper
            Piola transformation); ``False`` is pure transport.

        Returns
        -------
        Tensor2
            Transported tensor(s), same type tag.
        """
        F = np.asarray(F, dtype=np.float64)
        self._require_variance()
        if self._basis is not None:
            return self._transported(F, metric, True)
        if not self._kernel_variance():
            return self._lab_transport(F, metric, True)
        if self.single:
            cpp_result = self._to_cpp().push_forward(F, metric)
            return Tensor2._create(np.array(cpp_result.voigt).ravel(),
                                   self._type_str)
        data_2d = self._ensure_batch()
        F_batch = F[np.newaxis] if F.ndim == 2 else F
        result = _batch_push_forward(
            data_2d, _VTYPE_MAP[self._type_str], _to_f_cube(F_batch), metric)
        return self._rewrap(result)

    def pull_back(self, F, metric=True):
        """Pull-back (current to reference configuration) via the deformation gradient.

        Inverse of `push_forward`; same type-dependent transport and
        ``metric`` semantics.

        Parameters
        ----------
        F : array_like
            ``(3,3)`` deformation gradient, or ``(N,3,3)`` batch.
        metric : bool, optional
            ``True`` (default) includes the ``J = det(F)`` factor.

        Returns
        -------
        Tensor2
            Transported tensor(s), same type tag.
        """
        F = np.asarray(F, dtype=np.float64)
        self._require_variance()
        if self._basis is not None:
            return self._transported(F, metric, False)
        if not self._kernel_variance():
            return self._lab_transport(F, metric, False)
        if self.single:
            cpp_result = self._to_cpp().pull_back(F, metric)
            return Tensor2._create(np.array(cpp_result.voigt).ravel(),
                                   self._type_str)
        data_2d = self._ensure_batch()
        F_batch = F[np.newaxis] if F.ndim == 2 else F
        result = _batch_pull_back(
            data_2d, _VTYPE_MAP[self._type_str], _to_f_cube(F_batch), metric)
        return self._rewrap(result)

    def mises(self):
        """Von Mises equivalent (type-aware).

        Uses the stress definition ``sqrt(3/2 s_dev:s_dev)`` for stress/symmetric
        and the strain definition ``sqrt(2/3 e_dev:e_dev)`` for strain.

        Returns
        -------
        float or numpy.ndarray
            Scalar for a single tensor, ``(N,)`` for a batch.
        """
        metrics = self._metrics()
        if metrics is not None:
            # with T = G.m (mixed components): s:s = tr(T.T) - tr(T)^2 / 3
            T = metrics[0] @ self.mat
            tr = np.trace(T, axis1=-2, axis2=-1)
            s2 = np.sum(T * np.swapaxes(T, -1, -2), axis=(-2, -1)) - tr * tr / 3.0
            return np.sqrt((1.5 if self._type_str == "stress" else 2.0 / 3.0) * s2)
        if self._type_str == "none":
            raise ValueError("Mises is not defined for type 'none' (any 3x3, no Voigt convention)")
        d = self._data.copy()
        tr = d[..., 0] + d[..., 1] + d[..., 2]
        d[..., 0] -= tr / 3.0
        d[..., 1] -= tr / 3.0
        d[..., 2] -= tr / 3.0
        diag2 = d[..., 0]**2 + d[..., 1]**2 + d[..., 2]**2
        shear2 = d[..., 3]**2 + d[..., 4]**2 + d[..., 5]**2
        if self._type_str == "strain":
            return np.sqrt((2.0 / 3.0) * (diag2 + 0.5 * shear2))
        return np.sqrt(1.5 * (diag2 + 2.0 * shear2))

    def trace(self):
        """Trace ``t_kk``.

        Returns
        -------
        float or numpy.ndarray
            Scalar for a single tensor, ``(N,)`` for a batch.
        """
        if self._type_str == "none":
            return self._data[..., 0] + self._data[..., 4] + self._data[..., 8]
        metrics = self._metrics()
        if metrics is not None:
            return np.sum(metrics[0] * self.mat, axis=(-2, -1))
        return self._data[..., 0] + self._data[..., 1] + self._data[..., 2]

    def dev(self):
        """Deviatoric part ``t - tr(t)/3 I``.

        Returns
        -------
        Tensor2
            Deviatoric tensor(s), same type tag.
        """
        if self._type_str == "none":
            d = self._data.copy()
            tr = (d[..., 0] + d[..., 4] + d[..., 8]) / 3.0
            d[..., 0] -= tr
            d[..., 4] -= tr
            d[..., 8] -= tr
            return self._like(d)
        metrics = self._metrics()
        if metrics is not None:
            G, identity = metrics
            m = self.mat
            tr = np.sum(G * m, axis=(-2, -1))
            m = m - np.asarray(tr / 3.0)[..., np.newaxis, np.newaxis] * identity
            return self._like(_mat_to_voigt(m, self._type_str))
        d = self._data.copy()
        tr = d[..., 0] + d[..., 1] + d[..., 2]
        d[..., 0] -= tr / 3.0
        d[..., 1] -= tr / 3.0
        d[..., 2] -= tr / 3.0
        return self._like(d)

    def norm(self):
        """Frobenius norm ``sqrt(t_ij t_ij)`` (type-aware shear factors).

        Returns
        -------
        float or numpy.ndarray
            Scalar for a single tensor, ``(N,)`` for a batch.
        """
        if self._type_str == "none":
            return np.sqrt(np.sum(self._data * self._data, axis=-1))
        metrics = self._metrics()
        if metrics is not None:
            T = metrics[0] @ self.mat           # mixed components
            return np.sqrt(np.sum(T * np.swapaxes(T, -1, -2), axis=(-2, -1)))
        d = self._data
        diag2 = d[..., 0]**2 + d[..., 1]**2 + d[..., 2]**2
        shear2 = d[..., 3]**2 + d[..., 4]**2 + d[..., 5]**2
        if self._type_str == "strain":
            return np.sqrt(diag2 + 0.5 * shear2)
        return np.sqrt(diag2 + 2.0 * shear2)

    def det(self):
        r"""Determinant of the tensor (of its mixed components :math:`T^i{}_j`).

        In a natural basis :math:`\det \mathbf{T} = \det \hat{\mathbf{T}}\,
        \det \mathbf{g}` for contravariant components and
        :math:`\det \hat{\mathbf{T}} / \det \mathbf{g}` for covariant ones.

        Returns
        -------
        float or numpy.ndarray
            Scalar for a single tensor, ``(N,)`` for a batch.
        """
        d = np.linalg.det(self.mat)
        metrics = self._metrics()
        if metrics is not None:
            d = d * np.linalg.det(metrics[0])
        return d

    def eigvals(self):
        r"""Eigenvalues (principal values), in ascending order.

        They are those of the mixed components: in a natural basis the
        generalized problem :math:`\det(\hat{\mathbf{T}} - \lambda\,
        \mathbf{g}^{-1}) = 0` (contravariant) or
        :math:`\det(\hat{\mathbf{T}} - \lambda\,\mathbf{g}) = 0` (covariant),
        solved in symmetric form through a Cholesky factor of the metric.

        Returns
        -------
        numpy.ndarray
            ``(3,)`` for a single tensor, ``(N,3)`` for a batch.
        """
        m = self.mat
        if self._type_str == "none":
            w = np.linalg.eigvals(m)                     # possibly complex
            return np.sort(w.real, axis=-1) if np.all(np.abs(w.imag) < 1e-12 * (1 + np.abs(w.real))) else w
        metrics = self._metrics()
        if metrics is not None:
            L = np.linalg.cholesky(metrics[0])          # G = L L^T
            m = np.swapaxes(L, -1, -2) @ m @ L           # same spectrum as T.G
        return np.linalg.eigvalsh(m)

    def __mod__(self, other):
        """Double contraction ``A_ij B_ij`` (``a % b``): float for single, ``(N,)`` for batch."""
        if not isinstance(other, Tensor2):
            return NotImplemented
        ma, mb = _contraction_mats(self, other)
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
        type_str = _check_t2_type(type_str)
        n = _t2_ncomp(type_str)
        data = np.asarray(data, dtype=np.float64)
        if data.shape == (3, 3):
            return cls._create(_mat_to_voigt(data, type_str).ravel(), type_str)
        if data.ndim == 1 and data.size == n:
            return cls._create(data.copy(), type_str)
        if data.ndim == 2 and data.shape[1] == n:
            return cls._create(data.copy(), type_str)
        if data.ndim == 3 and data.shape[1:] == (3, 3):
            return cls._create(_mat_to_voigt(data, type_str), type_str)
        raise ValueError(
            f"Expected ({n},), (3,3), (N,{n}), or (N,3,3), got {data.shape}"
        )


# ======================================================================
# Tensor4
# ======================================================================

class Tensor4(_TensorBase):
    """A 4th-order tensor with a type tag driving rotation dispatch and Voigt factors.

    A single object transparently represents either one tensor or a batch:
    the stored numpy array is ``(6,6)`` (engineering Voigt) for a single
    tensor and ``(N,6,6)`` for a batch. The type is a string: ``"stiffness"``,
    ``"compliance"``, ``"strain_concentration"``, ``"stress_concentration"``,
    or ``"generic"``; it selects the Kelvin-Mandel congruence factors and the
    rotation/transport rules, and is propagated through operations
    (e.g. ``stiffness.inverse()`` is a compliance).

    The C++ backend works internally in the Kelvin-Mandel convention, so
    contraction, inverse and composition are plain 6x6 linear algebra; the
    engineering Voigt form is what `mat` / `voigt` expose.

    Construct through the typed factories (`stiffness`, `compliance`,
    `strain_concentration`, `stress_concentration`, `from_mat`, `from_mandel`),
    never through ``Tensor4(array)``.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        import simcoon as sim

        L = sim.Tensor4.stiffness(sim.L_iso([70000, 0.3], 'Enu'))
        eps = sim.Tensor2.strain(np.array([0.01, -0.003, -0.003, 0.005, 0., 0.]))
        sigma = L @ eps                  # contraction -> stress Tensor2
        M = L.inverse()                  # compliance Tensor4

    .. note::
        The underlying C++ ``tensor4`` class uses a lazy mutable Fastor cache
        that is **not** thread-safe.  However, the Python ``Tensor4`` class
        creates fresh C++ objects for each operation (via ``_to_cpp()``), so
        Python-level ``Tensor4`` instances are safe to use from multiple
        threads (subject to the GIL).
    """

    _single_ndim = 2
    # power of J = det(F) applied by a push-forward with metric=True
    _PIOLA_EXPONENT = {"stiffness": -1, "generic": -1, "compliance": 1,
                       "strain_concentration": 0, "stress_concentration": 0}
    # variance of the (output, input) index pairs unless declared otherwise
    _DEFAULT_VARIANCE = {"stiffness": (CONTRAVARIANT, CONTRAVARIANT), "generic": (CONTRAVARIANT, CONTRAVARIANT),
                         "compliance": (COVARIANT, COVARIANT),
                         "strain_concentration": (COVARIANT, CONTRAVARIANT),
                         "stress_concentration": (CONTRAVARIANT, COVARIANT)}
    # Voigt shear factors of the (output, input) pairs: True = tensorial (stress-like)
    _PAIR_TENSORIAL = {"stiffness": (True, True), "generic": (True, True),
                       "compliance": (False, False),
                       "strain_concentration": (False, True),
                       "stress_concentration": (True, False)}
    _KERNEL_TYPES = ("stiffness", "generic", "compliance")

    @staticmethod
    def _check_variance_arg(variance):
        if not (isinstance(variance, (tuple, list)) and len(variance) == 2):
            raise ValueError("A Tensor4 variance is a pair (output, input)")
        return (_check_variance(variance[0]), _check_variance(variance[1]))

    def _regauged(self, variance):
        """6x6 components after raising/lowering the pairs whose variance changes."""
        basis = self._basis
        tens = self._PAIR_TENSORIAL[self._type_str]
        data = self._data
        ops = []
        for k in range(2):
            if variance[k] == self._variance[k]:
                ops.append(None)
                continue
            B = basis._metric() if variance[k] == COVARIANT else basis._inverse_metric()
            ops.append(_congruence_operator(B, tens[k]))
        if ops[0] is not None:
            data = ops[0] @ data
        if ops[1] is not None:
            data = data @ np.swapaxes(ops[1], -1, -2)
        return data

    def __init__(self, data):
        if isinstance(data, list) and data and isinstance(data[0], Tensor4):
            if not all(t.single for t in data):
                raise ValueError("Cannot nest batches")
            ref = data[0]._type_str
            v = np.empty((len(data), 6, 6), dtype=np.float64)
            for i, t in enumerate(data):
                if t._type_str != ref:
                    raise ValueError(f"Mixed type: {ref} vs {t._type_str} at index {i}")
                if t._variance != data[0]._variance:
                    raise ValueError(f"Mixed variance: {data[0]._variance} vs {t._variance} at index {i}")
                v[i] = t._data
            self._data = v
            self._type_str = ref
            self._basis = _merge_bases([t._basis for t in data], [1] * len(data))
            self._variance = data[0]._variance
        else:
            raise TypeError("Use Tensor4.stiffness(), Tensor4.compliance(), etc.")

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _to_cpp(self):
        """Create a temporary _CppTensor4 for single-point C++ operations."""
        return _CppTensor4.from_mat(self._data, _T4TYPE_MAP[self._type_str])

    def _rotated_components(self, R, active):
        """6x6 components after the rotation congruence (active) or its inverse."""
        if self.single:
            return np.array(self._to_cpp().rotate(_to_cpp_rotation(R), active).mat)
        mats = _get_rotation_matrices(R, self._data.shape[0])
        return _from_f_cube(_batch_rotate(
            _to_f_cube(self._data), _T4TYPE_MAP[self._type_str], mats, active))

    def _changed_components(self, M):
        """6x6 components after the change of basis of matrix M (old vectors in the new basis)."""
        if self._kernel_variance():
            # one variance on the four indices: the transport kernel with F = M, no weight
            if self.single:
                return np.array(self._to_cpp().push_forward(M, False).mat)
            M_batch = M[np.newaxis] if M.ndim == 2 else M
            return _from_f_cube(_batch_push_forward(
                _to_f_cube(self._data), _T4TYPE_MAP[self._type_str],
                _to_f_cube(M_batch), False))
        tens = self._PAIR_TENSORIAL[self._type_str]
        left = _pair_operator(self._variance[0], tens[0], M)
        right = _pair_operator(self._variance[1], tens[1], M)
        return left @ self._data @ np.swapaxes(right, -1, -2)

    def _require_variance(self):
        pass    # every Tensor4 type has a variance

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
        """Create stiffness tensor(s) from an engineering Voigt 6x6.

        Parameters
        ----------
        data : array_like
            ``(6,6)`` matrix (e.g. from ``sim.L_iso``) or ``(N,6,6)`` batch.

        Returns
        -------
        Tensor4
        """
        return cls._typed_factory(data, "stiffness")

    @classmethod
    def compliance(cls, data):
        """Create compliance tensor(s) from an engineering Voigt 6x6.

        Parameters
        ----------
        data : array_like
            ``(6,6)`` matrix (e.g. from ``sim.M_iso``) or ``(N,6,6)`` batch.

        Returns
        -------
        Tensor4
        """
        return cls._typed_factory(data, "compliance")

    @classmethod
    def strain_concentration(cls, data):
        """Create strain concentration tensor(s) (``De = A : D`` maps, e.g. ``sim.A_R(F)``).

        Uses the engineering Voigt convention: the identity concentration is
        ``eye(6)`` and ``A.contract(strain)`` returns a strain.

        Parameters
        ----------
        data : array_like
            ``(6,6)`` matrix or ``(N,6,6)`` batch.

        Returns
        -------
        Tensor4
        """
        return cls._typed_factory(data, "strain_concentration")

    @classmethod
    def stress_concentration(cls, data):
        """Create stress concentration tensor(s) (``s_local = B : s`` maps).

        Parameters
        ----------
        data : array_like
            ``(6,6)`` matrix or ``(N,6,6)`` batch.

        Returns
        -------
        Tensor4
        """
        return cls._typed_factory(data, "stress_concentration")

    @classmethod
    def from_mat(cls, m, type_str):
        """Create from an engineering Voigt 6x6 with an explicit type.

        Parameters
        ----------
        m : array_like
            ``(6,6)`` matrix or ``(N,6,6)`` batch.
        type_str : str
            One of ``"stiffness"``, ``"compliance"``, ``"strain_concentration"``,
            ``"stress_concentration"``, ``"generic"``.

        Returns
        -------
        Tensor4
        """
        return cls._typed_factory(m, type_str)

    @classmethod
    def from_voigt(cls, v, type_str):
        """Alias for `from_mat` (explicit counterpart of `from_mandel`)."""
        return cls.from_mat(v, type_str)

    @classmethod
    def from_mandel(cls, m, type_str):
        """Create from a Kelvin-Mandel 6x6 (inverse of the `mandel` property).

        Parameters
        ----------
        m : array_like
            ``(6,6)`` Kelvin-Mandel matrix (per-type ``sqrt(2)`` congruence
            factors) or ``(N,6,6)`` batch.
        type_str : str
            One of ``"stiffness"``, ``"compliance"``, ``"strain_concentration"``,
            ``"stress_concentration"``, ``"generic"``.

        Returns
        -------
        Tensor4
        """
        _check_t4_type(type_str)
        m = np.asarray(m, dtype=np.float64)
        if m.shape != (6, 6) and not (m.ndim == 3 and m.shape[1:] == (6, 6)):
            raise ValueError(f"Expected (6,6) or (N,6,6), got {m.shape}")
        r, c = _t4_mandel_factors(type_str)
        eng = m.copy()
        eng[..., 3:, :] /= r
        eng[..., :, 3:] /= c
        return cls._create(eng, type_str)

    @classmethod
    def identity(cls, type_str="stiffness", basis=None):
        r"""Create the identity tensor of the given type (default ``"stiffness"``).

        For every type the identity contracts to the unchanged field
        (``I : t == t``); in the lab or an orthonormal basis its Kelvin-Mandel
        form is ``eye(6)``. In a natural basis (``basis=``) the stiffness-type
        identity is :math:`\tfrac12 (g^{ik} g^{jl} + g^{il} g^{jk})` and the
        compliance-type one the same with :math:`g_{ij}`; only the
        concentration (mixed) types keep ``eye(6)``.
        """
        _check_t4_type(type_str)
        cpp = _CppTensor4.identity(_T4TYPE_MAP[type_str])
        return _isotropic_in_basis(cls._create(np.array(cpp.mat), type_str), basis)

    @classmethod
    def volumetric(cls, type_str="stiffness", basis=None):
        r"""Create the volumetric (spherical) projector ``J = 1/3 I⊗I`` of the given type.

        In a natural basis (``basis=``) the stiffness-type projector is
        :math:`\tfrac13 g^{ij} g^{kl}`.
        """
        _check_t4_type(type_str)
        cpp = _CppTensor4.volumetric(_T4TYPE_MAP[type_str])
        return _isotropic_in_basis(cls._create(np.array(cpp.mat), type_str), basis)

    @classmethod
    def deviatoric(cls, type_str="stiffness", basis=None):
        """Create the deviatoric projector ``K = I - J`` of the given type (in ``basis``, default lab)."""
        _check_t4_type(type_str)
        cpp = _CppTensor4.deviatoric(_T4TYPE_MAP[type_str])
        return _isotropic_in_basis(cls._create(np.array(cpp.mat), type_str), basis)

    @classmethod
    def zeros(cls, type_str="stiffness"):
        """Create a single zero tensor of the given type (default ``"stiffness"``)."""
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

    @property
    def mandel(self):
        """Kelvin-Mandel 6x6 (per-type sqrt2 congruence; identity = eye(6) for every
        type): (6,6) single, (N,6,6) batch. Returns a copy."""
        r, c = _t4_mandel_factors(self._type_str)
        m = self._data.copy()
        m[..., 3:, :] *= r
        m[..., :, 3:] *= c
        return m

    # ------------------------------------------------------------------
    # Domain methods
    # ------------------------------------------------------------------

    @staticmethod
    def _infer_contraction_vtype(type_str):
        if type_str in ("stiffness", "stress_concentration", "generic"):
            return "stress"
        return "strain"

    def contract(self, t):
        """Double contraction with a Tensor2 (also available as ``L @ t``).

        The output type follows the tensor4 type: stiffness and
        stress_concentration produce a stress, compliance and
        strain_concentration produce a strain.

        Parameters
        ----------
        t : Tensor2
            Single or batch. Single ⊗ batch combinations broadcast; a single
            Tensor4 contracted with a batch Tensor2 collapses to one BLAS
            matrix-matrix product (the fast shared-tangent path).

        Returns
        -------
        Tensor2
            Contracted tensor(s) with the inferred type.
        """
        out_ts = self._infer_contraction_vtype(self._type_str)
        if t._type_str == "none":
            raise ValueError("A Tensor4 contracts a symmetric Tensor2 (Voigt), not type 'none'")
        basis = _common_basis(self._basis, t._basis)
        out_variance = self._variance[0]
        if basis is not None and basis._rotation is None:
            # natural basis: the contracted indices must be dual
            t._require_variance()
            expected = _DUAL[self._variance[1]]
            if t._variance != expected:
                raise ValueError(
                    f"In a natural basis a Tensor4 with input indices {self._variance[1]} "
                    f"contracts a {expected} Tensor2, got {t._variance} (use to_variance())")

        if self.single:
            result = (self._data @ t._data.T).T
            if t.single:
                return Tensor2._create(result.ravel(), out_ts, basis, out_variance)
            return Tensor2._create(result, out_ts, basis, out_variance)

        t2 = t._data[np.newaxis] if t.single else t._data
        # np.matmul broadcasts (1,6,1) -> (N,6,1) natively
        result = np.matmul(self._data, t2[..., np.newaxis]).squeeze(-1)
        return Tensor2._create(result, out_ts, basis, out_variance)

    def rotate(self, R, active=True):
        """Rotate the tensor(s) (active) or the frame they are written in (passive).

        Same rules as `Tensor2.rotate`, with the type-dependent Voigt congruence.

        Parameters
        ----------
        R : simcoon.Rotation or scipy.spatial.transform.Rotation
            Rotation(s) to apply. On a batch, a single rotation is broadcast;
            a batch of N rotations is applied slice-wise.
        active : bool, optional
            ``True`` (default) rotates the tensor; ``False`` the frame, and the
            result carries the turned `basis`.

        Returns
        -------
        Tensor4
            Same type tag.
        """
        return self._rotate(R, active)

    def push_forward(self, F, metric=True):
        """Push-forward (reference to current configuration) via the deformation gradient.

        Type-dependent transport: stiffness is fully contravariant
        (``F⊗F : L : F^T⊗F^T`` with ``1/J``), compliance fully covariant.
        For a lab tensor, concentration types raise (mixed variance: no
        transport implemented on lab components). A tensor with its own
        `basis` (any type) keeps its components and gets the convected basis
        ``F.A``, as for `Tensor2.push_forward`.

        Parameters
        ----------
        F : array_like
            ``(3,3)`` deformation gradient, or ``(N,3,3)`` batch (a single F
            is broadcast over a tensor batch).
        metric : bool, optional
            ``True`` (default) includes the ``J = det(F)`` factor;
            ``False`` is pure transport.

        Returns
        -------
        Tensor4
            Transported tensor(s), same type tag.
        """
        F = np.asarray(F, dtype=np.float64)
        if self._basis is not None:
            return self._transported(F, metric, True)
        if not self._kernel_variance():
            return self._lab_transport(F, metric, True)
        if self.single:
            cpp_result = self._to_cpp().push_forward(F, metric)
            return Tensor4._create(np.array(cpp_result.mat), self._type_str)
        F_batch = F[np.newaxis] if F.ndim == 2 else F
        result = _batch_push_forward(
            _to_f_cube(self._data), _T4TYPE_MAP[self._type_str],
            _to_f_cube(F_batch), metric)
        return self._rewrap(_from_f_cube(result))

    def pull_back(self, F, metric=True):
        """Pull-back (current to reference configuration) via the deformation gradient.

        Inverse of `push_forward`; same type-dependent transport and
        ``metric`` semantics.

        Parameters
        ----------
        F : array_like
            ``(3,3)`` deformation gradient, or ``(N,3,3)`` batch.
        metric : bool, optional
            ``True`` (default) includes the ``J = det(F)`` factor.

        Returns
        -------
        Tensor4
            Transported tensor(s), same type tag.
        """
        F = np.asarray(F, dtype=np.float64)
        if self._basis is not None:
            return self._transported(F, metric, False)
        if not self._kernel_variance():
            return self._lab_transport(F, metric, False)
        if self.single:
            cpp_result = self._to_cpp().pull_back(F, metric)
            return Tensor4._create(np.array(cpp_result.mat), self._type_str)
        F_batch = F[np.newaxis] if F.ndim == 2 else F
        result = _batch_pull_back(
            _to_f_cube(self._data), _T4TYPE_MAP[self._type_str],
            _to_f_cube(F_batch), metric)
        return self._rewrap(_from_f_cube(result))

    def inverse(self):
        """Invert the tensor (plain 6x6 inverse in Kelvin-Mandel).

        The type follows the algebra: stiffness ↔ compliance,
        concentration types keep their type, generic stays generic.

        Returns
        -------
        Tensor4
            Inverse tensor(s) with the inferred type.
        """
        if self.single:
            cpp_result = self._to_cpp().inverse()
            inv_ts = _T4TYPE_RMAP.get(cpp_result.type, self._type_str)
            return Tensor4._create(np.array(cpp_result.mat), inv_ts, self._basis, self._inverse_variance())
        result, inv_type_enum = _batch_inverse(
            _to_f_cube(self._data), _T4TYPE_MAP[self._type_str])
        inv_ts = _T4TYPE_RMAP.get(inv_type_enum, self._type_str)
        return Tensor4._create(_from_f_cube(result), inv_ts, self._basis, self._inverse_variance())

    def _inverse_variance(self):
        """The inverse map has the dual variances, swapped: (output, input) -> (dual input, dual output)."""
        return (_DUAL[self._variance[1]], _DUAL[self._variance[0]])

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

def _dyadic_variance(a, b):
    """(variance of a, variance of b) when both are set, else the default of a stiffness."""
    if a._variance is None or b._variance is None:
        return _DEFAULT
    return (a._variance, b._variance)


def _cpp_symmetric(t):
    """Temporary C++ tensor2 of a symmetric (Voigt) Tensor2 for the dyadic products."""
    _require_voigt_type(t._type_str)
    return _CppTensor2.from_voigt(t._data, _VTYPE_MAP[t._type_str])


def dyadic(a, b):
    """Dyadic (outer) product ``C_ijkl = a_ij b_kl`` of two single Tensor2.

    Parameters
    ----------
    a, b : Tensor2
        Single tensors.

    Returns
    -------
    Tensor4
        Stiffness-typed tensor.
    """
    if not (isinstance(a, Tensor2) and a.single and isinstance(b, Tensor2) and b.single):
        raise ValueError("dyadic requires two single Tensor2 arguments")
    cpp_a = _cpp_symmetric(a)
    cpp_b = _cpp_symmetric(b)
    result = _dyadic(cpp_a, cpp_b)
    return Tensor4._create(np.array(result.mat),
                           _T4TYPE_RMAP.get(result.type, "stiffness"),
                           _common_basis(a._basis, b._basis), _dyadic_variance(a, b))


def auto_dyadic(a):
    """Dyadic product of a single Tensor2 with itself: ``C_ijkl = a_ij a_kl`` -> stiffness Tensor4."""
    if not (isinstance(a, Tensor2) and a.single):
        raise ValueError("auto_dyadic requires a single Tensor2")
    cpp_a = _cpp_symmetric(a)
    result = _auto_dyadic(cpp_a)
    return Tensor4._create(np.array(result.mat),
                           _T4TYPE_RMAP.get(result.type, "stiffness"), a._basis, _dyadic_variance(a, a))


def sym_dyadic(a, b):
    """Symmetric Voigt dyadic product ``C = v(a) v(b)^T`` of two single Tensor2 -> stiffness Tensor4."""
    if not (isinstance(a, Tensor2) and a.single and isinstance(b, Tensor2) and b.single):
        raise ValueError("sym_dyadic requires two single Tensor2 arguments")
    cpp_a = _cpp_symmetric(a)
    cpp_b = _cpp_symmetric(b)
    result = _sym_dyadic(cpp_a, cpp_b)
    return Tensor4._create(np.array(result.mat),
                           _T4TYPE_RMAP.get(result.type, "stiffness"),
                           _common_basis(a._basis, b._basis), _dyadic_variance(a, b))


def auto_sym_dyadic(a):
    """Symmetric Voigt dyadic product of a single Tensor2 with itself -> stiffness Tensor4."""
    if not (isinstance(a, Tensor2) and a.single):
        raise ValueError("auto_sym_dyadic requires a single Tensor2")
    cpp_a = _cpp_symmetric(a)
    result = _auto_sym_dyadic(cpp_a)
    return Tensor4._create(np.array(result.mat),
                           _T4TYPE_RMAP.get(result.type, "stiffness"), a._basis, _dyadic_variance(a, a))


def double_contract(a, b):
    """Double contraction ``A_ij B_ij`` of two Tensor2 (single or batch).

    Parameters
    ----------
    a, b : Tensor2
        Single or batch; singles broadcast against batches.

    Returns
    -------
    numpy.ndarray
        ``(N,)`` values (``(1,)`` if both are single).
    """
    ma, mb = _contraction_mats(a, b)
    if ma.ndim == 2:
        ma = ma[np.newaxis]
    if mb.ndim == 2:
        mb = mb[np.newaxis]
    return np.sum(ma * mb, axis=(-2, -1))
