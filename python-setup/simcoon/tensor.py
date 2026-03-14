"""
Tensor2 and Tensor4 classes with type-tagged Voigt convention and rotation dispatch.

This module provides ``Tensor2`` and ``Tensor4`` classes that wrap simcoon's C++
tensor objects, providing:

- Automatic Voigt convention handling (stress vs strain factors)
- Type-dispatched rotation (stiffness, compliance, concentration tensors)
- Push-forward / pull-back operations with correct Voigt factor handling
- Integration with ``simcoon.Rotation``

Examples
--------
>>> import simcoon as smc
>>> import numpy as np

>>> # Build isotropic stiffness and its compliance inverse
>>> L = smc.Tensor4.stiffness(smc.L_iso([70000, 0.3], 'Enu'))
>>> M = L.inverse()  # automatically tagged as compliance

>>> # Build a strain state
>>> eps = smc.Tensor2.strain(np.array([0.01, -0.003, -0.003, 0.005, 0.002, 0.001]))

>>> # Contract: sigma = L : epsilon
>>> sigma = L @ eps
>>> sigma.vtype  # automatically inferred as stress
VoigtType.stress

>>> # Dyadic product
>>> C = smc.dyadic(sigma, sigma)  # returns Tensor4(stiffness)
"""

import numpy as np

from simcoon._core import (
    _CppTensor2,
    _CppTensor4,
    _dyadic,
    _auto_dyadic,
    VoigtType,
    Tensor4Type,
)


class Tensor2:
    """A 2nd-order tensor with type tag for Voigt convention and rotation dispatch.

    The ``VoigtType`` tag determines:

    - **Voigt vector factors**: stress uses ``[s11, s22, s33, s12, s13, s23]``,
      strain uses ``[e11, e22, e33, 2*e12, 2*e13, 2*e23]`` (doubled shear).
    - **Rotation rule**: stress rotates via ``QS``, strain via ``QE``.
    - **Push-forward/pull-back**: stress (contravariant) uses ``F·T·F^T``,
      strain (covariant) uses ``F^{-T}·T·F^{-1}``.

    Do not construct directly --- use the class methods ``stress()``, ``strain()``,
    ``from_mat()``, ``from_voigt()``, ``zeros()``, or ``identity()``.

    Parameters
    ----------
    cpp_obj : _CppTensor2
        Internal C++ backend object (not for direct use).
    """

    __slots__ = ("_cpp",)

    def __init__(self, cpp_obj):
        if not isinstance(cpp_obj, _CppTensor2):
            raise TypeError("Use Tensor2.stress(), Tensor2.strain(), etc.")
        self._cpp = cpp_obj

    # ------------------------------------------------------------------
    # Factory methods
    # ------------------------------------------------------------------

    @classmethod
    def stress(cls, data):
        """Create a stress tensor from a 3x3 matrix or 6-element Voigt vector.

        Parameters
        ----------
        data : array_like
            A (3, 3) symmetric matrix or (6,) Voigt vector
            ``[s11, s22, s33, s12, s13, s23]``.

        Returns
        -------
        Tensor2
            With ``VoigtType.stress``.
        """
        return cls._from_data(data, VoigtType.stress)

    @classmethod
    def strain(cls, data):
        """Create a strain tensor from a 3x3 matrix or 6-element Voigt vector.

        Parameters
        ----------
        data : array_like
            A (3, 3) symmetric matrix or (6,) Voigt vector
            ``[e11, e22, e33, 2*e12, 2*e13, 2*e23]`` (doubled shear).

        Returns
        -------
        Tensor2
            With ``VoigtType.strain``.
        """
        return cls._from_data(data, VoigtType.strain)

    @classmethod
    def from_mat(cls, m, vtype):
        """Create from a 3x3 numpy array with explicit VoigtType.

        Parameters
        ----------
        m : array_like, shape (3, 3)
            The 3x3 matrix representation.
        vtype : VoigtType
            The tensor type tag.

        Returns
        -------
        Tensor2
        """
        m = np.asarray(m, dtype=np.float64)
        return cls(_CppTensor2.from_mat(m, vtype))

    @classmethod
    def from_voigt(cls, v, vtype):
        """Create from a 6-element Voigt vector with explicit VoigtType.

        Parameters
        ----------
        v : array_like, shape (6,)
            The Voigt vector.
        vtype : VoigtType
            The tensor type tag.

        Returns
        -------
        Tensor2
        """
        v = np.asarray(v, dtype=np.float64).ravel()
        return cls(_CppTensor2.from_voigt(v, vtype))

    @classmethod
    def zeros(cls, vtype=VoigtType.stress):
        """Create a zero tensor.

        Parameters
        ----------
        vtype : VoigtType, optional
            Default ``VoigtType.stress``.

        Returns
        -------
        Tensor2
        """
        return cls(_CppTensor2.zeros(vtype))

    @classmethod
    def identity(cls, vtype=VoigtType.stress):
        """Create an identity tensor (3x3 eye).

        Parameters
        ----------
        vtype : VoigtType, optional
            Default ``VoigtType.stress``.

        Returns
        -------
        Tensor2
        """
        return cls(_CppTensor2.identity(vtype))

    # ------------------------------------------------------------------
    # Properties
    # ------------------------------------------------------------------

    @property
    def mat(self):
        """3x3 matrix representation as numpy array."""
        return self._cpp.mat

    @property
    def voigt(self):
        """6-element Voigt vector as numpy array (shape ``(6,)``).

        For stress: ``[s11, s22, s33, s12, s13, s23]``.
        For strain: ``[e11, e22, e33, 2*e12, 2*e13, 2*e23]``.
        """
        return self._cpp.voigt.ravel()

    @property
    def vtype(self):
        """VoigtType tag."""
        return self._cpp.vtype

    # ------------------------------------------------------------------
    # Methods
    # ------------------------------------------------------------------

    def is_symmetric(self, tol=1e-12):
        """Check if the tensor is symmetric.

        Parameters
        ----------
        tol : float, optional
            Tolerance for symmetry check (default 1e-12).

        Returns
        -------
        bool
        """
        return self._cpp.is_symmetric(tol)

    def rotate(self, R, active=True):
        """Rotate the tensor using a Rotation object.

        Dispatches based on VoigtType: stress uses ``QS``, strain uses ``QE``,
        generic/none uses direct 3x3 matrix rotation.

        Parameters
        ----------
        R : simcoon.Rotation
            The rotation to apply.
        active : bool, optional
            If True (default), apply active rotation (R * T * R^T).

        Returns
        -------
        Tensor2
            Rotated tensor with same VoigtType.
        """
        from simcoon.rotation import Rotation
        if isinstance(R, Rotation):
            cpp_rot = R._to_cpp()
        else:
            cpp_rot = R
        return Tensor2(self._cpp.rotate(cpp_rot, active))

    def push_forward(self, F):
        """Push-forward via deformation gradient F.

        For stress (contravariant): ``F * T * F^T``.
        For strain (covariant): ``F^{-T} * T * F^{-1}``.

        Parameters
        ----------
        F : array_like, shape (3, 3)
            Deformation gradient.

        Returns
        -------
        Tensor2
            Pushed-forward tensor with same VoigtType.
        """
        F = np.asarray(F, dtype=np.float64)
        return Tensor2(self._cpp.push_forward(F))

    def pull_back(self, F):
        """Pull-back via deformation gradient F.

        For stress (contravariant): ``F^{-1} * T * F^{-T}``.
        For strain (covariant): ``F^T * T * F``.

        Parameters
        ----------
        F : array_like, shape (3, 3)
            Deformation gradient.

        Returns
        -------
        Tensor2
            Pulled-back tensor with same VoigtType.
        """
        F = np.asarray(F, dtype=np.float64)
        return Tensor2(self._cpp.pull_back(F))

    def mises(self):
        """Compute von Mises equivalent.

        For stress: ``sqrt(3/2 * s_ij * s_ij)`` where ``s`` is deviatoric.
        For strain: ``sqrt(2/3 * e_ij * e_ij)`` where ``e`` is deviatoric.

        Returns
        -------
        float
        """
        return _mises_impl(self)

    def trace(self):
        """Compute trace (sum of diagonal elements).

        Returns
        -------
        float
        """
        m = self.mat
        return m[0, 0] + m[1, 1] + m[2, 2]

    # ------------------------------------------------------------------
    # Arithmetic
    # ------------------------------------------------------------------

    def __add__(self, other):
        if isinstance(other, Tensor2):
            return Tensor2(self._cpp + other._cpp)
        return NotImplemented

    def __sub__(self, other):
        if isinstance(other, Tensor2):
            return Tensor2(self._cpp - other._cpp)
        return NotImplemented

    def __mul__(self, other):
        if isinstance(other, (int, float)):
            return Tensor2(self._cpp * float(other))
        return NotImplemented

    def __rmul__(self, other):
        if isinstance(other, (int, float)):
            return Tensor2(self._cpp * float(other))
        return NotImplemented

    def __eq__(self, other):
        if isinstance(other, Tensor2):
            return self._cpp == other._cpp
        return NotImplemented

    def __repr__(self):
        return f"Tensor2(vtype={self.vtype.name})"

    # ------------------------------------------------------------------
    # Internal helper
    # ------------------------------------------------------------------

    @classmethod
    def _from_data(cls, data, vtype):
        data = np.asarray(data, dtype=np.float64)
        if data.shape == (3, 3):
            return cls(_CppTensor2.from_mat(data, vtype))
        elif data.shape == (6,) or (data.ndim == 1 and data.size == 6):
            return cls(_CppTensor2.from_voigt(data.ravel(), vtype))
        else:
            raise ValueError(
                f"Expected (3,3) matrix or (6,) vector, got shape {data.shape}"
            )


class Tensor4:
    """A 4th-order tensor with type tag for rotation dispatch and Voigt storage.

    The ``Tensor4Type`` tag determines:

    - **Rotation rule**: stiffness ``QS·L·QS^T``, compliance ``QE·M·QE^T``,
      strain_concentration ``QE·A·QS^T``, stress_concentration ``QS·B·QE^T``.
    - **Voigt factor convention** for Fastor full-index conversion.
    - **Contraction output type**: stiffness/stress_conc -> stress,
      compliance/strain_conc -> strain.

    Do not construct directly --- use the class methods ``stiffness()``,
    ``compliance()``, ``from_mat()``, etc.

    Parameters
    ----------
    cpp_obj : _CppTensor4
        Internal C++ backend object (not for direct use).
    """

    __slots__ = ("_cpp",)

    def __init__(self, cpp_obj):
        if not isinstance(cpp_obj, _CppTensor4):
            raise TypeError("Use Tensor4.stiffness(), Tensor4.compliance(), etc.")
        self._cpp = cpp_obj

    # ------------------------------------------------------------------
    # Factory methods
    # ------------------------------------------------------------------

    @classmethod
    def stiffness(cls, data):
        """Create a stiffness tensor from a 6x6 Voigt matrix.

        Parameters
        ----------
        data : array_like, shape (6, 6)
            Voigt stiffness matrix (no factor corrections on shear terms).

        Returns
        -------
        Tensor4
            With ``Tensor4Type.stiffness``.
        """
        data = np.asarray(data, dtype=np.float64)
        return cls(_CppTensor4.from_mat(data, Tensor4Type.stiffness))

    @classmethod
    def compliance(cls, data):
        """Create a compliance tensor from a 6x6 Voigt matrix.

        Parameters
        ----------
        data : array_like, shape (6, 6)
            Voigt compliance matrix (factor-2/4 on shear terms).

        Returns
        -------
        Tensor4
            With ``Tensor4Type.compliance``.
        """
        data = np.asarray(data, dtype=np.float64)
        return cls(_CppTensor4.from_mat(data, Tensor4Type.compliance))

    @classmethod
    def strain_concentration(cls, data):
        """Create a strain concentration tensor from a 6x6 Voigt matrix.

        Parameters
        ----------
        data : array_like, shape (6, 6)
            Voigt strain concentration matrix ``A`` where ``eps = A : eps_0``.

        Returns
        -------
        Tensor4
            With ``Tensor4Type.strain_concentration``.
        """
        data = np.asarray(data, dtype=np.float64)
        return cls(_CppTensor4.from_mat(data, Tensor4Type.strain_concentration))

    @classmethod
    def stress_concentration(cls, data):
        """Create a stress concentration tensor from a 6x6 Voigt matrix.

        Parameters
        ----------
        data : array_like, shape (6, 6)
            Voigt stress concentration matrix ``B`` where ``sigma = B : sigma_0``.

        Returns
        -------
        Tensor4
            With ``Tensor4Type.stress_concentration``.
        """
        data = np.asarray(data, dtype=np.float64)
        return cls(_CppTensor4.from_mat(data, Tensor4Type.stress_concentration))

    @classmethod
    def from_mat(cls, m, tensor_type):
        """Create from a 6x6 numpy array with explicit Tensor4Type.

        Parameters
        ----------
        m : array_like, shape (6, 6)
            The 6x6 Voigt matrix.
        tensor_type : Tensor4Type
            The tensor type tag.

        Returns
        -------
        Tensor4
        """
        m = np.asarray(m, dtype=np.float64)
        return cls(_CppTensor4.from_mat(m, tensor_type))

    @classmethod
    def identity(cls, tensor_type=Tensor4Type.stiffness):
        """Create symmetric 4th-order identity (Ireal): diag ``[1,1,1,0.5,0.5,0.5]``.

        Parameters
        ----------
        tensor_type : Tensor4Type, optional
            Default ``Tensor4Type.stiffness``.

        Returns
        -------
        Tensor4
        """
        return cls(_CppTensor4.identity(tensor_type))

    @classmethod
    def identity2(cls, tensor_type=Tensor4Type.stiffness):
        """Create 4th-order identity2 (Ireal2): diag ``[1,1,1,2,2,2]``.

        Parameters
        ----------
        tensor_type : Tensor4Type, optional
            Default ``Tensor4Type.stiffness``.

        Returns
        -------
        Tensor4
        """
        return cls(_CppTensor4.identity2(tensor_type))

    @classmethod
    def volumetric(cls, tensor_type=Tensor4Type.stiffness):
        """Create volumetric projector (Ivol).

        Parameters
        ----------
        tensor_type : Tensor4Type, optional
            Default ``Tensor4Type.stiffness``.

        Returns
        -------
        Tensor4
        """
        return cls(_CppTensor4.volumetric(tensor_type))

    @classmethod
    def deviatoric(cls, tensor_type=Tensor4Type.stiffness):
        """Create deviatoric projector (Idev).

        Parameters
        ----------
        tensor_type : Tensor4Type, optional
            Default ``Tensor4Type.stiffness``.

        Returns
        -------
        Tensor4
        """
        return cls(_CppTensor4.deviatoric(tensor_type))

    @classmethod
    def deviatoric2(cls, tensor_type=Tensor4Type.stiffness):
        """Create deviatoric projector2 (Idev2).

        Parameters
        ----------
        tensor_type : Tensor4Type, optional
            Default ``Tensor4Type.stiffness``.

        Returns
        -------
        Tensor4
        """
        return cls(_CppTensor4.deviatoric2(tensor_type))

    @classmethod
    def zeros(cls, tensor_type=Tensor4Type.stiffness):
        """Create zero tensor.

        Parameters
        ----------
        tensor_type : Tensor4Type, optional
            Default ``Tensor4Type.stiffness``.

        Returns
        -------
        Tensor4
        """
        return cls(_CppTensor4.zeros(tensor_type))

    # ------------------------------------------------------------------
    # Properties
    # ------------------------------------------------------------------

    @property
    def mat(self):
        """6x6 Voigt matrix as numpy array."""
        return self._cpp.mat

    @property
    def type(self):
        """Tensor4Type tag."""
        return self._cpp.type

    # ------------------------------------------------------------------
    # Methods
    # ------------------------------------------------------------------

    def contract(self, t):
        """Contract with a Tensor2: ``result = mat @ t.voigt()``.

        Output VoigtType is inferred from Tensor4Type:

        - stiffness / stress_concentration -> stress
        - compliance / strain_concentration -> strain

        Parameters
        ----------
        t : Tensor2
            The tensor to contract with.

        Returns
        -------
        Tensor2
            Result with VoigtType inferred from Tensor4Type.
        """
        return Tensor2(self._cpp.contract(t._cpp))

    def rotate(self, R, active=True):
        """Rotate the tensor using a Rotation object.

        Dispatches based on Tensor4Type to apply the correct rotation rule.

        Parameters
        ----------
        R : simcoon.Rotation
            The rotation to apply.
        active : bool, optional
            If True (default), apply active rotation.

        Returns
        -------
        Tensor4
            Rotated tensor with same Tensor4Type.
        """
        from simcoon.rotation import Rotation
        if isinstance(R, Rotation):
            cpp_rot = R._to_cpp()
        else:
            cpp_rot = R
        return Tensor4(self._cpp.rotate(cpp_rot, active))

    def push_forward(self, F):
        """Push-forward via deformation gradient F.

        Transforms the full-index tensor:
        ``C'_isrp = F_iL F_sJ F_rM F_pN C_LJMN``.

        The Voigt factor convention is respected for the tensor type,
        so compliance and concentration tensors are handled correctly.

        Parameters
        ----------
        F : array_like, shape (3, 3)
            Deformation gradient.

        Returns
        -------
        Tensor4
            Pushed-forward tensor with same Tensor4Type.
        """
        F = np.asarray(F, dtype=np.float64)
        return Tensor4(self._cpp.push_forward(F))

    def pull_back(self, F):
        """Pull-back via deformation gradient F.

        Inverse of push_forward: uses F^{-1} in the transformation.

        Parameters
        ----------
        F : array_like, shape (3, 3)
            Deformation gradient.

        Returns
        -------
        Tensor4
            Pulled-back tensor with same Tensor4Type.
        """
        F = np.asarray(F, dtype=np.float64)
        return Tensor4(self._cpp.pull_back(F))

    def inverse(self):
        """Invert the 6x6 Voigt matrix.

        Type inference: stiffness <-> compliance. Other types stay the same.

        Returns
        -------
        Tensor4
            Inverted tensor with inferred type.
        """
        return Tensor4(self._cpp.inverse())

    # ------------------------------------------------------------------
    # Arithmetic
    # ------------------------------------------------------------------

    def __add__(self, other):
        if isinstance(other, Tensor4):
            return Tensor4(self._cpp + other._cpp)
        return NotImplemented

    def __sub__(self, other):
        if isinstance(other, Tensor4):
            return Tensor4(self._cpp - other._cpp)
        return NotImplemented

    def __mul__(self, other):
        if isinstance(other, (int, float)):
            return Tensor4(self._cpp * float(other))
        return NotImplemented

    def __rmul__(self, other):
        if isinstance(other, (int, float)):
            return Tensor4(self._cpp * float(other))
        return NotImplemented

    def __matmul__(self, other):
        """Contract with a Tensor2 using the ``@`` operator."""
        if isinstance(other, Tensor2):
            return self.contract(other)
        return NotImplemented

    def __eq__(self, other):
        if isinstance(other, Tensor4):
            return self._cpp == other._cpp
        return NotImplemented

    def __repr__(self):
        return f"Tensor4(type={self.type.name})"


# ======================================================================
# Module-level free functions
# ======================================================================

def dyadic(a, b):
    """Dyadic (outer) product of two Tensor2 objects.

    Computes the 4th-order tensor ``C_IJ = a_I * b_J`` in Voigt notation,
    where the Voigt vectors use stress convention (symmetric average).

    Parameters
    ----------
    a : Tensor2
        First tensor.
    b : Tensor2
        Second tensor.

    Returns
    -------
    Tensor4
        Outer product with ``Tensor4Type.stiffness``.
    """
    return Tensor4(_dyadic(a._cpp, b._cpp))


def auto_dyadic(a):
    """Dyadic (outer) product of a Tensor2 with itself.

    Equivalent to ``dyadic(a, a)``.

    Parameters
    ----------
    a : Tensor2
        The tensor.

    Returns
    -------
    Tensor4
        Outer product with ``Tensor4Type.stiffness``.
    """
    return Tensor4(_auto_dyadic(a._cpp))


def _mises_impl(t):
    """Compute von Mises equivalent for a Tensor2."""
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
