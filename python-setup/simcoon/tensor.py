"""
Tensor2 and Tensor4 classes with type-tagged Voigt convention and rotation dispatch.

This module provides ``Tensor2`` and ``Tensor4`` classes that wrap simcoon's C++
tensor objects, providing:

- Automatic Voigt convention handling (stress vs strain factors)
- Type-dispatched rotation (stiffness, compliance, etc.)
- Push-forward / pull-back operations
- Integration with ``simcoon.Rotation``

Examples
--------
>>> import simcoon as smc
>>> import numpy as np

>>> # Build isotropic stiffness
>>> L = smc.Tensor4.stiffness(smc.L_iso([70000, 0.3], 'Enu'))

>>> # Build a strain state
>>> eps = smc.Tensor2.strain(np.array([0.01, -0.003, -0.003, 0.005, 0.002, 0.001]))

>>> # Contract: sigma = L : epsilon
>>> sigma = L @ eps
>>> sigma.vtype  # automatically inferred as stress
VoigtType.stress
"""

import numpy as np

from simcoon._core import (
    _CppTensor2,
    _CppTensor4,
    VoigtType,
    Tensor4Type,
)


class Tensor2:
    """A 2nd-order tensor with type tag for Voigt convention and rotation dispatch.

    Do not construct directly --- use the class methods ``stress()``, ``strain()``,
    ``from_mat()``, ``from_voigt()``, ``zeros()``, or ``identity()``.
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
            A (3, 3) matrix or (6,) Voigt vector.
        """
        return cls._from_data(data, VoigtType.stress)

    @classmethod
    def strain(cls, data):
        """Create a strain tensor from a 3x3 matrix or 6-element Voigt vector.

        Parameters
        ----------
        data : array_like
            A (3, 3) matrix or (6,) Voigt vector (with doubled shear).
        """
        return cls._from_data(data, VoigtType.strain)

    @classmethod
    def from_mat(cls, m, vtype):
        """Create from a 3x3 numpy array with explicit VoigtType."""
        m = np.asarray(m, dtype=np.float64)
        return cls(_CppTensor2.from_mat(m, vtype))

    @classmethod
    def from_voigt(cls, v, vtype):
        """Create from a 6-element Voigt vector with explicit VoigtType."""
        v = np.asarray(v, dtype=np.float64).ravel()
        return cls(_CppTensor2.from_voigt(v, vtype))

    @classmethod
    def zeros(cls, vtype=VoigtType.stress):
        """Create a zero tensor."""
        return cls(_CppTensor2.zeros(vtype))

    @classmethod
    def identity(cls, vtype=VoigtType.stress):
        """Create an identity tensor."""
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
        """6-element Voigt vector as numpy array (shape (6,))."""
        return self._cpp.voigt.ravel()

    @property
    def vtype(self):
        """VoigtType tag."""
        return self._cpp.vtype

    # ------------------------------------------------------------------
    # Methods
    # ------------------------------------------------------------------

    def is_symmetric(self, tol=1e-12):
        """Check if the tensor is symmetric."""
        return self._cpp.is_symmetric(tol)

    def rotate(self, R, active=True):
        """Rotate the tensor using a Rotation object.

        Parameters
        ----------
        R : simcoon.Rotation
            The rotation to apply.
        active : bool, optional
            If True (default), apply active rotation.
        """
        from simcoon.rotation import Rotation
        if isinstance(R, Rotation):
            cpp_rot = R._to_cpp()
        else:
            cpp_rot = R
        return Tensor2(self._cpp.rotate(cpp_rot, active))

    def push_forward(self, F):
        """Push-forward via deformation gradient F."""
        F = np.asarray(F, dtype=np.float64)
        return Tensor2(self._cpp.push_forward(F))

    def pull_back(self, F):
        """Pull-back via deformation gradient F."""
        F = np.asarray(F, dtype=np.float64)
        return Tensor2(self._cpp.pull_back(F))

    def mises(self):
        """Compute von Mises equivalent."""
        from simcoon._core import _CppTensor2
        # Use the C++ free function via the voigt approach
        from simcoon.tensor import _mises_impl
        return _mises_impl(self)

    def trace(self):
        """Compute trace."""
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

    Do not construct directly --- use the class methods ``stiffness()``,
    ``compliance()``, ``from_mat()``, etc.
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
        data : array_like
            A (6, 6) Voigt stiffness matrix.
        """
        data = np.asarray(data, dtype=np.float64)
        return cls(_CppTensor4.from_mat(data, Tensor4Type.stiffness))

    @classmethod
    def compliance(cls, data):
        """Create a compliance tensor from a 6x6 Voigt matrix.

        Parameters
        ----------
        data : array_like
            A (6, 6) Voigt compliance matrix.
        """
        data = np.asarray(data, dtype=np.float64)
        return cls(_CppTensor4.from_mat(data, Tensor4Type.compliance))

    @classmethod
    def from_mat(cls, m, tensor_type):
        """Create from a 6x6 numpy array with explicit Tensor4Type."""
        m = np.asarray(m, dtype=np.float64)
        return cls(_CppTensor4.from_mat(m, tensor_type))

    @classmethod
    def identity(cls, tensor_type=Tensor4Type.stiffness):
        """Create symmetric identity (Ireal)."""
        return cls(_CppTensor4.identity(tensor_type))

    @classmethod
    def volumetric(cls, tensor_type=Tensor4Type.stiffness):
        """Create volumetric projector (Ivol)."""
        return cls(_CppTensor4.volumetric(tensor_type))

    @classmethod
    def deviatoric(cls, tensor_type=Tensor4Type.stiffness):
        """Create deviatoric projector (Idev)."""
        return cls(_CppTensor4.deviatoric(tensor_type))

    @classmethod
    def zeros(cls, tensor_type=Tensor4Type.stiffness):
        """Create zero tensor."""
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
        """Contract with a Tensor2: result = mat @ t.voigt().

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
        """Rotate the tensor using a Rotation object."""
        from simcoon.rotation import Rotation
        if isinstance(R, Rotation):
            cpp_rot = R._to_cpp()
        else:
            cpp_rot = R
        return Tensor4(self._cpp.rotate(cpp_rot, active))

    def push_forward(self, F):
        """Push-forward via deformation gradient F."""
        F = np.asarray(F, dtype=np.float64)
        return Tensor4(self._cpp.push_forward(F))

    def pull_back(self, F):
        """Pull-back via deformation gradient F."""
        F = np.asarray(F, dtype=np.float64)
        return Tensor4(self._cpp.pull_back(F))

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
        """Contract with a Tensor2 using the @ operator."""
        if isinstance(other, Tensor2):
            return self.contract(other)
        return NotImplemented

    def __eq__(self, other):
        if isinstance(other, Tensor4):
            return self._cpp == other._cpp
        return NotImplemented

    def __repr__(self):
        return f"Tensor4(type={self.type.name})"


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
