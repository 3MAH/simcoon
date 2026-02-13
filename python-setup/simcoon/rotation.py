"""
Rotation class extending scipy.spatial.transform.Rotation with mechanics operations.

This module provides a ``Rotation`` class that inherits from scipy's ``Rotation``,
giving users access to all scipy rotation features (batch operations, ``mean()``,
``Slerp``, ``RotationSpline``, etc.) while adding simcoon's continuum-mechanics
methods (``apply_stress``, ``apply_stiffness``, etc.).
"""

import numpy as np
from scipy.spatial.transform import Rotation as ScipyRotation

from simcoon._core import _CppRotation


class Rotation(ScipyRotation):
    """Rotation class extending scipy with continuum-mechanics operations.

    Inherits all scipy.spatial.transform.Rotation functionality (batch rotations,
    ``mean()``, ``Slerp``, ``RotationSpline``, ``from_euler``, ``from_quat``, etc.)
    and adds methods for rotating Voigt-notation stress, strain, stiffness, and
    compliance tensors via the C++ backend.

    Examples
    --------
    >>> import simcoon as smc
    >>> import numpy as np

    >>> # All scipy factory methods return a simcoon Rotation
    >>> r = smc.Rotation.from_rotvec([0, 0, np.pi/4])
    >>> isinstance(r, smc.Rotation)
    True

    >>> # Rotate a stiffness matrix
    >>> L = smc.L_iso([70000, 0.3], 'Enu')
    >>> L_rot = r.apply_stiffness(L)

    >>> # Simcoon-specific: create rotation around a principal axis
    >>> r = smc.Rotation.from_axis_angle(np.pi/4, 3)  # 45 deg around z
    """

    # ------------------------------------------------------------------
    # Simcoon-specific factory methods
    # ------------------------------------------------------------------

    @classmethod
    def from_axis_angle(cls, angle, axis, degrees=False):
        """Create a rotation around a principal axis.

        Parameters
        ----------
        angle : float
            Rotation angle.
        axis : int
            Axis of rotation: 1 = x, 2 = y, 3 = z.
        degrees : bool, optional
            If True, *angle* is in degrees. Default is False (radians).

        Returns
        -------
        Rotation
        """
        if axis not in (1, 2, 3):
            raise ValueError(f"axis must be 1, 2, or 3, got {axis}")
        direction = np.zeros(3)
        direction[axis - 1] = 1.0
        return cls.from_rotvec(direction * (angle if not degrees else np.radians(angle)))

    @classmethod
    def from_scipy(cls, scipy_rot):
        """Create a simcoon Rotation from a plain scipy Rotation.

        This is a convenience method for upgrading a
        ``scipy.spatial.transform.Rotation`` to a ``simcoon.Rotation`` so
        that the mechanics methods (``apply_stress``, ``apply_stiffness``,
        etc.) become available.

        Parameters
        ----------
        scipy_rot : scipy.spatial.transform.Rotation
            A scipy Rotation object (single or batch).

        Returns
        -------
        Rotation
        """
        return cls.from_quat(scipy_rot.as_quat())

    # ------------------------------------------------------------------
    # Internal helper
    # ------------------------------------------------------------------

    def _to_cpp(self):
        """Convert to a _CppRotation for C++ method dispatch."""
        return _CppRotation.from_quat(self.as_quat())

    # ------------------------------------------------------------------
    # Mechanics methods (delegate to _CppRotation)
    # ------------------------------------------------------------------

    def apply_stress(self, sigma, active=True):
        """Apply rotation to a stress vector in Voigt notation.

        Parameters
        ----------
        sigma : array_like
            6-component stress vector [s11, s22, s33, s12, s13, s23].
        active : bool, optional
            If True (default), active rotation.

        Returns
        -------
        numpy.ndarray
            Rotated stress vector.
        """
        return self._to_cpp().apply_stress(np.asarray(sigma, dtype=float), active).ravel()

    def apply_strain(self, epsilon, active=True):
        """Apply rotation to a strain vector in Voigt notation.

        Parameters
        ----------
        epsilon : array_like
            6-component strain vector [e11, e22, e33, 2*e12, 2*e13, 2*e23].
        active : bool, optional
            If True (default), active rotation.

        Returns
        -------
        numpy.ndarray
            Rotated strain vector.
        """
        return self._to_cpp().apply_strain(np.asarray(epsilon, dtype=float), active).ravel()

    def apply_stiffness(self, L, active=True):
        """Apply rotation to a 6x6 stiffness matrix.

        Parameters
        ----------
        L : array_like
            6x6 stiffness matrix in Voigt notation.
        active : bool, optional
            If True (default), active rotation.

        Returns
        -------
        numpy.ndarray
            Rotated 6x6 stiffness matrix.
        """
        return self._to_cpp().apply_stiffness(np.asarray(L, dtype=float), active)

    def apply_compliance(self, M, active=True):
        """Apply rotation to a 6x6 compliance matrix.

        Parameters
        ----------
        M : array_like
            6x6 compliance matrix in Voigt notation.
        active : bool, optional
            If True (default), active rotation.

        Returns
        -------
        numpy.ndarray
            Rotated 6x6 compliance matrix.
        """
        return self._to_cpp().apply_compliance(np.asarray(M, dtype=float), active)

    def apply_strain_concentration(self, A, active=True):
        """Apply rotation to a 6x6 strain concentration tensor.

        Parameters
        ----------
        A : array_like
            6x6 strain concentration tensor in Voigt notation.
        active : bool, optional
            If True (default), active rotation.

        Returns
        -------
        numpy.ndarray
            Rotated strain concentration tensor: QE * A * QS^T.
        """
        return self._to_cpp().apply_strain_concentration(np.asarray(A, dtype=float), active)

    def apply_stress_concentration(self, B, active=True):
        """Apply rotation to a 6x6 stress concentration tensor.

        Parameters
        ----------
        B : array_like
            6x6 stress concentration tensor in Voigt notation.
        active : bool, optional
            If True (default), active rotation.

        Returns
        -------
        numpy.ndarray
            Rotated stress concentration tensor: QS * B * QE^T.
        """
        return self._to_cpp().apply_stress_concentration(np.asarray(B, dtype=float), active)

    def apply_tensor(self, m, inverse=False):
        """Apply rotation to a 3x3 tensor (matrix).

        Parameters
        ----------
        m : array_like
            3x3 tensor to rotate.
        inverse : bool, optional
            If True, apply inverse rotation. Default is False.

        Returns
        -------
        numpy.ndarray
            Rotated tensor: R * m * R^T (or R^T * m * R for inverse).
        """
        return self._to_cpp().apply_tensor(np.asarray(m, dtype=float), inverse)

    def as_voigt_stress_rotation(self, active=True):
        """Get 6x6 rotation matrix for stress tensors in Voigt notation.

        Parameters
        ----------
        active : bool, optional
            If True (default), active rotation; if False, passive.

        Returns
        -------
        numpy.ndarray
            6x6 stress rotation matrix (QS).
        """
        return self._to_cpp().as_voigt_stress_rotation(active)

    def as_voigt_strain_rotation(self, active=True):
        """Get 6x6 rotation matrix for strain tensors in Voigt notation.

        Parameters
        ----------
        active : bool, optional
            If True (default), active rotation; if False, passive.

        Returns
        -------
        numpy.ndarray
            6x6 strain rotation matrix (QE).
        """
        return self._to_cpp().as_voigt_strain_rotation(active)

    # ------------------------------------------------------------------
    # Compatibility helpers
    # ------------------------------------------------------------------

    def equals(self, other, tol=1e-12):
        """Check if this rotation equals another within tolerance.

        Accounts for quaternion antipodal equivalence (q and -q represent
        the same rotation).

        Parameters
        ----------
        other : Rotation
            Rotation to compare with.
        tol : float, optional
            Tolerance for comparison. Default is 1e-12.

        Returns
        -------
        bool
        """
        q1 = self.as_quat()
        q2 = other.as_quat()
        return np.allclose(q1, q2, atol=tol) or np.allclose(q1, -q2, atol=tol)

    def is_identity(self, tol=1e-12):
        """Check if this rotation is identity (no rotation).

        Parameters
        ----------
        tol : float, optional
            Tolerance for comparison. Default is 1e-12.

        Returns
        -------
        bool
        """
        return self.magnitude() < tol

    def slerp_to(self, other, t):
        """Spherical linear interpolation between this rotation and *other*.

        For full-featured interpolation over multiple keyframes, use
        ``scipy.spatial.transform.Slerp`` instead.

        Parameters
        ----------
        other : Rotation
            Target rotation.
        t : float
            Interpolation parameter in [0, 1] (0 = self, 1 = other).

        Returns
        -------
        Rotation
            Interpolated rotation.
        """
        from scipy.spatial.transform import Slerp
        key_rots = type(self).concatenate([self, other])
        interp = Slerp([0.0, 1.0], key_rots)
        return interp(t)
