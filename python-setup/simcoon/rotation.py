"""
Rotation class extending scipy.spatial.transform.Rotation with mechanics operations.

This module provides a ``Rotation`` class that inherits from scipy's ``Rotation``,
giving users access to all scipy rotation features (batch operations, ``mean()``,
``Slerp``, ``RotationSpline``, etc.) while adding simcoon's continuum-mechanics
methods (``apply_stress``, ``apply_stiffness``, etc.).
"""

import numpy as np
from scipy.spatial.transform import Rotation as ScipyRotation

from simcoon._core import (_CppRotation,
                           _batch_voigt_stress_rotation,
                           _batch_voigt_strain_rotation)


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
    # Ensure subclass type is preserved (scipy >= 1.17 may not propagate cls)
    # ------------------------------------------------------------------

    @classmethod
    def _ensure_cls(cls, result):
        """Wrap a scipy Rotation in this class if needed."""
        if type(result) is cls:
            return result
        return cls.from_quat(result.as_quat())

    @classmethod
    def from_quat(cls, quat):
        r = super().from_quat(quat)
        if type(r) is not cls:
            r.__class__ = cls
        return r

    @classmethod
    def from_matrix(cls, matrix):
        return cls._ensure_cls(super().from_matrix(matrix))

    @classmethod
    def from_rotvec(cls, rotvec, degrees=False):
        return cls._ensure_cls(super().from_rotvec(rotvec, degrees=degrees))

    @classmethod
    def from_euler(cls, seq, angles, degrees=False):
        return cls._ensure_cls(super().from_euler(seq, angles, degrees=degrees))

    @classmethod
    def from_mrp(cls, mrp):
        return cls._ensure_cls(super().from_mrp(mrp))

    @classmethod
    def identity(cls, num=None):
        return cls._ensure_cls(super().identity(num))

    @classmethod
    def random(cls, num=None, random_state=None):
        return cls._ensure_cls(super().random(num, random_state))

    @classmethod
    def concatenate(cls, rotations):
        return cls._ensure_cls(super().concatenate(rotations))

    def inv(self):
        return type(self)._ensure_cls(super().inv())

    def __mul__(self, other):
        return type(self)._ensure_cls(super().__mul__(other))

    def __getitem__(self, indexer):
        return type(self)._ensure_cls(super().__getitem__(indexer))

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
    # Internal helpers
    # ------------------------------------------------------------------

    @property
    def _is_batch(self):
        """True if this object stores more than one rotation."""
        return self.as_quat().ndim == 2

    def _to_cpp(self):
        """Convert to a _CppRotation for C++ method dispatch (single only)."""
        q = self.as_quat()
        if q.ndim == 2:
            raise ValueError(
                "Cannot convert a batch Rotation to a single _CppRotation. "
                "Batch operations are handled automatically by the Python methods."
            )
        return _CppRotation.from_quat(q)

    def _voigt_stress_matrices(self, active=True):
        """Return QS matrices: (6,6) for single, (N,6,6) for batch."""
        q = self.as_quat()
        if q.ndim == 1:
            return _CppRotation.from_quat(q).as_voigt_stress_rotation(active)
        return _batch_voigt_stress_rotation(q, active)

    def _voigt_strain_matrices(self, active=True):
        """Return QE matrices: (6,6) for single, (N,6,6) for batch."""
        q = self.as_quat()
        if q.ndim == 1:
            return _CppRotation.from_quat(q).as_voigt_strain_rotation(active)
        return _batch_voigt_strain_rotation(q, active)

    # ------------------------------------------------------------------
    # Mechanics methods — support single and batch (Gauss-point) operations
    #
    # Single rotation:
    #   sigma  (6,)   → (6,)
    #   L      (6,6)  → (6,6)
    #
    # Batch of N rotations:
    #   sigma  (6, N) → (6, N)       one stress per rotation
    #   L      (6, 6, N) → (6, 6, N) one stiffness per rotation
    # ------------------------------------------------------------------

    def apply_stress(self, sigma, active=True):
        """Apply rotation to stress vector(s) in Voigt notation.

        Parameters
        ----------
        sigma : array_like
            Single (6,) stress vector, or (6, N) array for batch
            (one column per Gauss point).
        active : bool, optional
            If True (default), active rotation.

        Returns
        -------
        numpy.ndarray
            Rotated stress: (6,) or (6, N).
        """
        sigma = np.asarray(sigma, dtype=float)
        if not self._is_batch:
            return self._to_cpp().apply_stress(sigma.ravel(), active).ravel()
        QS = self._voigt_stress_matrices(active)  # (N, 6, 6)
        return np.einsum("nij,jn->in", QS, sigma)

    def apply_strain(self, epsilon, active=True):
        """Apply rotation to strain vector(s) in Voigt notation.

        Parameters
        ----------
        epsilon : array_like
            Single (6,) strain vector, or (6, N) for batch.
        active : bool, optional
            If True (default), active rotation.

        Returns
        -------
        numpy.ndarray
            Rotated strain: (6,) or (6, N).
        """
        epsilon = np.asarray(epsilon, dtype=float)
        if not self._is_batch:
            return self._to_cpp().apply_strain(epsilon.ravel(), active).ravel()
        QE = self._voigt_strain_matrices(active)  # (N, 6, 6)
        return np.einsum("nij,jn->in", QE, epsilon)

    def apply_stiffness(self, L, active=True):
        """Apply rotation to 6x6 stiffness matrix/matrices.

        Parameters
        ----------
        L : array_like
            Single (6, 6) stiffness matrix, or (6, 6, N) for batch.
        active : bool, optional
            If True (default), active rotation.

        Returns
        -------
        numpy.ndarray
            Rotated stiffness: (6, 6) or (6, 6, N).
        """
        L = np.asarray(L, dtype=float)
        if not self._is_batch:
            return self._to_cpp().apply_stiffness(L, active)
        QS = self._voigt_stress_matrices(active)  # (N, 6, 6)
        # L_rot = QS @ L @ QS^T  for each n
        return np.einsum("nij,jkn,nlk->iln", QS, L, QS)

    def apply_compliance(self, M, active=True):
        """Apply rotation to 6x6 compliance matrix/matrices.

        Parameters
        ----------
        M : array_like
            Single (6, 6) compliance matrix, or (6, 6, N) for batch.
        active : bool, optional
            If True (default), active rotation.

        Returns
        -------
        numpy.ndarray
            Rotated compliance: (6, 6) or (6, 6, N).
        """
        M = np.asarray(M, dtype=float)
        if not self._is_batch:
            return self._to_cpp().apply_compliance(M, active)
        QE = self._voigt_strain_matrices(active)  # (N, 6, 6)
        # M_rot = QE @ M @ QE^T  for each n
        return np.einsum("nij,jkn,nlk->iln", QE, M, QE)

    def apply_strain_concentration(self, A, active=True):
        """Apply rotation to 6x6 strain concentration tensor(s).

        Parameters
        ----------
        A : array_like
            Single (6, 6) tensor, or (6, 6, N) for batch.
        active : bool, optional
            If True (default), active rotation.

        Returns
        -------
        numpy.ndarray
            Rotated tensor: QE * A * QS^T. Shape (6, 6) or (6, 6, N).
        """
        A = np.asarray(A, dtype=float)
        if not self._is_batch:
            return self._to_cpp().apply_strain_concentration(A, active)
        QE = self._voigt_strain_matrices(active)  # (N, 6, 6)
        QS = self._voigt_stress_matrices(active)  # (N, 6, 6)
        return np.einsum("nij,jkn,nlk->iln", QE, A, QS)

    def apply_stress_concentration(self, B, active=True):
        """Apply rotation to 6x6 stress concentration tensor(s).

        Parameters
        ----------
        B : array_like
            Single (6, 6) tensor, or (6, 6, N) for batch.
        active : bool, optional
            If True (default), active rotation.

        Returns
        -------
        numpy.ndarray
            Rotated tensor: QS * B * QE^T. Shape (6, 6) or (6, 6, N).
        """
        B = np.asarray(B, dtype=float)
        if not self._is_batch:
            return self._to_cpp().apply_stress_concentration(B, active)
        QS = self._voigt_stress_matrices(active)  # (N, 6, 6)
        QE = self._voigt_strain_matrices(active)  # (N, 6, 6)
        return np.einsum("nij,jkn,nlk->iln", QS, B, QE)

    def apply_tensor(self, m, inverse=False):
        """Apply rotation to 3x3 tensor(s).

        Parameters
        ----------
        m : array_like
            Single (3, 3) tensor, or (3, 3, N) for batch.
        inverse : bool, optional
            If True, apply inverse rotation. Default is False.

        Returns
        -------
        numpy.ndarray
            Rotated tensor: (3, 3) or (3, 3, N).
        """
        m = np.asarray(m, dtype=float)
        if not self._is_batch:
            return self._to_cpp().apply_tensor(m, inverse)
        # R matrices: (N, 3, 3)
        R = self.as_matrix()
        if inverse:
            # R^T @ m @ R for each n
            return np.einsum("nji,jkn,nkl->iln", R, m, R)
        # R @ m @ R^T for each n
        return np.einsum("nij,jkn,nlk->iln", R, m, R)

    def as_voigt_stress_rotation(self, active=True):
        """Get 6x6 rotation matrix for stress tensors in Voigt notation.

        Parameters
        ----------
        active : bool, optional
            If True (default), active rotation; if False, passive.

        Returns
        -------
        numpy.ndarray
            Single (6, 6) or batch (N, 6, 6) stress rotation matrix (QS).
        """
        return self._voigt_stress_matrices(active)

    def as_voigt_strain_rotation(self, active=True):
        """Get 6x6 rotation matrix for strain tensors in Voigt notation.

        Parameters
        ----------
        active : bool, optional
            If True (default), active rotation; if False, passive.

        Returns
        -------
        numpy.ndarray
            Single (6, 6) or batch (N, 6, 6) strain rotation matrix (QE).
        """
        return self._voigt_strain_matrices(active)

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
