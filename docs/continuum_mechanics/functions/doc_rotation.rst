The Rotation Library
====================

The rotation library provides comprehensive tools for 3D rotations in continuum mechanics.
The ``Rotation`` class inherits from ``scipy.spatial.transform.Rotation``, so all scipy
features (batch operations, ``mean()``, ``Slerp``, ``RotationSpline``, etc.) are available
directly, while simcoon adds continuum-mechanics methods for Voigt-notation tensors.

.. contents:: Contents
   :local:
   :depth: 2

The Rotation Class
------------------

The ``Rotation`` class extends ``scipy.spatial.transform.Rotation`` with methods for
rotating stress, strain, stiffness, and compliance tensors.  It uses unit quaternions
internally for numerical stability and efficient composition.

Creating Rotations
~~~~~~~~~~~~~~~~~~

All scipy factory methods are available and return a ``simcoon.Rotation`` instance:

.. code-block:: python

   import simcoon as smc
   import numpy as np

   # Identity rotation (no rotation)
   r = smc.Rotation.identity()

   # From Euler angles (scipy convention: uppercase = intrinsic)
   r = smc.Rotation.from_euler('ZXZ', [psi, theta, phi])

   # From Euler angles in degrees, extrinsic
   r = smc.Rotation.from_euler('xyz', [roll, pitch, yaw], degrees=True)

   # From rotation matrix
   R = np.array([[1, 0, 0], [0, 0, -1], [0, 1, 0]])  # 90° around x
   r = smc.Rotation.from_matrix(R)

   # From quaternion [qx, qy, qz, qw] (scalar-last)
   q = np.array([0, 0, np.sin(np.pi/4), np.cos(np.pi/4)])  # 90° around z
   r = smc.Rotation.from_quat(q)

   # From rotation vector (axis × angle)
   rotvec = np.array([0, 0, np.pi/2])  # 90° around z
   r = smc.Rotation.from_rotvec(rotvec)

   # From axis and angle (simcoon-specific)
   r = smc.Rotation.from_axis_angle(np.pi/2, 3)  # 90° around z (axis 3)

   # Random rotation (uniform distribution)
   r = smc.Rotation.random()

Converting Between Representations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   r = smc.Rotation.from_euler('ZXZ', [0.5, 0.3, 0.7])

   # To rotation matrix
   R = r.as_matrix()  # 3×3 numpy array

   # To quaternion
   q = r.as_quat()  # [qx, qy, qz, qw]

   # To Euler angles (scipy convention)
   angles = r.as_euler('ZXZ')  # [psi, theta, phi]
   angles_deg = r.as_euler('ZXZ', degrees=True)

   # To rotation vector
   rotvec = r.as_rotvec()
   rotvec_deg = r.as_rotvec(degrees=True)

   # To Voigt rotation matrices (simcoon-specific)
   QS = r.as_voigt_stress_rotation()  # 6×6 for stress
   QE = r.as_voigt_strain_rotation()  # 6×6 for strain

Applying Rotations
~~~~~~~~~~~~~~~~~~

**To 3D vectors:**

.. code-block:: python

   r = smc.Rotation.from_axis_angle(np.pi/2, 3)

   v = np.array([1.0, 0.0, 0.0])
   v_rot = r.apply(v)  # [0, 1, 0]

   # Inverse rotation
   v_back = r.apply(v_rot, inverse=True)  # [1, 0, 0]

**To 3×3 tensors:**

.. code-block:: python

   # Rotate a tensor: R · T · R^T
   T = np.eye(3)
   T_rot = r.apply_tensor(T)

**To Voigt notation tensors:**

.. code-block:: python

   # Stress vector [σ11, σ22, σ33, σ12, σ13, σ23]
   sigma = np.array([100.0, 50.0, 25.0, 10.0, 5.0, 2.0])
   sigma_rot = r.apply_stress(sigma)

   # Strain vector [ε11, ε22, ε33, 2ε12, 2ε13, 2ε23]
   epsilon = np.array([0.01, -0.005, -0.005, 0.002, 0.001, 0.0])
   epsilon_rot = r.apply_strain(epsilon)

   # Stiffness matrix (6×6)
   L = smc.L_iso([210e9, 0.3], "Enu")
   L_rot = r.apply_stiffness(L)

   # Compliance matrix (6×6)
   M = smc.M_iso([210e9, 0.3], "Enu")
   M_rot = r.apply_compliance(M)

   # Strain concentration tensor (6×6): A' = QE * A * QS^T
   A_global = r.apply_strain_concentration(A_local)

   # Stress concentration tensor (6×6): B' = QS * B * QE^T
   B_global = r.apply_stress_concentration(B_local)

Composing Rotations
~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   r1 = smc.Rotation.from_axis_angle(np.pi/4, 3)  # 45° around z
   r2 = smc.Rotation.from_axis_angle(np.pi/6, 1)  # 30° around x

   # Compose: apply r1 first, then r2
   r_combined = r2 * r1

   # This is equivalent to:
   v = np.array([1.0, 0.0, 0.0])
   v1 = r2.apply(r1.apply(v))
   v2 = r_combined.apply(v)
   # v1 == v2

   # Inverse
   r_inv = r1.inv()
   r_identity = r1 * r_inv  # Should be identity

Interpolating Rotations
~~~~~~~~~~~~~~~~~~~~~~~

SLERP (Spherical Linear Interpolation) is available via scipy's ``Slerp`` class:

.. code-block:: python

   from scipy.spatial.transform import Slerp

   r_start = smc.Rotation.identity()
   r_end = smc.Rotation.from_euler('ZXZ', [np.pi/2, np.pi/4, 0])

   key_rots = smc.Rotation.concatenate([r_start, r_end])
   slerp = Slerp([0, 1], key_rots)

   # Interpolate along path
   for t in np.linspace(0, 1, 11):
       r_t = slerp(t)
       print(f"t={t:.1f}: angle={np.degrees(r_t.magnitude()):.1f}°")

   # Or use the convenience method
   r_mid = r_start.slerp_to(r_end, 0.5)

Utility Methods
~~~~~~~~~~~~~~~

.. code-block:: python

   r = smc.Rotation.from_euler('ZXZ', [0.5, 0.3, 0.7])

   # Get rotation magnitude (angle)
   angle = r.magnitude()  # radians
   angle_deg = np.degrees(r.magnitude())

   # Check if identity
   is_id = r.is_identity()

   # Check equality (accounts for quaternion sign ambiguity)
   r2 = smc.Rotation.from_euler('ZXZ', [0.5, 0.3, 0.7])
   are_equal = r.equals(r2, tol=1e-12)

Scipy Integration
~~~~~~~~~~~~~~~~~

Since ``simcoon.Rotation`` inherits from ``scipy.spatial.transform.Rotation``,
all scipy features work directly:

.. code-block:: python

   from scipy.spatial.transform import Rotation as R

   # simcoon.Rotation IS a scipy Rotation
   r = smc.Rotation.from_axis_angle(np.pi/4, 3)
   isinstance(r, R)  # True

   # All scipy methods work directly
   r.as_mrp()           # Modified Rodrigues parameters
   r.as_davenport(...)  # Davenport angles

   # Batch operations
   rots = smc.Rotation.random(100)
   rots.mean()

   # RotationSpline
   from scipy.spatial.transform import RotationSpline
   times = [0, 1, 2, 3]
   key_rots = smc.Rotation.random(4)
   spline = RotationSpline(times, key_rots)
   r_interp = spline(1.5)

**Passing plain scipy Rotations to simcoon:**

If you already have a ``scipy.spatial.transform.Rotation`` (e.g. from another
library), you can upgrade it to a ``simcoon.Rotation`` with ``from_scipy()``:

.. code-block:: python

   from scipy.spatial.transform import Rotation as R

   # Created somewhere else, e.g. by a third-party library
   scipy_rot = R.from_euler('z', 45, degrees=True)

   # Upgrade to simcoon.Rotation to unlock mechanics methods
   r = smc.Rotation.from_scipy(scipy_rot)
   L_rot = r.apply_stiffness(L)   # now available

The C++ layer also accepts plain scipy ``Rotation`` objects directly wherever
a rotation parameter is expected (composition, ``equals``, ``slerp``), so
mixing scipy and simcoon rotations in the same expression works seamlessly.

Active vs Passive Rotations
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``active`` parameter on Voigt methods (``apply_stress``, ``apply_strain``,
``apply_stiffness``, ``apply_compliance``, ``as_voigt_stress_rotation``, ``as_voigt_strain_rotation``) controls the
rotation convention:

- **active=True** (default): **Alibi** rotation — rotates the physical object
  while the coordinate system stays fixed.
- **active=False**: **Alias** rotation — rotates the coordinate system while the
  object stays fixed. This is equivalent to the inverse active rotation.

.. code-block:: python

   r = smc.Rotation.from_axis_angle(np.pi/4, 3)

   # Active: rotate the stress tensor itself
   sigma_active = r.apply_stress(sigma, active=True)

   # Passive: express the same stress in a rotated coordinate system
   sigma_passive = r.apply_stress(sigma, active=False)

Input Validation
~~~~~~~~~~~~~~~~

All mechanics methods that accept NumPy arrays validate the input shape and raise
``ValueError`` with a descriptive message if the dimensions are wrong:

.. code-block:: python

   r = smc.Rotation.from_axis_angle(np.pi/4, 3)

   r.apply_stress(np.array([1.0, 2.0, 3.0]))
   # ValueError: sigma must have 6 elements, got 3

   r.apply_stiffness(np.eye(3))
   # ValueError: L must have shape (6, 6), got (3, 3)

Expected input shapes:

- ``apply``: 1D array, 3 elements (or Nx3 array for batch)
- ``apply_tensor``: 2D array, shape (3, 3)
- ``apply_stress``, ``apply_strain``: 1D array, 6 elements
- ``apply_stiffness``, ``apply_compliance``: 2D array, shape (6, 6)

.. note::

   Quaternions ``q`` and ``-q`` represent the same rotation. The ``equals()``
   method accounts for this antipodal equivalence.

Euler Angle Conventions
-----------------------

Simcoon uses scipy's Euler angle conventions:

- **Uppercase** letters (``'ZXZ'``, ``'XYZ'``, etc.) → intrinsic rotations (rotating frame)
- **Lowercase** letters (``'zxz'``, ``'xyz'``, etc.) → extrinsic rotations (fixed frame)

**Proper Euler angles** (axis sequence where first = last):

- ``ZXZ``, ``ZYZ``, ``XYX``, ``XZX``, ``YXY``, ``YZY``

**Tait-Bryan angles** (all axes different):

- ``XYZ`` (roll-pitch-yaw), ``XZY``, ``YXZ``, ``YZX``, ``ZXY``, ``ZYX``

.. code-block:: python

   # Intrinsic ZXZ
   r1 = smc.Rotation.from_euler('ZXZ', [psi, theta, phi])

   # Extrinsic xyz (same physical rotation as intrinsic ZYX)
   r2 = smc.Rotation.from_euler('xyz', [roll, pitch, yaw])

Practical Examples
------------------

Example 1: Rotating Orthotropic Material Properties
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import simcoon as smc
   import numpy as np

   # Define orthotropic material in local coordinates
   # [E1, E2, E3, nu12, nu13, nu23, G12, G13, G23]
   props = [150e9, 10e9, 10e9, 0.3, 0.3, 0.45, 5e9, 5e9, 3.5e9]
   L_local = smc.L_ortho(props, "EnuG")

   # Fiber orientation: 45° rotation around z-axis
   r = smc.Rotation.from_axis_angle(np.pi/4, 3)

   # Rotate stiffness to global coordinates
   L_global = r.apply_stiffness(L_local)

   print("Local stiffness L11:", L_local[0, 0] / 1e9, "GPa")
   print("Global stiffness L11:", L_global[0, 0] / 1e9, "GPa")

Example 2: Stress Transformation in a Rotated Element
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import simcoon as smc
   import numpy as np

   # Stress in global coordinates (Voigt notation)
   sigma_global = np.array([100.0, 50.0, 0.0, 25.0, 0.0, 0.0])  # MPa

   # Element orientation (Euler angles)
   psi, theta, phi = np.radians([30, 45, 60])

   # Create rotation (scipy convention: uppercase = intrinsic)
   r = smc.Rotation.from_euler('ZXZ', [psi, theta, phi])

   # Transform stress to local element coordinates
   sigma_local = r.inv().apply_stress(sigma_global)

   print("Global stress:", sigma_global)
   print("Local stress:", sigma_local)

Example 3: Averaging Orientations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import simcoon as smc
   import numpy as np

   # Generate random orientations
   rotations = smc.Rotation.random(10)

   # Use scipy's built-in mean
   r_mean = rotations.mean()

   print("Mean rotation angle:", np.degrees(r_mean.magnitude()), "degrees")

See Also
--------

- :doc:`doc_constitutive` - Constitutive matrices (L_iso, L_ortho, etc.)
- :doc:`doc_transfer` - Voigt notation conversions
- :doc:`doc_kinematics` - Finite strain kinematics
- :doc:`/cpp_api/simulation/rotation` - Full C++ API reference
