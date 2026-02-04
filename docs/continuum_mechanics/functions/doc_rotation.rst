The Rotation Library
====================

The rotation library provides comprehensive tools for 3D rotations in continuum mechanics.
Simcoon 2.0 introduces a powerful ``Rotation`` class alongside the existing rotation functions.

.. contents:: Contents
   :local:
   :depth: 2

The Rotation Class
------------------

The ``Rotation`` class provides an object-oriented interface for working with 3D rotations.
It uses unit quaternions internally for numerical stability and efficient composition.

Creating Rotations
~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import simcoon as smc
   import numpy as np

   # Identity rotation (no rotation)
   r = smc.Rotation.identity()

   # From Euler angles (default: zxz convention, intrinsic, radians)
   r = smc.Rotation.from_euler(psi, theta, phi)

   # From Euler angles with options
   r = smc.Rotation.from_euler(psi, theta, phi, conv="xyz", intrinsic=False, degrees=True)

   # From rotation matrix
   R = np.array([[1, 0, 0], [0, 0, -1], [0, 1, 0]])  # 90° around x
   r = smc.Rotation.from_matrix(R)

   # From quaternion [qx, qy, qz, qw] (scalar-last)
   q = np.array([0, 0, np.sin(np.pi/4), np.cos(np.pi/4)])  # 90° around z
   r = smc.Rotation.from_quat(q)

   # From rotation vector (axis × angle)
   rotvec = np.array([0, 0, np.pi/2])  # 90° around z
   r = smc.Rotation.from_rotvec(rotvec)

   # From axis and angle
   r = smc.Rotation.from_axis_angle(np.pi/2, 3)  # 90° around z (axis 3)

   # Random rotation (uniform distribution)
   r = smc.Rotation.random()

Converting Between Representations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   r = smc.Rotation.from_euler(0.5, 0.3, 0.7, "zxz")

   # To rotation matrix
   R = r.as_matrix()  # 3×3 numpy array

   # To quaternion
   q = r.as_quat()  # [qx, qy, qz, qw]

   # To Euler angles
   angles = r.as_euler("zxz")  # [psi, theta, phi]
   angles_deg = r.as_euler("zxz", degrees=True)

   # To rotation vector
   rotvec = r.as_rotvec()
   rotvec_deg = r.as_rotvec(degrees=True)

   # To Voigt rotation matrices
   QS = r.as_QS()  # 6×6 for stress
   QE = r.as_QE()  # 6×6 for strain

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

SLERP (Spherical Linear Interpolation) provides smooth interpolation:

.. code-block:: python

   r_start = smc.Rotation.identity()
   r_end = smc.Rotation.from_euler(np.pi/2, np.pi/4, 0, "zxz")

   # Interpolate at t=0.5 (halfway)
   r_mid = r_start.slerp(r_end, 0.5)

   # Interpolate along path
   for t in np.linspace(0, 1, 11):
       r_t = r_start.slerp(r_end, t)
       print(f"t={t:.1f}: angle={r_t.magnitude(degrees=True):.1f}°")

Utility Methods
~~~~~~~~~~~~~~~

.. code-block:: python

   r = smc.Rotation.from_euler(0.5, 0.3, 0.7, "zxz")

   # Get rotation magnitude (angle)
   angle = r.magnitude()  # radians
   angle_deg = r.magnitude(degrees=True)

   # Check if identity
   is_id = r.is_identity()

   # Check equality
   r2 = smc.Rotation.from_euler(0.5, 0.3, 0.7, "zxz")
   are_equal = r.equals(r2, tol=1e-12)

Euler Angle Conventions
-----------------------

Simcoon supports multiple Euler angle conventions:

**Proper Euler angles** (axis sequence where first = last):

- ``zxz`` (default in simcoon)
- ``zyz``
- ``xyx``, ``xzx``
- ``yxy``, ``yzy``

**Tait-Bryan angles** (all axes different):

- ``xyz`` (roll-pitch-yaw)
- ``xzy``, ``yxz``, ``yzx``, ``zxy``, ``zyx``

**Intrinsic vs Extrinsic:**

- **Intrinsic** (default): rotations about axes of the rotating frame
- **Extrinsic**: rotations about fixed world axes

.. code-block:: python

   # Intrinsic ZXZ (simcoon default)
   r1 = smc.Rotation.from_euler(psi, theta, phi, "zxz", intrinsic=True)

   # Extrinsic XYZ (aerospace convention)
   r2 = smc.Rotation.from_euler(roll, pitch, yaw, "xyz", intrinsic=False)

Rotation Functions
------------------

The following functions provide direct rotation operations without creating a ``Rotation`` object.

Rotation Matrix Generation
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: simcoon.fillR_angle

.. autofunction:: simcoon.fillR_euler

Voigt Rotation Matrices
~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: simcoon.fillQS_angle

.. autofunction:: simcoon.fillQS_R

.. autofunction:: simcoon.fillQE_angle

.. autofunction:: simcoon.fillQE_R

Vector and Matrix Rotation
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: simcoon.rotate_vec_R

.. autofunction:: simcoon.rotate_vec_angle

.. autofunction:: simcoon.rotate_mat_R

.. autofunction:: simcoon.rotate_mat_angle

Stress and Strain Rotation
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: simcoon.rotate_stress_angle

.. autofunction:: simcoon.rotate_stress_R

.. autofunction:: simcoon.rotate_strain_angle

.. autofunction:: simcoon.rotate_strain_R

Stiffness and Compliance Rotation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: simcoon.rotateL_angle

.. autofunction:: simcoon.rotateL_R

.. autofunction:: simcoon.rotateM_angle

.. autofunction:: simcoon.rotateM_R

Local-Global Transformations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: simcoon.rotate_l2g_L

.. autofunction:: simcoon.rotate_g2l_L

.. autofunction:: simcoon.rotate_l2g_M

.. autofunction:: simcoon.rotate_g2l_M

.. autofunction:: simcoon.rotate_l2g_strain

.. autofunction:: simcoon.rotate_g2l_strain

.. autofunction:: simcoon.rotate_l2g_stress

.. autofunction:: simcoon.rotate_g2l_stress

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

   # Create rotation
   r = smc.Rotation.from_euler(psi, theta, phi, "zxz")

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
   rotations = [smc.Rotation.random() for _ in range(10)]

   # Approximate mean using iterative SLERP
   r_mean = rotations[0]
   for i, r in enumerate(rotations[1:], 2):
       t = 1.0 / i  # Weight for new rotation
       r_mean = r_mean.slerp(r, t)

   print("Mean rotation angle:", r_mean.magnitude(degrees=True), "degrees")

See Also
--------

- :doc:`doc_constitutive` - Constitutive matrices (L_iso, L_ortho, etc.)
- :doc:`doc_transfer` - Voigt notation conversions
- :doc:`doc_kinematics` - Finite strain kinematics
- :doc:`/cpp_api/simulation/rotation` - Full C++ API reference
