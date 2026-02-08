The Rotation Library
====================

The rotation library provides comprehensive tools for 3D rotations in continuum mechanics.
The ``Rotation`` class provides a powerful object-oriented interface for all rotation operations.

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

   # From Euler angles (convention must be specified)
   r = smc.Rotation.from_euler(psi, theta, phi, "zxz")

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

   # From a scipy.spatial.transform.Rotation (requires scipy)
   from scipy.spatial.transform import Rotation as R
   scipy_rot = R.from_euler('z', 45, degrees=True)
   r = smc.Rotation.from_scipy(scipy_rot)

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
   QS = r.as_voigt_stress_rotation()  # 6×6 for stress
   QE = r.as_voigt_strain_rotation()  # 6×6 for strain

   # To scipy.spatial.transform.Rotation (requires scipy)
   scipy_rot = r.to_scipy()

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

Scipy Interoperability
~~~~~~~~~~~~~~~~~~~~~

Both simcoon and scipy use the scalar-last quaternion convention ``[qx, qy, qz, qw]``,
so conversion between the two is done via quaternion transfer with no trigonometric
or matrix computation overhead. Scipy is an **optional** dependency — it is only
imported when ``to_scipy()`` is called.

.. code-block:: python

   from scipy.spatial.transform import Rotation as R

   # scipy → simcoon
   scipy_rot = R.from_euler('z', 45, degrees=True)
   r = smc.Rotation.from_scipy(scipy_rot)

   # simcoon → scipy
   r = smc.Rotation.from_axis_angle(np.pi/4, 3)
   scipy_rot = r.to_scipy()

   # Round-trip is lossless
   q_original = scipy_rot.as_quat()
   q_roundtrip = smc.Rotation.from_scipy(scipy_rot).to_scipy().as_quat()
   # q_original == q_roundtrip

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

All methods that accept NumPy arrays validate the input shape and raise
``ValueError`` with a descriptive message if the dimensions are wrong:

.. code-block:: python

   r = smc.Rotation.from_axis_angle(np.pi/4, 3)

   r.apply(np.array([1.0, 2.0]))
   # ValueError: v must have 3 elements, got 2

   r.apply_stress(np.array([1.0, 2.0, 3.0]))
   # ValueError: sigma must have 6 elements, got 3

   smc.Rotation.from_matrix(np.eye(4))
   # ValueError: R must have shape (3, 3), got (4, 4)

   r.apply_stiffness(np.eye(3))
   # ValueError: L must have shape (6, 6), got (3, 3)

Expected input shapes:

- ``from_quat``: 1D array, 4 elements
- ``from_matrix``: 2D array, shape (3, 3)
- ``from_rotvec``: 1D array, 3 elements
- ``apply``: 1D array, 3 elements
- ``apply_tensor``: 2D array, shape (3, 3)
- ``apply_stress``, ``apply_strain``: 1D array, 6 elements
- ``apply_stiffness``, ``apply_compliance``: 2D array, shape (6, 6)

.. note::

   Quaternions ``q`` and ``-q`` represent the same rotation. The ``equals()``
   method accounts for this antipodal equivalence.

Euler Angle Conventions
-----------------------

Simcoon supports multiple Euler angle conventions:

**Proper Euler angles** (axis sequence where first = last):

- ``zxz``
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

   # Intrinsic ZXZ
   r1 = smc.Rotation.from_euler(psi, theta, phi, "zxz", intrinsic=True)

   # Extrinsic XYZ (aerospace convention)
   r2 = smc.Rotation.from_euler(roll, pitch, yaw, "xyz", intrinsic=False)

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
