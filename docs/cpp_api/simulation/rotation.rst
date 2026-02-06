==================
The Rotation Module
==================

The Rotation module provides comprehensive tools for 3D rotations in continuum mechanics
applications. It includes both a modern **Rotation class** (inspired by ``scipy.spatial.transform.Rotation``)
and **free functions** for tensor rotation in Voigt notation.

.. contents:: Table of Contents
   :local:
   :depth: 2

Overview
========

Simcoon 2.0 introduces a new ``Rotation`` class that encapsulates rotation operations using
unit quaternions as the internal representation. This design provides:

- **Numerical stability**: Avoids gimbal lock issues present in Euler angle representations
- **Efficient composition**: Quaternion multiplication is faster than matrix multiplication
- **Multiple representations**: Easy conversion between quaternions, matrices, Euler angles, and rotation vectors
- **Voigt notation support**: Direct application to stress/strain tensors and stiffness/compliance matrices

The existing free functions are preserved for backward compatibility and provide direct
operations on matrices and vectors.

The Rotation Class
==================

The ``Rotation`` class uses unit quaternions in **scalar-last convention** ``[qx, qy, qz, qw]``
as its internal representation.

C++ API
-------

.. code-block:: cpp

   #include <simcoon/Simulation/Maths/rotation.hpp>

   using namespace simcoon;

   // Create rotations
   Rotation r1 = Rotation::identity();
   Rotation r2 = Rotation::from_euler(0.5, 0.3, 0.7, "zxz");
   Rotation r3 = Rotation::from_axis_angle(M_PI/4, 3);  // 45° around z-axis
   Rotation r4 = Rotation::from_matrix(R);
   Rotation r5 = Rotation::random();

   // Convert to different representations
   arma::mat::fixed<3,3> R = r2.as_matrix();
   arma::vec::fixed<4> q = r2.as_quat();
   arma::vec::fixed<3> euler = r2.as_euler("zxz");
   arma::vec::fixed<3> rotvec = r2.as_rotvec();

   // Apply to vectors and tensors
   arma::vec::fixed<3> v_rot = r2.apply(v);
   arma::vec::fixed<6> sigma_rot = r2.apply_stress(sigma);
   arma::mat::fixed<6,6> L_rot = r2.apply_stiffness(L);

   // Compose rotations
   Rotation r_combined = r2 * r3;
   Rotation r_inv = r2.inv();

   // Interpolation
   Rotation r_mid = r2.slerp(r3, 0.5);

Python API
----------

.. code-block:: python

   import simcoon as smc
   import numpy as np

   # Create rotations
   r1 = smc.Rotation.identity()
   r2 = smc.Rotation.from_euler(0.5, 0.3, 0.7, "zxz")
   r3 = smc.Rotation.from_axis_angle(np.pi/4, 3)
   r4 = smc.Rotation.from_matrix(R)
   r5 = smc.Rotation.random()

   # Convert to different representations
   R = r2.as_matrix()
   q = r2.as_quat()
   euler = r2.as_euler("zxz")
   rotvec = r2.as_rotvec()

   # Apply to vectors and tensors
   v_rot = r2.apply(v)
   sigma_rot = r2.apply_stress(sigma)
   L_rot = r2.apply_stiffness(L)

   # Compose rotations
   r_combined = r2 * r3
   r_inv = r2.inv()

   # Interpolation
   r_mid = r2.slerp(r3, 0.5)

Factory Methods
---------------

.. list-table:: Factory Methods
   :widths: 30 70
   :header-rows: 1

   * - Method
     - Description
   * - ``identity()``
     - Create an identity rotation (no rotation)
   * - ``from_quat(quat)``
     - Create from quaternion [qx, qy, qz, qw] (scalar-last)
   * - ``from_matrix(R)``
     - Create from 3x3 rotation matrix
   * - ``from_euler(psi, theta, phi, conv, intrinsic, degrees)``
     - Create from Euler angles with specified convention
   * - ``from_rotvec(rotvec, degrees)``
     - Create from rotation vector (axis × angle)
   * - ``from_axis_angle(angle, axis, degrees)``
     - Create from angle and axis (1=x, 2=y, 3=z)
   * - ``random()``
     - Create uniformly distributed random rotation

Euler Angle Conventions
-----------------------

The ``from_euler`` and ``as_euler`` methods support multiple conventions:

**Proper Euler angles** (first and last axis are the same):

- ``zxz``, ``zyz``, ``xyx``, ``xzx``, ``yxy``, ``yzy``

**Tait-Bryan angles** (all three axes are different):

- ``xyz``, ``xzy``, ``yxz``, ``yzx``, ``zxy``, ``zyx``

**Parameters:**

- ``intrinsic`` (bool): If ``true`` (default), rotations are about the axes of the rotating
  coordinate system. If ``false``, rotations are about the fixed coordinate system (extrinsic).
- ``degrees`` (bool): If ``true``, angles are in degrees. Default is ``false`` (radians).

.. code-block:: cpp

   // Intrinsic ZXZ Euler angles
   Rotation r1 = Rotation::from_euler(psi, theta, phi, "zxz", true, false);

   // Extrinsic XYZ (aerospace convention)
   Rotation r2 = Rotation::from_euler(roll, pitch, yaw, "xyz", false, true);

Conversion Methods
------------------

.. list-table:: Conversion Methods
   :widths: 30 70
   :header-rows: 1

   * - Method
     - Description
   * - ``as_quat()``
     - Return quaternion [qx, qy, qz, qw]
   * - ``as_matrix()``
     - Return 3×3 rotation matrix
   * - ``as_euler(conv, intrinsic, degrees)``
     - Return Euler angles [psi, theta, phi]
   * - ``as_rotvec(degrees)``
     - Return rotation vector (axis × angle)
   * - ``as_voigt_stress_rotation(active)``
     - Return 6×6 stress rotation matrix
   * - ``as_voigt_strain_rotation(active)``
     - Return 6×6 strain rotation matrix

Apply Methods
-------------

**3D Vectors and Tensors:**

.. list-table:: Apply Methods for 3D Objects
   :widths: 30 70
   :header-rows: 1

   * - Method
     - Description
   * - ``apply(v, inverse)``
     - Rotate a 3D vector
   * - ``apply_tensor(m, inverse)``
     - Rotate a 3×3 tensor: R·m·R\ :sup:`T`

**Voigt Notation (Continuum Mechanics):**

.. list-table:: Apply Methods for Voigt Notation
   :widths: 30 70
   :header-rows: 1

   * - Method
     - Description
   * - ``apply_stress(sigma, active)``
     - Rotate 6-component stress vector
   * - ``apply_strain(epsilon, active)``
     - Rotate 6-component strain vector
   * - ``apply_stiffness(L, active)``
     - Rotate 6×6 stiffness matrix: QS·L·QS\ :sup:`T`
   * - ``apply_compliance(M, active)``
     - Rotate 6×6 compliance matrix: QE·M·QE\ :sup:`T`
   * - ``apply_strain_concentration(A, active)``
     - Rotate 6×6 strain concentration tensor: QE·A·QS\ :sup:`T`
   * - ``apply_stress_concentration(B, active)``
     - Rotate 6×6 stress concentration tensor: QS·B·QE\ :sup:`T`

The ``active`` parameter controls whether the rotation is **active** (alibi, rotating the object)
or **passive** (alias, rotating the coordinate system).

Operations
----------

.. list-table:: Rotation Operations
   :widths: 30 70
   :header-rows: 1

   * - Method
     - Description
   * - ``inv()``
     - Return inverse rotation
   * - ``magnitude(degrees)``
     - Return rotation angle
   * - ``r1 * r2``
     - Compose rotations (apply r2 first, then r1)
   * - ``r1 *= r2``
     - Compose in-place
   * - ``slerp(other, t)``
     - Spherical linear interpolation (t ∈ [0, 1])
   * - ``equals(other, tol)``
     - Check equality within tolerance
   * - ``is_identity(tol)``
     - Check if identity rotation

Rotation Free Functions
=======================

The following free functions provide direct rotation operations without creating a ``Rotation``
object. They are preserved for backward compatibility and direct matrix operations.

Vector and Matrix Rotation
--------------------------

.. doxygenfunction:: simcoon::rotate_vec(const arma::vec&, const arma::mat&)
   :project: simcoon

.. doxygenfunction:: simcoon::rotate_vec(const arma::vec&, const double&, const int&)
   :project: simcoon

.. doxygenfunction:: simcoon::rotate_mat(const arma::mat&, const arma::mat&)
   :project: simcoon

.. doxygenfunction:: simcoon::rotate_mat(const arma::mat&, const double&, const int&)
   :project: simcoon

Stress and Strain Rotation
--------------------------

.. doxygenfunction:: simcoon::rotate_stress(const arma::vec&, const double&, const int&, const bool&)
   :project: simcoon

.. doxygenfunction:: simcoon::rotate_stress(const arma::vec&, const arma::mat&, const bool&)
   :project: simcoon

.. doxygenfunction:: simcoon::rotate_strain(const arma::vec&, const double&, const int&, const bool&)
   :project: simcoon

.. doxygenfunction:: simcoon::rotate_strain(const arma::vec&, const arma::mat&, const bool&)
   :project: simcoon

Stiffness and Compliance Rotation
---------------------------------

.. doxygenfunction:: simcoon::rotate_stiffness(const arma::mat&, const double&, const int&, const bool&)
   :project: simcoon

.. doxygenfunction:: simcoon::rotate_stiffness(const arma::mat&, const arma::mat&, const bool&)
   :project: simcoon

.. doxygenfunction:: simcoon::rotate_compliance(const arma::mat&, const double&, const int&, const bool&)
   :project: simcoon

.. doxygenfunction:: simcoon::rotate_compliance(const arma::mat&, const arma::mat&, const bool&)
   :project: simcoon

Strain and Stress Concentration Tensor Rotation
------------------------------------------------

.. doxygenfunction:: simcoon::rotate_strain_concentration(const arma::mat&, const double&, const int&, const bool&)
   :project: simcoon

.. doxygenfunction:: simcoon::rotate_strain_concentration(const arma::mat&, const arma::mat&, const bool&)
   :project: simcoon

.. doxygenfunction:: simcoon::rotate_stress_concentration(const arma::mat&, const double&, const int&, const bool&)
   :project: simcoon

.. doxygenfunction:: simcoon::rotate_stress_concentration(const arma::mat&, const arma::mat&, const bool&)
   :project: simcoon

Python Functions
================

The following functions are available in Python for direct rotation operations:

.. autofunction:: simcoon.rotate_vec_R

.. autofunction:: simcoon.rotate_vec_angle

.. autofunction:: simcoon.rotate_mat_R

.. autofunction:: simcoon.rotate_mat_angle

.. autofunction:: simcoon.rotate_stress_angle

.. autofunction:: simcoon.rotate_stress_R

.. autofunction:: simcoon.rotate_strain_angle

.. autofunction:: simcoon.rotate_strain_R

.. autofunction:: simcoon.rotate_stiffness_angle

.. autofunction:: simcoon.rotate_stiffness_R

.. autofunction:: simcoon.rotate_compliance_angle

.. autofunction:: simcoon.rotate_compliance_R

.. autofunction:: simcoon.rotate_strain_concentration_angle

.. autofunction:: simcoon.rotate_strain_concentration_R

.. autofunction:: simcoon.rotate_stress_concentration_angle

.. autofunction:: simcoon.rotate_stress_concentration_R

Examples
========

Example 1: Material Orientation
-------------------------------

Rotating material properties from local to global coordinates:

.. code-block:: python

   import simcoon as smc
   import numpy as np

   # Define orthotropic stiffness in material coordinates
   L_local = smc.L_ortho([E1, E2, E3, nu12, nu13, nu23, G12, G13, G23], "EnuG")

   # Material orientation: Euler angles (radians)
   psi, theta, phi = 0.5, 0.3, 0.7

   # Method 1: Using Rotation class
   r = smc.Rotation.from_euler(psi, theta, phi, "zxz")
   L_global = r.apply_stiffness(L_local)

   # Method 2: Using free function with rotation matrix
   R = smc.Rotation.from_euler(psi, theta, phi, "zxz").as_matrix()
   L_global = smc.rotate_stiffness_R(L_local, R)

Example 2: Stress Transformation
--------------------------------

Transforming stress from global to local coordinates:

.. code-block:: python

   import simcoon as smc
   import numpy as np

   # Stress in global coordinates (Voigt notation)
   sigma_global = np.array([100.0, 50.0, 25.0, 10.0, 5.0, 2.0])

   # Create rotation from orientation
   r = smc.Rotation.from_euler(0.5, 0.3, 0.7, "zxz")

   # Transform to local coordinates (inverse rotation)
   sigma_local = r.inv().apply_stress(sigma_global)

Example 3: Rotation Composition
-------------------------------

Combining multiple rotations:

.. code-block:: cpp

   #include <simcoon/Simulation/Maths/rotation.hpp>

   using namespace simcoon;

   // First rotation: 45° around z-axis
   Rotation r1 = Rotation::from_axis_angle(M_PI/4, 3);

   // Second rotation: 30° around x-axis
   Rotation r2 = Rotation::from_axis_angle(M_PI/6, 1);

   // Combined rotation: apply r1 first, then r2
   Rotation r_combined = r2 * r1;

   // Apply to vector
   arma::vec::fixed<3> v = {1, 0, 0};
   arma::vec::fixed<3> v_rot = r_combined.apply(v);

   // This is equivalent to:
   arma::vec::fixed<3> v_rot2 = r2.apply(r1.apply(v));

Example 4: Interpolating Orientations
-------------------------------------

Smooth interpolation between two orientations using SLERP:

.. code-block:: python

   import simcoon as smc
   import numpy as np

   # Start and end orientations
   r_start = smc.Rotation.identity()
   r_end = smc.Rotation.from_euler(np.pi/2, np.pi/4, 0, "zxz")

   # Interpolate at 10 steps
   for t in np.linspace(0, 1, 10):
       r_t = r_start.slerp(r_end, t)
       R = r_t.as_matrix()
       print(f"t={t:.1f}: angle={r_t.magnitude(degrees=True):.1f}°")

Theory
======

Quaternion Representation
-------------------------

The ``Rotation`` class uses unit quaternions in scalar-last convention:

.. math::

   q = [q_x, q_y, q_z, q_w] = [\mathbf{v}, w]

where :math:`\|\mathbf{q}\| = 1`.

A rotation by angle :math:`\theta` around unit axis :math:`\mathbf{u}` is represented as:

.. math::

   q = \left[\mathbf{u} \sin\frac{\theta}{2}, \cos\frac{\theta}{2}\right]

Quaternion to Rotation Matrix
-----------------------------

The rotation matrix is computed from the quaternion as:

.. math::

   R = \begin{bmatrix}
   1 - 2(q_y^2 + q_z^2) & 2(q_x q_y - q_z q_w) & 2(q_x q_z + q_y q_w) \\
   2(q_x q_y + q_z q_w) & 1 - 2(q_x^2 + q_z^2) & 2(q_y q_z - q_x q_w) \\
   2(q_x q_z - q_y q_w) & 2(q_y q_z + q_x q_w) & 1 - 2(q_x^2 + q_y^2)
   \end{bmatrix}

Voigt Notation
--------------

For stress and strain tensors, simcoon uses Voigt notation:

.. math::

   \boldsymbol{\sigma} = [\sigma_{11}, \sigma_{22}, \sigma_{33}, \sigma_{12}, \sigma_{13}, \sigma_{23}]^T

.. math::

   \boldsymbol{\varepsilon} = [\varepsilon_{11}, \varepsilon_{22}, \varepsilon_{33}, 2\varepsilon_{12}, 2\varepsilon_{13}, 2\varepsilon_{23}]^T

The 6×6 rotation matrices :math:`Q_\sigma` (for stress) and :math:`Q_\varepsilon` (for strain)
are related to the 3×3 rotation matrix :math:`R` and satisfy:

.. math::

   \boldsymbol{\sigma}' = Q_\sigma \boldsymbol{\sigma}

.. math::

   \boldsymbol{\varepsilon}' = Q_\varepsilon \boldsymbol{\varepsilon}

.. math::

   L' = Q_\sigma L Q_\sigma^T

.. math::

   M' = Q_\varepsilon M Q_\varepsilon^T

SLERP Interpolation
-------------------

Spherical Linear Interpolation (SLERP) provides smooth interpolation between rotations:

.. math::

   \text{slerp}(q_0, q_1, t) = \frac{\sin((1-t)\Omega)}{\sin\Omega} q_0 + \frac{\sin(t\Omega)}{\sin\Omega} q_1

where :math:`\Omega = \arccos(q_0 \cdot q_1)` is the angle between the quaternions.

See Also
========

- :doc:`/continuum_mechanics/functions/doc_kinematics` - Finite strain kinematics
- :doc:`/continuum_mechanics/functions/doc_constitutive` - Constitutive matrices
- :doc:`/continuum_mechanics/functions/doc_transfer` - Voigt/tensor conversions
- :doc:`/continuum_mechanics/functions/doc_rotation` - Python user documentation
