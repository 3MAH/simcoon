"""
Rotation library Examples
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The simcoon ``Rotation`` class inherits from ``scipy.spatial.transform.Rotation``,
so all scipy factory methods and operations are available directly.  Plain scipy
``Rotation`` objects can be upgraded with ``from_scipy()`` or passed directly to
C++ methods that accept rotations.
"""

import numpy as np
from scipy.spatial.transform import Rotation as ScipyRotation
import simcoon as sim


# %%
# Rotation API examples
#
# These examples demonstrate the main rotation functions available in simcoon.

v = np.array([1.0, 0.0, 0.0])
m = np.eye(3)
angle = np.pi / 4  # 45 degrees
axis = 2  # y-axis (1=x, 2=y, 3=z)
active = True


# %%
# 1. Rotate a vector using the Rotation class
#
# This example shows how to build a rotation from an angle and axis
# and rotate a vector using :meth:`simcoon.Rotation.apply`.

rot = sim.Rotation.from_axis_angle(angle, axis)
Rmat = rot.as_matrix()
v_rot1 = rot.apply(v)
print("Rotated vector:", v_rot1)
print("Rotation matrix:\n", Rmat)


# %%
# 2. Rotate a vector using Euler angles (scipy syntax)
#
# Uppercase letters indicate intrinsic rotations (rotating frame).
# Lowercase letters indicate extrinsic rotations (fixed frame).

rot2 = sim.Rotation.from_euler("ZXZ", [np.pi / 6, np.pi / 4, np.pi / 3])
v_rot2 = rot2.apply(v)
print("Rotated vector:", v_rot2)


# %%
# 3. Rotate a matrix (3x3 tensor)
#
# This example shows how to rotate a 3x3 tensor using
# :meth:`simcoon.Rotation.apply_tensor`.

m_rot1 = rot.apply_tensor(m)
print("Rotated matrix:\n", m_rot1)


# %%
# 4. Rotate a matrix using a rotation from Euler angles

rot_euler = sim.Rotation.from_euler("ZXZ", [np.pi / 6, np.pi / 4, np.pi / 3])
m_rot2 = rot_euler.apply_tensor(m)
print("Rotated matrix:\n", m_rot2)


# %%
# 5. Build a rotation matrix from Euler angles

psi, theta, phi = np.pi / 6, np.pi / 4, np.pi / 3
r_euler = sim.Rotation.from_euler("ZXZ", [psi, theta, phi])
R_euler = r_euler.as_matrix()
print("Rotation matrix from Euler angles (ZXZ):\n", R_euler)


# %%
# Rotation of mechanical quantities
#
# These examples demonstrate how to rotate mechanical quantities such as stress,
# strain, and stiffness matrices.


# %%
# 6. Rotate a stress vector
#
# This example uses :meth:`simcoon.Rotation.apply_stress`.

stress = np.array([1, 2, 3, 4, 5, 6], dtype=float)
rot_stress1 = rot.apply_stress(stress, active)
print("Rotated stress:", rot_stress1)


# %%
# 7. Rotate a strain vector
#
# This example uses :meth:`simcoon.Rotation.apply_strain`.

strain = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6], dtype=float)
rot_strain1 = rot.apply_strain(strain, active)
print("Rotated strain:", rot_strain1)


# %%
# 8. Rotate a stiffness matrix (L)
#
# This example uses :meth:`simcoon.Rotation.apply_stiffness`.

L6 = np.eye(6)
rotL1 = rot.apply_stiffness(L6, active)
print("Rotated L:\n", rotL1)


# %%
# 9. Rotate a compliance matrix (M)
#
# This example uses :meth:`simcoon.Rotation.apply_compliance`.

M6 = np.eye(6)
rotM1 = rot.apply_compliance(M6, active)
print("Rotated M:\n", rotM1)


# %%
# 10. Rotate a strain concentration tensor (A)
#
# This example uses :meth:`simcoon.Rotation.apply_strain_concentration`.

A6 = np.eye(6)
rotA1 = rot.apply_strain_concentration(A6, active)
print("Rotated A:\n", rotA1)


# %%
# 11. Rotate a stress concentration tensor (B)
#
# This example uses :meth:`simcoon.Rotation.apply_stress_concentration`.

B6 = np.eye(6)
rotB1 = rot.apply_stress_concentration(B6, active)
print("Rotated B:\n", rotB1)


# %%
# 12. Rotate a cubic stiffness tensor
#
# Provide the elastic stiffness tensor for a cubic material and rotate it.

E = 70000.0
nu = 0.3
G = 23000.0
L = sim.L_cubic([E, nu, G], "EnuG")
print(np.array_str(L, precision=2, suppress_small=True))

d = sim.check_symetries(L, 1.0e-2)
print(d["umat_type"])
print(d["props"])

x = sim.L_cubic_props(L)
print(x)

alpha = np.pi / 4.0
rot_z = sim.Rotation.from_axis_angle(alpha, 3)
L_rotate = rot_z.apply_stiffness(L)
print(np.array_str(L_rotate, suppress_small=True))


# %%
# Scipy interoperability
#
# ``simcoon.Rotation`` inherits from ``scipy.spatial.transform.Rotation``,
# so every scipy feature is available.  Plain scipy ``Rotation`` objects can
# be upgraded to ``simcoon.Rotation`` via ``from_scipy()`` to unlock the
# mechanics methods.


# %%
# 13. simcoon.Rotation IS a scipy Rotation
#
# Type checks and scipy methods work directly.

r = sim.Rotation.from_axis_angle(np.pi / 4, 3)
print("isinstance(r, ScipyRotation):", isinstance(r, ScipyRotation))
print("Quaternion:", r.as_quat())
print("Rotation vector:", r.as_rotvec())
print("Euler ZXZ:", r.as_euler("ZXZ", degrees=True))


# %%
# 14. Upgrade a plain scipy Rotation to simcoon
#
# Use ``from_scipy()`` to add mechanics capabilities to an existing scipy
# rotation (e.g. one returned by a third-party library).

scipy_rot = ScipyRotation.from_euler("z", 45, degrees=True)
print("Type before:", type(scipy_rot))

smc_rot = sim.Rotation.from_scipy(scipy_rot)
print("Type after:", type(smc_rot))

# The quaternion is preserved exactly
print("Quaternions match:", np.allclose(scipy_rot.as_quat(), smc_rot.as_quat()))

# Now mechanics methods are available
L_iso = sim.L_iso([70000, 0.3], "Enu")
L_rotated = smc_rot.apply_stiffness(L_iso)
print("Rotated stiffness L[0,0]:", L_rotated[0, 0])


# %%
# 15. Verify scipy and simcoon produce the same rotation
#
# A simcoon rotation and a scipy rotation built from the same parameters
# give identical stiffness results.

angle_test = np.pi / 6
r_sim = sim.Rotation.from_axis_angle(angle_test, 3)
r_scipy = sim.Rotation.from_scipy(ScipyRotation.from_rotvec([0, 0, angle_test]))

L_sim = r_sim.apply_stiffness(L_iso)
L_scipy = r_scipy.apply_stiffness(L_iso)
print("Stiffness matrices match:", np.allclose(L_sim, L_scipy))


# %%
# 16. Compose simcoon and scipy rotations
#
# Because simcoon.Rotation inherits ``__mul__`` from scipy, composition
# with either type works directly.

r1 = sim.Rotation.from_axis_angle(np.pi / 4, 3)  # simcoon
r2 = ScipyRotation.from_euler("x", 30, degrees=True)  # plain scipy

r_composed = r1 * r2  # returns simcoon.Rotation
print("Composed type:", type(r_composed))
print("Composed is simcoon Rotation:", isinstance(r_composed, sim.Rotation))


# %%
# 17. Batch rotations and mean (scipy features)
#
# Generate random orientations and compute their mean — all scipy batch
# operations work out of the box.

rots = sim.Rotation.random(100)
r_mean = rots.mean()
print("Mean rotation type:", type(r_mean))
print("Mean rotation angle:", np.degrees(r_mean.magnitude()), "degrees")


# %%
# 18. Slerp interpolation
#
# Use scipy's ``Slerp`` class for multi-keyframe interpolation, or the
# convenience ``slerp_to()`` method for simple two-rotation interpolation.

from scipy.spatial.transform import Slerp

r_start = sim.Rotation.identity()
r_end = sim.Rotation.from_axis_angle(np.pi / 2, 3)

# Via scipy Slerp class
key_rots = sim.Rotation.concatenate([r_start, r_end])
slerp = Slerp([0, 1], key_rots)
for t in [0.0, 0.25, 0.5, 0.75, 1.0]:
    r_t = slerp(t)
    print(f"  t={t:.2f}: angle = {np.degrees(r_t.magnitude()):.1f} deg")

# Via convenience method
r_mid = r_start.slerp_to(r_end, 0.5)
print("Midpoint angle:", np.degrees(r_mid.magnitude()), "degrees")
