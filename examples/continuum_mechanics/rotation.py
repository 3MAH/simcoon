"""
Rotation library Examples
~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation as R
import simcoon as sim
import os


# %%
# Rotation API examples
#
# These examples demonstrate the main rotation functions available in simcoon.

v = np.array([1.0, 0.0, 0.0])
m = np.eye(3)
angle = np.pi / 4  # 45 degrees
axis = 2  # y-axis (1=x, 2=y, 3=z)
copy = True
active = True


# %%
# 1. Rotate a vector using a rotation matrix
#
# This example shows how to build a rotation matrix and rotate a vector.
# It uses :func:`simcoon.fillR_angle` and :func:`simcoon.rotate_vec_R`.


Rmat = sim.fillR_angle(angle, axis, active, copy)
v_rot1 = sim.rotate_vec_R(v, Rmat, copy)
print("Rotated vector (using R):", v_rot1)
print("Rotation matrix (using R):", Rmat)


# %%
# 2. Rotate a vector using angle/axis
#
# This example shows how to rotate a vector using an angle and an axis directly.
# It uses :func:`simcoon.rotate_vec_angle`.


v_rot2 = sim.rotate_vec_angle(v, angle, axis, copy)
print("Rotated vector (using angle/axis):", v_rot2)


# %%
# 3. Rotate a matrix using a rotation matrix
#
# This example shows how to rotate a matrix using a rotation matrix.
# It uses :func:`simcoon.rotate_mat_R`.


m_rot1 = sim.rotate_mat_R(m, Rmat, copy)
print("Rotated matrix (using R):\n", m_rot1)


# %%
# 4. Rotate a matrix using angle/axis
#
# This example shows how to rotate a matrix using an angle and an axis directly.
# It uses :func:`simcoon.rotate_mat_angle`.


m_rot2 = sim.rotate_mat_angle(m, angle, axis, copy)
print("Rotated matrix (using angle/axis):\n", m_rot2)


# %%
# 5. Build a rotation matrix from Euler angles
#
# This example shows how to create a rotation matrix from Euler angles.
# It uses :func:`simcoon.fillR_euler`.

psi, theta, phi = np.pi / 6, np.pi / 4, np.pi / 3
R_euler = sim.fillR_euler(psi, theta, phi, active, "zxz", copy)
print("Rotation matrix from Euler angles (zxz):\n", R_euler)


# %%
# Rotation of mechanical quantities
#
# These examples demonstrate how to rotate mechanical quantities such as stress,
# strain, and stiffness matrices.


# %%
# 6. Rotate a stress vector (single and batch)
#
# This example uses :func:`simcoon.rotate_stress_angle`.

stress = np.array([1, 2, 3, 4, 5, 6], dtype=float)
stress_batch = np.stack([stress, stress * 2], axis=1)
rot_stress1 = sim.rotate_stress_angle(stress, angle, axis, active, copy)
rot_stress2 = sim.rotate_stress_angle(stress_batch, angle, axis, active, copy)
print("Rotated stress (single):", rot_stress1)
print("Rotated stress (batch):\n", rot_stress2)


# %%
# 7. Rotate a strain vector (single and batch)
#
# This example uses :func:`simcoon.rotate_strain_angle`.

strain = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6], dtype=float)
strain_batch = np.stack([strain, strain * 2], axis=1)
rot_strain1 = sim.rotate_strain_angle(strain, angle, axis, active, copy)
rot_strain2 = sim.rotate_strain_angle(strain_batch, angle, axis, active, copy)
print("Rotated strain (single):", rot_strain1)
print("Rotated strain (batch):\n", rot_strain2)


# %%
# 8. Rotate a stiffness matrix (L)
#
# This example uses both :func:`simcoon.rotateL_angle` and
# :func:`simcoon.rotateL_R`.

L6 = np.eye(6)
rotL1 = sim.rotateL_angle(L6, angle, axis, active, copy)
rotL2 = sim.rotateL_R(L6, Rmat, active, copy)
print("Rotated L (angle):\n", rotL1)
print("Rotated L (R):\n", rotL2)


# %%
# 9. Rotate a compliance matrix (M)
#
# This example uses both :func:`simcoon.rotateM_angle` and
# :func:`simcoon.rotateM_R`.

M6 = np.eye(6)
rotM1 = sim.rotateM_angle(M6, angle, axis, active, copy)
rotM2 = sim.rotateM_R(M6, Rmat, active, copy)
print("Rotated M (angle):\n", rotM1)
print("Rotated M (R):\n", rotM2)


# %%
# 10. Rotate a strain concentration matrix (A)
#
# This example uses both :func:`simcoon.rotateA_angle` and
# :func:`simcoon.rotateA_R`.

A6 = np.eye(6)
rotA1 = sim.rotateA_angle(A6, angle, axis, active, copy)
rotA2 = sim.rotateA_R(A6, Rmat, active, copy)
print("Rotated A (angle):\n", rotA1)
print("Rotated A (R):\n", rotA2)


# %%
# 11. Rotate a stress concentration matrix (B)
#
# This example uses both :func:`simcoon.rotateB_angle` and
# :func:`simcoon.rotateB_R`.

B6 = np.eye(6)
rotB1 = sim.rotateB_angle(B6, angle, axis, active, copy)
rotB2 = sim.rotateB_R(B6, Rmat, active, copy)
print("Rotated B (angle):\n", rotB1)
print("Rotated B (R):\n", rotB2)


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

rot_matrix = np.array(
    [[np.cos(alpha), -np.sin(alpha), 0.0], [np.sin(alpha), np.cos(alpha), 0], [0, 0, 1]]
)
# print(R)

L_rotate = sim.rotateL_R(L, rot_matrix)
L_rotate_angle = sim.rotateL_angle(L, alpha, axis=3)

print(np.array_str(L_rotate, suppress_small=True))
print(np.array_str(L_rotate_angle, suppress_small=True))
