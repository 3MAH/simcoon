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
print("Rotation matrix:", Rmat)


# %%
# 2. Rotate a vector using a different constructor
#
# This example shows how to create a rotation from Euler angles
# and rotate a vector using :meth:`simcoon.Rotation.apply`.


rot2 = sim.Rotation.from_axis_angle(angle, axis)
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
#
# This example shows how to rotate a matrix using a rotation
# created from Euler angles.


rot_euler = sim.Rotation.from_axis_angle(angle, axis)
m_rot2 = rot_euler.apply_tensor(m)
print("Rotated matrix:\n", m_rot2)


# %%
# 5. Build a rotation matrix from Euler angles
#
# This example shows how to create a rotation from Euler angles using the Rotation class.

psi, theta, phi = np.pi / 6, np.pi / 4, np.pi / 3
r_euler = sim.Rotation.from_euler(psi, theta, phi, "zxz")
R_euler = r_euler.as_matrix()
print("Rotation matrix from Euler angles (zxz):\n", R_euler)


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
