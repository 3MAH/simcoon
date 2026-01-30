"""
Orthotropic Elasticity Example
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This example demonstrates the orthotropic elastic UMAT using the new Python Solver API.
"""

import numpy as np
import matplotlib.pyplot as plt
from simcoon.solver import Solver, Block, StepMeca

###################################################################################
# In orthotropic elastic materials, there are three mutually perpendicular planes of symmetry.
# The material has different mechanical properties in each of the three directions.
# Twelve parameters are required:
#
# 1. The Young modulus in direction 1: :math:`E_1`
# 2. The Young modulus in direction 2: :math:`E_2`
# 3. The Young modulus in direction 3: :math:`E_3`
# 4. The Poisson ratio :math:`\nu_{12}`
# 5. The Poisson ratio :math:`\nu_{13}`
# 6. The Poisson ratio :math:`\nu_{23}`
# 7. The shear modulus :math:`G_{12}`
# 8. The shear modulus :math:`G_{13}`
# 9. The shear modulus :math:`G_{23}`
# 10. The coefficient of thermal expansion :math:`\alpha_1`
# 11. The coefficient of thermal expansion :math:`\alpha_2`
# 12. The coefficient of thermal expansion :math:`\alpha_3`
#
# The elastic stiffness tensor for an orthotropic material is written in the Voigt notation as:
#
# .. math::
#
#    \mathbf{L} = \begin{pmatrix}
#        L_{11} & L_{12} & L_{13} & 0 & 0 & 0 \\
#        L_{12} & L_{22} & L_{23} & 0 & 0 & 0 \\
#        L_{13} & L_{23} & L_{33} & 0 & 0 & 0 \\
#        0 & 0 & 0 & G_{23} & 0 & 0 \\
#        0 & 0 & 0 & 0 & G_{13} & 0 \\
#        0 & 0 & 0 & 0 & 0 & G_{12}
#    \end{pmatrix}
#
# The thermal expansion tensor is:
#
# .. math::
#
#    \boldsymbol{\alpha} = \begin{pmatrix}
#        \alpha_1 & 0 & 0 \\
#        0 & \alpha_2 & 0 \\
#        0 & 0 & \alpha_3
#    \end{pmatrix}

umat_name = "ELORT"  # 5 character code for orthotropic elastic subroutine
nstatev = 1  # Number of internal variables

# Material parameters
E_1 = 4500.0     # Young's modulus in direction 1 (MPa)
E_2 = 2300.0     # Young's modulus in direction 2 (MPa)
E_3 = 2700.0     # Young's modulus in direction 3 (MPa)
nu_12 = 0.06     # Poisson ratio 12
nu_13 = 0.08     # Poisson ratio 13
nu_23 = 0.30     # Poisson ratio 23
G_12 = 2200.0    # Shear modulus 12 (MPa)
G_13 = 2100.0    # Shear modulus 13 (MPa)
G_23 = 2400.0    # Shear modulus 23 (MPa)
alpha_1 = 1.0e-5   # Thermal expansion in direction 1
alpha_2 = 2.5e-5   # Thermal expansion in direction 2
alpha_3 = 2.2e-5   # Thermal expansion in direction 3

props = np.array(
    [E_1, E_2, E_3, nu_12, nu_13, nu_23, G_12, G_13, G_23, alpha_1, alpha_2, alpha_3]
)

###################################################################################
# Create loading path using the new Python Solver API
# ---------------------------------------------------
# We define a uniaxial tension test along direction 1.

step = StepMeca(
    DEtot_end=np.array([0.01, 0, 0, 0, 0, 0]),  # 1% strain in direction 1
    Dsigma_end=np.array([0, 0, 0, 0, 0, 0]),
    control=['strain', 'stress', 'stress', 'stress', 'stress', 'stress'],
    Dn_init=50,
    Dn_mini=10,
    Dn_inc=100,
    time=1.0
)

block = Block(
    steps=[step],
    umat_name=umat_name,
    props=props,
    nstatev=nstatev,
    control_type='small_strain',
    corate_type='green_naghdi'
)

# Run the simulation
solver = Solver(blocks=[block])
history = solver.solve()

###################################################################################
# Extract results from history
# ----------------------------

e11 = np.array([h.Etot[0] for h in history])
e22 = np.array([h.Etot[1] for h in history])
e33 = np.array([h.Etot[2] for h in history])
s11 = np.array([h.sigma[0] for h in history])
s22 = np.array([h.sigma[1] for h in history])
s33 = np.array([h.sigma[2] for h in history])

###################################################################################
# Plotting the results
# ----------------------
#
# We plot the stress-strain curve in the loading direction (direction 1).

fig = plt.figure()

plt.grid(True)
plt.xlabel(r"Strain $\varepsilon_{11}$")
plt.ylabel(r"Stress $\sigma_{11}$ (MPa)")
plt.plot(e11, s11, c="blue", label="Loading direction 1")
plt.title("ELORT - Orthotropic Elasticity")
plt.legend(loc="best")

plt.show()

###################################################################################
# Verify orthotropic behavior
# ---------------------------

print("\nVerification of orthotropic behavior:")
print(f"Applied axial strain: {e11[-1]:.6f}")
print(f"Computed axial stress: {s11[-1]:.2f} MPa")
print(f"Expected stress (E_1 * epsilon): {E_1 * e11[-1]:.2f} MPa")
print(f"Transverse strain e22: {e22[-1]:.6f}")
print(f"Transverse strain e33: {e33[-1]:.6f}")
print(f"Orthotropy check (e22 != e33 for orthotropic materials): {not np.isclose(e22[-1], e33[-1])}")
