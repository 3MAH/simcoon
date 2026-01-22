"""
Orthotropic Elasticity Example
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""

import numpy as np
import simcoon as sim
import matplotlib.pyplot as plt
import os

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
E_1 = 4500.0  # Young's modulus in direction 1 (MPa)
E_2 = 2300.0  # Young's modulus in direction 2 (MPa)
E_3 = 2700.0  # Young's modulus in direction 3 (MPa)
nu_12 = 0.06  # Poisson ratio 12
nu_13 = 0.08  # Poisson ratio 13
nu_23 = 0.30  # Poisson ratio 23
G_12 = 2200.0  # Shear modulus 12 (MPa)
G_13 = 2100.0  # Shear modulus 13 (MPa)
G_23 = 2400.0  # Shear modulus 23 (MPa)
alpha_1 = 1.0e-5  # Thermal expansion in direction 1
alpha_2 = 2.5e-5  # Thermal expansion in direction 2
alpha_3 = 2.2e-5  # Thermal expansion in direction 3

psi_rve = 0.0
theta_rve = 0.0
phi_rve = 0.0
solver_type = 0
corate_type = 1

props = np.array(
    [E_1, E_2, E_3, nu_12, nu_13, nu_23, G_12, G_13, G_23, alpha_1, alpha_2, alpha_3]
)

path_data = "data"
path_results = "results"
pathfile = "ELORT_path.txt"
outputfile = "results_ELORT.txt"

sim.solver(
    umat_name,
    props,
    nstatev,
    psi_rve,
    theta_rve,
    phi_rve,
    solver_type,
    corate_type,
    path_data,
    path_results,
    pathfile,
    outputfile,
)

###################################################################################
# Plotting the results
# ----------------------
#
# We plot the stress-strain curve in the loading direction (direction 1).

outputfile_macro = os.path.join(path_results, "results_ELORT_global-0.txt")

fig = plt.figure()

e11, e22, e33, e12, e13, e23, s11, s22, s33, s12, s13, s23 = np.loadtxt(
    outputfile_macro,
    usecols=(8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19),
    unpack=True,
)

plt.grid(True)
plt.xlabel(r"Strain $\varepsilon_{11}$")
plt.ylabel(r"Stress $\sigma_{11}$ (MPa)")
plt.plot(e11, s11, c="blue", label="Loading direction 1")
plt.legend(loc="best")

plt.show()
