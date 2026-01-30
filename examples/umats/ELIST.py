"""
Transversely Isotropic Elasticity Example
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This example demonstrates the transversely isotropic elastic UMAT using the new Python Solver API.
"""

import numpy as np
import matplotlib.pyplot as plt
from simcoon.solver import Solver, Block, StepMeca

###################################################################################
# In transversely isotropic elastic materials, there is a single axis of symmetry.
# The material behaves isotropically in the plane perpendicular to this axis (the transverse plane).
# Eight parameters are required:
#
# 1. The axis of symmetry (1, 2, or 3)
# 2. The longitudinal Young modulus :math:`E_L`
# 3. The transverse Young modulus :math:`E_T`
# 4. The Poisson ratio in the transverse-longitudinal plane :math:`\nu_{TL}`
# 5. The Poisson ratio in the transverse-transverse plane :math:`\nu_{TT}`
# 6. The shear modulus in the longitudinal-transverse plane :math:`G_{LT}`
# 7. The coefficient of thermal expansion in the longitudinal direction :math:`\alpha_L`
# 8. The coefficient of thermal expansion in the transverse direction :math:`\alpha_T`
#
# The elastic stiffness tensor for a transversely isotropic material with axis 1 as the symmetry axis
# is written in the Voigt notation formalism as:
#
# .. math::
#
#    \mathbf{L} = \begin{pmatrix}
#        L_{11} & L_{12} & L_{12} & 0 & 0 & 0 \\
#        L_{12} & L_{22} & L_{23} & 0 & 0 & 0 \\
#        L_{12} & L_{23} & L_{22} & 0 & 0 & 0 \\
#        0 & 0 & 0 & G_{TT} & 0 & 0 \\
#        0 & 0 & 0 & 0 & G_{LT} & 0 \\
#        0 & 0 & 0 & 0 & 0 & G_{LT}
#    \end{pmatrix}
#
# where :math:`G_{TT} = E_T / (2(1+\nu_{TT}))` is the shear modulus in the transverse plane.
#
# The thermal expansion tensor is:
#
# .. math::
#
#    \boldsymbol{\alpha} = \begin{pmatrix}
#        \alpha_L & 0 & 0 \\
#        0 & \alpha_T & 0 \\
#        0 & 0 & \alpha_T
#    \end{pmatrix}

umat_name = "ELIST"  # 5 character code for transversely isotropic elastic subroutine
nstatev = 1  # Number of internal variables

# Material parameters
axis = 1       # Symmetry axis
E_L = 4500.0   # Longitudinal Young's modulus (MPa)
E_T = 2300.0   # Transverse Young's modulus (MPa)
nu_TL = 0.05   # Poisson ratio (transverse-longitudinal)
nu_TT = 0.3    # Poisson ratio (transverse-transverse)
G_LT = 2700.0  # Shear modulus (longitudinal-transverse)
alpha_L = 1.0e-5  # Thermal expansion (longitudinal)
alpha_T = 2.5e-5  # Thermal expansion (transverse)

props = np.array([axis, E_L, E_T, nu_TL, nu_TT, G_LT, alpha_L, alpha_T])

###################################################################################
# Create loading path using the new Python Solver API
# ---------------------------------------------------
# We define a uniaxial tension test along the longitudinal direction (direction 1).

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
plt.title("ELIST - Transversely Isotropic Elasticity")
plt.legend(loc="best")

plt.show()

###################################################################################
# Verify transverse isotropy
# --------------------------

print("\nVerification of transversely isotropic behavior:")
print(f"Applied axial strain: {e11[-1]:.6f}")
print(f"Computed axial stress: {s11[-1]:.2f} MPa")
print(f"Expected stress (E_L * epsilon): {E_L * e11[-1]:.2f} MPa")
print(f"Transverse strain e22: {e22[-1]:.6f}")
print(f"Transverse strain e33: {e33[-1]:.6f}")
print(f"Poisson effect check (e22 ~ e33 for transverse isotropy): {np.isclose(e22[-1], e33[-1])}")
