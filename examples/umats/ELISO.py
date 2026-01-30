"""
Isotropic elasticity examples
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This example demonstrates the isotropic elastic UMAT using the new Python Solver API.
"""

import numpy as np
import matplotlib.pyplot as plt
from simcoon.solver import Solver, Block, StepMeca

###################################################################################
# In thermoelastic isotropic materials three parameters are required:
#
# 1. The Young modulus :math:`E`,
# 2. The Poisson ratio :math:`\nu`,
# 3. The coefficient of thermal expansion :math:`\alpha`.
#
# The elastic stiffness tensor and the thermal expansion coefficients tensor are written in the Voigt notation formalism as
#
# .. math::
#
#    \mathbf{L} = \begin{pmatrix}
#        L_{1111} & L_{1122} & L_{1122} & 0 & 0 & 0 \\
#        L_{1122} & L_{1111} & L_{1122} & 0 & 0 & 0 \\
#        L_{1122} & L_{1122} & L_{1111} & 0 & 0 & 0 \\
#        0 & 0 & 0 & L_{1212} & 0 & 0 \\
#        0 & 0 & 0 & 0 & L_{1212} & 0 \\
#        0 & 0 & 0 & 0 & 0 & L_{1212}
#    \end{pmatrix},
#    \quad
#    \boldsymbol{\alpha} = \begin{pmatrix}
#        \alpha & 0 & 0 \\
#        0 & \alpha & 0 \\
#        0 & 0 & \alpha
#    \end{pmatrix}
#
# with
#
# .. math::
#
#    L_{1111} = \frac{E(1-\nu)}{(1+\nu)(1-2\nu)}, \quad
#    L_{1122} = \frac{E\nu}{(1+\nu)(1-2\nu)}, \quad
#    L_{1212} = \frac{E}{2(1+\nu)}.
#
# The tangent stiffness tensor in this case is :math:`\mathbf{L}^t = \mathbf{L}`.
# Moreover, the increment of the elastic strain is given by
#
# .. math::
#
#    \Delta\varepsilon^{\mathrm{el}}_{ij} = \Delta\varepsilon^{\mathrm{tot}}_{ij} - \alpha \Delta T \delta_{ij},
#
# where :math:`\delta_{ij}` implies the Kronecker delta operator.
# In the 1D case only one component of stress is computed, through the relation
#
# .. math::
#
#    \sigma^{\mathrm{fin}}_{11} = \sigma^{\mathrm{init}}_{11} + E \Delta\varepsilon^{\mathrm{el}}_{11}.
#
# In the plane stress case only three components of stress are computed, through the relations
#
# .. math::
#
#    \begin{pmatrix}
#        \sigma^{\mathrm{fin}}_{11} \\
#        \sigma^{\mathrm{fin}}_{22} \\
#        \sigma^{\mathrm{fin}}_{12}
#    \end{pmatrix}
#    =
#    \begin{pmatrix}
#        \sigma^{\mathrm{init}}_{11} \\
#        \sigma^{\mathrm{init}}_{22} \\
#        \sigma^{\mathrm{init}}_{12}
#    \end{pmatrix}
#    +
#    \frac{E}{1-\nu^2}
#    \begin{pmatrix}
#        1 & \nu & 0 \\
#        \nu & 1 & 0 \\
#        0 & 0 & \frac{1-\nu}{2}
#    \end{pmatrix}
#    \begin{pmatrix}
#        \Delta\varepsilon^{\mathrm{el}}_{11} \\
#        \Delta\varepsilon^{\mathrm{el}}_{22} \\
#        2\Delta\varepsilon^{\mathrm{el}}_{12}
#    \end{pmatrix}
#
# In the generalized plane strain/3D analysis case the stress tensor is computed through the relation
#
# .. math::
#
#    \sigma^{\mathrm{fin}}_{ij} = \sigma^{\mathrm{init}}_{ij} + L_{ijkl}~\Delta\varepsilon^{\mathrm{el}}_{kl}.


###################################################################################
# Define material properties
# --------------------------
# ELISO is the 5 character code for the elastic-isotropic subroutine

E = 700000.0      # Young's modulus (MPa)
nu = 0.2          # Poisson ratio
alpha = 1.0e-5    # Thermal expansion coefficient

props = np.array([E, nu, alpha])
nstatev = 1  # Number of internal state variables (only initial temperature stored)

###################################################################################
# Create loading path using the new Python Solver API
# ---------------------------------------------------
# We define a uniaxial tension test with strain control in direction 1
# and stress-free boundary conditions in the transverse directions.

# Uniaxial tension: strain in direction 1, stress-free in other directions
step = StepMeca(
    DEtot_end=np.array([0.01, 0, 0, 0, 0, 0]),  # 1% strain in direction 1
    Dsigma_end=np.array([0, 0, 0, 0, 0, 0]),    # Target stress increment (for stress-controlled)
    control=['strain', 'stress', 'stress', 'stress', 'stress', 'stress'],
    Dn_init=50,   # Number of increments
    Dn_mini=10,   # Minimum increments
    Dn_inc=100,   # Maximum increments
    time=1.0      # Time for this step
)

# Create block with material properties
block = Block(
    steps=[step],
    umat_name="ELISO",
    props=props,
    nstatev=nstatev,
    control_type='small_strain',
    corate_type='logarithmic'
)

###################################################################################
# Run the simulation
# ------------------

solver = Solver(blocks=[block])
history = solver.solve()

###################################################################################
# Extract results from history
# ----------------------------
# The history contains StateVariables objects at each converged increment

e11 = np.array([h.Etot[0] for h in history])
e22 = np.array([h.Etot[1] for h in history])
e33 = np.array([h.Etot[2] for h in history])
s11 = np.array([h.sigma[0] for h in history])
s22 = np.array([h.sigma[1] for h in history])
s33 = np.array([h.sigma[2] for h in history])

###################################################################################
# Plotting the results
# --------------------

fig = plt.figure()

plt.grid(True)
plt.plot(e11, s11, c="blue", label="Stress-strain response")
plt.xlabel("Strain")
plt.ylabel("Stress (MPa)")
plt.title("ELISO - Isotropic Elasticity (Uniaxial Tension)")
plt.legend()

plt.show()

###################################################################################
# Verify analytical solution
# --------------------------
# For uniaxial tension with isotropic elasticity:
# sigma_11 = E * epsilon_11

print(f"\nVerification:")
print(f"Applied strain: {e11[-1]:.6f}")
print(f"Computed stress: {s11[-1]:.2f} MPa")
print(f"Expected stress (E * epsilon): {E * e11[-1]:.2f} MPa")
print(f"Relative error: {abs(s11[-1] - E * e11[-1]) / (E * e11[-1]) * 100:.4f}%")
