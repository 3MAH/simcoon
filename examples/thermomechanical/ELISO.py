"""
Isotropic elasticity (thermomechanical)
=======================================
"""

import numpy as np
import simcoon as sim
import matplotlib.pyplot as plt
import os

plt.rcParams["figure.figsize"] = (18, 10)

###################################################################################
# In thermoelastic isotropic materials three mechanical parameters and two thermal
# parameters are required:
#
# 1. The density :math:`\rho`
# 2. The specific heat :math:`c_p`
# 3. The Young modulus :math:`E`
# 4. The Poisson ratio :math:`\nu`
# 5. The coefficient of thermal expansion :math:`\alpha`
#
# The elastic stiffness tensor is written in the Voigt notation formalism as
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
# The increment of the elastic strain is given by
#
# .. math::
#
#    \Delta\varepsilon^{\mathrm{el}}_{ij} = \Delta\varepsilon^{\mathrm{tot}}_{ij} - \alpha \Delta T \delta_{ij}
#
# In the thermomechanical framework, the thermal work terms :math:`W_t`, :math:`W_t^r`
# and :math:`W_t^{ir}` are also computed alongside the mechanical work terms.

umat_name = "ELISO"  # 5 character code for the elastic-isotropic subroutine
nstatev = 1  # Number of internal variables

# Material parameters
rho = 4.4  # Density
c_p = 0.656  # Specific heat capacity
E = 70000.0  # Young's modulus (MPa)
nu = 0.2  # Poisson ratio
alpha = 1.0e-5  # Thermal expansion coefficient

psi_rve = 0.0
theta_rve = 0.0
phi_rve = 0.0
solver_type = 0
corate_type = 2

props = np.array([rho, c_p, E, nu, alpha])

path_data = "../data"
path_results = "results"
pathfile = "THERM_ELISO_path.txt"
outputfile = "results_THERM_ELISO.txt"

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
# We plot the stress-strain curve, the temperature evolution, the mechanical work
# terms and the thermal work terms.

outputfile_macro = os.path.join(path_results, "results_THERM_ELISO_global-0.txt")

fig = plt.figure()

# Get the data
e11, e22, e33, e12, e13, e23, s11, s22, s33, s12, s13, s23 = np.loadtxt(
    outputfile_macro,
    usecols=(8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19),
    unpack=True,
)
time, T, Q, r = np.loadtxt(outputfile_macro, usecols=(4, 5, 6, 7), unpack=True)
Wm, Wm_r, Wm_ir, Wm_d, Wt, Wt_r, Wt_ir = np.loadtxt(
    outputfile_macro, usecols=(20, 21, 22, 23, 24, 25, 26), unpack=True
)

# Stress vs Strain
ax = fig.add_subplot(2, 2, 1)
plt.grid(True)
plt.tick_params(axis="both", which="major", labelsize=15)
plt.xlabel(r"Strain $\varepsilon_{11}$", size=15)
plt.ylabel(r"Stress $\sigma_{11}$ (MPa)", size=15)
plt.plot(e11, s11, c="black", label="direction 1")
plt.legend(loc="best")

# Temperature vs Time
ax = fig.add_subplot(2, 2, 2)
plt.grid(True)
plt.tick_params(axis="both", which="major", labelsize=15)
plt.xlabel("time (s)", size=15)
plt.ylabel(r"Temperature $\theta$ (K)", size=15)
plt.plot(time, T, c="black", label="temperature")
plt.legend(loc="best")

# Mechanical work vs Time
ax = fig.add_subplot(2, 2, 3)
plt.grid(True)
plt.tick_params(axis="both", which="major", labelsize=15)
plt.xlabel("time (s)", size=15)
plt.ylabel(r"$W_m$", size=15)
plt.plot(time, Wm, c="black", label=r"$W_m$")
plt.plot(time, Wm_r, c="green", label=r"$W_m^r$")
plt.plot(time, Wm_ir, c="blue", label=r"$W_m^{ir}$")
plt.plot(time, Wm_d, c="red", label=r"$W_m^d$")
plt.legend(loc="best")

# Thermal work vs Time
ax = fig.add_subplot(2, 2, 4)
plt.grid(True)
plt.tick_params(axis="both", which="major", labelsize=15)
plt.xlabel("time (s)", size=15)
plt.ylabel(r"$W_t$", size=15)
plt.plot(time, Wt, c="black", label=r"$W_t$")
plt.plot(time, Wt_r, c="green", label=r"$W_t^r$")
plt.plot(time, Wt_ir, c="blue", label=r"$W_t^{ir}$")
plt.legend(loc="best")

plt.show()
