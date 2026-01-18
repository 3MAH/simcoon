"""
Plasticity with Chaboche Hardening Example
=============================================
"""

import numpy as np
import matplotlib.pyplot as plt
from simcoon import simmit as sim
import os

plt.rcParams["figure.figsize"] = (18, 10)

###################################################################################
# The Chaboche plasticity model combines isotropic and kinematic hardening.
# This model is particularly suited for cyclic loading applications where the
# Bauschinger effect is important.
#
# Ten parameters are required:
#
# 1. The Young modulus :math:`E`
# 2. The Poisson ratio :math:`\nu`
# 3. The coefficient of thermal expansion :math:`\alpha`
# 4. The initial yield stress :math:`\sigma_Y`
# 5. The isotropic hardening saturation value :math:`Q`
# 6. The isotropic hardening rate :math:`b`
# 7. The first kinematic hardening modulus :math:`C_1`
# 8. The first kinematic hardening rate :math:`D_1`
# 9. The second kinematic hardening modulus :math:`C_2`
# 10. The second kinematic hardening rate :math:`D_2`
#
# The constitutive law is given by:
#
# .. math::
#
#   {\sigma}_{ij} & = L_{ijkl}\left({\varepsilon}^{\textrm{tot}}_{kl}-\alpha_{kl}\left(T-T^{\textrm{ref}}\right)-{\varepsilon}^{\textrm{p}}_{kl}\right) \\
#   \dot{\varepsilon}^{\textrm{p}}_{ij} & =\dot{p}\Lambda_{ij}, \quad \Lambda_{ij}=\frac{3}{2}\frac{\sigma'_{ij} - X_{ij}}{\overline{\sigma} - X} \\
#   \dot{X}_{ij} & = \sum_{k} \frac{2}{3} C_k \dot{\varepsilon}^{\textrm{p}}_{ij} - D_k X^{(k)}_{ij} \dot{p} \\
#   \dot{R} & = b(Q - R)\dot{p} \\
#   \Phi & =\overline{\sigma - X}-\sigma_{Y}-R\leq 0
#
# where :math:`X_{ij}` is the kinematic hardening (back stress) tensor and :math:`R` is the
# isotropic hardening variable.

umat_name = "EPCHA"  # 5 character code for Chaboche plasticity
nstatev = 33  # Number of internal variables

# Material parameters
E = 140000.0  # Young's modulus (MPa)
nu = 0.3  # Poisson ratio
alpha = 1.0e-6  # Thermal expansion coefficient
sigma_Y = 62.859017  # Initial yield stress (MPa)
Q = 416.004456  # Isotropic hardening saturation
b = 4.788635  # Isotropic hardening rate
C_1 = 30382.293921  # First kinematic hardening modulus
D_1 = 172.425687  # First kinematic hardening rate
C_2 = 195142.490843  # Second kinematic hardening modulus
D_2 = 3012.614659  # Second kinematic hardening rate

psi_rve = 0.0
theta_rve = 0.0
phi_rve = 0.0
solver_type = 0
corate_type = 1

props = np.array([E, nu, alpha, sigma_Y, Q, b, C_1, D_1, C_2, D_2])

path_data = "data"
path_results = "results"
pathfile = "EPCHA_path.txt"
outputfile = "results_EPCHA.txt"

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
# We plot the stress-strain hysteresis loop which shows the cyclic behavior
# including the Bauschinger effect from kinematic hardening.

outputfile_macro = os.path.join(path_results, "results_EPCHA_global-0.txt")

fig = plt.figure()

e11, e22, e33, e12, e13, e23, s11, s22, s33, s12, s13, s23 = np.loadtxt(
    outputfile_macro,
    usecols=(8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19),
    unpack=True,
)
time, T, Q_out, r = np.loadtxt(outputfile_macro, usecols=(4, 5, 6, 7), unpack=True)
Wm, Wm_r, Wm_ir, Wm_d = np.loadtxt(
    outputfile_macro, usecols=(20, 21, 22, 23), unpack=True
)

# First subplot: Stress vs Strain (hysteresis loop)
ax1 = fig.add_subplot(1, 2, 1)
plt.grid(True)
plt.tick_params(axis="both", which="major", labelsize=15)
plt.xlabel(r"Strain $\varepsilon_{11}$", size=15)
plt.ylabel(r"Stress $\sigma_{11}$ (MPa)", size=15)
plt.plot(e11, s11, c="blue", label="Chaboche model")
plt.legend(loc="best")

# Second subplot: Work terms vs Time
ax2 = fig.add_subplot(1, 2, 2)
plt.grid(True)
plt.tick_params(axis="both", which="major", labelsize=15)
plt.xlabel("time (s)", size=15)
plt.ylabel(r"$W_m$", size=15)
plt.plot(time, Wm, c="black", label=r"$W_m$")
plt.plot(time, Wm_r, c="green", label=r"$W_m^r$")
plt.plot(time, Wm_ir, c="blue", label=r"$W_m^{ir}$")
plt.plot(time, Wm_d, c="red", label=r"$W_m^d$")
plt.legend(loc="best")

plt.show()
