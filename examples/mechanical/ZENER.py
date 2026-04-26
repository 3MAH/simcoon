"""
Zener viscoelastic model
=========================
"""

import numpy as np
import matplotlib.pyplot as plt
import simcoon as sim
import os

plt.rcParams["figure.figsize"] = (18, 10)

###################################################################################
# The Poynting-Thomson (Zener) constitutive law is a rate-dependent, isotropic,
# linear viscoelastic model that accounts for thermal strains. It consists of
# an elastic spring in parallel with a Maxwell element (spring + dashpot in series).
#
# Seven parameters are required:
#
# 1. The thermoelastic Young's modulus :math:`E_0`
# 2. The thermoelastic Poisson's ratio :math:`\nu_0`
# 3. The coefficient of thermal expansion :math:`\alpha`
# 4. The viscoelastic Young's modulus of the Zener branch :math:`E_1`
# 5. The viscoelastic Poisson's ratio of the Zener branch :math:`\nu_1`
# 6. The bulk viscosity of the Zener branch :math:`\eta_B`
# 7. The shear viscosity of the Zener branch :math:`\eta_S`
#
# The viscoelastic material constitutive law is implemented using a
# *fast scalar updating method*. The updated stress is provided for 1D,
# plane stress, and generalized plane strain/3D analysis.

umat_name = "ZENER"  # 5 character code for the Zener model
nstatev = 8  # Number of internal variables

# Material parameters
E_0 = 3000.0  # Thermoelastic Young's modulus (MPa)
nu_0 = 0.4  # Thermoelastic Poisson's ratio
alpha = 0.0  # Thermal expansion coefficient
E_1 = 100.0  # Viscoelastic Young's modulus (MPa)
nu_1 = 0.3  # Viscoelastic Poisson's ratio
eta_S = 4000.0  # Shear viscosity
eta_B = eta_S / 4.0  # Bulk viscosity

psi_rve = 0.0
theta_rve = 0.0
phi_rve = 0.0
solver_type = 0
corate_type = 1

props = np.array([E_0, nu_0, alpha, E_1, nu_1, eta_B, eta_S])

path_data = "../data"
path_results = "results"
pathfile = "ZENER_path.txt"
outputfile = "results_ZENER.txt"

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
# We plot the stress-strain response which exhibits the characteristic
# rate-dependent behavior of the Zener viscoelastic model.

outputfile_macro = os.path.join(path_results, "results_ZENER_global-0.txt")

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

# First subplot: Stress vs Strain
ax1 = fig.add_subplot(1, 2, 1)
plt.grid(True)
plt.tick_params(axis="both", which="major", labelsize=15)
plt.xlabel(r"Strain $\varepsilon_{11}$", size=15)
plt.ylabel(r"Stress $\sigma_{11}$ (MPa)", size=15)
plt.plot(e11, s11, c="blue", label="Zener model")
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
