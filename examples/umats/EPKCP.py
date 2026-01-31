"""
Plasticity with Isotropic and Kinematic Hardening Example
============================================================

This example demonstrates the combined isotropic-kinematic hardening UMAT
using the new Python Solver API.
"""

import numpy as np
import matplotlib.pyplot as plt
from simcoon.solver import Solver, Block, StepMeca

plt.rcParams["figure.figsize"] = (18, 10)

###################################################################################
# This plasticity model combines both isotropic and kinematic hardening with a
# power-law isotropic hardening and linear kinematic hardening.
#
# Seven parameters are required:
#
# 1. The Young modulus :math:`E`
# 2. The Poisson ratio :math:`\nu`
# 3. The coefficient of thermal expansion :math:`\alpha`
# 4. The initial yield stress :math:`\sigma_Y`
# 5. The isotropic hardening parameter :math:`k`
# 6. The isotropic hardening exponent :math:`m`
# 7. The kinematic hardening modulus :math:`k_X`
#
# The constitutive law is given by:
#
# .. math::
#
#   {\sigma}_{ij} & = L_{ijkl}\left({\varepsilon}^{\textrm{tot}}_{kl}-\alpha_{kl}\left(T-T^{\textrm{ref}}\right)-{\varepsilon}^{\textrm{p}}_{kl}\right) \\
#   \dot{\varepsilon}^{\textrm{p}}_{ij} & =\dot{p}\Lambda_{ij}, \quad \Lambda_{ij}=\frac{3}{2}\frac{\sigma'_{ij} - X_{ij}}{\overline{\sigma - X}} \\
#   \dot{X}_{ij} & = \frac{2}{3} k_X \dot{\varepsilon}^{\textrm{p}}_{ij} \\
#   \Phi & =\overline{\sigma - X}-\sigma_{Y}-kp^m\leq 0
#
# where :math:`X_{ij}` is the kinematic hardening (back stress) tensor resulting from
# linear kinematic hardening, and the isotropic hardening follows a power-law
# evolution :math:`kp^m`.

umat_name = "EPKCP"  # 5 character code for combined isotropic-kinematic hardening
nstatev = 14  # Number of internal variables

# Material parameters
E = 67538.0        # Young's modulus (MPa)
nu = 0.349         # Poisson ratio
alpha = 1.0e-6     # Thermal expansion coefficient
sigma_Y = 300.0    # Initial yield stress (MPa)
k = 1500.0         # Isotropic hardening parameter
m = 0.3            # Isotropic hardening exponent
k_X = 2000.0       # Kinematic hardening modulus

props = np.array([E, nu, alpha, sigma_Y, k, m, k_X])

###################################################################################
# Create loading path using the new Python Solver API
# ---------------------------------------------------
# Define a uniaxial loading path.

step = StepMeca(
    DEtot_end=np.array([0.03, 0, 0, 0, 0, 0]),  # 3% strain
    Dsigma_end=np.array([0, 0, 0, 0, 0, 0]),
    control=['strain', 'stress', 'stress', 'stress', 'stress', 'stress'],
    Dn_init=300,
    Dn_mini=75,
    Dn_inc=600,
    time=1.0
)

block = Block(
    steps=[step],
    umat_name=umat_name,
    props=props,
    nstatev=nstatev,
    control_type='small_strain',
    corate_type='jaumann'
)

# Run the simulation
solver = Solver(blocks=[block])
history = solver.solve()

###################################################################################
# Extract results from history
# ----------------------------

e11 = np.array([h.Etot[0] for h in history])
s11 = np.array([h.sigma[0] for h in history])
time_arr = np.linspace(0, 1, len(history))
Wm = np.array([h.Wm[0] for h in history])
Wm_r = np.array([h.Wm[1] for h in history])
Wm_ir = np.array([h.Wm[2] for h in history])
Wm_d = np.array([h.Wm[3] for h in history])

###################################################################################
# Plotting the results
# ----------------------
#
# We plot the stress-strain curve showing both isotropic and kinematic hardening.

fig = plt.figure()

# First subplot: Stress vs Strain
ax1 = fig.add_subplot(1, 2, 1)
plt.grid(True)
plt.tick_params(axis="both", which="major", labelsize=15)
plt.xlabel(r"Strain $\varepsilon_{11}$", size=15)
plt.ylabel(r"Stress $\sigma_{11}$ (MPa)", size=15)
plt.plot(e11, s11, c="blue", label="EPKCP model")
plt.title("Stress-Strain Response")
plt.legend(loc="best")

# Second subplot: Work terms vs Time
ax2 = fig.add_subplot(1, 2, 2)
plt.grid(True)
plt.tick_params(axis="both", which="major", labelsize=15)
plt.xlabel("time (s)", size=15)
plt.ylabel(r"$W_m$", size=15)
plt.plot(time_arr, Wm, c="black", label=r"$W_m$")
plt.plot(time_arr, Wm_r, c="green", label=r"$W_m^r$")
plt.plot(time_arr, Wm_ir, c="blue", label=r"$W_m^{ir}$")
plt.plot(time_arr, Wm_d, c="red", label=r"$W_m^d$")
plt.title("Work Terms")
plt.legend(loc="best")

plt.suptitle("EPKCP - Combined Isotropic and Kinematic Hardening")
plt.tight_layout()
plt.show()

###################################################################################
# Verify plastic behavior
# -----------------------

print("\nEPKCP Model Results:")
print(f"Maximum strain: {max(e11):.4f}")
print(f"Maximum stress: {max(s11):.2f} MPa")
print(f"Yield stress: {sigma_Y:.2f} MPa")
print(f"Hardening contribution: {max(s11) - sigma_Y:.2f} MPa")
