"""
Plasticity with Chaboche Hardening Example
=============================================

This example demonstrates the Chaboche plasticity UMAT using the new Python Solver API.
"""

import numpy as np
import matplotlib.pyplot as plt
from simcoon.solver import Solver, Block, StepMeca

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
E = 140000.0         # Young's modulus (MPa)
nu = 0.3             # Poisson ratio
alpha = 1.0e-6       # Thermal expansion coefficient
sigma_Y = 62.859017  # Initial yield stress (MPa)
Q = 416.004456       # Isotropic hardening saturation
b = 4.788635         # Isotropic hardening rate
C_1 = 30382.293921   # First kinematic hardening modulus
D_1 = 172.425687     # First kinematic hardening rate
C_2 = 195142.490843  # Second kinematic hardening modulus
D_2 = 3012.614659    # Second kinematic hardening rate

props = np.array([E, nu, alpha, sigma_Y, Q, b, C_1, D_1, C_2, D_2])

###################################################################################
# Create loading path using the new Python Solver API
# ---------------------------------------------------
# Define a cyclic uniaxial loading to demonstrate the Bauschinger effect.

# Step 1: Tension to 1% strain
step1 = StepMeca(
    DEtot_end=np.array([0.01, 0, 0, 0, 0, 0]),
    Dsigma_end=np.array([0, 0, 0, 0, 0, 0]),
    control=['strain', 'stress', 'stress', 'stress', 'stress', 'stress'],
    Dn_init=100,
    Dn_mini=20,
    Dn_inc=200,
    time=1.0
)

# Step 2: Compression to -1% strain
step2 = StepMeca(
    DEtot_end=np.array([-0.02, 0, 0, 0, 0, 0]),  # -2% increment from +1%
    Dsigma_end=np.array([0, 0, 0, 0, 0, 0]),
    control=['strain', 'stress', 'stress', 'stress', 'stress', 'stress'],
    Dn_init=200,
    Dn_mini=40,
    Dn_inc=400,
    time=2.0
)

# Step 3: Tension back to +1% strain
step3 = StepMeca(
    DEtot_end=np.array([0.02, 0, 0, 0, 0, 0]),
    Dsigma_end=np.array([0, 0, 0, 0, 0, 0]),
    control=['strain', 'stress', 'stress', 'stress', 'stress', 'stress'],
    Dn_init=200,
    Dn_mini=40,
    Dn_inc=400,
    time=2.0
)

# Create block with material properties
block = Block(
    steps=[step1, step2, step3],
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
s11 = np.array([h.sigma[0] for h in history])
time_arr = np.array([i for i in range(len(history))])  # Increment counter as proxy for time
Wm = np.array([h.Wm[0] for h in history])
Wm_r = np.array([h.Wm[1] for h in history])
Wm_ir = np.array([h.Wm[2] for h in history])
Wm_d = np.array([h.Wm[3] for h in history])

###################################################################################
# Plotting the results
# ----------------------
#
# We plot the stress-strain hysteresis loop which shows the cyclic behavior
# including the Bauschinger effect from kinematic hardening.

fig = plt.figure()

# First subplot: Stress vs Strain (hysteresis loop)
ax1 = fig.add_subplot(1, 2, 1)
plt.grid(True)
plt.tick_params(axis="both", which="major", labelsize=15)
plt.xlabel(r"Strain $\varepsilon_{11}$", size=15)
plt.ylabel(r"Stress $\sigma_{11}$ (MPa)", size=15)
plt.plot(e11, s11, c="blue", label="Chaboche model")
plt.title("Stress-Strain Hysteresis Loop")
plt.legend(loc="best")

# Second subplot: Work terms vs Increment
ax2 = fig.add_subplot(1, 2, 2)
plt.grid(True)
plt.tick_params(axis="both", which="major", labelsize=15)
plt.xlabel("Increment", size=15)
plt.ylabel(r"$W_m$", size=15)
plt.plot(time_arr, Wm, c="black", label=r"$W_m$")
plt.plot(time_arr, Wm_r, c="green", label=r"$W_m^r$")
plt.plot(time_arr, Wm_ir, c="blue", label=r"$W_m^{ir}$")
plt.plot(time_arr, Wm_d, c="red", label=r"$W_m^d$")
plt.title("Work Terms")
plt.legend(loc="best")

plt.suptitle("EPCHA - Chaboche Plasticity with Kinematic Hardening")
plt.tight_layout()
plt.show()

###################################################################################
# Note on the Bauschinger effect
# ------------------------------
# The hysteresis loop shows the Bauschinger effect: upon load reversal, the
# material yields at a stress lower than the original yield stress due to
# the kinematic hardening (back stress) accumulation.

print("\nChaboche Model Results:")
print(f"Maximum tensile stress: {max(s11):.2f} MPa")
print(f"Maximum compressive stress: {min(s11):.2f} MPa")
print(f"Yield asymmetry (Bauschinger effect): {abs(max(s11)) - abs(min(s11)):.2f} MPa")
