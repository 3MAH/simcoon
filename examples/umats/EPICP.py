"""
Plasticity with isotropic hardening example
=============================================

This example demonstrates the elastic-plastic UMAT with isotropic hardening
using the new Python Solver API.
"""

import numpy as np
import matplotlib.pyplot as plt
from simcoon.solver import Solver, Block, StepMeca

plt.rcParams["figure.figsize"] = (18, 10)  # configure the figure output size

plt.rc("text", usetex=True)
plt.rc("font", family="serif")

# ###################################################################################
# The elastic-plastic (isotropic hardening) constitutive law implemented in Simcoon is a rate independent, isotropic, von Mises type material with power-law isotropic hardening.
# Eight parameters are required for the thermomechanical version:
# The parameters required are:
#
# 1. The density :math:`\rho`
# 2. The specific heat :math:`c_p`
# 3. The Young modulus :math:`E`
# 4. The Poisson ratio :math:`\nu`
# 5. The coefficient of thermal expansion :math:`\alpha`
# 6. The von Mises equivalent yield stress limit :math:`\sigma_{Y}`
# 7. The hardening parameter :math:`k`
# 8. The hardening exponent :math:`m`
#
# The constitutive law is given by the set of equations :
#
# .. math::
#
#   {\sigma}_{ij} & = L_{ijkl}\left({\varepsilon}^{\textrm{tot}}_{kl}-\alpha_{kl}\left(T-T^{\textrm{ref}}\right)-{\varepsilon}^{\textrm{p}}_{kl}\right) \\\\
#   \dot{\varepsilon}^{\textrm{p}}_{ij} & =\dot{p}\Lambda_{ij}, \quad \Lambda_{ij}=\frac{3}{2}\frac{\sigma'_{ij}}{\overline{\sigma}}, \quad \sigma'_{ij}=\sigma_{ij}-\frac{1}{3}\sigma_{kk}\delta_{ij}, \quad \overline{\sigma}=\sqrt{\frac{3}{2}\sigma'_{kl}\sigma'_{kl}}, \\\\
#   \Phi & =\overline{\sigma}-\sigma_{Y}-kp^m\leq 0, \quad \dot{p}\geq0,~~~ \dot{p}~\Phi=0
#
# where :math:`{\varepsilon}^{\textrm{p}}_{ij}` is the plastic strain tensor, :math:`p` is the plastic multiplier,
# :math:`\sigma'_{ij}` is the deviatoric part of the stress and :math:`\overline{\sigma}` is the von Mises equivalent
# stress (Lemaitre and Chaboche, 2002). Moreover, :math:`T^{\textrm{ref}}` is a reference temperature
# (usually the temperature at the beginning of the analysis).
#
# In Simcoon the elastoplastic material constitutive law is implemented using a *return mapping algorithm*,
# with use of the *convex cutting plane* algorithm (Simo and Hughes, 1998). The updated stress is provided for 1D,
# plane stress, and generalized plane strain/3D analysis according to the forms of elastic isotropic materials.
#
# The updated work, and internal heat production :math:`r` are determined with the algorithm presented in the *simcoon* documentation.
#
# As a start we should input the name of the UMAT as well as the list of parameters

umat_name = "EPICP"  # This is the 5 character code for the elastic-plastic subroutine
nstatev = 8  # The number of scalar variables required

E = 113800       # Young's modulus (MPa)
nu = 0.342       # Poisson ratio
alpha = 0.86e-5  # Thermal expansion coefficient
sigma_Y = 600    # Yield stress (MPa)
H = 1600         # Hardening parameter
beta = 0.25      # Hardening exponent

# Define the properties
props = np.array([E, nu, alpha, sigma_Y, H, beta])

# ###################################################################################
# Create loading path using the new Python Solver API
# ---------------------------------------------------
# Define a uniaxial tension-compression cycle

# Step 1: Tension to 2% strain
step1 = StepMeca(
    DEtot_end=np.array([0.02, 0, 0, 0, 0, 0]),
    Dsigma_end=np.array([0, 0, 0, 0, 0, 0]),
    control=['strain', 'stress', 'stress', 'stress', 'stress', 'stress'],
    Dn_init=200,
    Dn_mini=50,
    Dn_inc=400,
    time=1.0
)

# Create block with material properties
block = Block(
    steps=[step1],
    umat_name=umat_name,
    props=props,
    nstatev=nstatev,
    control_type='small_strain',
    corate_type='jaumann'
)

# Run the simulation
solver = Solver(blocks=[block])
history = solver.solve()

# ###################################################################################
# Plotting the results
# --------------------------------------
# This is it, now we just need to plot the results.
# In the left, we plot the stress vs strain curve, and in the right the different work terms vs time:
# - :meth:`Wm <simcoon.Wm>` the mechanical work,
# - :meth:`Wm_r <simcoon.Wm_r>` the recoverable mechanical work,
# - :meth:`Wm_ir <simcoon.Wm_ir>` the irrecoverable mechanical work,
# - :meth:`Wm_d <simcoon.Wm_d>` the dissipated mechanical work.
# ###################################################################################

# Extract data from history
e11 = np.array([h.Etot[0] for h in history])
s11 = np.array([h.sigma[0] for h in history])
time = np.linspace(0, 1, len(history))
Wm = np.array([h.Wm[0] for h in history])
Wm_r = np.array([h.Wm[1] for h in history])
Wm_ir = np.array([h.Wm[2] for h in history])
Wm_d = np.array([h.Wm[3] for h in history])

# Plot the results
fig = plt.figure()

ax = fig.add_subplot(1, 2, 1)
plt.grid(True)
plt.tick_params(axis="both", which="major", labelsize=15)
plt.xlabel(r"Strain $\varepsilon_{11}$", size=15)
plt.ylabel(r"Stress $\sigma_{11}$\,(MPa)", size=15)
plt.plot(e11, s11, c="black", label="direction 1")
plt.legend(loc=2)

ax = fig.add_subplot(1, 2, 2)
plt.grid(True)
plt.tick_params(axis="both", which="major", labelsize=15)
plt.xlabel("time (s)", size=15)
plt.ylabel(r"$W_m$", size=15)
plt.plot(time, Wm, c="black", label=r"$W_m$")
plt.plot(time, Wm_r, c="green", label=r"$W_m^r$")
plt.plot(time, Wm_ir, c="blue", label=r"$W_m^{ir}$")
plt.plot(time, Wm_d, c="red", label=r"$W_m^d$")
plt.legend(loc=2)

plt.suptitle("EPICP - Plasticity with Isotropic Hardening")
plt.tight_layout()
plt.show()

# ###################################################################################
# Here we test the increment size effect on the results
# ----------------------------------------------------------
# ###################################################################################

# Define different increment counts (starting at 50 to ensure convergence)
increments = [50, 100, 200, 500]
data = []

for inc in increments:
    step = StepMeca(
        DEtot_end=np.array([0.02, 0, 0, 0, 0, 0]),
        Dsigma_end=np.array([0, 0, 0, 0, 0, 0]),
        control=['strain', 'stress', 'stress', 'stress', 'stress', 'stress'],
        Dn_init=inc,
        Dn_mini=max(10, inc // 5),
        Dn_inc=inc * 2,
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

    solver = Solver(blocks=[block])
    history = solver.solve()

    e11 = np.array([h.Etot[0] for h in history])
    s11 = np.array([h.sigma[0] for h in history])
    time_arr = np.linspace(0, 1, len(history))
    Wm = np.array([h.Wm[0] for h in history])
    Wm_r = np.array([h.Wm[1] for h in history])
    Wm_ir = np.array([h.Wm[2] for h in history])
    Wm_d = np.array([h.Wm[3] for h in history])

    data.append({
        "e11": e11,
        "s11": s11,
        "time": time_arr,
        "Wm": Wm,
        "Wm_r": Wm_r,
        "Wm_ir": Wm_ir,
        "Wm_d": Wm_d,
    })

# ###################################################################################
# Plotting the results
# --------------------------------------
#
# In the left, we plot the stress vs strain curve, and in the right the different work terms vs time
# Note the ["D", "o", "x", None] markers used to differentiate the different increment sizes:
# ["1 increment", "10 increments", "100 increments", "1000 increments"]
#
# ###################################################################################

fig = plt.figure()

markers = ["D", "o", "x", None]
labels = ["50 increments", "100 increments", "200 increments", "500 increments"]
colors = ["black", "black", "black", "black"]

# First subplot: Stress vs Strain
ax1 = fig.add_subplot(1, 2, 1)
plt.grid(True)
plt.tick_params(axis="both", which="major", labelsize=15)
plt.xlabel(r"Strain $\varepsilon_{11}$", size=15)
plt.ylabel(r"Stress $\sigma_{11}$\,(MPa)", size=15)

for i, d in enumerate(data):
    if markers[i] is not None:
        plt.plot(
            d["e11"],
            d["s11"],
            linestyle="None",
            marker=markers[i],
            color=colors[i],
            markersize=10,
            label=labels[i],
        )
    else:
        plt.plot(d["e11"], d["s11"], c=colors[i], label=labels[i])
plt.legend(loc=2)

# Second subplot: Work terms vs Time
ax2 = fig.add_subplot(1, 2, 2)
plt.grid(True)
plt.tick_params(axis="both", which="major", labelsize=15)
plt.xlabel("time (s)", size=15)
plt.ylabel(r"$W_m$", size=15)

work_colors = ["black", "green", "blue", "red"]
work_keys = ["Wm", "Wm_r", "Wm_ir", "Wm_d"]
work_labels = [r"$W_m$", r"$W_m^r$", r"$W_m^{ir}$", r"$W_m^d$"]

for i, d in enumerate(data):
    for j, (wk, wc, wl) in enumerate(zip(work_keys, work_colors, work_labels)):
        if markers[i] is not None:
            plt.plot(
                d["time"],
                d[wk],
                linestyle="None",
                marker=markers[i],
                color=wc,
                markersize=10,
                label=wl if i == len(data) - 1 else None,  # Only label once
            )
        else:
            plt.plot(d["time"], d[wk], c=wc, label=wl)
plt.legend(loc=2)

plt.suptitle("Increment Size Effect on EPICP Results")
plt.tight_layout()
plt.show()
