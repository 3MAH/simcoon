"""
In-memory Python solver
=======================

Drive the simcoon material-point solver directly from Python: loading blocks
are built as objects, results come back as numpy arrays (fedoo-style layout),
and no path.txt / result files are involved.

Two demonstrations:

1. cyclic elastoplasticity (EPICP), stress-controlled;
2. superelastic SMA cycled under adiabatic conditions (coupled
   thermomechanical block): the transformation latent heat self-heats the
   material, which hardens the stress plateau.
"""

import matplotlib.pyplot as plt
import numpy as np

from simcoon import solver

###############################################################################
# Cyclic elastoplasticity (mechanical block)
# ------------------------------------------

props_epicp = np.array([70000.0, 0.3, 1.0e-5, 300.0, 1000.0, 0.3])
load = solver.StepMeca(control="stress", value=[450.0, 0, 0, 0, 0, 0], ninc=100)
unload = solver.StepMeca(control="stress", value=[0.0, 0, 0, 0, 0, 0], ninc=100)

res = solver.solve(
    solver.Block(steps=[load, unload], ncycle=3),
    "EPICP", props_epicp, nstatev=8,
)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4))
ax1.plot(100 * res["Strain"][0], res["Stress"][0], "-")
ax1.set_xlabel("strain $E_{11}$ (%)")
ax1.set_ylabel("stress $\\sigma_{11}$ (MPa)")
ax1.set_title("EPICP, 3 stress cycles")

###############################################################################
# Adiabatic superelastic SMA (thermomechanical block)
# ---------------------------------------------------
# Same material as examples/thermomechanical/SMA_T.py. The heat-flux control
# with Q = 0 makes the loading adiabatic: the exothermic forward
# transformation heats the sample, the reverse transformation cools it.

props_sma = np.array([
    5.5, 0.5, 0.5, 0,                 # rho, c_pA, c_pM, flagT
    67538.0, 67538.0, 0.349, 0.349,   # E_A, E_M, nu_A, nu_M
    1.0e-6, 1.0e-6,                   # alpha_A, alpha_M
    0.0, 0.0418, 0.021, 0.0,          # Hmin, Hmax, k1, sigmacrit
    10.0, 10.0,                       # C_A, C_M
    300.0, 290.0, 295.0, 305.0,       # Ms0, Mf0, As0, Af0
    0.2, 0.2, 0.2, 0.2,               # n1..n4
    300.0, 0.2, 2.0,                  # sigmacaliber, b_prager, n_prager
    1.0e-6, 1.0e-5, 1.0, 0.0,         # c_lambda, p0_lambda, n_lambda, alpha_lambda
])

loading = solver.StepThermomeca(
    control=["strain"] + ["stress"] * 5, value=[0.06, 0, 0, 0, 0, 0],
    ninc=300, time=30.0, thermal_control="heat_flux", Q=0.0,
)
unloading = solver.StepThermomeca(
    control=["strain"] + ["stress"] * 5, value=[0.0, 0, 0, 0, 0, 0],
    ninc=300, time=30.0, thermal_control="heat_flux", Q=0.0,
)

res_sma = solver.solve(
    solver.Block(steps=[loading, unloading]),
    "SMAUT", props_sma, nstatev=17, T_init=323.15,
)

e11 = 100 * res_sma["Strain"][0]
line = ax2.scatter(e11, res_sma["Stress"][0], c=res_sma["Temp"], s=4, cmap="coolwarm")
ax2.set_xlabel("strain $E_{11}$ (%)")
ax2.set_ylabel("stress $\\sigma_{11}$ (MPa)")
ax2.set_title("SMAUT, adiabatic superelastic cycle")
fig.colorbar(line, ax=ax2, label="temperature (K)")

fig.tight_layout()
plt.show()
