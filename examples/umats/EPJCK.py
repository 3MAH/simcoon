"""
Johnson-Cook plasticity example
================================

Demonstrates the Johnson-Cook constitutive law (EPJCK) with strain rate sensitivity.
The model is:

    sigma_Y = (A + B * p^n) * (1 + C * ln(edot/edot0)) * (1 - T*^m)

where T* = (T - T_ref) / (T_melt - T_ref).

This example runs three simulations at different strain rates to show rate sensitivity:
- Slow:   edot ~ 0.01 s^-1  (time = 10s for 10% strain)
- Medium: edot ~ 0.1 s^-1   (time = 1s for 10% strain)
- Fast:   edot ~ 100 s^-1   (time = 0.001s for 10% strain)
"""

import numpy as np
import matplotlib.pyplot as plt
import simcoon as sim
import os

plt.rcParams["figure.figsize"] = (18, 10)
dir = os.path.dirname(os.path.realpath("__file__"))

plt.rc("text", usetex=True)
plt.rc("font", family="serif")

# ###################################################################################
# Material parameters: AISI 4340 Steel (Johnson & Cook, 1985)
# ###################################################################################

umat_name = "EPJCK"
nstatev = 9

E = 200000        # Young's modulus (MPa)
nu = 0.33         # Poisson's ratio
alpha = 1.0e-5    # CTE
A_jc = 792        # JC initial yield stress (MPa)
B_jc = 510        # JC hardening coefficient (MPa)
n_jc = 0.26       # JC hardening exponent
C_jc = 0.014      # JC strain rate sensitivity
edot0 = 1.0       # JC reference strain rate (s^-1)
m_jc = 1.03       # JC thermal softening exponent
T_ref = 293.0     # Reference temperature (K)
T_melt = 1793.0   # Melting temperature (K)

psi_rve = 0.0
theta_rve = 0.0
phi_rve = 0.0
solver_type = 0
corate_type = 2

props = np.array([E, nu, alpha, A_jc, B_jc, n_jc, C_jc, edot0, m_jc, T_ref, T_melt])
path_data = "data"
path_results = "results"

# ###################################################################################
# Run simulations at different strain rates
# ###################################################################################

rate_configs = {
    "slow": {"pathfile": "EPJCK_path_slow.txt", "outputfile": "results_EPJCK_slow.txt", "label": r"$\dot{\varepsilon} \approx 0.01$ s$^{-1}$"},
    "medium": {"pathfile": "EPJCK_path.txt", "outputfile": "results_EPJCK_medium.txt", "label": r"$\dot{\varepsilon} \approx 0.1$ s$^{-1}$"},
    "fast": {"pathfile": "EPJCK_path_fast.txt", "outputfile": "results_EPJCK_fast.txt", "label": r"$\dot{\varepsilon} \approx 100$ s$^{-1}$"},
}

colors = {"slow": "blue", "medium": "black", "fast": "red"}

fig = plt.figure()
ax1 = fig.add_subplot(1, 2, 1)
ax2 = fig.add_subplot(1, 2, 2)

for key, cfg in rate_configs.items():
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
        cfg["pathfile"],
        cfg["outputfile"],
    )

    outputfile_global = cfg["outputfile"].replace(".txt", "_global-0.txt")
    P_global = os.path.join(dir, "results", outputfile_global)

    e11, s11 = np.loadtxt(P_global, usecols=(8, 14), unpack=True)
    time, = np.loadtxt(P_global, usecols=(4,), unpack=True)
    Wm, Wm_r, Wm_ir, Wm_d = np.loadtxt(P_global, usecols=(20, 21, 22, 23), unpack=True)

    ax1.plot(e11, s11, c=colors[key], label=cfg["label"])
    ax2.plot(time, Wm_d, c=colors[key], label=cfg["label"])

ax1.set_xlabel(r"Strain $\varepsilon_{11}$", size=15)
ax1.set_ylabel(r"Stress $\sigma_{11}$ (MPa)", size=15)
ax1.set_title("Strain rate sensitivity (Johnson-Cook)", size=15)
ax1.grid(True)
ax1.legend(loc="best", fontsize=12)
ax1.tick_params(axis="both", which="major", labelsize=15)

ax2.set_xlabel("Time (s)", size=15)
ax2.set_ylabel(r"Dissipated work $W_m^d$ (MPa)", size=15)
ax2.set_title("Dissipated work vs time", size=15)
ax2.grid(True)
ax2.legend(loc="best", fontsize=12)
ax2.tick_params(axis="both", which="major", labelsize=15)

plt.tight_layout()
plt.show()
