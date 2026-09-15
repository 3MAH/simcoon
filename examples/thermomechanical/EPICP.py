"""
Plasticity with isotropic hardening (thermomechanical)
========================================================
"""

import numpy as np
import matplotlib.pyplot as plt
import simcoon as sim
import os

plt.rcParams["figure.figsize"] = (18, 10)

###################################################################################
# The elastic-plastic (isotropic hardening) constitutive law implemented in simcoon
# is a rate independent, isotropic, von Mises type material with power-law isotropic
# hardening. Eight parameters are required for the thermomechanical version:
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
# The constitutive law is given by:
#
# .. math::
#
#   {\sigma}_{ij} & = L_{ijkl}\left({\varepsilon}^{\textrm{tot}}_{kl}-\alpha_{kl}\left(T-T^{\textrm{ref}}\right)-{\varepsilon}^{\textrm{p}}_{kl}\right) \\\\
#   \dot{\varepsilon}^{\textrm{p}}_{ij} & =\dot{p}\Lambda_{ij}, \quad \Lambda_{ij}=\frac{3}{2}\frac{\sigma'_{ij}}{\overline{\sigma}}, \quad \overline{\sigma}=\sqrt{\frac{3}{2}\sigma'_{kl}\sigma'_{kl}}, \\\\
#   \Phi & =\overline{\sigma}-\sigma_{Y}-kp^m\leq 0, \quad \dot{p}\geq0,~~~ \dot{p}~\Phi=0
#
# The updated work terms and internal heat production :math:`r` are determined
# with the thermomechanical algorithm.

umat_name = "EPICP"  # 5 character code for the elastic-plastic subroutine
nstatev = 8  # Number of internal variables

# Material parameters
rho = 4.4  # Density
c_p = 0.656  # Specific heat capacity
E = 113800.0  # Young's modulus (MPa)
nu = 0.342  # Poisson ratio
alpha = 0.86e-5  # Thermal expansion coefficient
sigma_Y = 500.0  # Yield stress (MPa)
H = 1600.0  # Hardening parameter
beta = 0.25  # Hardening exponent

psi_rve = 0.0
theta_rve = 0.0
phi_rve = 0.0
solver_type = 0
corate_type = 2

# Define the properties
props = np.array([rho, c_p, E, nu, alpha, sigma_Y, H, beta])

path_data = "../data"

# Run the simulation: the path file is parsed in Python and the case runs in
# memory, so no result file is written.
pathfile = "THERM_EPICP_path.json"

blocks, T_init, _ = sim.solver.load_path_json(os.path.join(path_data, pathfile))
res = sim.solver.solve(
    blocks,
    umat_name,
    props,
    nstatev,
    T_init=T_init,
    solver_type=solver_type,
    corate=corate_type,
    orientation=(psi_rve, theta_rve, phi_rve),
)

###################################################################################
# Plotting the results
# ----------------------
#
# We plot the stress-strain curve, the temperature evolution, the mechanical work
# terms (:math:`W_m`, :math:`W_m^r`, :math:`W_m^{ir}`, :math:`W_m^d`) and the
# thermal work terms (:math:`W_t`, :math:`W_t^r`, :math:`W_t^{ir}`).

fig = plt.figure()

# Get the data
e11, e22, e33, e12, e13, e23 = res["Strain"]
s11, s22, s33, s12, s13, s23 = res["Stress"]
time, T, Q, r = res["Time"], res["Temp"], res["Q"], res["r"]
Wm, Wm_r, Wm_ir, Wm_d = res["Wm"]
Wt, Wt_r, Wt_ir = res["Wt"]

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

###################################################################################
# Increment size effect
# -----------------------
#
# Here we test the effect of the increment size on the results. We run the same
# simulation with 1, 10, 100, and 1000 increments.

increments = [1, 10, 100, 1000]

# Run each case and collect its history: every path file is parsed in Python and
# the case runs in memory, so nothing is written to — or read back from — disk.
data = []
for inc in increments:
    pathfile = f"THERM_EPICP_path_{inc}.json"
    blocks, T_init, _ = sim.solver.load_path_json(os.path.join(path_data, pathfile))
    res_inc = sim.solver.solve(
        blocks,
        umat_name,
        props,
        nstatev,
        T_init=T_init,
        solver_type=solver_type,
        corate=corate_type,
        orientation=(psi_rve, theta_rve, phi_rve),
    )
    Wm_i, Wm_r_i, Wm_ir_i, Wm_d_i = res_inc["Wm"]
    Wt_i, Wt_r_i, Wt_ir_i = res_inc["Wt"]
    data.append(
        {
            "e11": res_inc["Strain"][0], "s11": res_inc["Stress"][0],
            "time": res_inc["Time"], "T": res_inc["Temp"],
            "Wm": Wm_i, "Wm_r": Wm_r_i, "Wm_ir": Wm_ir_i, "Wm_d": Wm_d_i,
            "Wt": Wt_i, "Wt_r": Wt_r_i, "Wt_ir": Wt_ir_i,
        }
    )

###################################################################################
# Plotting the increment size comparison
# -----------------------------------------
#
# We compare the stress-strain curves, the temperature evolution, the mechanical
# work terms and the thermal work terms for different increment sizes.

fig = plt.figure()

markers = ["D", "o", "x", None]
labels = ["1 increment", "10 increments", "100 increments", "1000 increments"]

# Stress vs Strain
ax = fig.add_subplot(2, 2, 1)
plt.grid(True)
plt.tick_params(axis="both", which="major", labelsize=15)
plt.xlabel(r"Strain $\varepsilon_{11}$", size=15)
plt.ylabel(r"Stress $\sigma_{11}$ (MPa)", size=15)
for i, d in enumerate(data):
    if markers[i] is not None:
        plt.plot(d["e11"], d["s11"], linestyle="None", marker=markers[i],
                 color="black", markersize=10, label=labels[i])
    else:
        plt.plot(d["e11"], d["s11"], c="black", label=labels[i])
plt.legend(loc="best")

# Temperature vs Time
ax = fig.add_subplot(2, 2, 2)
plt.grid(True)
plt.tick_params(axis="both", which="major", labelsize=15)
plt.xlabel("time (s)", size=15)
plt.ylabel(r"Temperature $\theta$ (K)", size=15)
for i, d in enumerate(data):
    if markers[i] is not None:
        plt.plot(d["time"], d["T"], linestyle="None", marker=markers[i],
                 color="black", markersize=10, label=labels[i])
    else:
        plt.plot(d["time"], d["T"], c="black", label=labels[i])
plt.legend(loc="best")

# Mechanical work vs Time
ax = fig.add_subplot(2, 2, 3)
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
            plt.plot(d["time"], d[wk], linestyle="None", marker=markers[i],
                     color=wc, markersize=10)
        else:
            plt.plot(d["time"], d[wk], c=wc, label=wl)
plt.legend(loc="best")

# Thermal work vs Time
ax = fig.add_subplot(2, 2, 4)
plt.grid(True)
plt.tick_params(axis="both", which="major", labelsize=15)
plt.xlabel("time (s)", size=15)
plt.ylabel(r"$W_t$", size=15)
therm_keys = ["Wt", "Wt_r", "Wt_ir"]
therm_labels = [r"$W_t$", r"$W_t^r$", r"$W_t^{ir}$"]
therm_colors = ["black", "green", "blue"]
for i, d in enumerate(data):
    for j, (wk, wc, wl) in enumerate(zip(therm_keys, therm_colors, therm_labels)):
        if markers[i] is not None:
            plt.plot(d["time"], d[wk], linestyle="None", marker=markers[i],
                     color=wc, markersize=10)
        else:
            plt.plot(d["time"], d[wk], c=wc, label=wl)
plt.legend(loc="best")

plt.show()
