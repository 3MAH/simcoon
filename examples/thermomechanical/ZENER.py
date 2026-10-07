"""
Zener viscoelastic model (thermomechanical)
=============================================
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import simcoon as sim

plt.rcParams["figure.figsize"] = (18, 10)

###################################################################################
# The Poynting-Thomson (Zener) constitutive law is a rate-dependent, isotropic,
# linear viscoelastic model that accounts for thermal strains. It consists of
# an elastic spring in parallel with a Maxwell element (spring + dashpot in series).
#
# Nine parameters are required for the thermomechanical version:
#
# 1. The density :math:`\rho`
# 2. The specific heat :math:`c_p`
# 3. The thermoelastic Young's modulus :math:`E_0`
# 4. The thermoelastic Poisson's ratio :math:`\nu_0`
# 5. The coefficient of thermal expansion :math:`\alpha`
# 6. The viscoelastic Young's modulus of the Zener branch :math:`E_1`
# 7. The viscoelastic Poisson's ratio of the Zener branch :math:`\nu_1`
# 8. The bulk viscosity of the Zener branch :math:`\eta_B`
# 9. The shear viscosity of the Zener branch :math:`\eta_S`
#
# The viscoelastic material constitutive law is implemented using a
# *fast scalar updating method*. The updated stress is provided for 1D,
# plane stress, and generalized plane strain/3D analysis.
# The updated work terms and internal heat production :math:`r` are
# determined with the thermomechanical algorithm.

umat_name = "ZENER"  # 5 character code for the Zener model
nstatev = 8  # Number of internal variables

# Material parameters
rho = 4.4  # Density
c_p = 0.656  # Specific heat capacity
E_0 = 3000.0  # Thermoelastic Young's modulus (MPa)
nu_0 = 0.4  # Thermoelastic Poisson's ratio
alpha = 0.86e-5  # Thermal expansion coefficient
E_1 = 1200.0  # Viscoelastic Young's modulus (MPa)
nu_1 = 0.3  # Viscoelastic Poisson's ratio
eta_B = 12500.0  # Bulk viscosity
eta_S = 400.0  # Shear viscosity

psi_rve = 0.0
theta_rve = 0.0
phi_rve = 0.0
solver_type = 0
corate_type = 0

# Define the properties
props = np.array([rho, c_p, E_0, nu_0, alpha, E_1, nu_1, eta_B, eta_S])

path_data = "../data"

# Run the simulation: the path file is parsed in Python and the case runs in
# memory, so no result file is written.
pathfile = "THERM_ZENER_path.json"

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
# We plot the stress-strain response, the temperature evolution, the mechanical
# work terms and the thermal work terms.

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
# Here we test the effect of the increment size on the results.

increments = [1, 10, 100, 1000]

# Run each case and collect its history: every path file is parsed in Python and
# the case runs in memory, so nothing is written to — or read back from — disk.
data = []
for inc in increments:
    pathfile = f"THERM_ZENER_path_{inc}.json"
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
