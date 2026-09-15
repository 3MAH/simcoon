"""
Orthotropic elasticity (thermomechanical)
=========================================
"""

import os
import numpy as np
import simcoon as sim
import matplotlib.pyplot as plt

plt.rcParams["figure.figsize"] = (18, 10)

###################################################################################
# In thermoelastic orthotropic materials, there are three mutually perpendicular
# planes of symmetry. The following parameters are required:
#
# 1. The density :math:`\rho`
# 2. The specific heat :math:`c_p`
# 3. The Young modulus in direction 1: :math:`E_1`
# 4. The Young modulus in direction 2: :math:`E_2`
# 5. The Young modulus in direction 3: :math:`E_3`
# 6. The Poisson ratio :math:`\nu_{12}`
# 7. The Poisson ratio :math:`\nu_{13}`
# 8. The Poisson ratio :math:`\nu_{23}`
# 9. The shear modulus :math:`G_{12}`
# 10. The shear modulus :math:`G_{13}`
# 11. The shear modulus :math:`G_{23}`
# 12. The coefficient of thermal expansion :math:`\alpha_1`
# 13. The coefficient of thermal expansion :math:`\alpha_2`
# 14. The coefficient of thermal expansion :math:`\alpha_3`
#
# The elastic stiffness tensor for an orthotropic material is written in the
# Voigt notation as:
#
# .. math::
#
#    \mathbf{L} = \begin{pmatrix}
#        L_{11} & L_{12} & L_{13} & 0 & 0 & 0 \\
#        L_{12} & L_{22} & L_{23} & 0 & 0 & 0 \\
#        L_{13} & L_{23} & L_{33} & 0 & 0 & 0 \\
#        0 & 0 & 0 & G_{12} & 0 & 0 \\
#        0 & 0 & 0 & 0 & G_{13} & 0 \\
#        0 & 0 & 0 & 0 & 0 & G_{23}
#    \end{pmatrix}
#
# The thermal expansion tensor is:
#
# .. math::
#
#    \boldsymbol{\alpha} = \begin{pmatrix}
#        \alpha_1 & 0 & 0 \\
#        0 & \alpha_2 & 0 \\
#        0 & 0 & \alpha_3
#    \end{pmatrix}

umat_name = "ELORT"  # 5 character code for orthotropic elastic subroutine
nstatev = 1  # Number of internal variables

# Material parameters
rho = 4.4  # Density
c_p = 0.656  # Specific heat capacity
E_1 = 4500.0  # Young's modulus in direction 1 (MPa)
E_2 = 2300.0  # Young's modulus in direction 2 (MPa)
E_3 = 2700.0  # Young's modulus in direction 3 (MPa)
nu_12 = 0.06  # Poisson ratio 12
nu_13 = 0.08  # Poisson ratio 13
nu_23 = 0.3  # Poisson ratio 23
G_12 = 2200.0  # Shear modulus 12 (MPa)
G_13 = 2100.0  # Shear modulus 13 (MPa)
G_23 = 2400.0  # Shear modulus 23 (MPa)
alpha_1 = 1.0e-5  # Thermal expansion in direction 1
alpha_2 = 2.5e-5  # Thermal expansion in direction 2
alpha_3 = 2.2e-5  # Thermal expansion in direction 3

psi_rve = 0.0
theta_rve = 0.0
phi_rve = 0.0
solver_type = 0
corate_type = 2

props = np.array(
    [rho, c_p, E_1, E_2, E_3, nu_12, nu_13, nu_23, G_12, G_13, G_23, alpha_1, alpha_2, alpha_3]
)

path_data = "../data"

###################################################################################
# Loading in direction 1
# ~~~~~~~~~~~~~~~~~~~~~~~~

pathfile = "THERM_ELISO_path_1.json"

blocks, T_init, _ = sim.solver.load_path_json(os.path.join(path_data, pathfile))
res_1 = sim.solver.solve(
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
# Loading in direction 2
# ~~~~~~~~~~~~~~~~~~~~~~~~

pathfile = "THERM_ELISO_path_2.json"

blocks, T_init, _ = sim.solver.load_path_json(os.path.join(path_data, pathfile))
res_2 = sim.solver.solve(
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
# Loading in direction 3
# ~~~~~~~~~~~~~~~~~~~~~~~~

pathfile = "THERM_ELISO_path_3.json"

blocks, T_init, _ = sim.solver.load_path_json(os.path.join(path_data, pathfile))
res_3 = sim.solver.solve(
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
# Plotting the results -- Loading direction 1
# -----------------------------------------------

fig = plt.figure()

e11, e22, e33, e12, e13, e23 = res_1["Strain"]
s11, s22, s33, s12, s13, s23 = res_1["Stress"]
time, T, Q, r = res_1["Time"], res_1["Temp"], res_1["Q"], res_1["r"]
Wm, Wm_r, Wm_ir, Wm_d = res_1["Wm"]
Wt, Wt_r, Wt_ir = res_1["Wt"]

ax = fig.add_subplot(2, 2, 1)
plt.grid(True)
plt.tick_params(axis="both", which="major", labelsize=15)
plt.xlabel(r"Strain $\varepsilon_{11}$", size=15)
plt.ylabel(r"Stress $\sigma_{11}$ (MPa)", size=15)
plt.plot(e11, s11, c="black", label="direction 1")
plt.legend(loc="best")

ax = fig.add_subplot(2, 2, 2)
plt.grid(True)
plt.tick_params(axis="both", which="major", labelsize=15)
plt.xlabel("time (s)", size=15)
plt.ylabel(r"Temperature $\theta$ (K)", size=15)
plt.plot(time, T, c="black", label="temperature")
plt.legend(loc="best")

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
# Plotting the results -- Loading direction 2
# -----------------------------------------------

fig = plt.figure()

e11, e22, e33, e12, e13, e23 = res_2["Strain"]
s11, s22, s33, s12, s13, s23 = res_2["Stress"]
time, T, Q, r = res_2["Time"], res_2["Temp"], res_2["Q"], res_2["r"]
Wm, Wm_r, Wm_ir, Wm_d = res_2["Wm"]
Wt, Wt_r, Wt_ir = res_2["Wt"]

ax = fig.add_subplot(2, 2, 1)
plt.grid(True)
plt.tick_params(axis="both", which="major", labelsize=15)
plt.xlabel(r"Strain $\varepsilon_{22}$", size=15)
plt.ylabel(r"Stress $\sigma_{22}$ (MPa)", size=15)
plt.plot(e22, s22, c="black", label="direction 2")
plt.legend(loc="best")

ax = fig.add_subplot(2, 2, 2)
plt.grid(True)
plt.tick_params(axis="both", which="major", labelsize=15)
plt.xlabel("time (s)", size=15)
plt.ylabel(r"Temperature $\theta$ (K)", size=15)
plt.plot(time, T, c="black", label="temperature")
plt.legend(loc="best")

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
# Plotting the results -- Loading direction 3
# -----------------------------------------------

fig = plt.figure()

e11, e22, e33, e12, e13, e23 = res_3["Strain"]
s11, s22, s33, s12, s13, s23 = res_3["Stress"]
time, T, Q, r = res_3["Time"], res_3["Temp"], res_3["Q"], res_3["r"]
Wm, Wm_r, Wm_ir, Wm_d = res_3["Wm"]
Wt, Wt_r, Wt_ir = res_3["Wt"]

ax = fig.add_subplot(2, 2, 1)
plt.grid(True)
plt.tick_params(axis="both", which="major", labelsize=15)
plt.xlabel(r"Strain $\varepsilon_{33}$", size=15)
plt.ylabel(r"Stress $\sigma_{33}$ (MPa)", size=15)
plt.plot(e33, s33, c="black", label="direction 3")
plt.legend(loc="best")

ax = fig.add_subplot(2, 2, 2)
plt.grid(True)
plt.tick_params(axis="both", which="major", labelsize=15)
plt.xlabel("time (s)", size=15)
plt.ylabel(r"Temperature $\theta$ (K)", size=15)
plt.plot(time, T, c="black", label="temperature")
plt.legend(loc="best")

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
