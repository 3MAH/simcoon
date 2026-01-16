"""
Micromechanics examples
~~~~~~~~~~~~~~~~~~~~~~~~~~~

This tutorial studies the evolution of the mechanical properties of a
composite material considering spherical reinforcement.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from simcoon import simmit as sim
from simcoon import parameter as par
import os

###################################################################################
# In this tutorial we will study the evolution of the mechanical properties of a 2-phase composite material considering spherical reinforcement
# Note that two micromechanical schemes will be utilized : Morio-Tanaka and Self-consistent
#
# First we define the number of state variables (if any) at the macroscopic level
# and the material properties.

dir = os.path.dirname(os.path.realpath("__file__"))
nstatev = 0

nphases = 2  # The number of phases
num_file = 0  # The num of the file that contains the subphases
int1 = 50
int2 = 50
n_matrix = 0

props = np.array([nphases, num_file, int1, int2, n_matrix], dtype="float")

###################################################################################
# There is a possibility to consider a misorientation between the test frame
# of reference and the composite material orientation, using Euler angles
# with the z–x–z convention. The stiffness tensor will be returned in the
# test frame of reference.
# We also specify which micromechanical scheme will be used:
# - `MIMTN` → Mori-Tanaka scheme
# - `MISCN` → Self-consistent scheme

path_data = dir + "/data"
path_keys = dir + "/keys"
pathfile = "path.txt"

param_list = par.read_parameters()

psi_rve = 0.0
theta_rve = 0.0
phi_rve = 0.0

###################################################################################
# Run the simulation to obtain the effective isotropic elastic properties
# for different volume fractions of reinforcement up to 50%.

concentration = np.arange(0.0, 0.51, 0.01)

E_MT = np.zeros(len(concentration))
umat_name = "MIMTN"
for i, x in enumerate(concentration):
    param_list[1].value = x
    param_list[0].value = 1.0 - x

    par.copy_parameters(param_list, path_keys, path_data)
    par.apply_parameters(param_list, path_data)

    L = sim.L_eff(umat_name, props, nstatev, psi_rve, theta_rve, phi_rve)
    p = sim.L_iso_props(L).flatten()
    E_MT[i] = p[0]


E_SC = np.zeros(len(concentration))
umat_name = "MISCN"
for i, x in enumerate(concentration):
    param_list[1].value = x
    param_list[0].value = 1.0 - x

    par.copy_parameters(param_list, path_keys, path_data)
    par.apply_parameters(param_list, path_data)

    L = sim.L_eff(umat_name, props, nstatev, psi_rve, theta_rve, phi_rve)
    p = sim.L_iso_props(L).flatten()
    E_SC[i] = p[0]


###################################################################################
# Plotting the results
# --------------------------------------
# This is it, now we just need to plot the results.
# In blue, we plot the effective Young's modulus vs volume fraction curve using the Mori-Tanaka scheme,
# in red using the Self-Consistent scheme, and with black crosses the experimental data from (Wang 2003)

fig = plt.figure()
plt.xlabel(r"Volume fraction (c)", size=15)
plt.ylabel(r"Effective Young\'s modulus ($E$, MPa)", size=15)
plt.plot(concentration, E_MT, c="blue", label="Mori-Tanaka")
plt.plot(concentration, E_SC, c="red", label="Self-Consistent")
expfile = path_data + "/" + "E_exp.txt"
c, E = np.loadtxt(expfile, usecols=(0, 1), unpack=True)
plt.plot(c, E, linestyle="None", marker="x", color="black", markersize=10)
plt.show()

###############################################################################
# Effect of aspect ratio on effective stiffness
# ---------------------------------------------
# The shape of the reinforcement significantly affects the effective properties.
# Let's study how the aspect ratio of prolate (fiber-like) and oblate (disc-like)
# inclusions influences the effective Young's modulus using the Mori-Tanaka scheme.
#
# The aspect ratio is defined as :math:`ar = a_1/a_3`, where:
#
# - :math:`ar > 1`: prolate ellipsoid (fiber-like)
# - :math:`ar = 1`: sphere
# - :math:`ar < 1`: oblate ellipsoid (disc-like)

# Fixed volume fraction of 20%
c_reinf = 0.20
param_list[1].value = c_reinf
param_list[0].value = 1.0 - c_reinf
par.copy_parameters(param_list, path_keys, path_data)
par.apply_parameters(param_list, path_data)

# Aspect ratios from oblate (0.01) to prolate (100)
aspect_ratios = np.logspace(-2, 2, 80)

E_eff_ar = np.zeros(len(aspect_ratios))
nu_eff_ar = np.zeros(len(aspect_ratios))

umat_name = "MIMTN"

# Save original phase file
phase_file = path_data + "/Nellipsoids0.dat"
with open(phase_file, "r") as f:
    original_content = f.read()

for i, ar in enumerate(aspect_ratios):
    # Read the phase file
    with open(phase_file, "r") as f:
        lines = f.readlines()

    # Modify the reinforcement phase (line 2, which is the inclusion phase)
    # The columns are: Number Coatingof umat save c psi_mat theta_mat phi_mat a1 a2 a3 ...
    # We need to update a1, a2, a3 (columns 8, 9, 10 in 0-indexed split)
    parts = lines[2].split()

    # Update semi-axes: a1/a3 = ar, a2 = a3 = 1
    parts[8] = str(ar)  # a1
    parts[9] = "1"  # a2
    parts[10] = "1"  # a3

    lines[2] = "\t".join(parts) + "\n"

    with open(phase_file, "w") as f:
        f.writelines(lines)

    L = sim.L_eff(umat_name, props, nstatev, psi_rve, theta_rve, phi_rve)
    p = sim.L_iso_props(L).flatten()
    E_eff_ar[i] = p[0]
    nu_eff_ar[i] = p[1]

# Restore original phase file
with open(phase_file, "w") as f:
    f.write(original_content)

###############################################################################
# Plot: Effective Young's modulus vs aspect ratio

fig, ax = plt.subplots(figsize=(10, 6))
ax.semilogx(aspect_ratios, E_eff_ar / 1000, "b-", linewidth=2, label="Mori-Tanaka")
ax.axvline(x=1, color="k", linestyle=":", alpha=0.5, label="Sphere (ar=1)")
ax.axhline(
    y=E_MT[20] / 1000,
    color="r",
    linestyle="--",
    alpha=0.7,
    label=f"Sphere ref (E={E_MT[20] / 1000:.1f} GPa)",
)

# Add shaded regions for oblate/prolate
ymin, ymax = ax.get_ylim()
ax.axvspan(0.01, 1, alpha=0.1, color="orange")
ax.axvspan(1, 100, alpha=0.1, color="green")
ax.text(
    0.05, ymax - 0.1 * (ymax - ymin), "Oblate\n(disc-like)", fontsize=9, ha="center"
)
ax.text(
    20, ymax - 0.1 * (ymax - ymin), "Prolate\n(fiber-like)", fontsize=9, ha="center"
)

ax.set_xlabel("Aspect ratio $a_1/a_3$", fontsize=12)
ax.set_ylabel("Effective Young's modulus $E_{eff}$ (GPa)", fontsize=12)
ax.set_title(
    f"Effect of inclusion shape on composite stiffness (c={c_reinf * 100:.0f}%, Mori-Tanaka)",
    fontsize=14,
)
ax.legend(fontsize=10, loc="lower right")
ax.grid(True, alpha=0.3)
ax.set_xlim([0.01, 100])
plt.tight_layout()
plt.show()
