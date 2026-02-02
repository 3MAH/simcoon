"""
Comparison of objective rates
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This example shows the well-known spurious oscillations that can occur in the
Zaremba-Jaumann rate when simulating elastic responses under large transformations.
The results are compared with the Green-Naghdi and logarithmic Xiao-Meyers-Bruhns rates,
which do not exhibit such oscillations.

This example uses the new Python Solver API.
"""

import numpy as np
import matplotlib.pyplot as plt
from simcoon.solver import Solver, Block, StepMeca

plt.rcParams["figure.figsize"] = (18, 10)  # configure the figure output size

plt.rc("text", usetex=True)
plt.rc("font", family="serif")


###################################################################################
# We first consider a material with isotropic elastic behavior defined by its Young modulus
# and Poisson ratio. The material is subjected to a large simple shear deformation.
# Note that this example is only illustrative since for large deformations
# elastic materials are not physically meaningful.

E = 70000.0      # Young's modulus (MPa)
nu = 0.3         # Poisson ratio
alpha = 1.0e-5   # Thermal expansion coefficient

props = np.array([E, nu, alpha])
nstatev = 1

###############################################################################
# In this example we compare three objective rates:
# - Jaumann (corate_type='jaumann')
# - Green-Naghdi (corate_type='green_naghdi')
# - Logarithmic (corate_type='logarithmic')
#
# The simulation consists of a simple shear deformation where we apply
# shear strain e12 up to a large value of 5.0

rate_configs = [
    ('Jaumann', 'jaumann'),
    ('Green-Naghdi', 'green_naghdi'),
    ('Logarithmic', 'logarithmic'),
]

colors = ["blue", "red", "green"]

###############################################################################
# Create loading path: Simple shear test
# For simple shear, we apply strain in the 12 (shear) component
# Note: In Voigt notation, index 3 is e23, index 4 is e13, index 5 is e12

# Simple shear: apply shear strain (engineering strain = 2 * tensor strain)
# For large simple shear up to gamma = 5.0
max_shear = 5.0  # Maximum engineering shear strain
n_increments = 500

###############################################################################
# Run simulations for each objective rate

results = {}

for rate_name, corate in rate_configs:
    print(f"Running simulation with {rate_name} rate...")

    # Create step with shear loading
    # e12 in Voigt notation is index 5 (0-indexed: e11, e22, e33, e23, e13, e12)
    step = StepMeca(
        DEtot_end=np.array([0, 0, 0, 0, 0, max_shear]),  # Shear strain e12
        Dsigma_end=np.array([0, 0, 0, 0, 0, 0]),
        control=['stress', 'stress', 'stress', 'stress', 'stress', 'strain'],
        Dn_init=n_increments,
        Dn_mini=n_increments // 5,
        Dn_inc=n_increments * 2,
        time=max_shear  # Time matches shear value for convenience
    )

    block = Block(
        steps=[step],
        umat_name="ELISO",
        props=props,
        nstatev=nstatev,
        control_type='small_strain',
        corate_type=corate
    )

    solver = Solver(blocks=[block])
    history = solver.solve()

    # Extract results
    e11 = np.array([h.Etot[0] for h in history])
    e22 = np.array([h.Etot[1] for h in history])
    e12 = np.array([h.Etot[5] for h in history])  # Shear strain
    s11 = np.array([h.sigma[0] for h in history])
    s22 = np.array([h.sigma[1] for h in history])
    s12 = np.array([h.sigma[5] for h in history])  # Shear stress
    time = np.linspace(0, max_shear, len(history))

    # Extract rotation from deformation gradient (approximation)
    # For simple shear, rotation angle can be estimated
    # R11 component from the rotation tensor
    R11 = np.array([h.R[0, 0] for h in history])
    rotation_angle = np.arccos(np.minimum(R11, 1.0))

    results[rate_name] = {
        'e11': e11,
        'e22': e22,
        'e12': e12,
        's11': s11,
        's22': s22,
        's12': s12,
        'time': time,
        'rotation': rotation_angle,
    }

###############################################################################
# Plotting the results
# Compare strain and rotation evolution for all three rates

fig, axes = plt.subplots(2, 2, figsize=(18, 10))

plot_info = [
    (0, 0, 'e11', r"Strain ($\varepsilon_{11}$)"),
    (0, 1, 'e12', r"Strain ($\varepsilon_{12}$)"),
    (1, 0, 'e22', r"Strain ($\varepsilon_{22}$)"),
    (1, 1, 'rotation', r"Rotation angle (rad)"),
]

for i, (rate_name, corate) in enumerate(rate_configs):
    data = results[rate_name]
    values = [data['e11'], data['e12'], data['e22'], data['rotation']]

    for ax_idx, (row, col, key, ylabel) in enumerate(plot_info):
        axes[row, col].plot(data['time'], values[ax_idx], c=colors[i], label=rate_name)

for row, col, key, ylabel in plot_info:
    axes[row, col].set_xlabel(r"Time (s)", size=15)
    axes[row, col].set_ylabel(ylabel, size=15)
    axes[row, col].legend(loc=2)
    axes[row, col].grid(True)
    axes[row, col].tick_params(axis="both", which="major", labelsize=12)

plt.suptitle("Comparison of Objective Rates under Simple Shear", fontsize=16)
plt.tight_layout()
plt.show()

###############################################################################
# Additional plot: Stress response comparison

fig2, axes2 = plt.subplots(1, 3, figsize=(18, 5))

stress_plot_info = [
    (0, 's11', r"Stress $\sigma_{11}$ (MPa)"),
    (1, 's22', r"Stress $\sigma_{22}$ (MPa)"),
    (2, 's12', r"Stress $\sigma_{12}$ (MPa)"),
]

for i, (rate_name, corate) in enumerate(rate_configs):
    data = results[rate_name]
    for ax_idx, (idx, key, ylabel) in enumerate(stress_plot_info):
        axes2[idx].plot(data['time'], data[key], c=colors[i], label=rate_name)

for idx, key, ylabel in stress_plot_info:
    axes2[idx].set_xlabel(r"Time (s)", size=15)
    axes2[idx].set_ylabel(ylabel, size=15)
    axes2[idx].legend(loc='best')
    axes2[idx].grid(True)
    axes2[idx].tick_params(axis="both", which="major", labelsize=12)

plt.suptitle("Stress Response Comparison", fontsize=16)
plt.tight_layout()
plt.show()

###############################################################################
# Note that the Jaumann rate exhibits spurious oscillations in the stress and strain response,
# while the Green-Naghdi and Logarithmic rates provide smooth responses.
# This is a well-known issue with the Jaumann rate when dealing with large simple shear transformation.
# The Green-Naghdi and Logarithmic rates do not suffer from this problem, making them more suitable
# for simulations involving large deformations and rotations.
#
# While logarithmic rates are often considered the most accurate for large deformations,
# please note that the induced rotation is however not correct. Only the Green-Naghdi rate provides the exact rotation
# for rigid body motions corresponding to the RU (or VR) decomposition.

print("\n" + "=" * 70)
print("Summary of Objective Rate Comparison")
print("=" * 70)
for rate_name in results:
    data = results[rate_name]
    print(f"\n{rate_name}:")
    print(f"  Final s12: {data['s12'][-1]:.2f} MPa")
    print(f"  Max s11:   {max(abs(data['s11'])):.2f} MPa")
    print(f"  Final rotation: {data['rotation'][-1]:.4f} rad")
