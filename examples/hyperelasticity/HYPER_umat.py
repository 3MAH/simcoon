"""
=========================================
Hyperelastic models - use the umat
=========================================

In this example, we compare five hyperelastic constitutive laws.
Each model is run for **uniaxial tension (UT)**,
**equibiaxial tension (ET)**, and **pure shear (PS)**.

We present one section per model.

This example uses the new Python Solver API.
"""

# sphinx_gallery_thumbnail_number = 1

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import NamedTuple, List, Tuple
from dataclasses import dataclass
from simcoon.solver import Solver, Block, StepMeca

###################################################################################
# Several hyperelastic isotropic materials are tested.
# They are compared to the well-known Treloar experimental data. The following
# models are tested:
#
# - The Neo-Hookean model
# - The Mooney-Rivlin model
# - The Isihara model
# - The Gent-Thomas model
# - The Swanson model
#
# Considering the Neo-Hookean model, the strain energy function is expressed as:
#
# .. math::
#    W = \frac{\mu}{2}\left(\bar{I}_1 - 3\right) + \kappa\left( J\ln J - J + 1\right)
#
# The parameters are:
# 1. The shear coefficient :math:`\mu`,
# 2. The bulk compressibility :math:`\kappa`,
#
# Considering the Mooney-Rivlin model, the strain energy function is expressed as:
#
# .. math::
#    W = C_{10} \left(\bar{I}_1 -3\right) + C_{01} \left(\bar{I}_2 -3\right) + \kappa \left( J \ln J - J + 1 \right)
#
# The parameters are:
# 1. The first governing parameter :math:`c_{10}`,
# 2. The second governing parameter :math:`c_{01}`,
# 3. The bulk compressibility :math:`\kappa`,


###############################################################################
# Data structures for material models and loading cases

@dataclass
class Umat:
    name: str
    parameters: List[float]
    colors: str


class loading_case(NamedTuple):
    name: str
    stretch_max: float
    comparison: List[Tuple[pd.Series, pd.Series]]


###############################################################################
# Reading experimental Treloar data

path_data = "comparison"
comparison_file_exp = "Treloar.txt"
comparison_exp = path_data + "/" + comparison_file_exp

df_exp = pd.read_csv(
    comparison_exp,
    sep=r"\s+",
    engine="python",
    names=["lambda_1", "P1_MPa", "lambda_2", "P2_MPa", "lambda_3", "P3_MPa"],
    header=0,
)

###############################################################################
# Defining the material models

Neo_Hookean_model = Umat(name="NEOHC", parameters=[0.5673, 1000.0], colors="blue")
Mooney_Rivlin_model = Umat(
    name="MOORI", parameters=[0.2588, -0.0449, 10000.0], colors="orange"
)
Isihara_model = Umat(
    name="ISHAH", parameters=[0.1161, 0.0136, 0.0114, 4000.0], colors="green"
)
Gent_Thomas_model = Umat(
    name="GETHH", parameters=[0.2837, 2.81e-11, 4000.0], colors="red"
)
Swanson_model = Umat(
    name="SWANH",
    parameters=[
        2,
        4000.0,
        2.83e-3,
        1.871e-13,
        1.684,
        -0.4302,
        2.82e-13,
        0.4643,
        9.141,
        0.7882,
    ],
    colors="purple",
)

list_umats = [
    Neo_Hookean_model,
    Mooney_Rivlin_model,
    Isihara_model,
    Gent_Thomas_model,
    Swanson_model,
]

###############################################################################
# Helper function to run a simulation using the new Solver API


def run_hyperelastic_simulation(umat_name, params, stretch_max, loading_type='UT'):
    """
    Run a hyperelastic simulation using the new Python Solver API.

    Parameters
    ----------
    umat_name : str
        Name of the UMAT (e.g., 'NEOHC', 'MOORI')
    params : list
        Material parameters
    stretch_max : float
        Maximum stretch ratio
    loading_type : str
        'UT' for uniaxial tension, 'PS' for pure shear, 'ET' for equibiaxial tension

    Returns
    -------
    tuple
        (stretch, stress) arrays
    """
    nstatev = 1
    props = np.array(params)

    # Calculate strain from stretch
    # For large deformations, use logarithmic strain: e = ln(lambda)
    strain_max = np.log(stretch_max)

    # Define loading path based on type
    # Using pure strain control for stability (approximate lateral strains for incompressibility)
    if loading_type == 'UT':
        # Uniaxial tension: for incompressible, lateral strain ~ -0.5 * axial strain
        lat_strain = -0.5 * strain_max  # Approximate for incompressible
        DEtot_end = np.array([strain_max, lat_strain, lat_strain, 0, 0, 0])
        control = ['strain', 'strain', 'strain', 'strain', 'strain', 'strain']
    elif loading_type == 'PS':
        # Pure shear: constrain direction 2, use incompressibility for direction 3
        lat_strain = -strain_max  # For pure shear, e33 ~ -e11 (incompressible)
        DEtot_end = np.array([strain_max, 0, lat_strain, 0, 0, 0])
        control = ['strain', 'strain', 'strain', 'strain', 'strain', 'strain']
    elif loading_type == 'ET':
        # Equibiaxial tension: for incompressible, e33 ~ -2*e11
        lat_strain = -2.0 * strain_max  # Approximate for incompressible
        DEtot_end = np.array([strain_max, strain_max, lat_strain, 0, 0, 0])
        control = ['strain', 'strain', 'strain', 'strain', 'strain', 'strain']
    else:
        raise ValueError(f"Unknown loading type: {loading_type}")

    step = StepMeca(
        DEtot_end=DEtot_end,
        Dsigma_end=np.array([0, 0, 0, 0, 0, 0]),
        control=control,
        Dn_init=200,
        Dn_mini=50,
        Dn_inc=400,
        time=1.0
    )

    block = Block(
        steps=[step],
        umat_name=umat_name,
        props=props,
        nstatev=nstatev,
        control_type='logarithmic',  # Use logarithmic strain for hyperelasticity
        corate_type='jaumann'
    )

    solver = Solver(blocks=[block])
    history = solver.solve()

    # Extract results
    # Stretch = exp(strain) for logarithmic strain
    e11 = np.array([h.Etot[0] for h in history])
    stretch = np.exp(e11)

    # For hyperelastic materials, convert Cauchy stress to 1st Piola-Kirchhoff
    # PK1 = J * sigma * F^{-T} -> for uniaxial: PK1_11 ≈ sigma_11 / lambda_1
    s11 = np.array([h.sigma[0] for h in history])

    # Approximate PK1 stress (nominal stress)
    # For incompressible/nearly incompressible: J ≈ 1
    PK1 = s11 / stretch

    return stretch, PK1


###############################################################################
# Plot hyperelastic models in uniaxial tension

model_colors = {
    "NEOHC": "blue",
    "MOORI": "orange",
    "ISHAH": "green",
    "GETHH": "red",
    "SWANH": "purple",
}

fig, axes = plt.subplots(2, 3, figsize=(12, 8))
axes = axes.flatten()

# Maximum stretch from experimental data
max_stretch_UT = df_exp["lambda_1"].max()

for i, umat in enumerate(list_umats):
    print(f"Running {umat.name} for uniaxial tension...")

    stretch, PK1 = run_hyperelastic_simulation(
        umat.name, umat.parameters, max_stretch_UT, 'UT'
    )

    # Plot model prediction
    axes[i].plot(
        stretch,
        PK1,
        color=model_colors[umat.name],
        label=f"{umat.name} prediction",
    )

    # Plot Treloar experimental data
    axes[i].plot(
        df_exp["lambda_1"],
        df_exp["P1_MPa"],
        linestyle="--",
        marker="o",
        color="black",
        lw=1,
        markerfacecolor="red",
        markeredgecolor="black",
        label="Treloar",
    )

    # Formatting
    axes[i].set_xlabel(r"$\lambda$")
    axes[i].set_ylabel("PK1 [MPa]")
    axes[i].set_title(umat.name)
    axes[i].grid(True)
    axes[i].legend()

fig.suptitle("Uniaxial tension", fontsize=14)
fig.delaxes(axes[5])
fig.tight_layout()
plt.show()

###############################################################################
# Plot hyperelastic models in pure shear
#
# Update parameters for pure shear fitting

Neo_Hookean_model.parameters = [0.3360, 1000.0]
Mooney_Rivlin_model.parameters = [0.2348, -0.065, 10000.0]
Isihara_model.parameters = [0.1601, 0.0037, 0.0031, 4000.0]
Gent_Thomas_model.parameters = [0.1629, 0.0376, 4000.0]
Swanson_model.parameters = [
    2,
    4000.0,
    0.0676,
    0.2861,
    0.2687,
    -0.4683,
    3.266e-11,
    0.0267,
    9.131,
    0.7157,
]

fig, axes = plt.subplots(2, 3, figsize=(12, 8))
axes = axes.flatten()

max_stretch_PS = df_exp["lambda_2"].max()

for i, umat in enumerate(list_umats):
    print(f"Running {umat.name} for pure shear...")

    stretch, PK1 = run_hyperelastic_simulation(
        umat.name, umat.parameters, max_stretch_PS, 'PS'
    )

    # Plot model prediction
    axes[i].plot(
        stretch,
        PK1,
        color=model_colors[umat.name],
        label=f"{umat.name} prediction",
    )

    # Plot Treloar experimental data
    axes[i].plot(
        df_exp["lambda_2"],
        df_exp["P2_MPa"],
        linestyle="--",
        marker="o",
        color="black",
        lw=1,
        markerfacecolor="red",
        markeredgecolor="black",
        label="Treloar",
    )

    # Formatting
    axes[i].set_xlabel(r"$\lambda$")
    axes[i].set_ylabel("PK1 [MPa]")
    axes[i].set_title(umat.name)
    axes[i].grid(True)
    axes[i].legend()

fig.suptitle("Pure shear", fontsize=14)
fig.delaxes(axes[5])
fig.tight_layout()
plt.show()

###############################################################################
# Plot hyperelastic models in equibiaxial tension
#
# Update parameters for equibiaxial tension fitting

Neo_Hookean_model.parameters = [0.4104, 1000.0]
Mooney_Rivlin_model.parameters = [0.1713, 0.0047, 10000.0]
Isihara_model.parameters = [0.1993, 0.0015, 0.0013, 4000.0]
Gent_Thomas_model.parameters = [0.2052, 2.22e-14, 4000.0]
Swanson_model.parameters = [
    2,
    4000.0,
    0.21,
    0.1036,
    -1833.0,
    -6.634,
    0.0074,
    0.266,
    1.429,
    -0.623,
]

fig, axes = plt.subplots(2, 3, figsize=(12, 8))
axes = axes.flatten()

max_stretch_ET = df_exp["lambda_3"].max()

for i, umat in enumerate(list_umats):
    print(f"Running {umat.name} for equibiaxial tension...")

    stretch, PK1 = run_hyperelastic_simulation(
        umat.name, umat.parameters, max_stretch_ET, 'ET'
    )

    # Plot model prediction
    axes[i].plot(
        stretch,
        PK1,
        color=model_colors[umat.name],
        label=f"{umat.name} prediction",
    )

    # Plot Treloar experimental data
    axes[i].plot(
        df_exp["lambda_3"],
        df_exp["P3_MPa"],
        linestyle="--",
        marker="o",
        color="black",
        lw=1,
        markerfacecolor="red",
        markeredgecolor="black",
        label="Treloar",
    )

    # Formatting
    axes[i].set_xlabel(r"$\lambda$")
    axes[i].set_ylabel("PK1 [MPa]")
    axes[i].set_title(umat.name)
    axes[i].grid(True)
    axes[i].legend()

fig.suptitle("Equibiaxial tension", fontsize=14)
fig.delaxes(axes[5])
fig.tight_layout()
plt.show()

print("\nHyperelastic model comparison complete.")
print("All simulations used the new Python Solver API.")
