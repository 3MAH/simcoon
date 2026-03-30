"""
Composite Parameter Identification using Mori-Tanaka and Self-Consistent
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This example demonstrates inverse identification of reinforcement elastic
properties in a two-phase composite material. Given experimental effective
Young's modulus data at various volume fractions (Wang 2003), we identify the
Young's modulus and Poisson's ratio of the reinforcement phase.

Two micromechanical schemes are compared:

- **Mori-Tanaka** (``MIMTN``): accurate at low-to-moderate volume fractions
- **Self-Consistent** (``MISCN``): better at high volume fractions (c > 0.3)

**Forward model**: simcoon mean-field homogenization (``L_eff``)
**Optimization**: ``scipy.optimize.differential_evolution`` (global optimizer)
**Key system**: simcoon ``Parameter`` for generic file templating
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import differential_evolution
import simcoon as sim
from simcoon import parameter as par
import os

###################################################################################
# Problem Setup
# -------------
#
# We have a glass-particle / epoxy composite. The matrix properties are known:
#
# - :math:`E_m = 2250` MPa, :math:`\\nu_m = 0.19`
#
# The reinforcement properties :math:`(E_f, \\nu_f)` are unknown and must be
# identified from experimental effective Young's modulus measurements at several
# volume fractions.
#
# **Key System**: simcoon's key system allows parameters to be injected into
# any input file via alphanumeric placeholders (e.g., ``@2p``). This makes the
# identification workflow generic and applicable to any simulation tool
# (simcoon solver, Mori-Tanaka, fedoo FE, etc.).


def load_experimental_data(filepath):
    """Load experimental E_eff vs volume fraction data."""
    data = np.loadtxt(filepath)
    return data[:, 0], data[:, 1]


###################################################################################
# Forward Model
# -------------
#
# For each set of candidate properties :math:`(E_f, \\nu_f)`, we compute the
# effective Young's modulus at all experimental volume fractions using simcoon's
# ``L_eff`` function.


def compute_E_eff(
    E_f, nu_f, concentrations, umat_name, param_list, path_keys, path_data,
    props_composite,
):
    """
    Compute effective Young's modulus at given volume fractions.

    Parameters
    ----------
    E_f : float
        Reinforcement Young's modulus [MPa]
    nu_f : float
        Reinforcement Poisson's ratio [-]
    concentrations : array-like
        Volume fractions of reinforcement
    umat_name : str
        Micromechanical scheme ("MIMTN" or "MISCN")
    param_list : list of Parameter
        Parameter objects with keys for file templating
    path_keys : str
        Path to template files
    path_data : str
        Path to working data directory (must be "data" relative to cwd)
    props_composite : numpy.ndarray
        Composite definition array [nphases, num_file, int1, int2, n_matrix]

    Returns
    -------
    numpy.ndarray
        Effective Young's modulus at each volume fraction [MPa]
    """
    E_eff = np.zeros(len(concentrations))

    for i, c in enumerate(concentrations):
        param_list[0].value = 1.0 - c
        param_list[1].value = c
        param_list[2].value = E_f
        param_list[3].value = nu_f

        par.copy_parameters(param_list, path_keys, path_data)
        par.apply_parameters(param_list, path_data)

        L = sim.L_eff(umat_name, props_composite, 0, 0.0, 0.0, 0.0)
        iso_props = sim.L_iso_props(L).flatten()
        E_eff[i] = iso_props[0]

    return E_eff


def cost_function(
    params_opt, c_exp, E_exp, umat_name, param_list, path_keys, path_data,
    props_composite,
):
    """MSE cost function for reinforcement property identification."""
    E_f, nu_f = params_opt
    try:
        E_pred = compute_E_eff(
            E_f, nu_f, c_exp, umat_name, param_list, path_keys, path_data,
            props_composite,
        )
        return np.mean((E_pred - E_exp) ** 2)
    except Exception:
        return 1e12


def identify_reinforcement(
    c_exp, E_exp, umat_name, param_list, path_keys, path_data,
    props_composite, bounds, verbose=True,
):
    """
    Identify reinforcement properties using differential evolution.

    Returns
    -------
    dict
        Identified parameters and optimization result
    """
    if verbose:
        print(f"\n{'=' * 60}")
        print(f"  Identification with {umat_name}")
        print(f"{'=' * 60}")

    result = differential_evolution(
        cost_function,
        bounds=bounds,
        args=(c_exp, E_exp, umat_name, param_list, path_keys, path_data,
              props_composite),
        strategy="best1bin",
        maxiter=200,
        popsize=15,
        tol=1e-8,
        mutation=(0.5, 1.0),
        recombination=0.7,
        seed=42,
        polish=True,
        disp=verbose,
    )

    E_f_opt, nu_f_opt = result.x
    if verbose:
        print(f"\n  E_f  = {E_f_opt:.0f} MPa")
        print(f"  nu_f = {nu_f_opt:.3f}")
        print(f"  MSE  = {result.fun:.2e} MPa^2")

    return {
        "E_f": E_f_opt,
        "nu_f": nu_f_opt,
        "mse": result.fun,
        "result": result,
        "umat_name": umat_name,
    }


###################################################################################
# Main Execution
# --------------
#
# We compare two identification strategies:
#
# 1. **Mori-Tanaka** (``MIMTN``): dilute approximation, best for c < 30%
# 2. **Self-Consistent** (``MISCN``): accounts for percolation, better at high c
#
# At high volume fractions (c = 0.5), Mori-Tanaka underestimates stiffness
# and compensates by overestimating :math:`E_f`. Self-Consistent gives more
# physically meaningful identified properties.

if __name__ == "__main__":

    # -----------------------------------------------------------------
    # Change to example directory (L_eff reads from data/ in cwd)
    # -----------------------------------------------------------------
    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)

    path_keys = "keys_ident"
    path_data = "data"

    # -----------------------------------------------------------------
    # Load experimental data (Wang 2003)
    # -----------------------------------------------------------------
    c_exp, E_exp = load_experimental_data("data/E_exp.txt")

    print("=" * 60)
    print(" COMPOSITE REINFORCEMENT IDENTIFICATION")
    print(" Glass/Epoxy — Mean-field + Differential Evolution")
    print("=" * 60)
    print(f"\nExperimental data (Wang 2003): {len(c_exp)} points")
    for c, E in zip(c_exp, E_exp):
        print(f"  c = {c:.1f}  ->  E_eff = {E:.0f} MPa")

    # -----------------------------------------------------------------
    # Load parameters (key system)
    # -----------------------------------------------------------------
    # Define parameters inline (no external file needed):
    #   @0p -> matrix concentration
    #   @1p -> reinforcement concentration
    #   @2p -> reinforcement E [MPa]
    #   @3p -> reinforcement nu
    param_list = [
        par.Parameter(0, bounds=(0.0, 1.0), key="@0p",
                      sim_input_files=["Nellipsoids0.dat"]),
        par.Parameter(1, bounds=(0.0, 1.0), key="@1p",
                      sim_input_files=["Nellipsoids0.dat"]),
        par.Parameter(2, bounds=(10000, 200000), key="@2p",
                      sim_input_files=["Nellipsoids0.dat"]),
        par.Parameter(3, bounds=(0.01, 0.45), key="@3p",
                      sim_input_files=["Nellipsoids0.dat"]),
    ]

    # Composite definition
    props_composite = np.array([2, 0, 50, 50, 0], dtype="float")

    # Bounds for identification (E_f, nu_f)
    bounds = [(10000, 200000), (0.01, 0.45)]

    # -----------------------------------------------------------------
    # Identification with Mori-Tanaka
    # -----------------------------------------------------------------
    result_MT = identify_reinforcement(
        c_exp, E_exp, "MIMTN", param_list, path_keys, path_data,
        props_composite, bounds,
    )

    # -----------------------------------------------------------------
    # Identification with Self-Consistent
    # -----------------------------------------------------------------
    result_SC = identify_reinforcement(
        c_exp, E_exp, "MISCN", param_list, path_keys, path_data,
        props_composite, bounds,
    )

    # -----------------------------------------------------------------
    # Summary
    # -----------------------------------------------------------------
    print(f"\n{'=' * 60}")
    print(f" SUMMARY")
    print(f"{'=' * 60}")
    print(f"  Reference glass: E_f ~ 73000 MPa, nu_f ~ 0.22")
    print(f"")
    print(f"  {'Scheme':<20} {'E_f [MPa]':>12} {'nu_f':>8} {'MSE':>14}")
    print(f"  {'-' * 56}")
    print(f"  {'Mori-Tanaka':<20} {result_MT['E_f']:>12.0f} "
          f"{result_MT['nu_f']:>8.3f} {result_MT['mse']:>14.2e}")
    print(f"  {'Self-Consistent':<20} {result_SC['E_f']:>12.0f} "
          f"{result_SC['nu_f']:>8.3f} {result_SC['mse']:>14.2e}")
    print()

    # Note on model limitations
    print("  Note: Mori-Tanaka underestimates stiffness at high volume")
    print("  fractions and compensates by overestimating E_f. Self-Consistent")
    print("  accounts for phase connectivity and gives more physically")
    print("  meaningful results for this dataset (c up to 50%).")

    # -----------------------------------------------------------------
    # Visualization
    # -----------------------------------------------------------------
    c_model = np.arange(0.0, 0.51, 0.01)

    E_MT = compute_E_eff(
        result_MT["E_f"], result_MT["nu_f"], c_model, "MIMTN",
        param_list, path_keys, path_data, props_composite,
    )
    E_SC = compute_E_eff(
        result_SC["E_f"], result_SC["nu_f"], c_model, "MISCN",
        param_list, path_keys, path_data, props_composite,
    )
    # Reference with handbook glass properties
    E_ref_MT = compute_E_eff(
        73000, 0.22, c_model, "MIMTN",
        param_list, path_keys, path_data, props_composite,
    )
    E_ref_SC = compute_E_eff(
        73000, 0.22, c_model, "MISCN",
        param_list, path_keys, path_data, props_composite,
    )

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # Left: Identified fits
    ax1.plot(c_exp, E_exp, "kx", markersize=10, markeredgewidth=2,
             label="Exp. (Wang 2003)")
    ax1.plot(c_model, E_MT, "b-", linewidth=2,
             label=f"MT identified (E_f={result_MT['E_f']:.0f})")
    ax1.plot(c_model, E_SC, "r-", linewidth=2,
             label=f"SC identified (E_f={result_SC['E_f']:.0f})")
    ax1.set_xlabel("Reinforcement volume fraction $c$", fontsize=12)
    ax1.set_ylabel("Effective Young's modulus $E_{eff}$ [MPa]", fontsize=12)
    ax1.set_title("Identified Reinforcement Properties", fontsize=13)
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)

    # Right: Reference vs Identified
    ax2.plot(c_exp, E_exp, "kx", markersize=10, markeredgewidth=2,
             label="Exp. (Wang 2003)")
    ax2.plot(c_model, E_ref_MT, "b--", linewidth=1.5, alpha=0.6,
             label="MT with E_f=73 GPa (handbook)")
    ax2.plot(c_model, E_ref_SC, "r--", linewidth=1.5, alpha=0.6,
             label="SC with E_f=73 GPa (handbook)")
    ax2.plot(c_model, E_MT, "b-", linewidth=2, alpha=0.8,
             label="MT identified")
    ax2.plot(c_model, E_SC, "r-", linewidth=2, alpha=0.8,
             label="SC identified")
    ax2.set_xlabel("Reinforcement volume fraction $c$", fontsize=12)
    ax2.set_ylabel("Effective Young's modulus $E_{eff}$ [MPa]", fontsize=12)
    ax2.set_title("Handbook vs Identified Properties", fontsize=13)
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.3)

    fig.suptitle(
        "Composite Reinforcement Identification — Glass/Epoxy\n"
        "Mori-Tanaka vs Self-Consistent + Differential Evolution",
        fontsize=14, fontweight="bold",
    )
    plt.tight_layout()
    plt.show()
