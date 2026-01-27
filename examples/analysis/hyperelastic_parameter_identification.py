"""
Parameter Identification for Mooney-Rivlin Hyperelastic Material
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This example demonstrates inverse parameter identification for hyperelastic
materials using experimental stress-strain data. We identify the Mooney-Rivlin
model parameters from Treloar's classical rubber elasticity experiments.

**Optimization**: scipy.optimize.differential_evolution (global optimizer)
**Cost function**: Mean Squared Error on stress (sklearn.metrics.mean_squared_error)
**Forward model**: simcoon hyperelastic stress functions
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import differential_evolution, fsolve
from sklearn.metrics import mean_squared_error
from simcoon import simmit as sim
import os

###################################################################################
# Methodology Overview
# --------------------
#
# **Inverse Problem Formulation**
#
# Parameter identification is an inverse problem: given experimental observations
# (stress-strain data), we seek the material parameters that minimize the
# discrepancy between model predictions and experiments.
#
# The optimization problem is:
#
# .. math::
#    \boldsymbol{\theta}^* = \arg\min_{\boldsymbol{\theta}} \mathcal{L}(\boldsymbol{\theta})
#
# where :math:`\boldsymbol{\theta}` are the material parameters and
# :math:`\mathcal{L}` is the cost function.
#
# **Cost Function: Stress-Based MSE**
#
# We use the Mean Squared Error (MSE) computed on the first Piola-Kirchhoff
# stress :math:`P_{11}`:
#
# .. math::
#    \text{MSE} = \frac{1}{N} \sum_{i=1}^{N} \left( P_{11}^{\text{model}}(\lambda_i; \boldsymbol{\theta}) - P_{11}^{\text{exp}}(\lambda_i) \right)^2
#
# This stress-based metric directly measures the model's ability to predict
# mechanical response, which is the primary quantity of interest in constitutive
# modeling.
#
# **Why Differential Evolution?**
#
# We use ``scipy.optimize.differential_evolution`` because:
#
# 1. **Global optimization**: Unlike gradient-based methods, it explores the
#    entire parameter space and avoids local minima
# 2. **Derivative-free**: No gradients required, robust for non-smooth objectives
# 3. **Bounded search**: Natural handling of physical parameter constraints
# 4. **Robustness**: Less sensitive to initial guess compared to local methods
#
# For hyperelastic models, the cost function landscape can be non-convex,
# making global optimization essential for reliable parameter identification.

###################################################################################
# Mooney-Rivlin Constitutive Model
# --------------------------------
#
# The Mooney-Rivlin model is a phenomenological hyperelastic model widely used
# for rubber-like materials. The strain energy function is:
#
# .. math::
#    W = C_{10}(\bar{I}_1 - 3) + C_{01}(\bar{I}_2 - 3) + U(J)
#
# where:
#
# - :math:`\bar{I}_1, \bar{I}_2` are the first and second isochoric invariants
#   of the left Cauchy-Green tensor :math:`\mathbf{b} = \mathbf{F}\mathbf{F}^T`
# - :math:`C_{10}, C_{01}` are the material parameters to identify
# - :math:`U(J) = \kappa(J \ln J - J + 1)` is the volumetric contribution
# - :math:`\kappa` is the bulk modulus (fixed, assuming near-incompressibility)
#
# The Cauchy stress is computed as:
#
# .. math::
#    \boldsymbol{\sigma} = \boldsymbol{\sigma}_{\text{iso}} + \boldsymbol{\sigma}_{\text{vol}}
#
# where simcoon provides ``sigma_iso_hyper_invariants`` and ``sigma_vol_hyper``
# to compute these contributions from the strain energy derivatives:
#
# .. math::
#    \frac{\partial W}{\partial \bar{I}_1} = C_{10}, \quad
#    \frac{\partial W}{\partial \bar{I}_2} = C_{01}, \quad
#    \frac{\partial U}{\partial J} = \kappa \ln J

###################################################################################
# Experimental Data: Treloar's Rubber Experiments
# -----------------------------------------------
#
# L.R.G. Treloar (1944) performed classical experiments on vulcanized rubber
# under three loading conditions:
#
# 1. **Uniaxial Tension (UT)**: :math:`\lambda_1 = \lambda, \lambda_2 = \lambda_3 = \lambda_t`
# 2. **Pure Shear (PS)**: :math:`\lambda_1 = \lambda, \lambda_2 = \lambda_t, \lambda_3 = 1`
# 3. **Equibiaxial Tension (ET)**: :math:`\lambda_1 = \lambda_2 = \lambda, \lambda_3 = \lambda_t`
#
# where :math:`\lambda_t` is the transverse stretch determined by equilibrium
# (zero transverse stress). The data provides stretch :math:`\lambda` and the
# corresponding first Piola-Kirchhoff stress :math:`P_{11}` [MPa].


def load_treloar_data(filepath: str) -> pd.DataFrame:
    """
    Load Treloar experimental data.

    The file contains columns for three loading cases:
    - lambda_1, P1_MPa: Uniaxial tension
    - lambda_2, P2_MPa: Pure shear
    - lambda_3, P3_MPa: Equibiaxial tension

    Parameters
    ----------
    filepath : str
        Path to Treloar.txt data file

    Returns
    -------
    pd.DataFrame
        Experimental data with columns for each loading case
    """
    df = pd.read_csv(
        filepath,
        sep=r"\s+",
        engine="python",
        names=["lambda_1", "P1_MPa", "lambda_2", "P2_MPa", "lambda_3", "P3_MPa"],
        header=0,
    )
    return df


###################################################################################
# Forward Model: Stress Prediction using simcoon
# ----------------------------------------------
#
# The forward model computes the first Piola-Kirchhoff stress for given
# material parameters and stretch values. This is the core function that
# will be called repeatedly during optimization.
#
# **Key steps:**
#
# 1. For each stretch :math:`\lambda_1`, solve for the transverse stretch
#    :math:`\lambda_t` that satisfies equilibrium (zero transverse stress)
# 2. Construct the deformation gradient :math:`\mathbf{F}`
# 3. Compute the left Cauchy-Green tensor :math:`\mathbf{b} = \mathbf{F}\mathbf{F}^T`
#    using ``sim.L_Cauchy_Green(F)``
# 4. Compute strain energy derivatives for Mooney-Rivlin
# 5. Compute Cauchy stress using ``sim.sigma_iso_hyper_invariants`` and
#    ``sim.sigma_vol_hyper``
# 6. Convert to first Piola-Kirchhoff stress: :math:`P_{11} = \sigma_{11}/\lambda_1`


def mooney_rivlin_pk1_stress(
    C10: float,
    C01: float,
    kappa: float,
    lambda_array: np.ndarray,
    loading_case: str = "UT",
) -> np.ndarray:
    """
    Compute first Piola-Kirchhoff stress for Mooney-Rivlin material.

    Uses simcoon's hyperelastic stress functions to compute the Cauchy stress
    from the strain energy derivatives, then converts to PK1 stress.

    Parameters
    ----------
    C10 : float
        First Mooney-Rivlin parameter [MPa]
    C01 : float
        Second Mooney-Rivlin parameter [MPa]
    kappa : float
        Bulk modulus for volumetric response [MPa]
    lambda_array : np.ndarray
        Array of principal stretch values in the loading direction
    loading_case : str
        Type of loading: "UT" (uniaxial tension), "PS" (pure shear),
        or "ET" (equibiaxial tension)

    Returns
    -------
    np.ndarray
        First Piola-Kirchhoff stress P_11 [MPa]
    """
    PK1_stress = []
    lambda_t_guess = 1.0  # Initial guess for transverse stretch

    for lam in lambda_array:
        # -----------------------------------------------------------------
        # Step 1: Find transverse stretch satisfying equilibrium
        # -----------------------------------------------------------------
        def equilibrium_residual(lambda_t_vec):
            """Residual function for finding equilibrium transverse stretch."""
            lt = lambda_t_vec[0]

            # Construct deformation gradient based on loading case
            if loading_case == "UT":
                # Uniaxial: F = diag(lambda, lambda_t, lambda_t)
                F = np.diag([lam, lt, lt])
                J = lam * lt**2
            elif loading_case == "PS":
                # Pure shear: F = diag(lambda, lambda_t, 1)
                F = np.diag([lam, lt, 1.0])
                J = lam * lt
            elif loading_case == "ET":
                # Equibiaxial: F = diag(lambda, lambda, lambda_t)
                F = np.diag([lam, lam, lt])
                J = lam**2 * lt
            else:
                raise ValueError(f"Unknown loading case: {loading_case}")

            # Compute left Cauchy-Green tensor using simcoon
            b = sim.L_Cauchy_Green(F)

            # Mooney-Rivlin strain energy derivatives
            dW_dI1_bar = C10
            dW_dI2_bar = C01
            dU_dJ = kappa * np.log(J)

            # Compute Cauchy stress using simcoon functions
            sigma_iso = sim.sigma_iso_hyper_invariants(
                float(dW_dI1_bar), float(dW_dI2_bar), b, J, False
            )
            sigma_vol = sim.sigma_vol_hyper(dU_dJ, b, J, False)
            sigma = sigma_iso + sigma_vol

            # Return stress component that should be zero at equilibrium
            if loading_case == "UT":
                # sigma_22 = sigma_33 = 0
                return 0.5 * (sigma[1, 1] + sigma[2, 2])
            elif loading_case == "PS":
                # sigma_22 = 0
                return sigma[1, 1]
            else:  # ET
                # sigma_33 = 0
                return sigma[2, 2]

        # Solve for equilibrium transverse stretch
        lambda_t = fsolve(equilibrium_residual, lambda_t_guess, full_output=False)[0]
        lambda_t_guess = lambda_t  # Use as starting point for next iteration

        # -----------------------------------------------------------------
        # Step 2: Compute final stress state at equilibrium
        # -----------------------------------------------------------------
        if loading_case == "UT":
            F = np.diag([lam, lambda_t, lambda_t])
            J = lam * lambda_t**2
        elif loading_case == "PS":
            F = np.diag([lam, lambda_t, 1.0])
            J = lam * lambda_t
        else:  # ET
            F = np.diag([lam, lam, lambda_t])
            J = lam**2 * lambda_t

        b = sim.L_Cauchy_Green(F)

        # Strain energy derivatives
        dW_dI1_bar = C10
        dW_dI2_bar = C01
        dU_dJ = kappa * np.log(J)

        # Compute Cauchy stress
        sigma_iso = sim.sigma_iso_hyper_invariants(
            float(dW_dI1_bar), float(dW_dI2_bar), b, J, False
        )
        sigma_vol = sim.sigma_vol_hyper(dU_dJ, b, J, False)
        sigma = sigma_iso + sigma_vol

        # Convert Cauchy stress to PK1: P = J * sigma * F^{-T}
        # For diagonal F: P_11 = sigma_11 / lambda_1
        PK1_stress.append(sigma[0, 0] / lam)

    return np.array(PK1_stress)


###################################################################################
# Cost Function: Stress-Based Mean Squared Error
# ----------------------------------------------
#
# The cost function quantifies the discrepancy between model predictions and
# experimental data. We use MSE computed on the PK1 stress values.
#
# **Normalized MSE (NMSE)**:
#
# When combining multiple loading cases with different stress magnitudes,
# we can normalize by the variance to ensure balanced weighting:
#
# .. math::
#    \text{NMSE} = \frac{\text{MSE}}{\text{Var}(P^{\text{exp}})}
#
# This ensures that loading cases with smaller stress magnitudes are not
# dominated by cases with larger stresses during combined fitting.


def cost_function_mse(
    params: np.ndarray,
    lambda_exp: np.ndarray,
    P_exp: np.ndarray,
    kappa: float,
    loading_case: str,
    normalize: bool = False,
) -> float:
    """
    Compute MSE cost between model prediction and experimental stress data.

    Parameters
    ----------
    params : np.ndarray
        Material parameters [C10, C01]
    lambda_exp : np.ndarray
        Experimental stretch values
    P_exp : np.ndarray
        Experimental PK1 stress values [MPa]
    kappa : float
        Fixed bulk modulus [MPa]
    loading_case : str
        Loading type ("UT", "PS", or "ET")
    normalize : bool
        If True, divide MSE by variance of experimental data

    Returns
    -------
    float
        MSE or NMSE value
    """
    C10, C01 = params

    try:
        # Predict stress using forward model
        P_model = mooney_rivlin_pk1_stress(C10, C01, kappa, lambda_exp, loading_case)

        # Compute MSE using sklearn
        mse = mean_squared_error(P_exp, P_model)

        if normalize:
            variance = np.var(P_exp)
            return mse / variance if variance > 1e-12 else mse
        return mse

    except Exception:
        # Return large cost if computation fails (e.g., numerical issues)
        return 1e10


def cost_function_combined(
    params: np.ndarray, data_dict: dict, kappa: float, normalize: bool = True
) -> float:
    """
    Combined cost function for simultaneous fitting to multiple loading cases.

    Parameters
    ----------
    params : np.ndarray
        Material parameters [C10, C01]
    data_dict : dict
        Dictionary mapping loading case names to (lambda, P) tuples
    kappa : float
        Fixed bulk modulus [MPa]
    normalize : bool
        If True, use NMSE for balanced weighting across cases

    Returns
    -------
    float
        Sum of (N)MSE values over all loading cases
    """
    total_cost = 0.0
    for case_name, (lambda_exp, P_exp) in data_dict.items():
        cost = cost_function_mse(params, lambda_exp, P_exp, kappa, case_name, normalize)
        total_cost += cost
    return total_cost


###################################################################################
# Parameter Identification using Differential Evolution
# -----------------------------------------------------
#
# ``scipy.optimize.differential_evolution`` is a stochastic global optimizer
# based on evolutionary algorithms. Key parameters:
#
# - **bounds**: Search space for each parameter
# - **strategy**: Mutation strategy ('best1bin' works well for most problems)
# - **popsize**: Population size (larger = more thorough search)
# - **mutation**: Mutation constant (controls exploration vs exploitation)
# - **recombination**: Crossover probability
# - **polish**: If True, refines result with L-BFGS-B (local optimizer)


def identify_parameters(
    lambda_exp: np.ndarray,
    P_exp: np.ndarray,
    kappa: float = 4000.0,
    loading_case: str = "UT",
    bounds: list = None,
    normalize: bool = False,
    seed: int = 42,
    verbose: bool = True,
) -> dict:
    """
    Identify Mooney-Rivlin parameters using differential evolution.

    Parameters
    ----------
    lambda_exp : np.ndarray
        Experimental stretch values
    P_exp : np.ndarray
        Experimental PK1 stress values [MPa]
    kappa : float
        Fixed bulk modulus [MPa] (default: 4000 for near-incompressibility)
    loading_case : str
        Loading type ("UT", "PS", or "ET")
    bounds : list
        Parameter bounds [(C10_min, C10_max), (C01_min, C01_max)]
    normalize : bool
        If True, use normalized MSE
    seed : int
        Random seed for reproducibility
    verbose : bool
        Print optimization progress

    Returns
    -------
    dict
        Optimized parameters and optimization results
    """
    if bounds is None:
        bounds = [(0.01, 2.0), (-1.0, 1.0)]

    if verbose:
        print(f"\n{'=' * 60}")
        print(f"Optimizing for {loading_case} loading case")
        print(f"{'=' * 60}")
        print(f"Parameter bounds: C10 in {bounds[0]}, C01 in {bounds[1]}")
        print(f"Bulk modulus (fixed): kappa = {kappa} MPa")

    result = differential_evolution(
        cost_function_mse,
        bounds=bounds,
        args=(lambda_exp, P_exp, kappa, loading_case, normalize),
        strategy="best1bin",
        maxiter=500,
        popsize=15,
        tol=1e-8,
        mutation=(0.5, 1.0),
        recombination=0.7,
        seed=seed,
        polish=True,  # Refine with L-BFGS-B
        disp=verbose,
    )

    C10_opt, C01_opt = result.x

    if verbose:
        print(f"\nOptimization completed!")
        print(f"  C10 = {C10_opt:.6f} MPa")
        print(f"  C01 = {C01_opt:.6f} MPa")
        print(f"  Final MSE = {result.fun:.6e} MPa^2")

    return {
        "C10": C10_opt,
        "C01": C01_opt,
        "kappa": kappa,
        "mse": result.fun,
        "success": result.success,
        "result": result,
    }


def identify_parameters_combined(
    data_dict: dict,
    kappa: float = 4000.0,
    bounds: list = None,
    normalize: bool = True,
    seed: int = 42,
    verbose: bool = True,
) -> dict:
    """
    Identify parameters by fitting to multiple loading cases simultaneously.

    Using NMSE (normalized MSE) ensures balanced contribution from each
    loading case, preventing cases with larger stress magnitudes from
    dominating the optimization.

    Parameters
    ----------
    data_dict : dict
        Dictionary {case_name: (lambda_exp, P_exp)} for each loading case
    kappa : float
        Fixed bulk modulus [MPa]
    bounds : list
        Parameter bounds [(C10_min, C10_max), (C01_min, C01_max)]
    normalize : bool
        If True, use NMSE (recommended for combined fitting)
    seed : int
        Random seed
    verbose : bool
        Print progress

    Returns
    -------
    dict
        Optimized parameters and results
    """
    if bounds is None:
        bounds = [(0.01, 2.0), (-1.0, 1.0)]

    if verbose:
        print(f"\n{'=' * 60}")
        print(f"Combined optimization for: {list(data_dict.keys())}")
        print(f"{'=' * 60}")
        print(f"Parameter bounds: C10 in {bounds[0]}, C01 in {bounds[1]}")
        print(f"Bulk modulus (fixed): kappa = {kappa} MPa")
        print(f"Using {'NMSE' if normalize else 'MSE'} for cost function")

    result = differential_evolution(
        cost_function_combined,
        bounds=bounds,
        args=(data_dict, kappa, normalize),
        strategy="best1bin",
        maxiter=500,
        popsize=15,
        tol=1e-8,
        mutation=(0.5, 1.0),
        recombination=0.7,
        seed=seed,
        polish=True,
        disp=verbose,
    )

    C10_opt, C01_opt = result.x

    if verbose:
        print(f"\nOptimization completed!")
        print(f"  C10 = {C10_opt:.6f} MPa")
        print(f"  C01 = {C01_opt:.6f} MPa")
        print(f"  Final combined cost = {result.fun:.6e}")

    return {
        "C10": C10_opt,
        "C01": C01_opt,
        "kappa": kappa,
        "cost": result.fun,
        "success": result.success,
        "result": result,
    }


###################################################################################
# Visualization Functions
# -----------------------


def plot_fit_comparison(
    lambda_exp: np.ndarray,
    P_exp: np.ndarray,
    params: dict,
    loading_case: str,
    title: str = None,
    ax=None,
):
    """
    Plot experimental data vs model prediction.

    Parameters
    ----------
    lambda_exp : np.ndarray
        Experimental stretch values
    P_exp : np.ndarray
        Experimental PK1 stress values [MPa]
    params : dict
        Dictionary with 'C10', 'C01', 'kappa' keys
    loading_case : str
        Loading type
    title : str
        Plot title (optional)
    ax : matplotlib.axes.Axes
        Axes to plot on (creates new figure if None)

    Returns
    -------
    matplotlib.axes.Axes
        The axes with the plot
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))

    # Model prediction on fine grid
    lambda_model = np.linspace(1.0, lambda_exp.max() * 1.02, 200)
    P_model = mooney_rivlin_pk1_stress(
        params["C10"], params["C01"], params["kappa"], lambda_model, loading_case
    )

    # Experimental data
    ax.plot(
        lambda_exp,
        P_exp,
        "o",
        markersize=8,
        markerfacecolor="red",
        markeredgecolor="black",
        label="Treloar experimental",
    )

    # Model prediction
    ax.plot(
        lambda_model,
        P_model,
        "-",
        linewidth=2,
        color="blue",
        label=f"Mooney-Rivlin (C10={params['C10']:.4f}, C01={params['C01']:.4f})",
    )

    ax.set_xlabel(r"Stretch $\lambda$ [-]", fontsize=12)
    ax.set_ylabel(r"PK1 Stress $P_{11}$ [MPa]", fontsize=12)
    ax.set_title(title or f"{loading_case} - Parameter Identification", fontsize=14)
    ax.legend(loc="upper left", fontsize=10)
    ax.grid(True, alpha=0.3)

    return ax


###################################################################################
# Main Execution
# --------------
#
# We demonstrate two identification strategies:
#
# 1. **Individual fitting**: Optimize parameters separately for each loading case
# 2. **Combined fitting**: Find a single parameter set that best fits all cases
#
# The individual fits provide loading-case-specific parameters (useful for
# understanding material behavior), while combined fitting gives a universal
# parameter set for general use.

if __name__ == "__main__":
    # Locate data file relative to script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    data_path = os.path.join(
        script_dir, "..", "hyperelasticity", "comparison", "Treloar.txt"
    )

    print("=" * 70)
    print(" MOONEY-RIVLIN PARAMETER IDENTIFICATION")
    print(" Using scipy.differential_evolution and sklearn.mean_squared_error")
    print("=" * 70)

    # -------------------------------------------------------------------------
    # Load and prepare experimental data
    # -------------------------------------------------------------------------
    df = load_treloar_data(data_path)
    print(f"\nLoaded Treloar data: {len(df)} data points")

    # Extract data for each loading case (remove NaN values)
    ut_mask = ~df["lambda_1"].isna() & ~df["P1_MPa"].isna()
    ps_mask = ~df["lambda_2"].isna() & ~df["P2_MPa"].isna()
    et_mask = ~df["lambda_3"].isna() & ~df["P3_MPa"].isna()

    lambda_ut = df.loc[ut_mask, "lambda_1"].values
    P_ut = df.loc[ut_mask, "P1_MPa"].values

    lambda_ps = df.loc[ps_mask, "lambda_2"].values
    P_ps = df.loc[ps_mask, "P2_MPa"].values

    lambda_et = df.loc[et_mask, "lambda_3"].values
    P_et = df.loc[et_mask, "P3_MPa"].values

    print(f"\nData summary:")
    print(
        f"  Uniaxial Tension (UT): {len(lambda_ut)} points, "
        f"lambda in [{lambda_ut.min():.2f}, {lambda_ut.max():.2f}]"
    )
    print(
        f"  Pure Shear (PS):       {len(lambda_ps)} points, "
        f"lambda in [{lambda_ps.min():.2f}, {lambda_ps.max():.2f}]"
    )
    print(
        f"  Equibiaxial (ET):      {len(lambda_et)} points, "
        f"lambda in [{lambda_et.min():.2f}, {lambda_et.max():.2f}]"
    )

    # Fixed bulk modulus (near-incompressible rubber)
    kappa = 4000.0

    # -------------------------------------------------------------------------
    # Individual fitting for each loading case
    # -------------------------------------------------------------------------
    print("\n" + "=" * 70)
    print(" INDIVIDUAL FITTING (one parameter set per loading case)")
    print("=" * 70)

    params_ut = identify_parameters(lambda_ut, P_ut, kappa, "UT", verbose=True)
    params_ps = identify_parameters(lambda_ps, P_ps, kappa, "PS", verbose=True)
    params_et = identify_parameters(lambda_et, P_et, kappa, "ET", verbose=True)

    # Reference values from literature (Steinmann et al., 2012)
    print("\n" + "-" * 70)
    print("Comparison with literature (Steinmann et al., 2012):")
    print("-" * 70)
    print(
        f"{'Case':<8} {'C10 (lit.)':<12} {'C10 (ident.)':<14} {'C01 (lit.)':<12} {'C01 (ident.)':<14}"
    )
    print("-" * 70)
    print(
        f"{'UT':<8} {0.2588:<12.4f} {params_ut['C10']:<14.4f} {-0.0449:<12.4f} {params_ut['C01']:<14.4f}"
    )
    print(
        f"{'PS':<8} {0.2348:<12.4f} {params_ps['C10']:<14.4f} {-0.065:<12.4f} {params_ps['C01']:<14.4f}"
    )
    print(
        f"{'ET':<8} {0.1713:<12.4f} {params_et['C10']:<14.4f} {0.0047:<12.4f} {params_et['C01']:<14.4f}"
    )

    # -------------------------------------------------------------------------
    # Combined fitting
    # -------------------------------------------------------------------------
    print("\n" + "=" * 70)
    print(" COMBINED FITTING (single parameter set for all loading cases)")
    print("=" * 70)

    data_combined = {
        "UT": (lambda_ut, P_ut),
        "PS": (lambda_ps, P_ps),
        "ET": (lambda_et, P_et),
    }

    params_combined = identify_parameters_combined(
        data_combined, kappa=kappa, normalize=True, verbose=True
    )

    # -------------------------------------------------------------------------
    # Visualization
    # -------------------------------------------------------------------------
    print("\n" + "=" * 70)
    print(" GENERATING PLOTS")
    print("=" * 70)

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))

    # Row 1: Individual fits
    plot_fit_comparison(
        lambda_ut, P_ut, params_ut, "UT", title="UT - Individual fit", ax=axes[0, 0]
    )
    plot_fit_comparison(
        lambda_ps, P_ps, params_ps, "PS", title="PS - Individual fit", ax=axes[0, 1]
    )
    plot_fit_comparison(
        lambda_et, P_et, params_et, "ET", title="ET - Individual fit", ax=axes[0, 2]
    )

    # Row 2: Combined fit applied to all cases
    plot_fit_comparison(
        lambda_ut, P_ut, params_combined, "UT", title="UT - Combined fit", ax=axes[1, 0]
    )
    plot_fit_comparison(
        lambda_ps, P_ps, params_combined, "PS", title="PS - Combined fit", ax=axes[1, 1]
    )
    plot_fit_comparison(
        lambda_et, P_et, params_combined, "ET", title="ET - Combined fit", ax=axes[1, 2]
    )

    fig.suptitle(
        "Mooney-Rivlin Parameter Identification using Treloar Data\n"
        "(Top: Individual fits, Bottom: Combined fit)",
        fontsize=14,
        fontweight="bold",
    )
    plt.tight_layout()

    # -------------------------------------------------------------------------
    # Summary
    # -------------------------------------------------------------------------
    print("\n" + "=" * 70)
    print(" SUMMARY OF IDENTIFIED PARAMETERS")
    print("=" * 70)
    print(f"{'Case':<12} {'C10 [MPa]':>12} {'C01 [MPa]':>12} {'MSE [MPa^2]':>14}")
    print("-" * 52)
    print(
        f"{'UT':<12} {params_ut['C10']:>12.4f} {params_ut['C01']:>12.4f} {params_ut['mse']:>14.2e}"
    )
    print(
        f"{'PS':<12} {params_ps['C10']:>12.4f} {params_ps['C01']:>12.4f} {params_ps['mse']:>14.2e}"
    )
    print(
        f"{'ET':<12} {params_et['C10']:>12.4f} {params_et['C01']:>12.4f} {params_et['mse']:>14.2e}"
    )
    print("-" * 52)
    print(
        f"{'Combined':<12} {params_combined['C10']:>12.4f} {params_combined['C01']:>12.4f} {params_combined['cost']:>14.2e}"
    )
    print("=" * 70)

    print("\nNote: The combined fit uses NMSE (normalized MSE) which explains")
    print("the different cost magnitude compared to individual MSE values.")

    plt.show()
