"""
Material Parameter Identification Example
=========================================

This example demonstrates how to use the simcoon identification module
to calibrate material parameters from experimental data.

We identify the Young's modulus and yield stress of an elastic-plastic
material from a simulated uniaxial tension test.
"""

import numpy as np
import matplotlib.pyplot as plt

from simcoon.solver import Solver, Block, StepMeca
from simcoon.identification import (
    IdentificationProblem,
    levenberg_marquardt,
    differential_evolution,
    hybrid_optimization,
    compute_sensitivity,
    identifiability_check,
)


###############################################################################
# Generate synthetic experimental data
# ------------------------------------
# First, we generate "experimental" data by running a simulation with known
# parameters. In practice, you would load actual experimental data.

# True material parameters (what we want to identify)
E_true = 200000.0      # Young's modulus (MPa)
nu_true = 0.3          # Poisson ratio (fixed)
sigma_Y_true = 350.0   # Yield stress (MPa)
H_true = 1000.0        # Hardening modulus (MPa)
n_true = 0.4           # Hardening exponent

# Generate "experimental" data
print("Generating synthetic experimental data...")

def run_simulation(params):
    """Run a plasticity simulation with given parameters."""
    E, sigma_Y, H, n = params
    nu = nu_true  # Fixed
    alpha = 0.0   # No thermal expansion

    props = np.array([E, nu, alpha, sigma_Y, H, n])

    # Uniaxial tension to 5% strain
    step = StepMeca(
        DEtot_end=np.array([0.05, 0, 0, 0, 0, 0]),
        Dsigma_end=np.array([0, 0, 0, 0, 0, 0]),
        control=['strain', 'stress', 'stress', 'stress', 'stress', 'stress'],
        Dn_init=100,
        Dn_mini=20,
        Dn_inc=200,
        time=1.0
    )

    block = Block(
        steps=[step],
        umat_name="EPICP",
        props=props,
        nstatev=8,
        control_type='small_strain',
        corate_type='logarithmic'
    )

    solver = Solver(blocks=[block])
    history = solver.solve()

    # Extract strain-stress curve
    strain = np.array([h.Etot[0] for h in history])
    stress = np.array([h.sigma[0] for h in history])

    return {'strain': strain, 'stress': stress}

# Generate experimental data with true parameters
true_params = np.array([E_true, sigma_Y_true, H_true, n_true])
exp_results = run_simulation(true_params)
exp_strain = exp_results['strain']
exp_stress = exp_results['stress']

# Add some noise to simulate experimental uncertainty
np.random.seed(42)
noise_level = 5.0  # MPa
exp_stress_noisy = exp_stress + np.random.normal(0, noise_level, len(exp_stress))

print(f"Generated {len(exp_strain)} data points")
print(f"Strain range: {exp_strain[0]:.4f} to {exp_strain[-1]:.4f}")
print(f"Stress range: {exp_stress_noisy[0]:.1f} to {exp_stress_noisy[-1]:.1f} MPa")


###############################################################################
# Define the identification problem
# ---------------------------------
# We'll try to identify E and sigma_Y, keeping other parameters fixed.

def simulation_wrapper(params):
    """Wrapper that uses only the parameters being identified."""
    E, sigma_Y = params
    # Use fixed values for other parameters
    full_params = np.array([E, sigma_Y, H_true, n_true])
    return run_simulation(full_params)

# Define parameters to identify
parameters = [
    {'name': 'E', 'bounds': (150000, 250000), 'initial': 180000},
    {'name': 'sigma_Y', 'bounds': (200, 500), 'initial': 300},
]

# Create identification problem
problem = IdentificationProblem(
    parameters=parameters,
    simulate=simulation_wrapper,
    exp_data={'stress': exp_stress_noisy},
    cost_type='mse',
)

print(f"\nIdentification problem:")
print(f"  Parameters to identify: {problem.parameter_names}")
print(f"  Initial guess: {problem.get_initial()}")


###############################################################################
# Run identification with Levenberg-Marquardt
# -------------------------------------------

print("\n" + "=" * 60)
print("Running Levenberg-Marquardt optimization...")
print("=" * 60)

result_lm = levenberg_marquardt(problem, verbose=1)

print(f"\nLevenberg-Marquardt Result:")
print(result_lm)
print(f"\nTrue values:       E = {E_true:.0f}, sigma_Y = {sigma_Y_true:.0f}")
print(f"Identified values: E = {result_lm.x[0]:.0f}, sigma_Y = {result_lm.x[1]:.0f}")


###############################################################################
# Run identification with Differential Evolution (global search)
# --------------------------------------------------------------

print("\n" + "=" * 60)
print("Running Differential Evolution (global optimization)...")
print("=" * 60)

result_de = differential_evolution(
    problem,
    maxiter=50,
    popsize=10,
    polish=True,
    verbose=False,
)

print(f"\nDifferential Evolution Result:")
print(result_de)


###############################################################################
# Sensitivity analysis
# --------------------
# Check how sensitive the stress response is to each parameter.

print("\n" + "=" * 60)
print("Sensitivity Analysis")
print("=" * 60)

sensitivities = compute_sensitivity(problem, result_lm.x)
stress_sens = sensitivities['stress']

print(f"\nRelative sensitivity of stress to parameters:")
print(f"  Mean |dStress/dE * E/Stress|:       {np.mean(np.abs(stress_sens[:, 0])):.4f}")
print(f"  Mean |dStress/dSigmaY * SigmaY/Stress|: {np.mean(np.abs(stress_sens[:, 1])):.4f}")


###############################################################################
# Identifiability check
# ---------------------

print("\n" + "=" * 60)
print("Identifiability Check")
print("=" * 60)

ident_check = identifiability_check(problem, result_lm.x)

print(f"\nIdentifiable: {ident_check['identifiable']}")
print(f"Condition number: {ident_check['condition_number']:.2e}")
print(f"\nRecommendations:")
for rec in ident_check['recommendations']:
    print(f"  - {rec}")


###############################################################################
# Plot results
# ------------

fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# Plot 1: Stress-strain comparison
ax1 = axes[0]
ax1.plot(exp_strain * 100, exp_stress_noisy, 'ko', markersize=3, alpha=0.5, label='Experimental')
ax1.plot(exp_strain * 100, exp_stress, 'b-', linewidth=2, label='True model')

# Run with identified parameters
identified_params = np.array([result_lm.x[0], result_lm.x[1], H_true, n_true])
id_results = run_simulation(identified_params)
ax1.plot(id_results['strain'] * 100, id_results['stress'], 'r--', linewidth=2, label='Identified')

ax1.set_xlabel('Strain (%)')
ax1.set_ylabel('Stress (MPa)')
ax1.set_title('Stress-Strain Comparison')
ax1.legend()
ax1.grid(True, alpha=0.3)

# Plot 2: Sensitivity visualization
ax2 = axes[1]
strain_percent = exp_strain * 100
ax2.plot(strain_percent, np.abs(stress_sens[:, 0]), 'b-', label=f'E ({problem.parameter_names[0]})')
ax2.plot(strain_percent, np.abs(stress_sens[:, 1]), 'r-', label=f'sigma_Y ({problem.parameter_names[1]})')
ax2.set_xlabel('Strain (%)')
ax2.set_ylabel('Relative Sensitivity')
ax2.set_title('Parameter Sensitivity')
ax2.legend()
ax2.grid(True, alpha=0.3)

# Plot 3: Correlation matrix
ax3 = axes[2]
corr = ident_check['correlation_matrix']
im = ax3.imshow(corr, cmap='RdBu', vmin=-1, vmax=1)
ax3.set_xticks([0, 1])
ax3.set_yticks([0, 1])
ax3.set_xticklabels(problem.parameter_names)
ax3.set_yticklabels(problem.parameter_names)
ax3.set_title('Parameter Correlation')
plt.colorbar(im, ax=ax3)
for i in range(2):
    for j in range(2):
        ax3.text(j, i, f'{corr[i, j]:.3f}', ha='center', va='center')

plt.tight_layout()
plt.savefig('identification_results.png', dpi=150)
plt.show()

print("\nResults saved to 'identification_results.png'")


###############################################################################
# Summary
# -------

print("\n" + "=" * 60)
print("SUMMARY")
print("=" * 60)
print(f"\nTrue parameters:")
print(f"  E = {E_true:.0f} MPa")
print(f"  sigma_Y = {sigma_Y_true:.0f} MPa")

print(f"\nIdentified parameters (Levenberg-Marquardt):")
print(f"  E = {result_lm.x[0]:.0f} MPa (error: {100 * abs(result_lm.x[0] - E_true) / E_true:.2f}%)")
print(f"  sigma_Y = {result_lm.x[1]:.0f} MPa (error: {100 * abs(result_lm.x[1] - sigma_Y_true) / sigma_Y_true:.2f}%)")

print(f"\nIdentified parameters (Differential Evolution):")
print(f"  E = {result_de.x[0]:.0f} MPa")
print(f"  sigma_Y = {result_de.x[1]:.0f} MPa")

print(f"\nOptimization statistics:")
print(f"  LM iterations: {result_lm.n_iterations}")
print(f"  LM function evaluations: {result_lm.n_function_evals}")
print(f"  DE function evaluations: {result_de.n_function_evals}")
