"""
Example: Isotropic Plasticity with Isotropic Hardening (EPICP)

This example demonstrates:
1. Elasto-plastic material behavior
2. Cyclic loading with hysteresis
3. Tracking internal state variables (plastic strain, etc.)

The EPICP model uses:
- Isotropic J2 plasticity with associated flow
- Isotropic hardening: sigma_Y = sigma_Y0 + k * p^m
"""

import numpy as np
import matplotlib.pyplot as plt

from simcoon.solver import (
    Solver, Block, StepMeca, StateVariablesM,
    load_simulation_json
)


def example_cyclic_plasticity():
    """Cyclic loading of elasto-plastic material."""
    print("=" * 60)
    print("EPICP: Cyclic Plasticity Example")
    print("=" * 60)

    # Material properties for EPICP
    # E, nu, alpha, sigmaY, k, m
    E = 67538.0
    nu = 0.349
    alpha = 1e-6
    sigmaY = 300.0  # Initial yield stress
    k = 1500.0      # Hardening modulus
    m = 0.3         # Hardening exponent

    props = np.array([E, nu, alpha, sigmaY, k, m])

    # Cyclic loading: tension -> compression -> tension
    max_strain = 0.08

    step_tension = StepMeca(
        DEtot_end=np.array([max_strain, 0, 0, 0, 0, 0]),
        control=['strain', 'stress', 'stress', 'stress', 'stress', 'stress'],
        Dn_init=100,
        Dn_mini=10,
        Dn_inc=1000
    )

    step_compression = StepMeca(
        DEtot_end=np.array([-2*max_strain, 0, 0, 0, 0, 0]),
        control=['strain', 'stress', 'stress', 'stress', 'stress', 'stress'],
        Dn_init=100,
        Dn_mini=10,
        Dn_inc=1000
    )

    step_return = StepMeca(
        DEtot_end=np.array([max_strain, 0, 0, 0, 0, 0]),
        control=['strain', 'stress', 'stress', 'stress', 'stress', 'stress'],
        Dn_init=100,
        Dn_mini=10,
        Dn_inc=1000
    )

    block = Block(
        steps=[step_tension, step_compression, step_return],
        umat_name="EPICP",
        props=props,
        nstatev=8,  # EPICP has 8 internal state variables
        control_type='small_strain'
    )

    # Initialize
    sv = StateVariablesM(nstatev=8)
    sv.T = 323.15  # Temperature

    # Solve
    solver = Solver(blocks=[block], max_iter=20, tol=1e-9)
    history = solver.solve(sv)

    # Extract results
    strains_11 = [h.Etot[0] for h in history]
    stresses_11 = [h.sigma[0] for h in history]

    # State variable 0 is typically equivalent plastic strain
    plastic_strains = [h.statev[0] if len(h.statev) > 0 else 0 for h in history]

    print(f"Number of converged increments: {len(history)}")
    print(f"Max plastic strain: {max(plastic_strains):.6f}")

    return strains_11, stresses_11, plastic_strains


def example_from_json():
    """Load from JSON files."""
    print("\n" + "=" * 60)
    print("EPICP: JSON File Example")
    print("=" * 60)

    try:
        sim = load_simulation_json('material.json', 'path.json')

        sv = StateVariablesM(nstatev=sim['material']['nstatev'])
        sv.T = sim['initial_temperature']

        solver = Solver(blocks=sim['blocks'], max_iter=20, tol=1e-9)
        history = solver.solve(sv)

        strains = [h.Etot[0] for h in history]
        stresses = [h.sigma[0] for h in history]

        print(f"Loaded from JSON successfully")
        print(f"Material: {sim['material']['name']}")
        print(f"Converged increments: {len(history)}")

        return strains, stresses

    except FileNotFoundError:
        print("  (JSON files not found in current directory)")
        return None, None


if __name__ == '__main__':
    # Run cyclic plasticity example
    strains, stresses, plastic = example_cyclic_plasticity()

    # Plot results
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Stress-strain curve
    axes[0].plot(strains, stresses, 'b-', linewidth=1.5)
    axes[0].set_xlabel('Total Strain')
    axes[0].set_ylabel('Stress (MPa)')
    axes[0].set_title('EPICP: Cyclic Stress-Strain Response')
    axes[0].grid(True)
    axes[0].axhline(y=0, color='k', linewidth=0.5)
    axes[0].axvline(x=0, color='k', linewidth=0.5)

    # Plastic strain evolution
    axes[1].plot(range(len(plastic)), plastic, 'r-', linewidth=1.5)
    axes[1].set_xlabel('Increment')
    axes[1].set_ylabel('Equivalent Plastic Strain')
    axes[1].set_title('Plastic Strain Accumulation')
    axes[1].grid(True)

    plt.tight_layout()
    plt.savefig('epicp_results.png', dpi=150)
    plt.show()
