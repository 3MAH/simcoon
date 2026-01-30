"""
Example: Isotropic Linear Elasticity (ELISO) with Python Solver

This example demonstrates:
1. Loading material/path from JSON files
2. Using the Python solver directly (programmatic API)
3. Plotting stress-strain curves

The Python solver replaces the deprecated C++ solver (sim.solver) with a more
flexible Python-native implementation.
"""

import numpy as np
import matplotlib.pyplot as plt

# Import from the new solver module
from simcoon.solver import (
    Solver, Block, StepMeca, StateVariablesM,
    load_material_json, load_path_json, load_simulation_json
)


def example_programmatic():
    """Example using programmatic API (no files needed)."""
    print("=" * 60)
    print("Example 1: Programmatic API")
    print("=" * 60)

    # Material properties for ELISO (E, nu, alpha)
    E = 70000.0
    nu = 0.3
    alpha = 1e-5
    props = np.array([E, nu, alpha])

    # Create uniaxial tension step
    # - Strain controlled in direction 11
    # - Stress-free in all other directions (uniaxial condition)
    step = StepMeca(
        DEtot_end=np.array([0.02, 0, 0, 0, 0, 0]),  # 2% axial strain
        Dsigma_end=np.zeros(6),
        control=['strain', 'stress', 'stress', 'stress', 'stress', 'stress'],
        Dn_init=10,
        Dn_inc=100,
        time=1.0
    )

    # Create block with material
    block = Block(
        steps=[step],
        umat_name="ELISO",
        props=props,
        nstatev=1,
        control_type='small_strain',
        corate_type='jaumann'
    )

    # Initialize state variables
    sv = StateVariablesM(nstatev=1)
    sv.T = 290.0  # Initial temperature

    # Create and run solver
    solver = Solver(blocks=[block], max_iter=10, tol=1e-9)
    history = solver.solve(sv)

    # Extract results
    strains_11 = [h.Etot[0] for h in history]
    stresses_11 = [h.sigma[0] for h in history]

    print(f"Number of converged increments: {len(history)}")
    print(f"Final strain: {strains_11[-1]:.6f}")
    print(f"Final stress: {stresses_11[-1]:.2f} MPa")
    print(f"Effective E (from stress/strain): {stresses_11[-1]/strains_11[-1]:.2f} MPa")

    return strains_11, stresses_11


def example_json_files():
    """Example loading from JSON files."""
    print("\n" + "=" * 60)
    print("Example 2: JSON File API")
    print("=" * 60)

    # Load simulation from JSON files
    sim_config = load_simulation_json('material.json', 'path.json')

    # Initialize state variables
    sv = StateVariablesM(nstatev=sim_config['material']['nstatev'])
    sv.T = sim_config['initial_temperature']

    # Run solver
    solver = Solver(blocks=sim_config['blocks'], max_iter=10, tol=1e-9)
    history = solver.solve(sv)

    # Extract results
    strains_11 = [h.Etot[0] for h in history]
    stresses_11 = [h.sigma[0] for h in history]

    print(f"Material: {sim_config['material']['name']}")
    print(f"Initial temperature: {sim_config['initial_temperature']} K")
    print(f"Number of blocks: {len(sim_config['blocks'])}")
    print(f"Number of converged increments: {len(history)}")

    return strains_11, stresses_11


def example_cyclic():
    """Example with cyclic loading."""
    print("\n" + "=" * 60)
    print("Example 3: Cyclic Loading")
    print("=" * 60)

    E = 70000.0
    nu = 0.3
    props = np.array([E, nu, 1e-5])

    # Loading step
    step_load = StepMeca(
        DEtot_end=np.array([0.02, 0, 0, 0, 0, 0]),
        control=['strain', 'stress', 'stress', 'stress', 'stress', 'stress'],
        Dn_init=10
    )

    # Unloading step
    step_unload = StepMeca(
        DEtot_end=np.array([-0.02, 0, 0, 0, 0, 0]),
        control=['strain', 'stress', 'stress', 'stress', 'stress', 'stress'],
        Dn_init=10
    )

    # Block with both steps, repeated 3 times
    block = Block(
        steps=[step_load, step_unload],
        umat_name="ELISO",
        props=props,
        nstatev=1,
        ncycle=3  # Repeat 3 cycles
    )

    sv = StateVariablesM(nstatev=1)
    solver = Solver(blocks=[block])
    history = solver.solve(sv)

    strains_11 = [h.Etot[0] for h in history]
    stresses_11 = [h.sigma[0] for h in history]

    print(f"Number of cycles: 3")
    print(f"Total increments: {len(history)}")

    return strains_11, stresses_11


if __name__ == '__main__':
    # Run examples
    strains1, stresses1 = example_programmatic()

    try:
        strains2, stresses2 = example_json_files()
    except FileNotFoundError:
        print("  (Skipping JSON example - files not in current directory)")
        strains2, stresses2 = None, None

    strains3, stresses3 = example_cyclic()

    # Plot results
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))

    axes[0].plot(strains1, stresses1, 'b-', linewidth=2)
    axes[0].set_xlabel('Strain')
    axes[0].set_ylabel('Stress (MPa)')
    axes[0].set_title('Programmatic API')
    axes[0].grid(True)

    if strains2:
        axes[1].plot(strains2, stresses2, 'r-', linewidth=2)
        axes[1].set_xlabel('Strain')
        axes[1].set_ylabel('Stress (MPa)')
        axes[1].set_title('JSON File API')
        axes[1].grid(True)

    axes[2].plot(strains3, stresses3, 'g-', linewidth=2)
    axes[2].set_xlabel('Strain')
    axes[2].set_ylabel('Stress (MPa)')
    axes[2].set_title('Cyclic Loading')
    axes[2].grid(True)

    plt.tight_layout()
    plt.savefig('eliso_results.png', dpi=150)
    plt.show()
