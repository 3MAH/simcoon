"""
Tutorial 01B: Parametric Study of Composite Effective Properties.

This example demonstrates computing effective stiffness as a function
of inclusion volume fraction, comparing Mori-Tanaka (MIMTN) and
Self-Consistent (MISCN) homogenization schemes.
"""

import numpy as np
import matplotlib.pyplot as plt
import simcoon as sim
from simcoon.solver import Ellipsoid, Phase, save_ellipsoids_json, save_phases_json
import os
import tempfile

# Create temporary directory for JSON files
with tempfile.TemporaryDirectory() as tmpdir:

    # Material properties
    E_matrix = 70000.0      # Matrix Young's modulus (MPa)
    nu_matrix = 0.3         # Matrix Poisson's ratio
    E_inclusion = 400000.0  # Inclusion Young's modulus (MPa)
    nu_inclusion = 0.2      # Inclusion Poisson's ratio

    # Homogenization parameters
    nphases = 2
    num_file = 0
    int1 = 50
    int2 = 50
    n_matrix = 0

    props = np.array([nphases, num_file, int1, int2, n_matrix], dtype='float')
    nstatev = 0

    # RVE orientation
    psi_rve = 0.0
    theta_rve = 0.0
    phi_rve = 0.0

    # Volume fraction range
    concentration = np.arange(0.0, 0.51, 0.05)

    E_MT = np.zeros(len(concentration))  # Mori-Tanaka results
    E_SC = np.zeros(len(concentration))  # Self-Consistent results

    for i, vf in enumerate(concentration):
        # Update phases with current volume fraction
        matrix = Phase(
            number=0,
            umat_name='ELISO',
            props=np.array([E_matrix, nu_matrix, 0.0]),
            volume_fraction=1.0 - vf
        )

        inclusion = Phase(
            number=1,
            umat_name='ELISO',
            props=np.array([E_inclusion, nu_inclusion, 0.0]),
            volume_fraction=vf
        )

        # Save phases
        phases_file = os.path.join(tmpdir, 'phases0.json')
        save_phases_json([matrix, inclusion], phases_file)

        # Update ellipsoid geometry
        ellipsoid = Ellipsoid(
            number=0,
            coatingof=0,
            umat_name='ELISO',
            save=0,
            volume_fraction=vf,
            a1=1.0, a2=1.0, a3=1.0,  # Sphere
            props=np.array([E_inclusion, nu_inclusion, 0.0]),
            nstatev=1
        )

        # Save ellipsoids
        ellipsoids_file = os.path.join(tmpdir, 'ellipsoids0.json')
        save_ellipsoids_json([ellipsoid], ellipsoids_file)

        # Compute Mori-Tanaka
        L_MT = sim.L_eff("MIMTN", props, nstatev, psi_rve, theta_rve, phi_rve,
                         tmpdir, tmpdir)
        E_MT[i] = sim.L_iso_props(L_MT)[0]

        # Compute Self-Consistent
        L_SC = sim.L_eff("MISCN", props, nstatev, psi_rve, theta_rve, phi_rve,
                         tmpdir, tmpdir)
        E_SC[i] = sim.L_iso_props(L_SC)[0]

        print(f"vf = {vf:.2f}: E_MT = {E_MT[i]:.1f}, E_SC = {E_SC[i]:.1f}")

    # Plot results
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(concentration, E_MT, 'b-', linewidth=2, label='Mori-Tanaka')
    ax.plot(concentration, E_SC, 'r--', linewidth=2, label='Self-Consistent')

    # Voigt and Reuss bounds
    E_voigt = (1 - concentration) * E_matrix + concentration * E_inclusion
    E_reuss = 1.0 / ((1 - concentration) / E_matrix + concentration / E_inclusion)
    ax.fill_between(concentration, E_reuss, E_voigt, alpha=0.2, color='gray',
                    label='Voigt-Reuss bounds')

    ax.set_xlabel('Inclusion Volume Fraction', fontsize=12)
    ax.set_ylabel('Effective Young\'s Modulus (MPa)', fontsize=12)
    ax.set_title('Composite Effective Properties: MT vs SC', fontsize=14)
    ax.legend(loc='upper left')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('composite_effective_properties.png', dpi=150)
    plt.show()

    print("\nPlot saved to composite_effective_properties.png")
