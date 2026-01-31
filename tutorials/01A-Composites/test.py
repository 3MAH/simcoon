"""
Tutorial 01A: Composite Effective Properties using Mori-Tanaka.

This example demonstrates computing effective stiffness of a two-phase
composite using the Mori-Tanaka homogenization scheme (MIMTN).
"""

import numpy as np
import matplotlib.pyplot as plt
import simcoon as sim
from simcoon.solver import Ellipsoid, Phase, save_ellipsoids_json, save_phases_json
import os
import tempfile

# Create temporary directory for JSON files
with tempfile.TemporaryDirectory() as tmpdir:

    # Define matrix phase (phase 0)
    matrix = Phase(
        number=0,
        umat_name='ELISO',
        props=np.array([70000.0, 0.3, 0.0]),  # E, nu, alpha
        volume_fraction=0.7
    )

    # Define inclusion phase (phase 1) - spherical inclusions
    inclusion = Phase(
        number=1,
        umat_name='ELISO',
        props=np.array([400000.0, 0.2, 0.0]),  # E, nu, alpha (stiffer)
        volume_fraction=0.3
    )

    # Save phases to JSON
    phases_file = os.path.join(tmpdir, 'phases0.json')
    save_phases_json([matrix, inclusion], phases_file)

    # Define ellipsoidal geometry (spherical: a1=a2=a3=1)
    ellipsoid = Ellipsoid(
        number=0,
        coatingof=0,
        umat_name='ELISO',
        save=0,
        volume_fraction=0.3,
        a1=1.0, a2=1.0, a3=1.0,  # Sphere
        props=np.array([400000.0, 0.2, 0.0]),
        nstatev=1
    )

    # Save ellipsoids to JSON
    ellipsoids_file = os.path.join(tmpdir, 'ellipsoids0.json')
    save_ellipsoids_json([ellipsoid], ellipsoids_file)

    # Homogenization parameters
    nphases = 2
    num_file = 0  # File index
    int1 = 50     # Integration points theta
    int2 = 50     # Integration points phi
    n_matrix = 0  # Matrix phase number

    props = np.array([nphases, num_file, int1, int2, n_matrix], dtype='float')
    nstatev = 0

    # RVE orientation (no rotation)
    psi_rve = 0.0
    theta_rve = 0.0
    phi_rve = 0.0

    # Compute effective stiffness using Mori-Tanaka
    umat_name = "MIMTN"
    L_eff = sim.L_eff(umat_name, props, nstatev, psi_rve, theta_rve, phi_rve,
                      tmpdir, tmpdir)

    # Extract isotropic properties from effective stiffness
    E_eff, nu_eff = sim.L_iso_props(L_eff)[:2]

    print("Effective Properties (Mori-Tanaka):")
    print(f"  E_eff  = {E_eff:.1f} MPa")
    print(f"  nu_eff = {nu_eff:.4f}")
    print(f"\nEffective Stiffness Tensor L:")
    print(L_eff)
