"""
Example: Micromechanics Homogenization with JSON Configuration

This example demonstrates:
1. Loading phase configurations from JSON files
2. Creating phases programmatically
3. Converting legacy .dat files to JSON

Micromechanics models in Simcoon:
- MIHEN: Mori-Tanaka with ellipsoidal inclusions
- MIMTN: Mori-Tanaka with layers (laminates)
- MISCN: Self-consistent with ellipsoidal inclusions
- MIPLN: Self-consistent with layers

The micromechanics module can be used without building simcoon._core:
    from simcoon.solver.micromechanics import (
        Layer, Ellipsoid,
        load_layers_json, save_layers_json,
        load_ellipsoids_json, save_ellipsoids_json,
    )
"""

import os
from pathlib import Path

import numpy as np

# Import from the standalone micromechanics module
# This works without building simcoon._core
from simcoon.solver.micromechanics import (
    MaterialOrientation,
    GeometryOrientation,
    Layer,
    Ellipsoid,
    load_layers_json,
    save_layers_json,
    load_ellipsoids_json,
    save_ellipsoids_json,
    convert_legacy_layers,
    convert_legacy_ellipsoids,
)


def example_programmatic_layers():
    """Create laminate layers programmatically and save to JSON."""
    print("=" * 60)
    print("Example 1: Programmatic Layer Definition")
    print("=" * 60)

    # Define a 3-layer laminate: [0/90/0] composite
    layers = [
        Layer(
            number=0,
            umat_name="ELISO",
            save=1,
            concentration=0.4,
            material_orientation=MaterialOrientation(psi=0, theta=0, phi=0),
            geometry_orientation=GeometryOrientation(psi=0, theta=90, phi=-90),
            nstatev=1,
            props=np.array([70000.0, 0.3, 1e-5])
        ),
        Layer(
            number=1,
            umat_name="ELISO",
            save=1,
            concentration=0.2,
            material_orientation=MaterialOrientation(psi=90, theta=0, phi=0),
            geometry_orientation=GeometryOrientation(psi=0, theta=90, phi=-90),
            nstatev=1,
            props=np.array([70000.0, 0.3, 1e-5])
        ),
        Layer(
            number=2,
            umat_name="ELISO",
            save=1,
            concentration=0.4,
            material_orientation=MaterialOrientation(psi=0, theta=0, phi=0),
            geometry_orientation=GeometryOrientation(psi=0, theta=90, phi=-90),
            nstatev=1,
            props=np.array([70000.0, 0.3, 1e-5])
        ),
    ]

    print(f"Created {len(layers)} layers for [0/90/0] laminate:")
    for lyr in layers:
        print(f"  Layer {lyr.number}: c={lyr.concentration}, "
              f"mat_ori=({lyr.material_orientation.psi}, {lyr.material_orientation.theta}, {lyr.material_orientation.phi})")

    save_layers_json('laminate_090.json', layers, prop_names=['E', 'nu', 'alpha'])
    print("\nSaved to laminate_090.json")

    return layers


def example_programmatic_ellipsoids():
    """Create ellipsoidal inclusions programmatically and save to JSON."""
    print("\n" + "=" * 60)
    print("Example 2: Programmatic Ellipsoid Definition")
    print("=" * 60)

    # Glass fiber / Epoxy composite
    ellipsoids = [
        # Matrix phase (isotropic epoxy)
        Ellipsoid(
            number=0,
            coatingof=0,
            umat_name="ELISO",
            save=1,
            concentration=0.65,
            material_orientation=MaterialOrientation(psi=0, theta=0, phi=0),
            a1=1.0,
            a2=1.0,
            a3=1.0,
            geometry_orientation=GeometryOrientation(psi=0, theta=0, phi=0),
            nstatev=1,
            props=np.array([3500.0, 0.35, 60e-6])
        ),
        # Fiber inclusions (prolate spheroids with aspect ratio 20)
        Ellipsoid(
            number=1,
            coatingof=0,
            umat_name="ELISO",
            save=1,
            concentration=0.35,
            material_orientation=MaterialOrientation(psi=0, theta=0, phi=0),
            a1=20.0,
            a2=1.0,
            a3=1.0,
            geometry_orientation=GeometryOrientation(psi=0, theta=0, phi=0),
            nstatev=1,
            props=np.array([72000.0, 0.22, 5e-6])
        ),
    ]

    print("Glass fiber / Epoxy composite configuration:")
    for ell in ellipsoids:
        print(f"  Phase {ell.number}: {ell.shape_type}, c={ell.concentration}, "
              f"E={ell.props[0]:.0f} MPa, axes=({ell.a1}, {ell.a2}, {ell.a3})")

    save_ellipsoids_json('glass_epoxy.json', ellipsoids, prop_names=['E', 'nu', 'alpha'])
    print("\nSaved to glass_epoxy.json")

    # Simple bounds on effective modulus
    matrix, fibers = ellipsoids
    E_voigt = fibers.concentration * fibers.props[0] + matrix.concentration * matrix.props[0]
    E_reuss = 1 / (fibers.concentration / fibers.props[0] + matrix.concentration / matrix.props[0])

    print(f"\nSimple bounds on effective modulus:")
    print(f"  Voigt (upper): E = {E_voigt:.0f} MPa")
    print(f"  Reuss (lower): E = {E_reuss:.0f} MPa")

    return ellipsoids


def example_convert_legacy():
    """Convert legacy .dat files to JSON format."""
    print("\n" + "=" * 60)
    print("Example 3: Convert Legacy Files to JSON")
    print("=" * 60)

    legacy_dir = Path(__file__).parent.parent.parent.parent / 'testBin' / 'Libraries' / 'Phase' / 'data'

    legacy_layers = legacy_dir / 'Nlayers0.dat'
    if legacy_layers.exists():
        layers = convert_legacy_layers(legacy_layers)
        save_layers_json('converted_layers.json', layers)
        print(f"Converted {legacy_layers.name} -> converted_layers.json ({len(layers)} layers)")
    else:
        print(f"Legacy file not found: {legacy_layers}")

    legacy_ellipsoids = legacy_dir / 'Nellipsoids0.dat'
    if legacy_ellipsoids.exists():
        ellipsoids = convert_legacy_ellipsoids(legacy_ellipsoids)
        save_ellipsoids_json('converted_ellipsoids.json', ellipsoids)
        print(f"Converted {legacy_ellipsoids.name} -> converted_ellipsoids.json ({len(ellipsoids)} ellipsoids)")
    else:
        print(f"Legacy file not found: {legacy_ellipsoids}")


if __name__ == '__main__':
    # Change to the script directory
    script_dir = Path(__file__).parent
    os.chdir(script_dir)

    # Run examples
    example_programmatic_layers()
    example_programmatic_ellipsoids()
    example_convert_legacy()

    print("\n" + "=" * 60)
    print("All examples completed!")
    print("=" * 60)
