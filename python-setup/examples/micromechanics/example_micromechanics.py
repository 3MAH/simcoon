"""
Example: Micromechanics Homogenization with JSON Configuration

This example demonstrates:
1. Creating phases programmatically
2. Saving/loading JSON configurations
3. Working with layers and ellipsoids

Micromechanics models in Simcoon:
- MIHEN: Mori-Tanaka with ellipsoidal inclusions
- MIMTN: Mori-Tanaka with layers (laminates)
- MISCN: Self-consistent with ellipsoidal inclusions
- MIPLN: Self-consistent with layers
"""

import os
from pathlib import Path

import numpy as np

from simcoon.solver.micromechanics import (
    MaterialOrientation,
    GeometryOrientation,
    Layer,
    Ellipsoid,
    load_layers_json,
    save_layers_json,
    load_ellipsoids_json,
    save_ellipsoids_json,
)


def example_laminate():
    """Create a [0/90/0] laminate and save to JSON."""
    print("=" * 60)
    print("Example 1: Laminate Definition")
    print("=" * 60)

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
              f"orientation=({lyr.material_orientation.psi}, {lyr.material_orientation.theta}, {lyr.material_orientation.phi})")

    save_layers_json('laminate_090.json', layers, prop_names=['E', 'nu', 'alpha'])
    print("\nSaved to laminate_090.json")

    # Reload and verify
    loaded = load_layers_json('laminate_090.json')
    print(f"Reloaded {len(loaded)} layers")

    return layers


def example_fiber_composite():
    """Create a glass fiber / epoxy composite."""
    print("\n" + "=" * 60)
    print("Example 2: Fiber Composite")
    print("=" * 60)

    ellipsoids = [
        # Matrix (epoxy)
        Ellipsoid(
            number=0,
            coatingof=0,
            umat_name="ELISO",
            save=1,
            concentration=0.65,
            material_orientation=MaterialOrientation(psi=0, theta=0, phi=0),
            a1=1.0, a2=1.0, a3=1.0,  # Sphere
            geometry_orientation=GeometryOrientation(psi=0, theta=0, phi=0),
            nstatev=1,
            props=np.array([3500.0, 0.35, 60e-6])
        ),
        # Fibers (prolate spheroids, aspect ratio 20)
        Ellipsoid(
            number=1,
            coatingof=0,
            umat_name="ELISO",
            save=1,
            concentration=0.35,
            material_orientation=MaterialOrientation(psi=0, theta=0, phi=0),
            a1=20.0, a2=1.0, a3=1.0,
            geometry_orientation=GeometryOrientation(psi=0, theta=0, phi=0),
            nstatev=1,
            props=np.array([72000.0, 0.22, 5e-6])
        ),
    ]

    print("Glass fiber / Epoxy composite:")
    for ell in ellipsoids:
        print(f"  Phase {ell.number}: {ell.shape_type}, c={ell.concentration}, "
              f"E={ell.props[0]:.0f} MPa")

    save_ellipsoids_json('glass_epoxy.json', ellipsoids, prop_names=['E', 'nu', 'alpha'])
    print("\nSaved to glass_epoxy.json")

    # Voigt-Reuss bounds
    matrix, fibers = ellipsoids
    E_voigt = fibers.concentration * fibers.props[0] + matrix.concentration * matrix.props[0]
    E_reuss = 1 / (fibers.concentration / fibers.props[0] + matrix.concentration / matrix.props[0])
    print(f"\nSimple bounds:")
    print(f"  Voigt: E = {E_voigt:.0f} MPa")
    print(f"  Reuss: E = {E_reuss:.0f} MPa")

    return ellipsoids


def example_coated_inclusions():
    """Create a composite with coated inclusions."""
    print("\n" + "=" * 60)
    print("Example 3: Coated Inclusions")
    print("=" * 60)

    ellipsoids = [
        # Matrix
        Ellipsoid(
            number=0,
            coatingof=0,
            umat_name="ELISO",
            concentration=0.6,
            a1=1.0, a2=1.0, a3=1.0,
            nstatev=1,
            props=np.array([3000.0, 0.4, 50e-6])
        ),
        # Core (stiff particle)
        Ellipsoid(
            number=1,
            coatingof=0,
            umat_name="ELISO",
            concentration=0.2,
            a1=1.0, a2=1.0, a3=1.0,
            nstatev=1,
            props=np.array([400000.0, 0.2, 5e-6])
        ),
        # Coating (interphase)
        Ellipsoid(
            number=2,
            coatingof=1,  # Coats phase 1
            umat_name="ELISO",
            concentration=0.2,
            a1=1.0, a2=1.0, a3=1.0,
            nstatev=1,
            props=np.array([10000.0, 0.35, 30e-6])
        ),
    ]

    print("Coated particle composite:")
    for ell in ellipsoids:
        coating_info = f" (coats phase {ell.coatingof})" if ell.coatingof > 0 else ""
        print(f"  Phase {ell.number}: c={ell.concentration}, E={ell.props[0]:.0f} MPa{coating_info}")

    save_ellipsoids_json('coated_particles.json', ellipsoids, prop_names=['E', 'nu', 'alpha'])
    print("\nSaved to coated_particles.json")

    return ellipsoids


if __name__ == '__main__':
    script_dir = Path(__file__).parent
    os.chdir(script_dir)

    example_laminate()
    example_fiber_composite()
    example_coated_inclusions()

    print("\n" + "=" * 60)
    print("All examples completed!")
    print("=" * 60)
