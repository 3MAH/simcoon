Heterogeneous Materials Examples
================================

This directory contains examples demonstrating Simcoon's capabilities for simulating
heterogeneous materials, including Eshelby tensor computations and micromechanical
homogenization schemes for predicting effective properties of composite materials.

Examples
--------

**effective_props.py** - Effective Properties vs Volume Fraction and Aspect Ratio
    Studies the evolution of the mechanical properties of a 2-phase composite material
    considering spherical and ellipsoidal reinforcements. Compares Mori-Tanaka and
    Self-Consistent schemes, and explores the effect of inclusion aspect ratio on
    effective stiffness.

    Key concepts:

    - Mori-Tanaka (MIMTN) and Self-Consistent (MISCN) schemes
    - Effect of volume fraction on effective properties
    - Effect of aspect ratio (prolate vs oblate inclusions)
    - Comparison with experimental data

**eshelby_tensors.py** - Eshelby Tensor Computations
    Demonstrates Eshelby tensor computations for various inclusion shapes including
    spheres, cylinders, prolate and oblate ellipsoids, and penny-shaped cracks.

    Key concepts:

    - Analytical Eshelby tensors for special cases
    - Numerical integration for general ellipsoids
    - Hill's interaction tensor
    - Convergence of oblate ellipsoid to penny-shaped crack limit

**homogenization.py** - Two-Phase Composite Homogenization
    Tutorial studying the mechanical properties of a 2-phase composite material,
    visualizing how Eshelby tensor components vary with aspect ratio.

    Key concepts:

    - Effective elastic properties computation
    - Eshelby tensor component variation with aspect ratio
    - Isotropic effective properties from isotropic phases

Phase Configuration
-------------------

Phase configurations (ellipsoids, layers, cylinders) are defined using the
Python ``simcoon.solver.micromechanics`` module with JSON format:

.. code-block:: python

    from simcoon.solver.micromechanics import (
        Ellipsoid, load_ellipsoids_json, save_ellipsoids_json,
    )

    # Load phases from JSON
    phases = load_ellipsoids_json('data/ellipsoids0.json')

    # Modify programmatically
    phases[1].a1 = 20.0  # Change aspect ratio
    phases[1].concentration = 0.35

    # Save back
    save_ellipsoids_json('data/ellipsoids0.json', phases)

Quick Start
-----------

**Computing effective properties:**

.. code-block:: python

    import numpy as np
    import simcoon as sim

    # Define micromechanics parameters
    nstatev = 0
    nphases = 2
    num_file = 0  # Uses ellipsoids0.json, phases0.json
    int1 = 50     # Integration points
    int2 = 50
    n_matrix = 0  # Matrix phase index

    props = np.array([nphases, num_file, int1, int2, n_matrix], dtype='float')

    # Euler angles for RVE orientation (z-x-z convention)
    psi_rve, theta_rve, phi_rve = 0.0, 0.0, 0.0

    # Compute effective stiffness using Mori-Tanaka scheme
    L_eff = sim.L_eff('MIMTN', props, nstatev, psi_rve, theta_rve, phi_rve)

    # Extract isotropic properties (E, nu) if applicable
    E_nu = sim.L_iso_props(L_eff).flatten()
    print(f"Effective E = {E_nu[0]:.1f} MPa, nu = {E_nu[1]:.4f}")

**Computing Eshelby tensors:**

.. code-block:: python

    import numpy as np
    import simcoon as sim

    nu = 0.3  # Matrix Poisson ratio

    # Analytical solutions for special shapes
    S_sphere = sim.Eshelby_sphere(nu)
    S_cylinder = sim.Eshelby_cylinder(nu)
    S_prolate = sim.Eshelby_prolate(nu, ar=5.0)   # ar = a1/a3 > 1
    S_oblate = sim.Eshelby_oblate(nu, ar=0.2)     # ar = a1/a3 < 1
    S_penny = sim.Eshelby_penny(nu)               # Penny-shaped crack

    # Numerical solution for general ellipsoid in anisotropic matrix
    E = 70000.0  # Young's modulus (MPa)
    L = sim.L_iso([E, nu], 'Enu')
    a1, a2, a3 = 2.0, 1.0, 0.5  # Semi-axes
    mp, np_int = 50, 50         # Integration points

    S_general = sim.Eshelby(L, a1, a2, a3, mp, np_int)

    # Hill's interaction tensor
    T_II = sim.T_II(L, a1, a2, a3, mp, np_int)

Running the Examples
--------------------

.. code-block:: bash

    # From the repository root
    cd examples/heterogeneous
    python effective_props.py
    python eshelby_tensors.py
    python homogenization.py

Data Files
----------

The ``data/`` directory contains:

- ``ellipsoids0.json`` - Ellipsoid phase definitions (geometry, orientation, material)
- ``phases0.json`` - Phase material properties
- ``E_exp.txt`` - Experimental data for validation (Wang 2003)

Micromechanical Schemes
-----------------------

Simcoon supports the following homogenization schemes:

- **MIMTN** - Mori-Tanaka scheme (dilute to moderate concentrations)
- **MISCN** - Self-Consistent scheme (higher concentrations, interpenetrating phases)

The Mori-Tanaka scheme is generally more accurate for matrix-inclusion composites
at low to moderate reinforcement volume fractions, while the Self-Consistent scheme
is appropriate for interpenetrating microstructures or higher concentrations.

See Also
--------

- :doc:`../analysis/README` - Eshelby tensor validation examples
- :doc:`../continuum_mechanics/README` - Tensor operations and rotations
- `API Documentation <https://3mah.github.io/simcoon-docs/>`_
