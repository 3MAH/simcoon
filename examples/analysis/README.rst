Analysis and Processing Examples
=================================

This directory contains examples illustrating Simcoon's capabilities for simulating
mechanical and thermomechanical responses and post-processing results.

Examples
--------

**objective_rates.py** - Objective Stress Rates
    Compares different objective stress rates (Jaumann, Green-Naghdi, Logarithmic)
    for large deformation simulations. Demonstrates how the choice of corotational
    rate affects stress predictions under simple shear loading.

    Key concepts:

    - Corotational stress rates
    - Large strain formulation
    - Simple shear deformation

**directional_stiffness.py** - Directional Elastic Properties
    Computes and visualizes the directional Young's modulus for anisotropic
    materials. Shows how elastic stiffness varies with loading direction.

    Key concepts:

    - Stiffness tensor manipulation
    - Directional properties
    - 3D visualization of elastic anisotropy

**eshelby_numerical_vs_analytical.py** - Eshelby Tensor Validation
    Compares numerical and analytical Eshelby tensor calculations for various
    inclusion shapes (sphere, prolate, oblate ellipsoids). Validates the
    implementation against known analytical solutions.

    Key concepts:

    - Eshelby tensor computation
    - Micromechanics fundamentals
    - Numerical validation

Quick Start
-----------

All examples use the new Python Solver API (v2.0):

.. code-block:: python

    import numpy as np
    from simcoon.solver import Solver, Block, StepMeca

    # Define material and loading
    props = np.array([210000.0, 0.3, 1e-5])  # E, nu, alpha

    step = StepMeca(
        DEtot_end=np.array([0.1, 0, 0, 0, 0, 0]),
        control=['strain', 'stress', 'stress', 'stress', 'stress', 'stress'],
        Dn_init=100
    )

    block = Block(
        steps=[step],
        umat_name="ELISO",
        props=props,
        nstatev=1,
        corate_type='logarithmic'  # Try 'jaumann' or 'green_naghdi'
    )

    solver = Solver(blocks=[block])
    history = solver.solve()

Running the Examples
--------------------

.. code-block:: bash

    # From the repository root
    cd examples/analysis
    python objective_rates.py
    python directional_stiffness.py
    python eshelby_numerical_vs_analytical.py

See Also
--------

- :doc:`../umats/README` - Constitutive law examples
- :doc:`../continuum_mechanics/README` - Tensor operations and stress measures
- `Migration Guide <../../docs/migration_guide.md>`_ - Upgrading from v1.x
