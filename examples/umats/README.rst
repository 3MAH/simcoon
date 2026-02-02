Constitutive Laws Examples
==========================

This directory contains examples demonstrating Simcoon's constitutive laws library
using the Python Solver API (v2.0).

Available Material Models
-------------------------

**Elastic Models:**

- **ELISO** - Isotropic linear elasticity
- **ELIST** - Transversely isotropic elasticity
- **ELORT** - Orthotropic elasticity

**Plasticity Models:**

- **EPICP** - Plasticity with isotropic hardening (power-law)
- **EPKCP** - Plasticity with combined isotropic and kinematic hardening
- **EPCHA** - Plasticity with Chaboche hardening (cyclic plasticity)

Examples
--------

**ELISO.py** - Isotropic Elasticity
    Demonstrates uniaxial tension simulation with isotropic elastic material.
    Shows basic usage of `Solver`, `Block`, and `StepMeca` classes.

    Material parameters: ``[E, nu, alpha]``

**ELIST.py** - Transversely Isotropic Elasticity
    Simulates uniaxial tension for a fiber-reinforced composite material with
    transverse isotropy. Compares loading parallel and perpendicular to fibers.

    Material parameters: ``[EL, ET, nuTL, nuTT, GLT, alphaL, alphaT]``

**ELORT.py** - Orthotropic Elasticity
    Demonstrates fully orthotropic elastic behavior typical of wood or layered
    composites. Shows directional dependence of elastic response.

    Material parameters: ``[E1, E2, E3, nu12, nu13, nu23, G12, G13, G23, alpha1, alpha2, alpha3]``

**EPICP.py** - Isotropic Plasticity
    Simulates plastic deformation with power-law isotropic hardening.
    Demonstrates elastic-plastic transition and hardening behavior.

    Material parameters: ``[E, nu, alpha, sigma_Y, H, n]``

**EPKCP.py** - Combined Hardening Plasticity
    Shows combined isotropic and kinematic hardening behavior. Demonstrates
    the Bauschinger effect under cyclic loading.

    Material parameters: ``[E, nu, alpha, sigma_Y, H, n, C, gamma]``

**EPCHA.py** - Chaboche Cyclic Plasticity
    Advanced cyclic plasticity with multiple backstress tensors. Suitable for
    predicting fatigue behavior and ratcheting.

    Material parameters: ``[E, nu, alpha, sigma_Y, H, n, C1, gamma1, C2, gamma2]``

Quick Start
-----------

.. code-block:: python

    import numpy as np
    from simcoon.solver import Solver, Block, StepMeca

    # Define isotropic elastic material
    props = np.array([210000.0, 0.3, 1e-5])  # E, nu, alpha

    # Create uniaxial tension step
    step = StepMeca(
        DEtot_end=np.array([0.01, 0, 0, 0, 0, 0]),  # 1% axial strain
        control=['strain', 'stress', 'stress', 'stress', 'stress', 'stress'],
        Dn_init=50
    )

    # Create simulation block
    block = Block(
        steps=[step],
        umat_name="ELISO",
        props=props,
        nstatev=1,
        control_type='small_strain',
        corate_type='logarithmic'
    )

    # Run simulation
    solver = Solver(blocks=[block])
    history = solver.solve()

    # Extract results
    strain = np.array([h.Etot[0] for h in history])
    stress = np.array([h.sigma[0] for h in history])

Cyclic Loading Example
----------------------

.. code-block:: python

    import numpy as np
    from simcoon.solver import Solver, Block, StepMeca

    # Plasticity with kinematic hardening
    props = np.array([200000, 0.3, 0, 350, 1000, 0.4, 50000, 200])

    # Tension step
    step1 = StepMeca(
        DEtot_end=np.array([0.02, 0, 0, 0, 0, 0]),
        control=['strain', 'stress', 'stress', 'stress', 'stress', 'stress'],
        Dn_init=100
    )

    # Compression step (shows Bauschinger effect)
    step2 = StepMeca(
        DEtot_end=np.array([-0.02, 0, 0, 0, 0, 0]),
        control=['strain', 'stress', 'stress', 'stress', 'stress', 'stress'],
        Dn_init=100
    )

    block = Block(
        steps=[step1, step2],
        umat_name="EPKCP",
        props=props,
        nstatev=10
    )

    solver = Solver(blocks=[block])
    history = solver.solve()

Running the Examples
--------------------

.. code-block:: bash

    # From the repository root
    cd examples/umats
    python ELISO.py
    python EPICP.py
    python EPCHA.py

State Variables
---------------

Each history entry contains:

- ``Etot`` - Total strain (Voigt: [e11, e22, e33, e12, e13, e23])
- ``sigma`` - Cauchy stress (Voigt notation)
- ``statev`` - Internal state variables (model-dependent)
- ``F0``, ``F1`` - Deformation gradient (start/end of increment)
- ``Lt`` - Tangent stiffness matrix (6x6)
- ``Wm`` - Work measures [Wm, Wm_r, Wm_ir, Wm_d]

See Also
--------

- :doc:`../analysis/README` - Analysis and post-processing examples
- :doc:`../continuum_mechanics/README` - Tensor operations
- `Migration Guide <../../docs/migration_guide.md>`_ - Upgrading from v1.x
