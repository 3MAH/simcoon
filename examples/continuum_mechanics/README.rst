Continuum Mechanics Examples
============================

This directory contains examples demonstrating Simcoon's continuum mechanics
operations, including tensor manipulations, stress measures, rotations, and
yield criteria.

Examples
--------

**constitutive_relations.py** - Stiffness and Compliance Tensors
    Demonstrates building stiffness (L) and compliance (M) tensors for various
    material symmetries: isotropic, transversely isotropic, orthotropic, and
    monoclinic.

    Key functions:

    - ``L_iso(E, nu)`` - Isotropic stiffness
    - ``M_iso(E, nu)`` - Isotropic compliance
    - ``L_isotrans(EL, ET, nuTL, nuTT, GLT, axis)`` - Transverse isotropy
    - ``L_ortho(...)`` - Orthotropic stiffness

**stress_measures.py** - Stress Tensor Conversions
    Shows conversions between different stress measures used in finite
    deformation mechanics:

    - Cauchy stress (sigma)
    - First Piola-Kirchhoff stress (PK1)
    - Second Piola-Kirchhoff stress (PK2)
    - Kirchhoff stress (tau)
    - Mandel stress

    Key functions:

    - ``Cauchy_to_PKII(sigma, F)``
    - ``PKII_to_Cauchy(PKII, F)``
    - ``Cauchy_to_PKI(sigma, F)``

**stress_transfer_helpers.py** - Push-Forward and Pull-Back Operations
    Demonstrates tensor transformation operations between reference and
    current configurations.

    Key concepts:

    - Covariant vs contravariant transformations
    - Push-forward and pull-back of stress tensors
    - Piola transformation

**rotation.py** - Tensor Rotations
    Shows how to rotate tensors using rotation matrices. Demonstrates rotation
    of stiffness tensors for computing effective properties of rotated materials.

    Key functions:

    - ``rotate_strain(E, R)`` - Rotate strain tensor
    - ``rotate_stress(sigma, R)`` - Rotate stress tensor
    - ``rotate_L(L, R)`` - Rotate 4th-order stiffness tensor
    - ``fillR_euler(psi, theta, phi)`` - Create rotation matrix from Euler angles

**yield_criteria.py** - Yield Function Evaluation
    Demonstrates evaluation of various yield criteria:

    - von Mises (J2 plasticity)
    - Tresca (maximum shear stress)
    - Drucker (pressure-dependent)
    - Hill (anisotropic)

    Key functions:

    - ``Mises_stress(sigma)`` - von Mises equivalent stress
    - ``tr(sigma)`` - Trace (hydrostatic component)
    - ``dev(sigma)`` - Deviatoric stress

Quick Start
-----------

**Building stiffness tensors:**

.. code-block:: python

    import simcoon as sim
    import numpy as np

    # Isotropic material
    E, nu = 210000.0, 0.3
    L = sim.L_iso(E, nu)  # 6x6 stiffness matrix
    M = sim.M_iso(E, nu)  # 6x6 compliance matrix

    # Verify L and M are inverses
    np.testing.assert_allclose(L @ M, np.eye(6), atol=1e-10)

    # Transversely isotropic (fiber along axis 1)
    EL, ET = 150000, 10000
    nuTL, nuTT = 0.3, 0.4
    GLT = 5000
    L_trans = sim.L_isotrans(EL, ET, nuTL, nuTT, GLT, axis=1)

**Stress conversions:**

.. code-block:: python

    import simcoon as sim
    import numpy as np

    # Cauchy stress (Voigt notation)
    sigma = np.array([100, 50, 0, 25, 0, 0])

    # Deformation gradient (simple extension)
    F = np.array([
        [1.1, 0, 0],
        [0, 0.95, 0],
        [0, 0, 0.95]
    ])

    # Convert to 2nd Piola-Kirchhoff
    PKII = sim.Cauchy_to_PKII(sigma, F)

    # Convert back
    sigma_back = sim.PKII_to_Cauchy(PKII, F)
    np.testing.assert_allclose(sigma, sigma_back, atol=1e-10)

**Rotating tensors:**

.. code-block:: python

    import simcoon as sim
    import numpy as np

    # Create rotation matrix (45 degrees about z-axis)
    psi, theta, phi = 45.0, 0.0, 0.0  # Euler angles in degrees
    R = sim.fillR_euler(np.radians(psi), np.radians(theta), np.radians(phi))

    # Rotate stiffness tensor
    L = sim.L_iso(210000, 0.3)
    L_rotated = sim.rotate_L(L, R)

    # For isotropic material, rotation has no effect
    np.testing.assert_allclose(L, L_rotated, atol=1e-10)

**Yield criteria:**

.. code-block:: python

    import simcoon as sim
    import numpy as np

    # Uniaxial stress state
    sigma = np.array([350, 0, 0, 0, 0, 0])  # MPa

    # von Mises equivalent stress
    sigma_eq = sim.Mises_stress(sigma)
    print(f"von Mises stress: {sigma_eq:.1f} MPa")  # 350 MPa

    # Deviatoric stress
    s = sim.dev(sigma)
    print(f"Deviatoric: {s}")

Running the Examples
--------------------

.. code-block:: bash

    # From the repository root
    cd examples/continuum_mechanics
    python constitutive_relations.py
    python stress_measures.py
    python rotation.py
    python yield_criteria.py

Voigt Notation Convention
-------------------------

Simcoon uses Voigt notation for symmetric tensors:

- **Strain:** ``[e11, e22, e33, 2*e12, 2*e13, 2*e23]``
- **Stress:** ``[s11, s22, s33, s12, s13, s23]``

The factor of 2 on shear strains ensures energy consistency: ``W = sigma . epsilon``.

See Also
--------

- :doc:`../umats/README` - Constitutive law examples
- :doc:`../analysis/README` - Analysis and post-processing
- `API Documentation <https://3mah.github.io/simcoon-docs/>`_
