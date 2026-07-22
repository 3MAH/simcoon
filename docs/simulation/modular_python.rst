==========================================
Building materials in Python (``modular``)
==========================================

.. module:: simcoon.modular

The :mod:`simcoon.modular` package is the Python builder for the composable
``MODUL`` engine: a constitutive model is assembled from an **elasticity**
definition plus any combination of **mechanisms** (plasticity,
viscoelasticity, damage), and serialized into the self-describing props
stream that the C++ engine consumes. The same object drives the in-memory
solver, the file solver, or any FEA coupling — see the
:doc:`UMAT catalog <umat_catalog>` for how ``MODUL`` fits among the other
material names. Units are MPa; Voigt order is [11, 22, 33, 12, 13, 23].

Quick start
-----------

.. code-block:: python

    from simcoon.modular import (
        ModularMaterial, IsotropicElasticity,
        Plasticity, VonMisesYield, VoceHardening, ArmstrongFrederickHardening,
    )
    from simcoon import solver

    mat = ModularMaterial(
        elasticity=IsotropicElasticity(C1=210000.0, C2=0.3, alpha=1.2e-5,
                                       convention="Enu"),
        mechanisms=[
            Plasticity(
                sigma_Y=300.0,
                yield_criterion=VonMisesYield(),
                isotropic_hardening=VoceHardening(Q=200.0, b=10.0),
                kinematic_hardening=ArmstrongFrederickHardening(C=20000.0, D=100.0),
            ),
        ],
    )
    print(mat.summary())

    step = solver.StepMeca(control=['strain'] + ['stress'] * 5,
                           value=[0.02, 0, 0, 0, 0, 0], ninc=100)
    res = solver.solve(step, "MODUL", mat.props, mat.nstatev)

``mat.props`` and ``mat.nstatev`` are all any caller needs — they work
identically with ``material.dat`` files, :func:`simcoon.umat` point
evaluation, micromechanics phase files (``Nellipsoids``/``Nlayers``) and FEA
couplings such as fedoo.

Elasticity
----------

Exactly one elasticity definition per material. The elastic constants are
**ordinal slots** ``C1..Cn`` whose meaning is fixed by the ``convention``
argument (enum, int, or string aliases):

.. list-table::
   :header-rows: 1
   :widths: 26 34 40

   * - Class
     - Parameters
     - Conventions (string aliases)
   * - :class:`IsotropicElasticity`
     - ``C1, C2, alpha``
     - ``"Enu"`` (C1=E, C2=nu, default), ``"nuE"``, ``"Kmu"``/``"KG"``,
       ``"muK"``, ``"lambdamu"``, ``"mulambda"``
   * - :class:`CubicElasticity`
     - ``C1, C2, C3, alpha``
     - ``"EnuG"`` (default), ``"Cii"`` (C11, C12, C44)
   * - :class:`TransverseIsotropicElasticity`
     - ``EL, ET, nuTL, nuTT, GLT, alpha_L, alpha_T, axis``
     - ``"EnuG"`` only (axis = isotropy axis, default 3)
   * - :class:`OrthotropicElasticity`
     - ``C1..C9, alpha1..alpha3``
     - ``"EnuG"`` (E1,E2,E3,nu12,nu13,nu23,G12,G13,G23; default), ``"Cii"``

``alpha`` are the thermal-expansion coefficients (per direction where
applicable).

Yield criteria
--------------

Used inside :class:`Plasticity`; the criterion defines the equivalent stress
:math:`\sigma_{eq}` of the yield function
:math:`f = \sigma_{eq}(\boldsymbol{\sigma} - \mathbf{X}) - (\sigma_Y + R(p))`:

.. list-table::
   :header-rows: 1
   :widths: 28 30 42

   * - Class
     - Parameters
     - Criterion
   * - :class:`VonMisesYield`
     - —
     - :math:`\sqrt{\tfrac{3}{2}\,\mathbf{s}:\mathbf{s}}`
   * - :class:`TrescaYield`
     - —
     - Maximum shear stress
   * - :class:`DruckerYield`
     - ``b, n``
     - Drucker :math:`J_2`–:math:`J_3` criterion
   * - :class:`HillYield`
     - ``F, G, H, L, M, N``
     - Hill 1948 quadratic anisotropy
   * - :class:`DFAYield`
     - ``F, G, H, L, M, N, K``
     - Deshpande–Fleck–Ashby (pressure-sensitive)
   * - :class:`AnisotropicYield`
     - ``P11, P22, P33, P12, P13, P23, P44, P55, P66``
     - Full quadratic form. **P must be an admissible (convex) form**:
       symmetric with zero row sums on the normal block — an indefinite P
       produces ``sqrt(<0) = NaN``.

Isotropic hardening
-------------------

The isotropic hardening law :math:`R(p)` of the accumulated plastic strain
:math:`p`:

.. list-table::
   :header-rows: 1
   :widths: 30 24 46

   * - Class
     - Parameters
     - Law
   * - :class:`NoIsotropicHardening`
     - —
     - :math:`R = 0` (perfect plasticity, default)
   * - :class:`LinearIsotropicHardening`
     - ``H``
     - :math:`R = H\,p`
   * - :class:`PowerLawHardening`
     - ``k, m``
     - :math:`R = k\,p^m` — for :math:`m < 1` the singular onset slope is
       C\ :sup:`1`-regularized below :math:`p = 10^{-6}` (exact above)
   * - :class:`VoceHardening`
     - ``Q, b``
     - :math:`R = Q\,(1 - e^{-b\,p})`
   * - :class:`CombinedVoceHardening`
     - ``terms = ((Q_1, b_1), ...)``
     - :math:`R = \sum_i Q_i\,(1 - e^{-b_i p})` (standard independent sum —
       note this differs from the removed legacy EPCHG coupling, see the
       :doc:`catalog <umat_catalog>`)

Kinematic hardening
-------------------

The back stress :math:`\mathbf{X}` shifting the yield surface. The stored
internal variable is the **back-strain** :math:`\boldsymbol{\alpha}`
(thermodynamic variable); :math:`\mathbf{X} = \tfrac{2}{3} C \boldsymbol{\alpha}`:

.. list-table::
   :header-rows: 1
   :widths: 32 26 42

   * - Class
     - Parameters
     - Law
   * - :class:`NoKinematicHardening`
     - —
     - :math:`\mathbf{X} = 0` (default)
   * - :class:`PragerHardening`
     - ``C``
     - Linear: :math:`\dot{\boldsymbol{\alpha}} = \dot{p}\,\mathbf{n}`
   * - :class:`ArmstrongFrederickHardening`
     - ``C, D``
     - :math:`\dot{\mathbf{X}} = \tfrac{2}{3}C\,\dot{\boldsymbol{\varepsilon}}^p - D\,\mathbf{X}\dot{p}`
   * - :class:`ChabocheHardening`
     - ``terms = ((C_1, D_1), ...)``
     - :math:`\mathbf{X} = \sum_i \mathbf{X}_i`, each branch
       Armstrong–Frederick

Mechanisms
----------

.. list-table::
   :header-rows: 1
   :widths: 26 40 34

   * - Class
     - Parameters
     - Physics
   * - :class:`Plasticity`
     - ``sigma_Y, yield_criterion, isotropic_hardening, kinematic_hardening``
     - Rate-independent plasticity (Fischer–Burmeister return mapping)
   * - :class:`Viscoelasticity`
     - ``terms = ((E_i, nu_i, etaB_i, etaS_i), ...)``
     - Generalized Maxwell (Prony) branches: per branch a spring
       (:math:`E_i, \nu_i`) in series with bulk/shear dashpots
       (:math:`\eta_B, \eta_S`) — same rheology and layout as the kept
       ``PRONK`` kernel
   * - :class:`Damage`
     - ``Y_0, Y_c, damage_type, A, n``
     - Scalar stiffness-degradation damage; ``damage_type`` selects the
       evolution law: ``LINEAR``, ``EXPONENTIAL``, ``POWER_LAW`` (uses
       ``A, n``) or ``WEIBULL``

Multiple mechanisms compose additively on the inelastic strain; the
registration order defines the statev layout (see the
:doc:`catalog <umat_catalog>` statev section).

Tangent operator and finite strain
----------------------------------

``MODUL`` honors the solver's ``tangent_mode`` (continuum or algorithmic,
algorithmic being the 2.0 default — :doc:`solver`). Under the finite-strain
control types the composition acts as a Hencky hyperelastic law on the
logarithmic strain and requires ``corate_type = 3`` (log_R).

See also
--------

- :doc:`umat_catalog` — where ``MODUL`` and the adapter-served legacy names
  meet, props streams and statev layouts
- :doc:`umat_tutorial` — writing a dedicated UMAT by hand instead
- ``examples/mechanical/MODUL.py`` — runnable gallery example
