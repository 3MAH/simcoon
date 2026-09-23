====================================
Constitutive model (UMAT) catalog
====================================

Every constitutive law is selected by its 5-character ``umat_name``. Since the
modular UMAT framework, four kinds of implementation coexist behind those
names — **the calling convention is identical for all of them** (same
``umat_name``, same props, same solver/FEA usage):

- **modular (native)**: the composable ``MODUL`` engine, configured from a
  props stream (see :mod:`simcoon.modular` for the Python builder).
- **modular (adapter)**: a legacy name whose dedicated kernel was removed
  after its equivalence with a ``MODUL`` configuration was proven against the
  retained reference kernels — bit-identical for the elastic, power-law and
  Hill families; within 2e-3 relative for the Voce/Chaboche families, whose
  legacy incremental ``Hp`` update differs from the modular closed form
  (tests in ``test_modular.py`` and ``Treference_umats``); a translator maps
  the legacy props to the modular configuration at each call.
- **legacy (kept)**: a dedicated, self-contained implementation kept either
  for pedagogy (readable single-file reference of the CCP return mapping) or
  because no modular equivalent exists.
- **Python (external)**: the ``PYEXT`` name serves a law written in Python
  (:class:`simcoon.PythonUMAT`) through a registered callback; see :doc:`python_umat`.

Small-strain mechanical models
==============================

.. list-table::
   :header-rows: 1
   :widths: 8 26 12 40 14

   * - Name
     - Physics
     - Engine
     - Props (in order)
     - Notes
   * - ELISO
     - Isotropic elasticity
     - modular (adapter)
     - E, nu, alpha
     -
   * - ELIST
     - Transversely isotropic elasticity
     - modular (adapter)
     - axis, EL, ET, nuTL, nuTT, GLT, alpha_L, alpha_T
     -
   * - ELORT
     - Orthotropic elasticity
     - modular (adapter)
     - E1, E2, E3, nu12, nu13, nu23, G12, G13, G23, alpha1, alpha2, alpha3
     -
   * - EPICP
     - Von Mises + power-law isotropic hardening
     - legacy (kept)
     - E, nu, alpha, sigmaY, k, m
     - Pedagogical reference of the CCP return mapping
   * - EPKCP
     - Von Mises + power-law isotropic + Prager kinematic
     - modular (adapter)
     - E, nu, alpha, sigmaY, k, m, kX
     - Legacy Prager writes X = kX·a; the modular twin uses X = (2/3)C·a,
       i.e. C = 1.5·kX (handled by the adapter)
   * - EPCHA
     - Von Mises + Voce + 2x Armstrong-Frederick
     - legacy (kept)
     - E, nu, alpha, sigmaY, Q, b, C1, D1, C2, D2
     - Pedagogical reference of kinematic hardening in CCP
   * - EPHIL / EPTRI
     - Hill yield + power-law isotropic hardening
     - modular (adapter)
     - E, nu, alpha, sigmaY, k, m, F, G, H, L, M, N
     - Bit-identical to the modular twin (machine precision)
   * - EPHAC
     - Cubic elasticity + Hill + Voce + 2x AF
     - modular (adapter)
     - E, nu, G, alpha, sigmaY, Q, b, C1, D1, C2, D2, F, G, H, L, M, N
     - statev columns beyond the modular layout are unused (see below)
   * - EPANI
     - Cubic elasticity + 9-parameter anisotropic yield + Voce + 2x AF
     - modular (adapter)
     - E, nu, G, alpha, sigmaY, Q, b, C1, D1, C2, D2, P11, P22, P33, P12,
       P13, P23, P44, P55, P66
     - P must be an admissible (convex) quadratic form: symmetric with zero
       row sums on the normal block; an indefinite P yields sqrt(<0) = NaN
   * - EPDFA
     - Cubic elasticity + Deshpande-Fleck-Ashby yield + Voce + 2x AF
     - modular (adapter)
     - E, nu, G, alpha, sigmaY, Q, b, C1, D1, C2, D2, F, G, H, L, M, N, K
     -
   * - EPCHG
     - Cubic elasticity + selectable yield + N-term "Voce" + N-term Chaboche
     - modular (adapter)
     - E, nu, G, alpha, sigmaY, N_iso, N_kin, criteria(0=Mises, 1=Hill,
       2=DFA, 3=anisotropic), (Q_i, b_i) x N_iso, (C_i, D_i) x N_kin,
       criterion parameters
     - The legacy N-term isotropic hardening couples all terms through a
       single Hp (dHp/dp = sum b_i (Q_i - Hp)): mathematically ONE effective
       Voce with b_eff = sum(b_i), Q_eff = sum(b_i Q_i)/sum(b_i) — not the
       standard combined-Voce sum. The adapter maps accordingly.
   * - EPHIN
     - N Hill yield surfaces, each with power-law isotropic hardening
     - modular (adapter)
     - E, nu, alpha, N, then per surface: sigmaY, k, m, F, G, H, L, M, N
     - The removed legacy kernel was defective for N >= 2 (NaN even for
       identical or inactive second surfaces); the modular engine handles
       multiple surfaces correctly, so N >= 2 is now functional.
   * - ZENER
     - Generalized KELVIN chain, 1 branch (standard solid)
     - legacy (kept)
     - E0, nu0, alpha, E1, nu1, etaB1, etaS1
     - No modular equivalent (Kelvin branches in series; the modular
       viscoelasticity is a generalized Maxwell/Prony model)
   * - ZENNK
     - Generalized KELVIN chain, N branches
     - legacy (kept)
     - E0, nu0, alpha, N, then per branch: E_i, nu_i, etaB_i, etaS_i
     - Same rheology note as ZENER — NOT equivalent to PRONK despite the
       identical props layout (measured 86% response difference)
   * - PRONK
     - Generalized Maxwell (Prony series), N branches
     - legacy (kept)
     - E0, nu0, alpha, N, then per branch: E_i, nu_i, etaB_i, etaS_i
     - Pedagogical reference; the modular Viscoelasticity mechanism is its
       proven twin (< 0.1%)
   * - LLDM0
     - Ductile damage (Lemaitre-Ladeveze-Dufailly)
     - legacy (kept)
     - see header
     - Modular equivalence not yet established (audit pending)
   * - MODUL
     - Composable modular UMAT (elasticity + N mechanisms)
     - modular (native)
     - self-describing stream — build it with
       :class:`simcoon.modular.ModularMaterial`
     - Also available under finite strain (NLGEOM control types 2-6), where
       the composition is a Hencky hyperelastic law on the logarithmic
       strain; requires ``corate_type = 3`` (log_R) — any other corate is
       rejected with a ``RuntimeError`` (hyper/hypo consistency)

Shape memory alloys, finite strain, multiscale, plugins
========================================================

Unchanged dedicated implementations (out of the modular scope):

- **SMA**: SMADI/SMADC/SMAAI/SMAAC (unified),
  SMRDI/SMRDC/SMRAI/SMRAC (unified with reorientation), SMAMO/SMAMC (monocrystal).
- **Finite strain**: HYPOO (hypoelastic orthotropic), SNTVE (Saint-Venant),
  NEOHI/NEOHC (Neo-Hookean), MOORI, YEOHH, ISHAH, GETHH, SWANH, HOLZA
  (invariant-based hyperelasticity); OGDEN (isochoric principal
  stretches, props = ``N, kappa, mu_1, alpha_1, ...``). The compressible
  ones take one optional trailing prop selecting the volumetric term
  :math:`U(J)`: absent or 0 for :math:`\kappa (J \ln J - J + 1)`, 1 for
  :math:`\frac{\kappa}{2} (J - 1)^2` (NEOHI has the latter form built in,
  with :math:`\kappa = 2 / D_1`).

  HOLZA is the Gasser-Ogden-Holzapfel model, the only anisotropic one of the
  set: an isotropic neo-Hookean matrix reinforced by :math:`n` families of
  dispersed fibres,

  .. math::

     W = C_{10} \left(\bar{I}_1 - 3\right)
       + \sum_i \frac{k_1}{2 k_2}
         \left[ \exp\left(k_2 \left(\bar{I}^{*}_{4,i} - 1\right)^2\right) - 1 \right]
       + U(J),

  with :math:`\bar{I}^{*}_{4,i} = \kappa_d \bar{I}_1 + (1 - 3 \kappa_d) \bar{I}_{4,i}`
  and :math:`\bar{I}_{4,i} = \mathbf{a}_{0,i} \cdot \bar{\mathbf{C}} \, \mathbf{a}_{0,i}`.
  The dispersion :math:`\kappa_d \in [0, 1/3]` interpolates between perfectly
  aligned fibres (:math:`\kappa_d = 0`, the Holzapfel-Gasser-Ogden 2000 model)
  and an isotropic distribution (:math:`\kappa_d = 1/3`). The fibre term is inactive
  wherever :math:`\bar{I}^{*}_{4,i} < 1`, which is the switch of the original
  Gasser-Ogden-Holzapfel formulation. Note it is a condition on the *generalized*
  invariant, not on fibre compression: once :math:`\kappa_d > 0` the term picks up a
  :math:`\kappa_d \bar{I}_1` contribution, so a fibre with :math:`\bar{I}_{4,i} < 1` can
  still be active (at :math:`\kappa_d = 0.2` a fibre shortened to 0.39 of its length
  gives :math:`\bar{I}^{*}_4 = 1.12`).

  props = ``C10, k1, k2, kappa_d, n_fam, a0x_1, a0y_1, a0z_1, ..., kappa``,
  where each :math:`\mathbf{a}_{0,i}` is a unit direction **in the local
  material frame** (the solver's material orientation places it globally, as
  for ELIST/ELORT). From Python the directions are given as a
  :class:`simcoon.Rotation` applied to :math:`\mathbf{e}_1`, one entry per
  family, which keeps Euler angles and their gimbal lock off the path to the
  kernel::

      sim.modular.HolzapfelElasticity(
          C10=0.0354, k1=0.0107, k2=7.48, kappa_d=0.0,
          fibres=sim.Rotation.from_euler('zxz', [[0, 0, 40], [0, 0, -40]],
                                         degrees=True),
          kappa=1000.)

  The same potential is available as a MODUL elasticity block, and may be composed
  with any mechanism. Composed with **damage** -- anisotropic tissue with softening --
  it is exact: damage subtracts no inelastic strain (it scales the stiffness instead),
  so the elastic stretch is still the total one, and its driving force is built from
  the current anisotropic tangent.

  .. warning::

     Composed with a mechanism that *does* subtract an inelastic strain
     (**plasticity**, **viscoelasticity**), the fibre convection is
     **approximate**. The block carries :math:`\mathbf{a}_{0,i}` from the reference
     configuration and pushes it forward with the elastic stretch, so the inelastic
     strain does not reorient the fibres. The composition is well posed and
     converges -- the return mapping is handed the anisotropic tangent and its
     consistency condition holds exactly -- but the response is only as good as that
     assumption: exact while the inelastic strain is small or leaves the fibre
     directions fixed, degrading as it reorients them. Representing the convection
     exactly would require the convected directions as state variables, which the
     additive corotational kinematics of the modular UMAT cannot express (there is
     no plastic deformation gradient to convect with).

     A second, separate caveat applies to **viscoelasticity** only: every Prony
     branch is built as an *isotropic* :math:`\mathbf{L}_i(E_i, \nu_i)`, so the viscous
     response carries none of the fibre anisotropy while the equilibrium response does.
     That is a property of the viscoelastic mechanism's :math:`(E_i, \nu_i)`
     parameterization rather than of HOLZA -- a Prony branch cannot follow an ELORT or
     ELIST block's symmetry either.

     Neither is rejected at run time; both are modelling choices left to the user.

  A single scalar damage variable, finally, degrades matrix and fibres at the same
  rate. The Holzapfel damage literature instead carries separate variables -- one on
  the isotropic term and one per fibre family -- since collagen and ground substance
  damage very differently. simcoon's damage mechanism is a single scalar, so the
  composition models uniform softening, not anisotropic damage.

  MUSCL is **activated skeletal muscle**, and the only law in simcoon driven by
  something that is neither strain nor temperature. A 5-parameter Mooney-Rivlin ground
  matrix whose stiffness rises with the activation :math:`a`, plus an along-fibre
  force law:

  .. math::

     W = s(a) \Big[ C_{10}(\bar{I}_1 - 3) + C_{01}(\bar{I}_2 - 3)
       + C_{20}(\bar{I}_1 - 3)^2 + C_{11}(\bar{I}_1 - 3)(\bar{I}_2 - 3)
       + C_{02}(\bar{I}_2 - 3)^2 \Big]
       + \sum_i \Phi(\bar{\lambda}_i; a) + s(a)\, U(J),

  with :math:`s(a) = 1 + (s_{max} - 1)a` and
  :math:`\bar{\lambda} = \sqrt{\bar{I}^{*}_{4}}` the isochoric fibre stretch. The
  matrix is Nazari et al.'s Eq. (1) and is a superset of neo-Hookean, Mooney-Rivlin
  and second-order Yeoh; setting only :math:`C_{10}` and :math:`C_{20}` with no fibre
  term reproduces those laws bit for bit.

  One potential covers a published family, selected by the leading prop
  ``fibre_law``, which changes the *interpretation* of the fibre parameters and never
  their number -- the same convention idiom the linear elasticity blocks use. With
  :math:`\hat{\lambda} = \bar{\lambda}/\lambda_{opt}` and
  :math:`f_d = \partial W / \partial \bar{\lambda}`:

  .. list-table::
     :header-rows: 1
     :widths: 12 46 42

     * - ``fibre_law``
       - :math:`f_d(\bar{\lambda})`
       - reproduces
     * - 0 ``NONE``
       - no fibre term
       - Nazari et al. (2010, 2011): activation raises the matrix stiffness only, the
         contractile force coming from elsewhere (1-D cable elements in their model)
     * - 1 ``SIMPLE``
       - :math:`a\,\sigma_{max}`
       - ArtiSynth ``SimpleForceMuscle``
     * - 2 ``GENERIC``
       - :math:`a\,\sigma_{max} + P_1\left(e^{P_2(\bar{\lambda}-1)}-1\right)/\bar{\lambda}`
       - ArtiSynth ``GenericMuscle``
     * - 3 ``BLEMKER``
       - :math:`\sigma_{max}\left(a f_a(\hat{\lambda}) + f_p(\hat{\lambda})\right)/\lambda_{opt}`
       - Blemker et al. (2005); ArtiSynth ``BlemkerMuscle``; FEBio

  :math:`f_p` is exponential between :math:`\lambda_{opt}` and :math:`\lambda^{*}`
  then linear, and :math:`f_a` is the three-piece Hill parabola, zero outside
  :math:`\hat{\lambda} \in [0.4, 1.6]`. As published both are only :math:`C^0` at
  :math:`\hat{\lambda} = 1, 0.6, 1.4`; simcoon continues them across those junctions so
  the tangent is not decided by round-off, exactly beyond a :math:`10^{-3}`
  neighbourhood.

  props = ``fibre_law, C10, C01, C20, C11, C02, s_max, act, sigma_max, lambda_opt,
  lambda_star, P1, P2, zero_below_opt, kappa_d, n_fam, a0x_1, a0y_1, a0z_1, ..., kappa``.
  From Python use the named constructors, which carry each source's published values::

      sim.modular.MuscleElasticity.blemker(
          fibres=sim.Rotation.from_euler('zxz', [[0, 0, 0]], degrees=True),
          kappa=1000., activation=0.5)
      sim.modular.MuscleElasticity.nazari()      # s_max = 10, no fibre term

  .. warning::

     ``P1`` means **different things** in ``GENERIC`` and ``BLEMKER``: a stress in the
     former, dimensionless in the latter. ArtiSynth ships ``0.05`` as the default for
     both, which in ``GenericMuscle`` at :math:`\sigma_{max} = 30` kPa gives 0.65 Pa of
     passive fibre stress at :math:`\bar{\lambda} = 1.4` against 30 kPa active -- five
     orders of magnitude apart. Likewise ``zero_below_opt`` is true in
     ``BlemkerMuscle`` and **false** in ``GenericMuscle``, which therefore produces a
     negative passive fibre stress in compression. Both are explicit props here rather
     than hidden constants; do not carry a number across from one law to the other.

  .. note::

     **The activation is a driven input.** It is a prop, and ``umat_modular`` re-parses
     its props on every call, so a caller supplies a time-varying, per-integration-point
     activation by rewriting one row of a ``(nprops, n_points)`` props array between
     increments -- nothing in the UMAT interface changes. ``props[0]`` of a
     :class:`~simcoon.modular.ModularMaterial` is the elasticity type, so the absolute
     index is :attr:`~simcoon.modular.ModularMaterial.activation_index`; read it from
     there rather than hard-coding it. The solver's own ``props`` are fixed for a run,
     so in-house a ramp is built as one block per activation level -- which is how
     Nazari's static analyses were run.

     Overlapping muscles are the driver's business, not the law's: Nazari takes the
     **maximum** activation over the muscles sharing an element, Buchaillard et al. the
     **sum**. One UMAT call sees one activation.

  .. warning::

     ``s_max > 1`` together with an active fibre law is **rejected**. Scaling the
     passive stiffness with activation is Nazari's surrogate for the transverse
     stress-stiffening a real contractile fibre produces -- his own later 3-D muscle
     element drops it for exactly that reason -- so composing the two double-counts the
     same physics. Use ``s_max = 1`` with a fibre law, or ``fibre_law = 0`` with
     ``s_max > 1``.

  .. note::

     **An activated muscle is pre-stressed at zero strain**, unlike every other law in
     simcoon: at the optimal length the Hill curve peaks at 1, so the fibre carries
     exactly :math:`a\,\sigma_{max}` of deviatoric stress with no deformation at all.
     The ground-state stiffness the UMAT reports as ``L`` is the *matrix* one and does
     not include the active fibre's contribution; the composed mechanisms are handed the
     current tangent, which does.

     The law is hyperelastic only at **frozen** activation. The activation is prescribed
     rather than governed, so the material point is thermodynamically open, drawing
     chemical energy the mechanical balance does not see. :math:`W_{m,d}` and
     :math:`W_{m,ir}` stay 0 -- booking the active term as dissipation would make
     :math:`W_{m,d} < 0`, which is forbidden -- but along a path where the activation
     varies :math:`W_{m,r}` is only the mechanical work residue, not the stored energy,
     and it **can go negative**. Over a closed cycle in strain and activation,
     :math:`\oint \boldsymbol{\tau} : \mathrm{d}\boldsymbol{\varepsilon}
     = -\oint (\partial W / \partial a)\, \mathrm{d}a`, the metabolic input.

  The fibre-convection caveats stated above for HOLZA apply verbatim. Blemker's full
  model adds Criscione shear terms in :math:`\bar{I}_5` and cross-derivatives in
  :math:`(\bar{I}_1, \bar{I}_4, \bar{I}_5)`, which the invariant framework's additive
  separability cannot express; those are not part of MUSCL.
- **Multiscale**: MIHEN, MIMTN, MISCN, MIPLN. Their sub-phases are passed in
  memory (``phases=``, see :doc:`python_solver`) and their ``props`` hold only the
  scheme's settings: ``[mp, np]`` for MIHEN, ``[mp, np, n_matrix]`` for MIMTN,
  ``[mp, np, n_matrix]`` or ``[mp, np, n_matrix, start]`` for MISCN, nothing for
  MIPLN (``mp``, ``np``: integration points of the Eshelby integrals; ``n_matrix``:
  index of the matrix phase in the list; ``start``: first guess of the
  self-consistent iteration, 1 Mori-Tanaka (default) or 0 homogeneous strain with
  ``n_matrix < 0``). The length is checked, so the pre-2.0 layout with its leading
  ``nphases`` and file-number slots is refused rather than misread.
- **Plugins**: UMEXT (external dylib), UMABA (Abaqus wrapper).
- **Python**: PYEXT (registered Python law; ``props``/``nstatev`` from the object).

State variable (statev) layout for adapter-served names
========================================================

The modular engine claims the FIRST ``required_nstatev`` slots of the caller's
statev array — always within the legacy allocation, so array sizes never
change. For ELISO/ELIST/ELORT (``T_init``) and EPHIL/EPTRI
(``T_init, p, EP(6)``) the column meaning is identical to the removed kernels.
For the Chaboche-family names and EPKCP the columns re-mean: the layout is
``T_init | p, EP(6) | back-strains a_i(6) ...`` in mechanism registration
order; trailing legacy slots are left untouched. Code that read specific
legacy statev columns (e.g. the stored X_i of EPHAC) must be updated to the
modular layout.

Tangent-operator mode
=====================

All models receive the solver's ``tangent_mode`` (named constants in
``parameter.hpp`` / ``sim.tangent_*``): 0 = none (Lt = elastic L, explicit
integration), 1 = continuum, 2 = algorithmic/Simo-Hughes (**default**),
3 = closest-point (reserved). Pre-2.0 numbering was 0 = continuum,
1 = algorithmic — see :doc:`solver` for the migration note. The
finite-strain hyperelastic models ignore the mode (their tangent is always
the exact one of the hyperelastic law).

Validation and performance
==========================

Each adapter-served name is validated at two levels:

- **Translator correctness**: a pytest equivalence test
  (``simcoon-python-builder/test/test_core/test_modular.py``, the
  ``*_matches_modul`` family) proves the legacy name and the explicit
  ``MODUL`` configuration are bit-identical through the solver.
- **Independent physics**: the removed legacy kernels are retained VERBATIM
  as test-only reference oracles under
  ``test/Libraries/Umat/reference_kernels/`` (compiled only into the
  ``Treference_umats`` gtest, never into ``libsimcoon``, not dispatchable by
  name). Every adapter is driven side by side with its reference kernel on a
  cyclic strain path each test run — machine precision for the
  elastic/power-law families, < 2e-3 for the Voce/Chaboche family (legacy
  incremental vs modular closed-form Voce integration).

A benchmark row per family lives in ``bench/bench_legacy_vs_modular.py``:
results show no measurable adapter overhead, and the modular
engine runs at 0.8-1.2x the speed of the removed kernels on all families.
