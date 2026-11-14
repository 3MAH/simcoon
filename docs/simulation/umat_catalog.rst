====================================
Constitutive model (UMAT) catalog
====================================

Every constitutive law is selected by its ``umat_name`` (five letters, or the ``_closed``
name of a kept dedicated kernel). The names live in one table, the registry of
``umat_registry.cpp`` (``umat_registry()``, one row per name): the kernel that serves the name on
each face of the dispatch (small-strain mechanical, finite-strain mechanical, thermomechanical),
its finite-strain conventions and its capabilities (condensed stress states of ``sim.umat``, the
elastic stiffness of the mean-field schemes). Every route reads it: the solver, ``sim.umat`` /
``sim.umat_T``, the Abaqus and Ansys bridges, the micromechanics. Since the
modular UMAT framework, four kinds of implementation coexist behind those
names — **the calling convention is identical for all of them** (same
``umat_name``, same props, same solver/FEA usage):

- **modular (native)**: the composable ``MODUL`` engine, configured from a
  props stream (see :mod:`simcoon.modular` for the Python builder).
- **modular (adapter)**: a legacy name whose dedicated kernel was removed
  after its equivalence with a ``MODUL`` configuration was proven against the
  retained reference kernels — bit-identical for the elastic, power-law and
  Hill families and, since the modular engine runs the same cutting plane, for
  the Voce/Chaboche families as well (tests in ``test_modular.py`` and
  ``Treference_umats``); a translator maps
  the legacy props to the modular configuration at each call.
- **legacy (kept)**: a dedicated, self-contained implementation kept because
  no modular equivalent exists.
- **modular (adapter); ``NAME_closed`` runs the dedicated kernel**: the name
  runs the modular engine; the dedicated closed-form kernel stays reachable
  under the ``_closed`` suffix (benchmarks, comparisons, the readable
  single-file reference of the CCP return mapping), on every route the plain
  name takes (solver, ``sim.umat`` / ``sim.umat_T``, Abaqus ``CMNAME``).
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
     - modular (adapter); ``EPICP_closed`` runs the dedicated kernel
     - E, nu, alpha, sigmaY, k, m
     - The dedicated kernel (``EPICP_closed``, mechanical and thermomechanical)
       is the pedagogical reference of the CCP return mapping with the
       closed-form radial-return tangent, 1.3x faster than the generic
       assembly of the modular engine; the hardening goes through the modular
       PowerLawHardening (C1 onset blend below p = 1e-6 for m < 1, exact
       above), as EPJCK
   * - EPKCP
     - Von Mises + power-law isotropic + Prager kinematic
     - modular (adapter)
     - E, nu, alpha, sigmaY, k, m, kX
     - The legacy kernel writes :math:`X = k_X a`, the modular twin
       :math:`X = \tfrac{2}{3} C a`, i.e. :math:`C = 1.5 k_X` (handled by the adapter)
   * - EPCHA
     - Von Mises + Voce + 2x Armstrong-Frederick
     - modular (adapter); ``EPCHA_closed`` runs the dedicated kernel
     - E, nu, alpha, sigmaY, Q, b, C1, D1, C2, D2
     - The dedicated kernel (``EPCHA_closed``) is the pedagogical reference of
       kinematic hardening in CCP, the modular composition to round-off on a
       cyclic path (same cutting plane); the modular statev are
       ``[T_init, p, EP, a_1, a_2]`` (20, the back-strains), the back stresses
       and ``H_p`` of the dedicated layout are not stored
   * - EPJCK
     - Von Mises + Johnson-Cook yield stress (power-law hardening x logarithmic
       strain-rate factor x thermal softening)
     - modular (adapter, ``JohnsonCookHardening`` block); ``EPJCK_closed`` runs the
       dedicated kernel; the thermomechanical face keeps the dedicated kernel under
       both names
     - E, nu, alpha, A, B, n, C, edot0, m, T_ref, T_melt
     - Rate-dependent: the plastic strain rate is Dp/Dt, fully implicit in the
       return mapping; the rate factor is clamped at 1 below edot0 and the cutting
       plane starts at the kink (``multiplier_predictor``). 9 statev
       (T_init, p, EP, edot_p), the same layout and meaning under the adapter.
       Thermomechanical twin with (rho, c_p) prepended, dPhi/dT from the thermal
       softening (dedicated kernel: the hardening entropy has no modular hook yet). tangent_mode 0 returns
       the elastic L in the mechanical kernel (explicit integration); the
       thermomechanical twin promotes it to the continuum operator, as EPICP_T.
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
     - Standard linear solid, Poynting-Thomson (Kelvin) form: spring in series with one Kelvin-Voigt unit
     - modular (adapter); ``ZENER_closed`` runs the dedicated kernel
     - E0, nu0, alpha, E1, nu1, etaB1, etaS1
     - The adapter runs the one-branch Maxwell composition with the Kelvin unit
       mapped per isotropic mode (:math:`X_M = X_0^2/(X_0 + X_1)`,
       :math:`\eta_M = \eta_K X_0^2/(X_0 + X_1)^2`, ``simcoon.modular.zener_as_maxwell``):
       same response to the round-off of the two inverses; the statev ``EV`` is then
       the Maxwell branch strain, not the Kelvin viscous strain of the dedicated
       kernel (``ZENER_closed``, mechanical and thermomechanical)
   * - ZENNK
     - Generalized KELVIN chain, N branches
     - legacy (kept)
     - E0, nu0, alpha, N, then per branch: E_i, nu_i, etaB_i, etaS_i
     - Kelvin units in series: NOT equivalent to PRONK despite the identical
       props layout (measured 86% response difference); the Kelvin-to-Maxwell
       map of ZENER has no closed form for N >= 2
   * - PRONK
     - Generalized Maxwell (Prony series), N branches
     - modular (adapter); ``PRONK_closed`` runs the dedicated kernel
     - E0, nu0, alpha, N, then per branch: E_i, nu_i, etaB_i, etaS_i
     - The dedicated kernel (``PRONK_closed``, mechanical and thermomechanical)
       is the pedagogical reference; the modular Viscoelasticity mechanism
       takes the same closed-form step and is identical to it (round-off); the
       modular statev are ``[T_init, (v_i, EV_i) x N]`` (1 + 7N), the total
       viscous strain of the dedicated layout is not stored
   * - LLDM0
     - Ladevèze-Le Dantec damage of the unidirectional ply (d22, d12), plasticity in effective stress
     - legacy (kept)
     - see header (18 props, 10 statev)
     - Modular equivalence not yet established (audit pending)
   * - MODUL
     - Composable modular UMAT (elasticity + N mechanisms)
     - modular (native)
     - self-describing stream — build it with
       :class:`simcoon.modular.ModularMaterial`
     - Also available under finite strain (NLGEOM control types 2-6), where
       the composition is a Hencky hyperelastic law on the logarithmic
       strain; requires ``corate_type = 3`` (log_R) — any other corate is
       rejected with a ``RuntimeError`` (hyper/hypo consistency). Thermomechanical
       twin with (rho, c_p) prepended: linear elasticity blocks, every
       mechanism (:doc:`modular_python`)

Shape memory alloys, finite strain, multiscale, plugins
========================================================

Unchanged dedicated implementations (out of the modular scope):

- **SMA**: SMADI/SMADC/SMAAI/SMAAC (unified),
  SMRDI/SMRDC/SMRAI/SMRAC (unified with reorientation), SMAMO/SMAMC (monocrystal).
- **Finite strain**: HYPOO (hypoelastic orthotropic), SNTVE (Saint-Venant),
  OGDEN (isochoric principal stretches, props = ``N, kappa, mu_1, alpha_1, ...``)
  are kernels built on :math:`\mathbf{F}`. The invariant-based hyperelastic names
  NEOHI/NEOHC (Neo-Hookean), MOORI, YEOHH, ISHAH, GETHH, SWANH and HOLZA run the
  ``HYPER_INVARIANTS`` block of the modular engine through a props adapter (the
  same potential, the same props): a potential of the logarithmic strain, served
  under corate 3 (log_R) and the total formulation only, as ``MODUL``; ``statev``
  is ``[T_init]``. ``NEOHC_closed`` and ``NEOHI_closed`` keep the kernels built on
  :math:`\mathbf{F}`, which serve every corate and the native material tangent
  (comparisons, ``test_finite_tangents``). The compressible potentials take one
  optional trailing prop selecting the volumetric term :math:`U(J)`: absent or 0
  for :math:`\kappa (J \ln J - J + 1)`, 1 for :math:`\frac{\kappa}{2} (J - 1)^2`
  (NEOHI is the Neo-Hookean potential with the latter form, :math:`\mu = E/2(1+\nu)`,
  :math:`\kappa = E/3(1-2\nu)`, i.e. :math:`C_{10} = \mu/2`, :math:`1/D_1 = \kappa/2`).

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

  **Parameters in the papers and in other codes.** The matrix coefficient is
  :math:`C_{10}` here and in Abaqus, the shear modulus :math:`c = 2 C_{10}` (written
  :math:`\mu` in 2000) in the papers, FEBio and matadi. The volumetric terms differ, but all
  are fixed by the same ground-state bulk modulus :math:`K = U''(1)`; the papers treat the
  material as incompressible.

  .. list-table::
     :header-rows: 1
     :widths: 18 20 22 20 20

     * - quantity
       - Gasser et al. 2006 (Holzapfel et al. 2000)
       - simcoon ``HolzapfelElasticity`` / ``HOLZA``
       - Abaqus ``*ANISOTROPIC HYPERELASTIC, HOLZAPFEL``
       - FEBio ``Holzapfel-Gasser-Ogden``
     * - matrix
       - :math:`c` (:math:`\mu`)
       - ``C10`` :math:`= c/2`
       - :math:`C_{10} = c/2`
       - ``c``
     * - fibre stiffness
       - :math:`k_1`
       - ``k1``
       - :math:`k_1`
       - ``k1``
     * - fibre exponent
       - :math:`k_2`
       - ``k2``
       - :math:`k_2`
       - ``k2``
     * - dispersion
       - :math:`\kappa` (0 in 2000)
       - ``kappa_d``
       - :math:`\kappa`
       - ``kappa``
     * - fibre families
       - two, at :math:`\pm\gamma` (:math:`\pm\beta` in 2000) from the circumferential
         direction, in the tangent plane
       - ``fibres``: one rotation of :math:`\mathbf{e}_1` per family in the local frame,
         ``Rotation.from_euler("zxz", [[0, 0, g], [0, 0, -g]], degrees=True)``;
         props ``n_fam`` and the unit directions
       - two local directions (``*ORIENTATION``, ``LOCAL DIRECTIONS``)
       - ``gamma``, angle in a plane of the local frame
     * - bulk modulus
       - incompressible
       - ``kappa`` :math:`= K`; ``volumetric="log"``:
         :math:`K (J \ln J - J + 1)`, ``"quadratic"``: :math:`K (J-1)^2/2`
       - :math:`D = 2/K`, :math:`U = \frac{1}{D}\left(\frac{J^2 - 1}{2} - \ln J\right)`
       - ``k`` (FEBio's own :math:`U(J)`, see its manual)

  The 2006 paper rewrites the 2000 model around a generalized structure tensor; the
  notations map as follows (the 2006 model with :math:`\kappa = 0` is the 2000 one).

  .. list-table::
     :header-rows: 1
     :widths: 18 28 32 22

     * -
       - Holzapfel, Gasser, Ogden 2000
       - Gasser, Ogden, Holzapfel 2006
       - simcoon
     * - structure tensor
       - :math:`\mathbf{A}_i = \mathbf{a}_{0i} \otimes \mathbf{a}_{0i}`
       - :math:`\mathbf{H}_i = \kappa \mathbf{I} + (1 - 3\kappa)\,\mathbf{a}_{0i} \otimes \mathbf{a}_{0i}`
       - built from ``kappa_d`` and ``fibres``
     * - fibre strain
       - :math:`\bar{I}_4 - 1`, :math:`\bar{I}_6 - 1` (:math:`\bar{I}_{4,6} = \bar{\mathbf{C}} : \mathbf{A}_{1,2}`)
       - :math:`\bar{E}_i = \mathbf{H}_i : \bar{\mathbf{C}} - 1
         = \kappa (\bar{I}_1 - 3) + (1 - 3\kappa)(\bar{I}_{4,i} - 1)`
       - the same quantity, written :math:`\bar{I}^{*}_{4,i} - 1`
         (:math:`\bar{I}^{*}_{4,i} = \kappa_d \bar{I}_1 + (1 - 3\kappa_d) \bar{I}_{4,i} = \mathbf{H}_i : \bar{\mathbf{C}}`)
     * - fibre energy
       - :math:`\frac{k_1}{2 k_2} \sum_{i = 4, 6} \left[ e^{k_2 (\bar{I}_i - 1)^2} - 1 \right]`
       - :math:`\bar{\Psi}_{f i} = \frac{k_1}{2 k_2} \left[ e^{k_2 \bar{E}_i^2} - 1 \right]`
       - same, summed over the families
     * - dispersion
       - none (aligned fibres)
       - :math:`\kappa = \frac{1}{4} \int_0^\pi \rho(\Theta) \sin^3\Theta \, d\Theta`, with
         :math:`\rho` a von Mises density of concentration :math:`b`
       - ``kappa_d`` (``kappa`` is the bulk modulus in simcoon)
     * - fibre angle
       - :math:`\beta`
       - :math:`\gamma`
       - ``fibres``
     * - wall thickness (tube)
       - :math:`H`
       - :math:`H` (also the letter of the structure tensor)
       - not a material parameter

  matadi (``holzapfel_gasser_ogden``) takes ``c``, ``k1``, ``k2``, ``kappa`` and the fibre
  ``angle`` about an ``axis``; ``axis=2, angle=g`` is the plane and the convention of the
  rotation above (checked to 1.6e-14 on the stress against simcoon, quadratic :math:`U(J)`).
  The example :ref:`sphx_glr_examples_hyperelasticity_GOH2006_tube_inflation.py` reproduces figure 7 of Gasser et al. (2006) with the
  parameters of the iliac adventitia, :math:`c = 7.64` kPa, :math:`k_1 = 996.6` kPa,
  :math:`k_2 = 524.6`, :math:`\gamma = 49.98^\circ`, :math:`\kappa = 0.226`.

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
  the kinematics are exact: damage subtracts no inelastic strain (it scales the
  stiffness instead), so the elastic stretch is the total one and the fibres follow it.
  Its driving force is the thermodynamic force :math:`Y = -\partial\Psi/\partial D = \Psi_0`, the
  undamaged energy, evaluated as
  :math:`\tfrac12\,\boldsymbol{\tau}_{eff} : \mathbf{M}_t : \boldsymbol{\tau}_{eff}` on the effective
  stress :math:`\boldsymbol{\tau}/(1-D)`: exact for a linear block, an approximation of
  :math:`\Psi_0` through the tangent compliance for this nonlinear potential.

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
- **Multiscale**: MIHEN, MIMTN, MISCN, MIPLN. Their sub-phases are passed in
  memory (``phases=``, see :doc:`python_solver`) and their ``props`` hold only the
  scheme's settings: ``[mp, np]`` for MIHEN, ``[mp, np, n_matrix]`` for MIMTN,
  ``[mp, np, n_matrix]`` or ``[mp, np, n_matrix, start]`` for MISCN, nothing for
  MIPLN (``mp``, ``np``: integration points of the Eshelby integrals; ``n_matrix``:
  index of the matrix phase in the list; ``start``: first guess of the
  self-consistent iteration, 1 Mori-Tanaka (default) or 0 homogeneous strain with
  ``n_matrix < 0``). The length and the values are checked, so the pre-2.0 layout with
  its leading ``nphases`` and file-number slots is refused rather than misread.

  Optional trailing slots set numerical options; all default to 0, the historical
  behaviour (giving a slot gives the ones before it). The same settings by name:
  ``MoriTanaka``, ``SelfConsistent``, ``Homogeneous`` and ``PeriodicLayers`` of
  ``simcoon.solver.micromechanics`` build these lists (see :doc:`python_solver`):
  ``[mp, np, n_matrix, closed_form, rule]`` for MIMTN and
  ``[mp, np, n_matrix, start, closed_form, rule, anderson, rtol]`` for MISCN.

  * ``closed_form = 1``: the Eshelby tensor is taken in closed form (sphere, prolate or
    oblate spheroid of axis 1) when the reference medium is isotropic, and integrated
    otherwise. Exact, and no integration: ``L_eff`` of a self-consistent composite of
    spheres is 7 to 10 times faster (20 times with ``anderson = 2``).
  * ``rule = 1``: half-sphere quadrature, ``(mp+1)/2`` Gauss-Legendre points in
    ``cos(theta)`` and ``np`` trapezoidal points in ``phi`` (the integrand is centrally
    symmetric and periodic in ``phi``): half the evaluations of the default
    Gauss-Legendre rule for the same accuracy or better.
  * ``anderson = m`` (MISCN, ``L_eff`` only): Anderson acceleration of depth ``m`` of the
    self-consistent fixed point (``m = 2`` is a good choice: 65 iterations become about
    11 at a 50 % concentration). The solver does one self-consistent sweep per
    localization iteration and ignores it.
  * ``rtol > 0`` (MISCN, ``L_eff`` only): relative stopping criterion
    ``||L - L_prev|| <= rtol ||L||`` of that fixed point instead of the absolute
    ``1e-6`` (in the stiffness unit), so the iteration count no longer depends on the unit.
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

.. _stress-measure-tangent-rate:

Stress measure and tangent rate
===============================

Every **native** simcoon kernel is Kirchhoff-native. In the logarithmic framework an
isotropic stored energy per *reference* volume gives
:math:`\boldsymbol{\tau} = \partial W / \partial \ln \mathbf{V}` (for an anisotropic one the
identity needs :math:`\boldsymbol{\tau}` coaxial with :math:`\mathbf{V}`), so :math:`\boldsymbol{\tau}`
is the route stress, and the Cauchy stress is the derived output
:math:`\boldsymbol{\sigma} = \boldsymbol{\tau}/J`, formed only at the boundaries -- the
solver's output sinks and the Python wrapper. Nothing on the route carries Cauchy.

**The stress measure** a kernel writes into ``sigma`` is declared, one entry per kernel, in
the registry (``umat_entry::convention``, ``umat_registry.cpp``; read through
``output_convention_of``). A name without a finite-strain kernel is rejected by the finite
route rather than given a default -- a missing declaration is an error of exactly
:math:`J`, and invisible near :math:`J = 1`:

- ``kirchhoff`` -- every native kernel.
- ``cauchy`` -- only the host plugins ``UMABA`` / ``UMANS``, whose contract belongs to the host
  code (Cauchy stress, a Jaumann or Green-Naghdi Jacobian). ``select_umat_M_finite`` hands them
  :math:`\boldsymbol{\sigma}_n = \boldsymbol{\tau}_n / J_n`, converts their stress to
  :math:`\boldsymbol{\tau}` on the way out and their Jacobian to the box of the solver's corate
  (``UMEXT`` is a Kirchhoff box, as ``PYEXT``; see :doc:`../external`).

``HYPOO`` is native too: it integrates the corotational *Kirchhoff* rate,
:math:`\boldsymbol{\tau}_{n+1} = \boldsymbol{\tau}_n + \mathbf{L} : \Delta\boldsymbol{\varepsilon}^{el}`,
and is the rate counterpart of ``ELORT`` (total form
:math:`\boldsymbol{\tau} = \mathbf{L} : \boldsymbol{\varepsilon}^{el}`). With the same
:math:`\mathbf{L}`, isotropic or anisotropic, the two are identical on any path for corates
0 to 3: with the material axes following the body, the transported strain and stress satisfy
the same recursion (a ~74 % orthotropic shear gap measured before 2.1 came from lab-fixed
axes, not from the rate form). They differ under corate 4 even for an isotropic
:math:`\mathbf{L}`: ``ELORT`` is then the total Almansi law
:math:`\boldsymbol{\tau} = \lambda\,\mathrm{tr}(\mathbf{e}_A)\mathbf{I} + 2\mu\,\mathbf{e}_A`, ``HYPOO``
integrates the Oldroyd rate of :math:`\boldsymbol{\tau}`, whose isochoric solution is
:math:`\boldsymbol{\tau} = \mu(\mathbf{b} - \mathbf{I})`. In simple shear both give
:math:`\tau_{12} = \mu\gamma`; the normal stresses differ, :math:`\tau_{11} = \mu\gamma^2`,
:math:`\tau_{22} = \tau_{33} = 0` for ``HYPOO`` against :math:`\tau_{11} = \tau_{33} = -\lambda\gamma^2/2`,
:math:`\tau_{22} = -(\lambda/2 + \mu)\gamma^2` for ``ELORT``. Under corate 5 they differ for an
anisotropic :math:`\mathbf{L}` only (about 6 % measured in orthotropic shear): the similarity
transport by a rotation-free stretch does not commute with it. For an
anisotropic :math:`\mathbf{L}` the common law :math:`\boldsymbol{\tau} = \mathbf{L}_R :
\ln\mathbf{V}` is Cauchy-elastic, not hyperelastic.

**Material frame.** On the finite route the kernels fed the logarithmic strain (ELISO, ELIST,
ELORT, HYPOO, the plasticity names, ``MODUL``, ``PYEXT``) run in a frame that follows the
material: the accumulated corate rotation for corates 0 to 3, the polar rotation of
:math:`\mathbf{F}` for 4 and 5. Their anisotropy axes (orthotropic and transversely isotropic
stiffness, Hill-type criteria, HOLZA fibres, in ``MODUL`` and through the hyperelastic names)
therefore rotate with the body, as with Abaqus ``*ORIENTATION`` under NLGEOM, and their tensorial
``Statev`` are material-frame components. The kernels built from :math:`\mathbf{F}` (SNTVE,
``NEOHC_closed``, ``NEOHI_closed``, OGDEN) are objective by construction and keep the lab frame.

**The tangent rate is not declared, it is deduced from the solver's ``corate_type``.** Every
kernel receives the solver's ``corate_type`` and must return :math:`\mathbf{L}_t` expressed in
it. A kernel that builds its tangent from :math:`\mathbf{F}` -- the finite hyperelastic family --
converts the spatial (Lie/Oldroyd) closed form of the potential in one step, with
``Dtau_LieDD_2_DtauDe_corate``; a kernel handed the solver's already-corotated strain
increment is in that rate for free. Per corate: 0 Jaumann and 1 Green-Naghdi are spin/rate
corrections, 2 (XBM) and 3 (log_R) share the exact spectral map, and 4 (Truesdell) is the
convected box, which *is* the spatial (Lie) tangent -- an identity: the Kirchhoff stress is
transported upper-convected and the strain lower-convected, so the strain is the Almansi strain
:math:`\mathbf{e}_A = \tfrac12(\mathbf{I} - \mathbf{b}^{-1})` (the ``Strain`` output; see
:doc:`output`). Logarithmic control (control type 3) still prescribes :math:`\ln\mathbf{V}`
under corate 4: the solver rebuilds it from the stored strain,
:math:`\ln\mathbf{V} = -\tfrac12\ln(\mathbf{I} - 2\mathbf{e}_A)`.
5 (log_F) is the chain rule through its own increment :math:`\mathbf{D}_e = \mathbb{A}^F :
\mathbf{D}\,\Delta t`, :math:`\mathbb{C}^J : (\mathbb{A}^F)^{-1}`, equal to the log box for
isotropic laws (both finite-difference verified).

.. note::

   The solver transports the total strain, the start stress and the tensorial internal
   variables with the variance each corate requires: rotation for 0-3 (in the material frame the
   frame-relative increment is the identity, so a kernel is called with
   :math:`\Delta\mathbf{R} = \mathbf{I}`); for 4, lower-convected strain-like and upper-convected
   stress-like quantities; for 5, similarity. For 4 and 5 the dispatcher applies it to the
   internal variables each kernel declares in the registry (``umat_registry.cpp``: EPICP,
   EPJCK, EPCHA, the elastic and hyperelastic kernels); a kernel that has not declared them (the
   plasticity names served by the modular engine, ``PYEXT``) is refused under corates 4 and 5,
   and ``MODUL`` under any corate but 3 (it runs under the total formulation). Under the total formulation nothing is transported
   (see :ref:`the total formulation <strain-space-route>` below). ``EPCHA`` is also refused under corate 4: it stores
   :math:`\mathbf{X}_i = \tfrac23 C_i\,\mathbf{a}_i`, and the Truesdell rate convects the
   back-strain and the back-stress differently.

**Mechanical work.** ``Wm`` is the work per reference volume,
:math:`\int \mathbf{P} : \mathrm{d}\mathbf{F}`. A kernel accumulates
:math:`\tfrac12(\hat{\boldsymbol{\tau}}_n + \boldsymbol{\tau}_{n+1}) : \Delta\mathbf{e}` on the strain
increment it is handed, :math:`\hat{\boldsymbol{\tau}}_n` being the start stress transported to the
end configuration. Per corate:

- 0, 1 (Jaumann, Green-Naghdi): :math:`\Delta\mathbf{e} = \mathbf{D}\,\Delta t`, the kernel work
  converges to the true work at second order.
- 2 (XBM): the rate is conjugate, :math:`\mathbf{D} = (\ln\mathbf{V})^{\circ\log}`, so there is no
  continuum gap.
- 3, 5 (log_R, log_F): :math:`\boldsymbol{\tau} : \dot{\mathbf{e}} \neq \boldsymbol{\tau} : \mathbf{D}` as
  soon as :math:`\boldsymbol{\tau}` is not coaxial with :math:`\mathbf{V}` (anisotropy, plasticity).
- 4 (Truesdell): the kernel work is exactly the trapezoid
  :math:`\tfrac12(\mathbf{S}_n + \mathbf{S}_{n+1}) : \Delta\mathbf{E}`.

Under 2, 3 and 5 the solver replaces the kernel work by the midpoint stress power with both
stresses in the lab frame, :math:`\tfrac12(\boldsymbol{\tau}_n + \boldsymbol{\tau}_{n+1}) :
\mathrm{sym}\big(2(\mathbf{F}_1 - \mathbf{F}_0)(\mathbf{F}_1 + \mathbf{F}_0)^{-1}\big)`, which is the
trapezoid of :math:`\mathbf{P} : \mathrm{d}\mathbf{F}` itself and zero under a rigid rotation (``Delta_work_conjugacy``).
The difference is added to :math:`W_m` and :math:`W_m^r`, so :math:`W_m^r` then holds the stored
energy plus the work the rate does not conjugate: under 3 and 5, in non-coaxial states, it is not
a state function. :func:`simcoon.umat` does the same when it is given
:math:`\mathbf{F}_0, \mathbf{F}_1`, recovering the lab start stress from the one it is passed
(transported by the caller) as :math:`\mathrm{sym}(\Delta\mathbf{R}^{-1}\hat{\boldsymbol{\tau}}_n\Delta\mathbf{R})`.

.. note::

   Under the rate formulation ``MODUL`` passes corate 3 unconditionally rather than the solver's value, because it
   evaluates its potential at :math:`\mathbf{V}^{el}` in the corotated frame
   (:math:`\mathbf{R} = \mathbf{I}`), where the log box with respect to
   :math:`\ln \mathbf{V}^{el}` **is**
   :math:`\partial \boldsymbol{\tau} / \partial \boldsymbol{\varepsilon}^{el}`. Any other
   corate is already refused upstream under NLGEOM, so passing one down would not be more
   general -- it would be wrong. Under the total formulation the block is told so
   (``KinematicForm::total``) and returns :math:`\mathbf{T}` and
   :math:`\partial\mathbf{T}/\partial\mathbf{H}^{el}` instead.

.. _strain-space-route:

**The total formulation (strain space,** ``kinematic_form="total"`` **).** Under the rate
formulation a law works on :math:`\ln\mathbf{V}` (with the logarithmic rates) and the Kirchhoff stress
:math:`\boldsymbol{\tau}`; under the total formulation it works on :math:`\ln\mathbf{U}` and the
stress :math:`\mathbf{T}` dual to it. It is not an objective rate: no rate and no spin are
involved, and ``corate`` must be left unset. This is the total-Lagrangian
framework of Miehe, Apel and Lambrecht (*Comput. Methods Appl. Mech. Engrg.* 191, 2002,
5383-5425). A log-strain box is fed the material Hencky strain
:math:`\mathbf{H} = \tfrac12\ln\mathbf{C} = \ln\mathbf{U}`, a state function of
:math:`\mathbf{F}` in reference components; the frame increment is the identity and the internal
variables are never transported. Its output is the stress :math:`\mathbf{T}` conjugate to
:math:`\mathbf{H}` per reference volume,
:math:`\mathbf{T}:\dot{\mathbf{H}} = \mathbf{S}:\dot{\mathbf{E}} = \boldsymbol{\tau}:\mathbf{D}`, and
the geometry is applied once, after the box:

.. math::

   \mathbf{S} = \mathbf{T}:\mathbb{P}, \qquad
   \frac{\mathrm{d}\mathbf{S}}{\mathrm{d}\mathbf{E}} = \mathbb{P}^T:\frac{\partial\mathbf{T}}{\partial\mathbf{H}}:\mathbb{P} + \mathbf{T}:\mathbb{L}, \qquad
   \mathbb{P} = 2\,\frac{\partial\mathbf{H}}{\partial\mathbf{C}}, \quad
   \mathbb{L} = 4\,\frac{\partial^2\mathbf{H}}{\partial\mathbf{C}\,\partial\mathbf{C}} .

This is the one exception to "the kernel's ``sigma`` is :math:`\boldsymbol{\tau}`": under the total
formulation the ``sigma`` of a box kernel is :math:`\mathbf{T}`. No kernel changes; the dispatcher converts at entry
(``Kirchoff2LogStressU``) and at exit (``LogStressU2Kirchoff``), and :math:`\boldsymbol{\tau}` stays the only
stress on the route. In the eigenbasis of :math:`\mathbf{C}`, with :math:`\lambda_a` the principal
stretches,

.. math::

   \boldsymbol{\tau} = \mathbf{R}\Big[\sum_{a,b}\kappa_{ab}\,\hat T_{ab}\,\mathbf{n}_a\otimes\mathbf{n}_b\Big]\mathbf{R}^T,
   \qquad \kappa_{ab} = \frac{t_{ab}}{\sinh t_{ab}}, \quad t_{ab} = \ln\frac{\lambda_a}{\lambda_b},

which in Voigt form is :math:`(\mathbf{A}^R)^T` applied to :math:`\mathbf{R}\,\mathbf{T}\,\mathbf{R}^T`.

*Relation to the rate formulation under corate 3.* Under log_R a box kernel already runs in the
rotated frame, on :math:`\mathbf{R}^T(\ln\mathbf{V})\mathbf{R} = \mathbf{H}`: the two formulations
share the box, its strain and its internal variables. They differ in two ways.

- The strain increment is the closed form :math:`\mathbf{H}(\mathbf{F}_1) - \mathbf{H}(\mathbf{F}_0)`
  here, against the integral of :math:`\mathbf{A}^R:\mathbf{D}` there, whose residue is first order
  in the step. An elastic law is therefore exact in a single increment.
- The box output is read as :math:`\mathbf{T}` here and as :math:`\mathbf{R}^T\boldsymbol{\tau}\,\mathbf{R}`
  there. Since :math:`\kappa_{aa} = 1` and :math:`\kappa_{ab} \le 1`, the two readings are identical
  when the stress is coaxial with :math:`\mathbf{U}` (isotropic elasticity, proportional
  plasticity), and differ otherwise. These are two constitutive statements, not two accuracies
  of one. A yield criterion holds on :math:`\mathbf{T}` here and on :math:`\boldsymbol{\tau}` there: in
  simple shear to :math:`\gamma = 1` with :math:`\sigma_Y = 250` MPa, the Cauchy shear and normal
  stresses are 140.9 and -15.7 MPa here against 144.4 and -8.4 MPa under log_R, and the difference
  does not vanish with the step. An anisotropic elastic law :math:`\mathbf{T} = \mathbf{L}:\mathbf{H}`
  derives from a stored energy with its axes on the reference frame, so a closed path stores no
  net work; read as :math:`\boldsymbol{\tau} = \mathbf{R}(\mathbf{L}:\mathbf{H})\mathbf{R}^T` it is only
  Cauchy-elastic.

*Reading one run in the measures of the other.* :math:`\ln\mathbf{V} = \mathbf{R}\,\mathbf{H}\,\mathbf{R}^T`
and :math:`\boldsymbol{\tau} \leftrightarrow \mathbf{T}` are changes of measure, valid for any state:
``sim.stress_convert`` has the keys ``"Kirchoff2LogStressU"`` and ``"LogStressU2Kirchoff"``, and the solver
results offer ``LogStrainU``, ``LogStrainV`` and ``LogStressU`` whenever the strain of the run is
a state function of :math:`\mathbf{F}` (see :doc:`output`). On a path of prescribed
:math:`\mathbf{F}`, and for a log-strain box kernel, a run under corate 3 and a run under the
total formulation then compare as follows: the box output is the same, :math:`\mathbf{R}^T\boldsymbol{\tau}_{\text{rate}}\mathbf{R} =
\mathbf{T}_{\text{total}}`, and the two Kirchhoff stresses differ by the map above,
:math:`\boldsymbol{\tau}_{\text{total}} = (\mathbf{A}^R)^T\boldsymbol{\tau}_{\text{rate}}`, up to the
first-order quadrature of the log_R strain. Under the other rates the box does not see
:math:`\mathbf{H}` in the polar frame and this comparison has no meaning.

*Scope.* Log-strain box kernels (``ELISO``, ``ELIST``, ``ELORT``, ``EPICP``, ``EPJCK``, ``EPCHA``,
the plasticity and hyperelastic names served by the modular engine, and ``PYEXT``, whose Python
law then receives :math:`\ln\mathbf{U}` and must return :math:`\mathbf{T}`), and the kernels built on
:math:`\mathbf{F}` (``NEOHI_closed``, ``NEOHC_closed``, ``OGDEN``, ``SNTVE``). The latter return the
true :math:`\boldsymbol{\tau}` in both
formulations, so their stress does not depend on it; the dispatcher forms their work on the dual
pair, :math:`\tfrac12(\mathbf{T}_n + \mathbf{T}_{n+1}):\Delta\mathbf{H}` with
:math:`\mathbf{T}` computed from :math:`\boldsymbol{\tau}`, and their box is
:math:`\partial\mathbf{T}/\partial\mathbf{H}`. The rate law ``HYPOO`` integrates from its start
stress and is refused. ``MODUL`` is a log-strain box and is served: its mechanisms work on
:math:`\mathbf{H}` and :math:`\mathbf{T}` as those of ``EPICP`` do, and its hyperelastic block
returns :math:`\mathbf{T} = \partial\Psi/\partial\mathbf{H}^{el}` with
:math:`\partial\mathbf{T}/\partial\mathbf{H}^{el}`. A fibre block (Holzapfel) is then the law built
on :math:`\mathbf{F}` at any step, where under corate 3 it converges to it at first order.
``"logarithmic"`` control is refused too (the stored strain is :math:`\ln\mathbf{U}`). The work of
a box kernel is its own :math:`\tfrac12(\mathbf{T}_n + \mathbf{T}_{n+1}):\Delta\mathbf{H}`: the
trapezoidal rule on a conjugate pair, second order in the step, with no correction to apply. Through :func:`simcoon.umat` the caller passes
``kinematic_form="total"``, :math:`\mathbf{H}(\mathbf{F}_0)`,
:math:`\mathbf{H}(\mathbf{F}_1) - \mathbf{H}(\mathbf{F}_0)`,
:math:`\Delta\mathbf{R} = \mathbf{I}`, and both deformation gradients; the stress stays Cauchy in
and out, ``tangent_output="box"`` is :math:`\partial\mathbf{T}/\partial\mathbf{H}` and its
``Lt_convert`` key to :math:`\mathrm{d}\mathbf{S}/\mathrm{d}\mathbf{E}` is ``DTDH_2_DSDE``.

**The Python contract.** ``sim.umat`` returns **Cauchy** stress and the **Kirchhoff box**
tangent :math:`\partial \hat{\boldsymbol{\tau}} / \partial \mathbf{D}_e` with no :math:`J`,
and, under the rate formulation, takes the objective rate as ``kinematic_form=`` (``corate=``, its
former name; default ``3``, log_R)
selecting the rate that tangent is expressed in. The default is contract-preserving: corates 2 and 3 resolve to the same exact map, so it
returns precisely the box the finite kernels used to produce unconditionally. Rescaling that
tangent by :math:`1/J` at the boundary would break ``Lt_convert`` by exactly :math:`J`.

Integrator and tangent operator
===============================

A run makes three numerical choices for the laws with a return mapping, given to ``solve``
and ``sim.umat`` by name:

.. list-table::
   :header-rows: 1

   * - Option
     - Chooses
     - Values
     - Changes the converged answer
   * - ``integrator``
     - the discrete update of the local step
     - ``"cutting_plane"`` (default), ``"closest_point"`` (backward Euler)
     - yes: two different schemes
   * - ``ncp``
     - how the Kuhn-Tucker conditions are solved inside it
     - ``"fischer_burmeister"`` (default), ``"min"``
     - no, to the local tolerance (see below)
   * - ``tangent_mode``
     - the operator handed to the global Newton
     - ``"none"``, ``"continuum"``, ``"algorithmic"`` (default)
     - no: global convergence only

``"algorithmic"`` is the consistent operator of the integrator that ran: Simo-Hughes after the
cutting plane, the exact operator of the closest-point map after the closest point. The closest
point is offered with that operator only (``integrator="closest_point"`` with ``"none"`` or
``"continuum"`` raises ``ValueError``). ``tangent_mode = 3`` (``"closest_point"``), the former way
to ask for the closest point, is no longer accepted: it raises ``ValueError``, pointing to
``integrator="closest_point"``. Inside, the kernels still read one code (``kernel_tangent_mode``,
``parameter.hpp``): 0 = none (Lt = elastic L, explicit integration), 1 = continuum,
2 = algorithmic/Simo-Hughes (**default**), 3 = the closest-point integrator with its exact
operator. Code 3 reaches only the laws that implement it (``check_integrator_serves``; closest-point
coverage below), so a
Python law (``PYEXT``) receives 0, 1 or 2. Pre-2.0 numbering was 0 = continuum,
1 = algorithmic — see :doc:`solver` for the migration note. The
finite-strain hyperelastic models ignore the mode (their tangent is always
the exact one of the hyperelastic law). So do the linear viscoelastic models (``ZENER``,
``ZENNK``, ``PRONK`` and their thermomechanical versions) except for mode 0: their
backward-Euler step is solved in closed form, so the tangent, and for the thermomechanical
versions ``dSdT`` and the heat-source derivatives ``drdE``, ``drdT``, are the exact
derivatives of the discrete update. For plasticity, the algorithmic operator is exact with a
von Mises criterion and isotropic hardening only; with kinematic hardening or a Hill, DFA,
Drucker or Tresca criterion the cutting-plane update is path dependent and the operator is
approximate (see the note on the scope of the algorithmic tangent in :doc:`solver`).

**Closest-point coverage (2.2).** The closest-point integrator with its exact operator runs for
``MODUL`` rows whose criterion has a flow Hessian (von Mises, Hill, DFA, Ani; with every
hardening law and a hyperelastic block) and for the legacy names the modular engine serves,
for the dedicated kernels ``EPICP_closed``, ``EPCHA_closed`` (``EPICP_closed`` / ``EPKCP_closed``
thermomechanical), and for the SMA kernels ``SMADI`` / ``SMADC`` /
``SMAAI`` / ``SMAAC`` and ``SMRDI`` / ``SMRDC`` / ``SMRAI`` / ``SMRAC`` (analytic total
derivatives through the mean transformation strain and the reorientation back-strain; about
2x the cost of mode 2 per increment on a superelastic loop, 3x with the reorientation row,
up to 5x on a non-proportional path, for an operator exact to :math:`10^{-8}`). For
``EPJCK`` (and the adapter-served von Mises + isotropic names) mode 3 is **identical to mode
2 by construction**: with a radial flow the cutting-plane update already is the closest-point
projection and the algorithmic operator is the exact Jacobian of the discrete map, including
the rate term. Any other law, the thermomechanical SMA kernel and ``SMA_mono`` included, refuses
``integrator="closest_point"`` with ``ValueError`` rather than running its cutting plane
silently; inside ``MODUL``, a Tresca or Drucker row keeps the cutting plane and its mode-2
operator. A condensed stress state (plane stress, beam,
truss) keeps mode 3: the law runs 3D and the state is solved around it.

Numerical options of the local Newton
=====================================

These options change *how* the local return mapping is solved, not the constitutive law. They
are settings of the run, not material parameters: they are given to ``solve`` or ``sim.umat``
by name, and the ``props`` stay those of the law.

``ncp``: a faster local solve for simple plasticity
---------------------------------------------------

``solve(..., ncp="min")`` or ``sim.umat(..., ncp="min")`` switches the complementarity function
of the cutting plane from Fischer-Burmeister (``ncp="fischer_burmeister"``, the default, the
historical results bit for bit) to the min function, the classic cutting-plane step. The modular
engine (``MODUL`` and the legacy names it serves, mechanical and thermomechanical) and
``EPICP_closed`` read it, in their cutting plane and in their closest point (the step of the
reduced Newton). Any other law or a mean-field model raises ``ValueError`` for ``ncp="min"``.

**Use it** for a single yield surface with isotropic or linear (Prager) kinematic hardening.
The return is the same to the local tolerance, and the local iterations drop from about 6 to 3.
In ``MODUL`` this includes a yield surface composed with viscoelastic (Prony) or damage
mechanisms: their rows are always satisfied and uncoupled, so the min step applies to the yield
surface next to them. With several coupled yield surfaces the option has no effect:
Fischer-Burmeister is kept.
If the min step makes the error grow, the solver switches back to Fischer-Burmeister for the
rest of the increment.

**Do not use it** with Armstrong-Frederick (Chaboche) back stresses. ``MODUL`` refuses the
combination, and the Chaboche kernels do not take the option. Those back stresses are integrated
along the cutting-plane iterates, so the full min steps land on a different answer, without
any growth in the residual that the fallback could detect. On coarse increments of the
Chaboche identification programmes, this answer is up to 46 % off the closest-point update,
where Fischer-Burmeister stays within 4 %. For an accurate Chaboche return, use
``integrator="closest_point"``.

How the two functions work
~~~~~~~~~~~~~~~~~~~~~~~~~~

At each increment the law seeks the plastic multiplier :math:`\Delta p` of each yield surface
that satisfies the Kuhn-Tucker conditions

.. math::

   \Phi \le 0, \qquad \Delta p \ge 0, \qquad \Phi\,\Delta p = 0,

with :math:`\Phi` the yield function (for EPICP,
:math:`\Phi = \sigma_{Mises} - R(p) - \sigma_Y`): either the point stays elastic
(:math:`\Delta p = 0`, :math:`\Phi \le 0`) or it lies on the surface (:math:`\Phi = 0`,
:math:`\Delta p > 0`). The cutting plane solves them by Newton on :math:`\Delta p`, with the
linearisation :math:`\Phi(\Delta p + \delta p) \approx \Phi + B\,\delta p`,
:math:`B = \partial\Phi/\partial\Delta p` (for von Mises with isotropic hardening
:math:`B \approx -(3G + H)`, :math:`H` the hardening slope). A complementarity function turns
the three conditions into one equation. Both variants use the normalised multiplier
:math:`\Delta p^* = |B|\,\Delta p`, a stress like :math:`\Phi`.

**Fischer-Burmeister** (``ncp="fischer_burmeister"``) solves

.. math::

   \sqrt{\Phi^2 + \Delta p^{*2}} + \Phi - \Delta p^* = 0 .

The equation is smooth away from :math:`(0, 0)` and decides by itself whether the surface is
active: Newton never chooses a branch. It extends as is to a coupled system of several surfaces
(multi-mechanism ``MODUL``, the SMA kernels). Its first step is short: from
:math:`\Delta p = 0` on a violated surface (:math:`\Phi > 0`), the residual is :math:`2\Phi`
and its derivative :math:`2B - |B| = 3B`, so the step is :math:`-2\Phi/(3B)`, two thirds of the
full cutting-plane step :math:`-\Phi/B`. The return then needs about 6 local iterations. The
stopping criterion is the residual divided by :math:`\sigma_Y`, summed over the surfaces.

**Min** (``ncp="min"``) solves

.. math::

   \min(\Delta p^*, -\Phi) = 0

by semi-smooth Newton (primal-dual active set, Hintermueller, Ito and Kunisch 2002), which
picks a branch at each iteration:

* active surface (:math:`\Phi > 0` or :math:`-\Phi < \Delta p^*`): it solves
  :math:`\Phi + B\,\delta p = 0`, the full step :math:`\delta p = -\Phi/B` of the classic cutting
  plane (Ortiz and Simo 1986);
* inactive surface: the multiplier goes back to zero, :math:`\delta p = -\Delta p`.

With linear hardening :math:`\Phi` is affine in :math:`\Delta p`, so the first step is the exact
return, and one more iteration confirms it. With a nonlinear hardening (power law) the return
takes about 3 iterations instead of 6. The stopping criterion is twice the residual of the chosen
branch over :math:`\sigma_Y`: on the first step this equals the Fischer-Burmeister residual, and
near the solution it is at least as strict.

**Fallback.** Within an increment, the min step switches to Fischer-Burmeister for the rest of
the local solve as soon as its error grows from one iteration to the next (the full step
overshoots, or the branch choice oscillates) or :math:`B` vanishes. Coupled surfaces stay on
Fischer-Burmeister from the start, because the plain semi-smooth Newton can cycle there, the set
of active surfaces changing from one iteration to the next without converging.

**Which rows count.** The min step applies row by row when the local Jacobian :math:`B` is
diagonal, i.e. when the rows are independent scalar problems. That is the case of a single yield
surface, and of a yield surface next to the rows of the viscoelastic and damage mechanisms of
``MODUL``: those are integrated outside the Newton (the Prony branches in closed form, the damage
by a fixed point), so their rows are always satisfied (:math:`\Phi = -Y^{crit}`), carry no
multiplier and no cross term, and stay at zero. Two yield surfaces share the stress, so their
rows are coupled (off-diagonal terms of :math:`B`) and keep Fischer-Burmeister.

.. list-table:: What ``ncp="min"`` does, case by case
   :header-rows: 1

   * - Case
     - Result
   * - EPICP, linear isotropic hardening
     - one exact step and one check; the Fischer-Burmeister return to :math:`10^{-10}`
   * - EPICP, power-law isotropic hardening
     - about 3 iterations instead of 6; the Fischer-Burmeister return to :math:`10^{-7}`
   * - EPKCP, ``MODUL`` with isotropic and Prager hardening
     - the Fischer-Burmeister return to the local tolerance (EPKCP and ``MODUL`` agree to
       :math:`10^{-12}`)
   * - ``MODUL``, one yield surface with viscoelastic (Prony) or damage mechanisms
     - the min step on the yield surface, the other rows untouched; the Fischer-Burmeister
       return to the local tolerance
   * - ``MODUL`` with several yield surfaces (coupled)
     - no effect: Fischer-Burmeister is kept
   * - ``MODUL`` with an Armstrong-Frederick back stress
     - refused (``ValueError``) under the cutting plane, see above; taken under the closest point
   * - any other law, any other thermomechanical block, a mean-field model
     - refused (``ValueError``); the sub-phases of a mean-field model keep Fischer-Burmeister
   * - ``integrator="closest_point"``
     - the min step on the reduced Newton system of the closest point, every composition
       (Armstrong-Frederick included): one backward-Euler system, the same root as
       Fischer-Burmeister (measured 1e-10 to 4e-10 on a coarse non-proportional path, the
       local tolerance), where the cutting plane with Fischer-Burmeister sits 1 to 4 % away
       for an Armstrong-Frederick back stress

On a batch of :math:`10^4` EPICP points (``sim.umat``, one thread), one increment takes 20 ms
with Fischer-Burmeister and 14 ms with ``ncp="min"``. On a cyclic uniaxial ``MODUL`` solve
(150 increments) the whole run drops by 17 % with one Prony branch next to the yield surface,
29 % with three, 23 % with a damage mechanism.

EPICP consistent tangent
------------------------

With ``tangent_mode = 2``, EPICP builds the consistent tangent of the radial return in closed
form (Simo and Taylor 1985),
:math:`3K\,\mathbb{I}_{vol} + 2G\theta\,\mathbb{I}_{dev} - 2G\bar\theta\,\mathbf{N}\otimes\mathbf{N}`.
It is the operator of the generic assembly to round-off at the converged return, without its
6x6 inverse.

Choosing the tangent and the integrator for EPCHA (0D uniaxial)
----------------------------------------------------------------

* ``tangent_mode = 1`` (continuum) gives the same iterates as the default mode 2 (to
  :math:`10^{-14}` on the Chaboche identification programmes) for about 27 % less time.
* ``integrator="closest_point"`` is a different update,
  backward Euler: 1 to 2 % off the cutting plane there, at about 2.3 times the cost.

Global Newton tolerance
-----------------------

The global Newton has its own option, ``solve(..., precision_rel=...)``
(``solver_params::precision_rel``). The residual tolerance becomes
:math:`\max(\mathrm{precision}, \mathrm{precision\_rel}\,|\boldsymbol\tau|)`, see
:doc:`python_solver`.

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
  cyclic strain path each test run — machine precision for every family (the
  Chaboche row of the gtest keeps a 2e-3 tolerance, measured at 4e-15).

A benchmark row per family lives in ``benchmark/bench_legacy_vs_modular.py``:
results show no measurable adapter overhead, and the modular
engine runs at 0.8-1.2x the speed of the removed kernels on all families.
