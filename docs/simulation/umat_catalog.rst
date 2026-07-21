====================================
Constitutive model (UMAT) catalog
====================================

Every constitutive law is selected by its 5-character ``umat_name``. Since the
modular UMAT framework, three kinds of implementation coexist behind those
names — **the calling convention is identical for all of them** (same
``umat_name``, same props, same solver/FEA usage):

- **modular (native)**: the composable ``MODUL`` engine, configured from a
  props stream (see :mod:`simcoon.modular` for the Python builder).
- **modular (adapter)**: a legacy name whose dedicated kernel was removed
  after its equivalence with a ``MODUL`` configuration was proven
  (bit-identical results, dedicated tests in ``test_modular.py``); a
  translator maps the legacy props to the modular configuration at each call.
- **legacy (kept)**: a dedicated, self-contained implementation kept either
  for pedagogy (readable single-file reference of the CCP return mapping) or
  because no modular equivalent exists.

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

- **SMA**: SMAUT/SMANI/SMADI/SMADC/SMAAI/SMAAC (unified), SMRDI/SMRDC/SMRAI/
  SMRAC (unified with reorientation), SMAMO/SMAMC (monocrystal).
- **Finite strain**: HYPOO (hypoelastic orthotropic), SNTVE (Saint-Venant),
  NEOHI/NEOHC (Neo-Hookean), MOORI, YEOHH, ISHAH, GETHH, SWANH
  (invariant-based hyperelasticity).
- **Multiscale**: MIHEN, MIMTN, MISCN, MIPLN.
- **Plugins**: UMEXT (external dylib), UMABA (Abaqus wrapper).

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
