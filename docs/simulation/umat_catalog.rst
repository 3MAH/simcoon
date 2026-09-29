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
     - Pedagogical reference; the modular Viscoelasticity mechanism takes the
       same closed-form step and is identical to it (round-off)
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
  the kinematics are exact: damage subtracts no inelastic strain (it scales the
  stiffness instead), so the elastic stretch is the total one and the fibres follow it.
  Its driving force is the thermodynamic force :math:`Y = -\partial\psi/\partial D = \psi_0`, the
  undamaged energy, evaluated as
  :math:`\tfrac12\,\boldsymbol{\tau}_{eff} : \mathbf{M}_t : \boldsymbol{\tau}_{eff}` on the effective
  stress :math:`\boldsymbol{\tau}/(1-D)`: exact for a linear block, an approximation of
  :math:`\psi_0` through the tangent compliance for this nonlinear potential.

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
``output_convention_of`` (``umat_smart.cpp``). A kernel that fails to declare it is
rejected rather than given a default -- a missing declaration is an error of exactly
:math:`J`, and invisible near :math:`J = 1`:

- ``kirchhoff`` -- every native kernel.
- ``cauchy`` -- only the ``UMEXT`` / ``UMABA`` plugin adapters, whose contract belongs
  to the host code (Abaqus ``DDSDDE`` is Cauchy-based). ``select_umat_M_finite`` converts
  those to :math:`\boldsymbol{\tau}` on the way out.

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
stiffness, Hill-type criteria, HOLZA fibres in ``MODUL``) therefore rotate with the body, as with
Abaqus ``*ORIENTATION`` under NLGEOM, and their tensorial ``Statev`` are material-frame
components. The kernels built from :math:`\mathbf{F}` (SNTVE, the Neo-Hookean and invariant
family, OGDEN, HOLZA) are objective by construction and keep the lab frame.

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
   internal variables each kernel declares in ``umat_conventions`` (``umat_smart.cpp``: EPICP,
   EPCHA, the elastic and hyperelastic kernels); a kernel that has not declared them (the
   plasticity names served by the modular engine, ``PYEXT``) is refused under corates 4 and 5,
   and ``MODUL`` under any corate but 3. ``EPCHA`` is also refused under corate 4: it stores
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

   ``MODUL`` passes corate 3 unconditionally rather than the solver's value, because it
   evaluates its potential at :math:`\mathbf{V}^{el}` in the corotated frame
   (:math:`\mathbf{R} = \mathbf{I}`), where the log box with respect to
   :math:`\ln \mathbf{V}^{el}` **is**
   :math:`\partial \boldsymbol{\tau} / \partial \boldsymbol{\varepsilon}^{el}`. Any other
   corate is already refused upstream under NLGEOM, so passing one down would not be more
   general -- it would be wrong.

**The Python contract.** ``sim.umat`` returns **Cauchy** stress and the **Kirchhoff box**
tangent :math:`\partial \hat{\boldsymbol{\tau}} / \partial \mathbf{D}_e` with no :math:`J`,
and takes ``corate=`` (default ``3``, log_R) selecting the rate that tangent is expressed
in. The default is contract-preserving: corates 2 and 3 resolve to the same exact map, so it
returns precisely the box the finite kernels used to produce unconditionally. Rescaling that
tangent by :math:`1/J` at the boundary would break ``Lt_convert`` by exactly :math:`J`.

Tangent-operator mode
=====================

All models receive the solver's ``tangent_mode`` (named constants in
``parameter.hpp`` / ``sim.tangent_*``): 0 = none (Lt = elastic L, explicit
integration), 1 = continuum, 2 = algorithmic/Simo-Hughes (**default**),
3 = closest-point (reserved). Pre-2.0 numbering was 0 = continuum,
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
