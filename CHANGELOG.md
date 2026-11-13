# Changelog

All notable changes to simcoon are listed here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and simcoon uses
[semantic versioning](https://semver.org/). Each pull request adds its entry
under **Unreleased**; the section takes the version number at the release.

## In progress

Known gaps and planned work, not yet in a release. An item moves to
**Unreleased** with the pull request that closes it.

- `integrator="closest_point"` is provided with its algorithmic tangent only
  (`tangent_mode="algorithmic"`): its continuum operator is not implemented.
- Other complementarity functions to measure against Fischer-Burmeister and
  `min`, under both integrators: active set, scaled Fischer-Burmeister.
- The numerical options of the mean-field schemes (MIMTN/MISCN closed form,
  rule, Anderson, tolerance) are still props slots: they are per phase (nested
  schemes), so they need a named builder rather than a run option.
- fedoo does not pass `ncp` and `integrator` through yet.
- Under the rounding policy of 2.4 (results move at the round-off level; the gate is
  `dump_results.py compare --rtol 1e-11`, the heat tables judged against the scale the dumper
  stores with them; the zero-tolerance pass stays available to pin a defect), still to do: the
  Fastor to Armadillo rewrite. `mul6`/`mul66` keep their fixed order: they replace the BLAS, they do
  not imitate it.
- Migration of the laws to the modular engine (MODUL), state in
  `benchmark/README.md`. Routed through the adapter (code 201): ELISO, ELIST,
  ELORT, EPKCP, EPHIL, EPTRI, EPHAC, EPANI, EPDFA, EPCHG, EPHIN. To do:
  - adapters for the hyperelastic laws whose MODUL equivalent exists (NEOHC,
    MOORI, YEOHH, ISHAH, GETHH, SWANH, HOLZA; NEOHI to confirm);
  - no MODUL equivalent yet: OGDEN, SNTVE, HYPOO, EPJCK, LLDM0, ZENNK (N >= 2: the
    Kelvin-to-Maxwell map has no closed form), the SMA family, the micromechanics schemes;
  - the thermomechanical dispatch (`umat_T`): `MODUL` has its face (every block, every
    mechanism) and serves EPICP, EPKCP and PRONK through adapters; under the cutting plane
    `MODUL` is 1.3x slower than EPICP_closed on plastic increments (J2: the generic algorithmic
    assembly where the dedicated kernel has the closed radial-return form) and within 10 % of
    EPKCP_closed (Prager) and PRONK_closed; under the closest point it is faster than the
    dedicated kernels (0.9x), both running the generic solver (`benchmark/README.md`).

## Planned for 3.0

Removals and renames scheduled for the next major version (breaking, kept in 2.x).

- The migration message for `tangent_mode=3` / `"closest_point"`: the check
  `tangent_mode in ("closest_point", 3)` of `solve` (`simcoon/solver/core.py`)
  and its C++ counterpart in `kernel_tangent_mode` (`umat_smart.cpp`). The
  value then gets the generic error of an unknown tangent mode; the internal
  kernel code 3 of the closest point is unchanged.
- Possibly rename `ZENER` / `ZENNK` after their rheology: they are the Poynting-Thomson
  (Kelvin) form of the standard linear solid and its generalized Kelvin chain, not the Zener
  (Maxwell) form; the 2.x names stay for compatibility.

## [Unreleased]

### Added

- The external plugins run under finite strain (NLGEOM), with no change to the plugin API: the
  finite dispatch serves UMEXT, UMABA and UMANS as log-strain boxes in the material frame. UMEXT
  follows the simcoon convention exactly like PYEXT and the built-in kernels (logarithmic strain in,
  Kirchhoff stress and box tangent out: a small-strain plugin runs unchanged; the elastic test plugin
  is ELISO to round-off). The host plugins keep their conventions (`StressMeasure::cauchy`, now
  complete on both sides of the box): the Cauchy start stress `tau_n / J_n` in, the Cauchy stress and
  the host Jacobian out, converted to the route (`tau = J sigma`; Abaqus DDSDDE in the Jaumann rate,
  Ansys dsdePl in the Green-Naghdi rate with its `sigma (x) I` term added, both mapped through dS/dE to
  the box of the solver's corate; the hosts' symmetric Jacobians give the symmetric part of the
  exact tangent). The host work trapezoid is replaced by the Kirchhoff one. Corates 0 to 3 are
  served, the total formulation for UMEXT; corates 4 and 5 are refused (no declared statev layout);
  `sim.umat` still refuses the plugins. The sample `UMAT_ABAQUS_ELASTIC.f` now returns the NLGEOM
  Jacobian (`sigma (x) I` term). Tests: the three plugins against ELISO (UMEXT bit-identical, host
  Cauchy stress = `C : e` to 1e-9), the converted host tangent against central differences, a rigid
  spin, the admission; `test_external_plugins_refused_under_finite_strain` is gone.

## [2.4.0] - 2026-11-12

The modular engine (`MODUL`) gains its thermomechanical face and the finite-strain
thermomechanical route (#155), and becomes the default kernel of EPICP, EPCHA, PRONK and ZENER
(the dedicated kernels stay under `NAME_closed`). One registry of the law names replaces the
hand-kept dispatch tables (#156), the kernels take their natural fixed-size form with results at
round-off (#154), and the closest point gets by-reference callbacks and `ncp="min"` (#157).
Breaking at the C++ API only (closest-point callbacks, three removed accessors, convention
types moved to `umat_registry.hpp`): see **Breaking changes**.

### Breaking changes

- C++ API of the closest-point solver (`return_mapping.hpp`): the row callbacks of
  `ReturnMechanism` (`dPhi_dsigma`, `Lambda`, `dLambda_dsigma`) and the hooks `K` and
  `flow_state_coupling` of `ReturnStateHooks` write into buffers the solver owns (a 6, a 6x6,
  the NxN block, the 6xN state chain, sized once per call) instead of returning a `vec`, a
  `mat` or a `std::vector<vec>` by value. `Phi`, `update_state`, `elastic_response` and the
  single-mechanism overload are unchanged.
- One registry of the law names (`umat_registry.hpp`, `umat_registry()`: for each name the kernel
  of each face of the dispatch as a `UmatKernel`, its finite-strain `umat_convention`, the
  condensed stress states `sim.umat` has validated and the elastic stiffness of `L_eff`)
  replaces the seven tables that were kept by hand: the three name-to-id maps of the dispatch,
  the convention table, the capability sets of the closest point, `ncp` and the total form, the
  `L_eff` map and the two of the Python wrapper. The dispatch switches read the kernel, never an
  integer code; a row without its `case` is a thrown `std::logic_error`. C++ API: `finite_umat_names()`,
  `small_strain_umat_names()` and `thermo_umat_id()` are removed (`umat_entry_of(name)` answers
  all three), `StressMeasure`, `StatevTensor` and `umat_convention` move to `umat_registry.hpp`
  (included by `umat_smart.hpp`), `sub_phase_shape` and `get_L_elastic` read the registry.
  Results are bit-identical; the names, props and ids of the Python API are unchanged.
- The modular engine is the default for every legacy name: EPICP, EPCHA and PRONK (mechanical)
  and EPICP, EPKCP and PRONK (thermomechanical) are served by the modular adapters (code 201,
  `legacy_adapters.cpp`, the thermomechanical face through `umat_legacy_modular_T`), as EPKCP
  and the Hill/Chaboche families already were. The dedicated closed-form kernels answer the
  `_closed` names (`EPICP_closed`, `EPCHA_closed`, `PRONK_closed`; `EPKCP_closed` on the
  thermomechanical face) on every route (solver, `sim.umat` / `sim.umat_T`, Abaqus `CMNAME`,
  which is now read whole up to its first blank, `_R` excluded, the `_closed` suffix in any
  case; a token outside the name tables that does not end with `_closed`, `EPICP_STEEL`, falls
  back to its first five letters with one warning per name). Results of EPICP and PRONK move at the solver tolerance (parity 1e-10 to 1e-11), those
  of EPCHA at round-off as well (same cutting plane; the `testBin/Umats/EPCHA` reference is
  replayed by both). The modular statev layouts replace the dedicated ones:
  EPICP slot for slot (`[T_init, p, EP]`), EPCHA `[T_init, p, EP, a_1, a_2]` (20, no stored
  back stress or `H_p`), PRONK `[T_init, (v_i, EV_i) x N]` (1 + 7N, no total viscous strain).
  EPICP and EPCHA keep the corate coverage their kernels had, 0 to 5 (layouts declared for
  4/5); EPCHA under Truesdell is no longer refused (the modular composition stores the
  back-strains only); PRONK joins the finite-strain route (corates 4/5 refused: no declared
  layout). The `MODUL` guard (corate 3 or the total form) stays on the `MODUL` name.
- The tangent assembly (`tangent_assembly.hpp`) takes its per-mechanism quantities as matrices:
  `kappa_j` and `dPhidsigma_l` are `6 x N` (one column per mechanism / criterion), the flow
  Hessians `6 x 6N` (column blocks), `HessianProvider` returns that matrix and
  `ContinuumTangent::P_epsilon` is `6 x N` (`P_epsilon.col(l)` where `P_epsilon[l]` was read);
  `ReturnMappingResult` carries the same forms. The single-mechanism overloads are unchanged (a
  `vec` is a `6 x 1` matrix for the N-form). `StrainMechanism::dPhi_dsigma` / `kappa` return a
  `const arma::mat&` (`6 x N`, engineering Voigt) instead of `std::vector<tensor2>`,
  `dLambda_dsigma` a `const arma::mat*` (`6 x 6N`) instead of `std::vector<tensor4>`,
  `inelastic_strain` a `const arma::vec&`, and `dissipation_sensitivity` takes `P_epsilon` as a
  `6 x k` matrix. `KinematicHardening::compute_backstresses` fills a `6 x n` matrix (one column
  per branch, stress Voigt) instead of `std::vector<tensor2>`, and `hardening_modulus`, `update`,
  `refresh_state` take the flow normal as an engineering strain Voigt `vec` instead of a
  strain-typed `tensor2` (the `tensor2` conveniences `total_backstress` and
  `hardening_modulus(n, ivc)` stay). A plugin or a test that wrote the vector or tensor forms
  adapts with `.col()` / `.cols()` / `.voigt()`.
- `sigma_in` and `sigma_in_start` are removed from `state_variables_M` and `state_variables_T`
  (members, the long constructors and `update`, which lose two arguments), from the legacy
  `umat_external` declaration and from the plugin example: no kernel wrote them and nothing read
  them (the solver block that built `sigma_in_red` without using it is removed too). The plugin
  entry points (`umat_plugin_ext_api`, UMABA, UMANS) never carried them.
- The thermal derivatives have one shape: `dSdT` and `drdE` are `vec::fixed<6>` members of
  `state_variables_T` and `arma::vec &` arguments of the ten `_T` kernels, of `heat_source_increment`,
  `condense_tangent_T`, `abaqus2smart_T` / `smart2abaqus_T` and `external_umat_T`; `drdT` is a
  `double`. They were `mat` of shape 1x6, 6x1 or 1x1 depending on the writer, vectorised and reshaped
  at every crossing. Results are bit-identical; `sim.umat_T` and the solver results keep their
  array shapes.
- The hyperelastic block of `MODUL` reads one more prop, the entropic fraction `omega`, after
  `alpha`: `[potential, n_params, params..., alpha, omega]` (`ElasticityModule::props_count`
  is `4 + n_params`). The Python blocks write it (`omega=0.` by default, results bit-identical);
  a hand-written stream must append the slot.

### Added

- `sim.umat_T` takes `ncp` (the thermomechanical face of `MODUL` and its adapters read it, as
  the solver's thermomechanical blocks now do); `benchmark/dump_results.py batch_T`.
- The Abaqus bridges (`umat_singleM`, `umat_singleT`) run the finite route for every `(NDI, NSHR)`
  under NLGEOM: a condensed state (plane stress, shells, beams, trusses) solves its free
  stretches from the strain it keeps in STATEV (`condensed_F`, moved from the Python wrapper
  into `condensed_state.hpp`; `abaqus2smart_finite` takes the condensed state,
  `finite_box_strain_end`), a Newton on the free components of the box strain increment around
  the 3D finite dispatch (`run_umat_M` / `run_umat_T`); the host's `F33` is a placeholder. It
  used to fall back to the corotational small-strain route silently. Checked against the 3D point
  run at the solved thickness stretch (stress-free there, in-plane stress and `RPL` to 1e-8).
- The thermomechanical face of `MODUL` (`umat_modular_T`, dispatch id 200, `ModularUMAT::run_T`):
  the name `MODUL` with the props `[rho, c_p, stream]` serves the thermomechanical blocks of
  the solver, `sim.umat_T` and the Abaqus `umat_singleT`, under small and finite strain (both
  forms), every tangent mode, the closest-point integrator and `ncp="min"`. The heat of the
  increment and its derivatives are written once on the mechanisms (the entropy
  `c_0 ln(T/T_0) + alpha : sigma`, the dissipation of the increment and
  `StrainMechanism::dissipation_sensitivity`), so a composition equivalent to `EPICP` or
  `EPKCP` is that kernel to round-off (stress, temperature, the four works, `Wt`), and a
  Maxwell composition is `PRONK` (`ViscoelasticMechanism::dissipation_sensitivity`, exact:
  linear branches). `ZENER` is the one-branch Maxwell composition with the mapped branch
  (`X_1^M = X_0^2/(X_0 + X_1)`, `eta^M = eta^K X_0^2/(X_0 + X_1)^2` per isotropic mode), to
  1e-10. A damage mechanism dissipates the trapezoid `1/2 (Y_n + Y_n+1) dD` on its driving
  force, its heat derivatives on the effective stress (`D'(Y)` of the committed state; exact
  against finite differences, alone and on top of plasticity). A hyperelastic block carries
  its temperature dependence through its new `omega` slot, `psi = (1 - omega + omega
  theta/theta_0) W_iso + U(J)`: `omega = 0` (default) is energetic thermoelasticity, the limit of
  the linear blocks; `omega = 1` the entropic elasticity of a rubber (Gough-Joule heating,
  thermoelastic inversion at `lambda - 1 = alpha theta_0`, reproduced), the volumetric term
  energetic, the heat capacity at fixed deformation `c_0` (`hyper_potential_energy` gives the
  value of every isochoric potential). Serves every block and every mechanism, `omega != 0`
  composed with mechanisms included: the entropy carries `-W_iso(H_e)/theta_0` on the elastic
  strain `H_e = eps - eps_in - alpha dtheta`, its derivatives through the chain `dH_e = L^-1
  (dsigma_eff - (omega/theta_0) T_iso dtheta)` on the tangents of the effective stress (exact
  against finite differences with plasticity, viscoelasticity and damage).

### Changed

- `ncp="min"` under `integrator="closest_point"`: the closest point takes the min step (with its
  Fischer-Burmeister fallback) on its reduced Newton system, for every composition the modular
  engine serves and for `EPICP_closed`, Armstrong-Frederick back stresses included: one
  backward-Euler system, the same root as Fischer-Burmeister. Measured on a coarse
  non-proportional path (tension then shear, 5 and 20 increments): 1.5e-10 to 4e-10 relative
  between the two steps for AF, Prager, J2 and EPICP_closed, where the cutting plane with
  Fischer-Burmeister sits 1.3e-2 to 3.8e-2 from the closest point for AF. The refusal of
  `ncp="min"` with an Armstrong-Frederick back stress stays under the cutting plane (it is
  raised by `ModularUMAT::run` now, not `set_ncp`); the closest-point refusal of `ncp="min"`
  is lifted (`check_ncp_serves` loses its `tangent_mode` argument). The Fischer-Burmeister
  path is bit for bit. Measured on a batch of 997 points under the closest point (`sim.umat`,
  one thread, six increments): MODUL with an Armstrong-Frederick back stress 34.0 to 25.5 ms
  (-25 %), EPICP_closed 23.6 to 18.7 ms (-21 %), the two returns 4e-11 apart.
- `benchmark/dump_results.py` stores a reference scale with the heat tables: `r_scale` for each
  thermomechanical batch increment and `Q_scale` / `r_scale` for each solver history,
  `|dWm_d| + |dWt| + |r|` per increment, the works the heat is a difference of; `compare`
  judges `Q` and `r` against it when both dumps carry it, so the gate of the rounding policy is
  `--rtol 1e-11` alone (`--field-rtol` stays for a table without a scale).
- `ZENER` runs the modular engine (adapter `translate_ZENER`, `legacy_adapters.cpp`) on every
  face: the Kelvin unit in series with the spring becomes one Maxwell branch, per isotropic mode
  `X_M = X_0^2/(X_0 + X_1)`, `eta_M = eta_K X_0^2/(X_0 + X_1)^2` (`simcoon.modular.zener_as_maxwell`
  gives the branch). The dedicated Kelvin kernel answers `ZENER_closed` (mechanical and
  thermomechanical). Same instantaneous and long-term moduli and relaxation time: the response
  moves at the round-off of the two inverses (measured 5e-15 relative on the relaxation path of
  `PRONK_path.json`, below 1e-10 on the thermomechanical face, works and heat included). The statev layout keeps its count
  (8) but `EV` (statev 2..7) is the Maxwell branch strain, not the Kelvin viscous strain of
  `ZENER_closed`; `ZENER` gains the finite-strain route (log_R and the total form, as PRONK).
  `bench_legacy_vs_modular.py` times the pair.
- Speed of the closest point, bit for bit: the solver keeps one fixed-size scratch per row
  (gradient, flow, `kappa`, state chain, `M^-1 kappa`, `dLambda/dsigma`) and the callbacks
  write into it; no temporary, no `std::vector<vec>`, no `zeros(N, N)` per iteration. Measured
  on the `sim.umat` / `sim.umat_T` loops of `benchmark/README.md` (2000 points, 1 thread,
  median per plastic increment): EPICP_closed 4.92 to 4.23 us (-14 %, thermomechanical -15 %), MODUL
  J2 4.85 to 3.71 us (-24 %, -16 %), EPKCP_closed thermomechanical -9 %, MODUL Prager -15 % /
  -11 %, MODUL hyperelastic + J2 -11 % / -8 %. The MODUL/EPICP_closed ratio under the closest
  point is now 0.88 (mechanical) and 0.96 (thermomechanical).
- The `I_2` terms of the isochoric hyperelastic response (`isochoric_invariants`,
  `tau_iso_hyper_invariants`, `L_iso_hyper_invariants`, hence every invariant-based kernel and
  the hyperelastic block of `MODUL`) are evaluated in their cofactor form,
  `I_1 dev(b_bar) - dev(b_bar^2) = -dev(cof b_bar)` and `I_2 = tr cof(b_bar)`, instead of the
  difference of two `O(lambda^4)` terms. Measured against an extended-precision reference at a
  stretch `lambda`: the difference form keeps `16 - 3 log10(lambda)` digits (1.5e-7 relative at
  `lambda = 1000`), the cofactor form is exact to round-off in the principal frame and 3 to 10
  times better in a rotated frame, where the representation of `b` itself sets the limit. Results
  of the laws with an `I_2` dependence (MOORI, ISHAH, GETHH, SWANH) move at the round-off level
  (measured on a uniaxial path to a log strain of 1.5: tangent below 6e-15 relative, work below
  4e-16, stress below 2e-11 through the Newton of the stress-free components); every other law is
  bit for bit.
- Speed of the thermomechanical faces and of the closest point, results at round-off
  (`dump_results.py compare --rtol 1e-11` passes; the mechanical cutting plane is bit for
  bit): the transposed products of the heat derivatives go through `mul6_At`
  (`small_dense.hpp`, fixed order, no BLAS: MODUL and the dedicated `_T` kernels, -6 %);
  `IsotropicHardening::R_and_dR` evaluates `R` and `dR/dp` with one power or exponential
  (`dR = m R / p` for the power law); the viscoelastic mechanism sums its weighted branch
  strain only when a branch moved and keeps its strain map (`I - dEV/dE`) from `predict()`,
  applied with `mul66` (MODUL viscoelastic thermomechanical -13 %); the entropy of an entropic
  hyperelastic block at the start of the increment comes from `ElasticityModule::iso_energy`
  (the potential alone, no stress or tangent) through a two-entry memo (the end state of the
  last call, or the same start again: no evaluation in a sequential solve; thermal face of
  the entropic block with plasticity -50 %); the generic closest-point solver keeps its
  per-iteration vectors, reduced system, active set and step across iterations (-6 % on every
  closest-point law, bit for bit).
- The plastic work of the modular engine has the thermodynamic split of the dedicated kernels:
  `Wm_ir` is the isotropic hardening work `1/2 (R_n + R_n+1) Dp`, the energy stored in the back
  stresses, the trapezoid `1/2 (X_i,n + X_i,n+1) : D a_i` (`KinematicHardening::stored_increment`),
  goes to `Wm_r` (a reverse loading recovers it), and `Wm_d` is the rest of
  `1/2 (sigma_n + sigma_n+1) : D EP`. `Wm` and the stress are unchanged; `Wm_r` / `Wm_ir` / `Wm_d`
  move for `MODUL` and the adapter-served names (EPKCP, EPCHG, EPHIL...): the isotropic hardening
  energy was dissipated and the back-stress one booked irrecoverable, at the end point.
- The thermomechanical face of `MODUL` shares the multiplier assembly of its mechanical tangent
  (`ContinuumTangent` of `compute_tangent`, kept in the instance) instead of rebuilding it at the
  committed state, as the dedicated kernels do: `P_epsilon`, `invBhat` and the flux vectors are
  those of `dSdE` (the closest-point ingredients, or the caches of the last cutting-plane
  iterate). Stress, state, works, `r`, `dSdE` and the closest-point results are bit-identical;
  `drdE` of the cutting-plane modes moves by the solver tolerance (7e-10 relative), `dSdT` at
  round-off. The plastic heat derivatives take the flow direction of the last constraint
  evaluation (the direction the update moved along, as EPICP_T) and the flow Hessian
  `dLambda_dsigma` evaluated for the tangent instead of re-evaluating both at the committed
  stress (1e-15 on the measured paths).
- The modular engine and the tangent assembly in natural form, as the mechanical kernels of
  2.3: no `tensor2` / `tensor4` round trip on the hot path (the mechanism caches are read as
  engineering-Voigt matrices, the flow Hessian without the Mandel scaling: tangent bits move at
  3e-16), the per-call `std::vector`s of `compute_tangent` and `run_T` are instance scratch,
  `sigma_eff_`, `sigma_start_`, `L_cur_`, `Lt_eff_` are fixed-size, `small_lu` runs its loops
  with compile-time bounds when `n == Nmax` (same operations, same order), the von Mises
  criterion evaluates its equivalent and its flow on one deviator (`flow_voigt`, the kernel of
  `eta_stress`, now in `tensor.hpp`) and the kinematic hardening works on the engineering
  components (back stresses, flow updates and backward-Euler refreshes entry by entry, the
  `<n, n>` of the hardening modulus as a weighted contraction: the stress moves at 1e-14 with
  back stresses). Every other output is bit-identical. `sim.umat_T` on 2000 points of a cyclic
  plastic path, 1 thread, mean over the plastic increments, before -> after the whole cycle:
  `MODUL` J2 power law 3.34 -> 1.72 us per point (EPICP_T 1.68 -> 1.39), `MODUL` Prager 3.94 ->
  2.07 (EPKCP_T 2.35 -> 2.07); `sim.umat` `MODUL` J2 1.83 -> 1.37, Prager 2.22 -> 1.63 (EPICP 1.0). `benchmark/dump_results.py batch_T` dumps
  `sim.umat_T` (legacy kernels and six MODUL compositions, both integrators) for the gate.
- The dissipation of the damage mechanism is the trapezoid `1/2 (Y_n + Y_n+1) dD` on the driving
  force (it was the end-point rule `Y_n+1 dD`): the quadrature of every other work of the engine,
  exact for the linear law. `Wm_d` of a damaged composition moves at the increment order; the
  stress does not. `D'(Y)`, infinite at `Y_0+` for the power-law and Weibull laws with `n < 1`,
  enters `drdE` as it enters the mechanical tangent: not regularised.
- The invariants of `contimech.hpp` (`dev`, `Mises_*`, `J2_*`, `J3_*`, `norm_*`, `eta_*`, `deta_stress`)
  and `batch_mises` share one scalar kernel on the engineering-Voigt 6-vector: `contract_voigt`,
  `dev_voigt` and `det_sym_voigt` (`small_dense.hpp`), with the shear weight of a contraction
  (`voigt_shear_weight(a, b)`), the Mises factor and the Mises kernel (`mises_voigt`) declared once
  per `Tensor2Type` next to `voigt_shear_in`, and the Mandel shear factors derived from one table of
  exponents (`mandel_shear_exponents`, `mandel_shear_product`). `J3` is the closed-form determinant
  of the deviator (no LAPACK dispatch). `eta_strain` and `eta_norm_strain` keep the strain
  convention of their input, as before, now stated in their contract. The pair of scaled copies the
  2.3 kernels carried, and `dot6`, are gone. `batch_mises` refuses an array without 6 rows.
- `iso_operator(a, b, Tensor4Type)` builds `a Ivol + b Idev` in closed form; `Ireal`, `Ivol`,
  `Idev`, `Ireal2`, `Idev2`, `L_iso`, `M_iso`, `H_iso` and the hyperelastic volumetric operator are
  instances of it and return `mat::fixed<6,6>`; `dev`, `eta_*` return `vec::fixed<6>`, `deta_stress`
  `mat::fixed<6,6>`; `auto_sym_dyadic_operator` takes a `mat::fixed<3,3>` and returns a
  `mat::fixed<6,6>` (a `mat` argument converts). The tensor members of `state_variables` (`Etot`,
  `sigma`, `F0`, `Lt`, `Wm`...) are fixed-size: assigning another size throws `std::logic_error`;
  `statev` stays dynamic. `small_lu::solve_columns_in_place` solves a matrix in place.
- Natural form of the perf kernels of 2.3: the hyperelastic invariant tangent sums its five
  operators in order, the Eshelby quadrature folds its two weights once (no compiler pragma left in
  a `.cpp`), the single-mechanism algorithmic tangent factorizes `M` with `small_lu` instead of a
  LAPACK inverse. Results move at the round-off level (below 2e-14 relative on every solver, batch
  and micromechanics table; the heat `Q`/`r` of the thermomechanical SMA laws, differences of
  large terms, up to 3e-6 of their own 1e-7 scale). `benchmark/dump_results.py compare` takes
  `--field-rtol` for such tables.

### Fixed

- The thermomechanical face of `MODUL` built `dS/dT` and `P_theta` on the raw elastic operator
  and fluxes (`-L alpha`, `kappa`) while its tangent went through the algorithmic `M^-1`: exact
  only where the flow Hessian annihilates them (isotropic blocks, J2 flows, a spherical thermal
  stress rate), 2e-4 off for an entropic hyperelastic block with kinematic hardening (its stress
  rate is deviatoric) and inexact for any anisotropic block. The thermal trial stress rate and
  the fluxes now go through `ContinuumTangent::L_tilde` and `kappa_tilde` (`M^-1 L` and
  `M^-1 kappa` of the algorithmic assembly, `L` and `kappa` of the continuum one, no solve on
  the thermal face): `dS/dT` is as exact as `dS/dE` for every composition (1e-8 against finite
  differences under the closest point; results of the isotropic compositions move at round-off).
- The heat derivatives of a plasticity mechanism composed with damage took their flow Hessian
  at the damaged stress: the turning term is now written on the effective stress
  (`StrainMechanism::DissipationInputs`), exact for a Prager back stress under damage.
- The thermomechanical face of `MODUL` reads the kinematic form (`select_umat_T`, `sim.umat_T`):
  under the total form a hyperelastic block returns `T` and `dT/dH` as the mechanical face does.
- The solver never passed `ncp` to its thermomechanical state (`MODUL` thermomechanical blocks
  ran Fischer-Burmeister whatever the run asked); `sim.umat_T` had no `ncp` argument.
- The complementarity function asked before a modular composition was configured was not
  checked against it: `ncp="min"` on an Armstrong-Frederick composition passed its first call.
- The plugin bridge copied the material name into a 5-byte `CMNAME` buffer (an overflow for
  every 5-letter name); it now fills the 80-character space-padded `CMNAME` Abaqus defines.
- The heat derivatives `drdE` / `drdT` of the plastic thermomechanical kernels (EPICP, EPKCP,
  EPJCK) held the flow direction fixed over the increment: on a non-proportional increment they
  were off by 1e-4 of the heat (the 1e-3 tolerance of the heat-contract test). The turning of the
  direction, `1/2 Dp (d eta/d sigma) sigma_n` (the end stress drops out by homogeneity), is in the
  derivative: EPICP and EPJCK are exact to
  1e-10, EPKCP to 1e-6 (its back-stress direction is still held fixed).
- The viscoelastic branches of `MODUL` were driven by the total strain: in a mechanical block with
  a temperature change (`T_final != T_init`) they relaxed the thermal strain too (free thermal
  expansion created viscous strain and a dissipation; 2 % of the stress and 10 % of the work on a
  30 K ramp), while `PRONK` drives them with `eps - alpha (T - T_init)`. The mechanisms now see the
  thermally corrected strain (`StrainMechanism::predict`, `compute_work`); isothermal results are
  bit-identical.
- The thermomechanical solver loop recomputes `DTime` from the accepted fraction after an inforced
  increment, as the mechanical loop does: the heat of an inforced increment was integrated over the
  time of the fraction the cut had refused.
- The rate (`logarithmic_R`) and total forms of a thermomechanical block are documented as one law
  on paths with fixed principal stretch directions and two laws otherwise, with the measured
  orders (tests on a coaxial stretch and on simple shear).

## [2.3.0] - 2026-11-09

Condensed stress states (#139, #141, #145), finite-strain and thermomechanical Abaqus bridges
(#143), LLDM0 ply damage (#144), `ncp` / `integrator` run options, least-squares identification,
and the performance campaign (#147 to #153: lock-free slices, fixed-size kernels, cached
homogenization, Armadillo local storage of 36 elements). An intermediate release: the
**In progress** section above lists the fixed-size and natural-form work scheduled for 2.4.

### Breaking changes

- `copy_parameters`, `apply_parameters`, `copy_constants` and `apply_constants`
  are removed. `KeyTemplate(src, params, consts).write(dst)` copies the
  templates and substitutes the keys in one call, about 5 times faster
  (templates read once, each working file written once).
- The C++ kernel `umat_damage_weibull` is removed: no UMAT name dispatched to it. The
  Weibull damage function `damage_weibull` (Functions/damage.hpp, `sim.damage_weibull`) stays.
- The C++ kernel `umat_prony_Nfast_plastic` (no header, no UMAT name) is removed: Prony
  viscoelasticity with J2 plasticity is a `ModularMaterial` (MODUL).
- `tangent_mode=3` / `"closest_point"` is refused with a `ValueError`: the local
  integrator is now chosen apart, with `integrator="closest_point"`, and
  `tangent_mode` only takes `"none"`, `"continuum"` or `"algorithmic"`.
- `integrator="closest_point"` is refused with a `ValueError` for a law without a
  closest-point update instead of running its cutting plane silently: it serves
  EPICP, EPCHA, EPJCK, MODUL and its adapters and the SMA kernels (mechanical),
  EPICP and EPKCP (thermomechanical). A Python law (PYEXT) receives
  `tangent_mode` 0, 1 or 2 only.
- The numerical option slots of the props are removed: EPICP `props[6]`
  (min-NCP) and `props[7]` (closed-form tangent), EPKCP `props[7]`, and
  `ModularMaterial(ncp=...)`. The complementarity function is the run option
  `ncp`, and the EPICP closed-form consistent tangent is always used.
- `condense="approximate"` / `"exact"` are renamed `"spatial_face"` /
  `"material_face"` (#145).
- The thermomechanical UMATs return the heat of the increment `R` (storage
  included), not a rate (#140).

### Added

- The Ansys USERMAT plugin (`UMANS`, `umat_plugin_ans_api`) is dispatched on the small-strain
  route, as UMABA, and tested by `TUMABA` on the same case. The example plugin is rebuilt on the
  `fea_transfer` component maps (now public: `fea_components`, `ansys_components`, `fe2vec`,
  `fe2mat`...): it read the tangent and wrote the rotation row-major (transposed), hid four Wm
  slots in the USERMAT's statev, and sent condensed states as 3D. `TUMANS` checks the data layout
  on a USERMAT with a non-symmetric tangent and a non-trivial rotation. The Fortran example did
  not compile in fixed form. `sim.umat` names the plugins in its error (`is_external_plugin`).
- `KeyTemplate.write(dst, test=i)`: a constant holds one value per test
  (`input_values`, e.g. the fibre angle of each specimen) and the value of test
  `i` is substituted. A constant with a single value warns: it belongs in the
  template, not in a key.
- `identification(..., method="least_squares")`: bounded trust-region
  reflective least squares, parameters scaled to [0, 1] (logarithmically for
  positive bounds), `n_starts` starting points, `n_threads` for the Jacobian
  columns or the generations of differential evolution. `calc_residuals`
  returns the residuals whose sum of squares is `calc_cost`. Differential
  evolution remains the default.
- Run options `ncp="fischer_burmeister" | "min"` and
  `integrator="cutting_plane" | "closest_point"` in `solve`, `sim.umat` and
  `sim.umat_T`. `min` applies to EPICP, EPKCP and MODUL (also next to Prony and
  damage mechanisms), not to Armstrong-Frederick back stresses.
- `simcoon.rotate_voigt`: batch rotation of Voigt stresses, strains and
  stiffness or concentration tensors by one rotation matrix per point, without
  the GIL.
- Native dS/dE of the finite-strain laws built on F, on demand.
- C++: `cubic_coefficients` and `is_isotropic(L, Tensor4Type)` (constitutive.hpp; the isotropy
  relation follows the Voigt convention of the type), `spectral_congruence` (eig_sym3.hpp),
  `add_linearop_eigsym` (contimech.hpp), replacing private copies.
- `isochoric_invariants` takes a batch, (3,3,N) tensors or (6,N) Voigt columns, as
  `isochoric_pstretch` does; both share one implementation (parallel, without the GIL).
- Condensed stress states (`nshr`, plane stress, shells, beams, trusses) in
  `sim.umat`, the solver and the Abaqus bridges, with `condensation_serves`
  (#139, #141, #145).
- Finite strain (NLGEOM) and thermomechanical finite strain in the Abaqus
  bridges and the solver, `tangent_output` for `umat_T` (#143).
- LLDM0 ply damage model (Ladevèze–Le Dantec) (#144).
- STATEV size checks and typed exceptions in the Abaqus UMATs: numerical
  failures request a step cut, input errors abort (#146).
- `slice_view`: lock-free aliases of cube slices in parallel loops (#147).
- GOH 2006 tube-inflation example for HOLZA (#142).
- `benchmark/`: batch benchmarks and bit-for-bit dumps of the solver outputs.

### Changed

- Armadillo keeps matrices of up to 36 elements in local storage (`ARMA_MAT_PREALLOC`, default 16;
  CMake cache `SIMCOON_ARMA_MAT_PREALLOC`): the 6x6 operators of the kernels no longer allocate.
  The definition changes the layout of `arma::Mat` and is exported with the `simcoon` CMake target,
  so everything linked to simcoon (bindings, UMEXT/UMABA/UMANS plugins) must be compiled with it:
  link the target through CMake, or pass `-DARMA_MAT_PREALLOC=36` by hand. Results bitwise
  unchanged.
- `kelvin_series_step` and `maxwell_parallel_step` take a `state_derivatives` flag (default
  `true`): ZENER, ZENNK and PRONK skip the thermomechanical derivatives they do not read. Results
  bitwise unchanged; the two exported symbols change.
- `pstretch`, `isochoric_pstretch` and `isochoric_pstretch_from_V/b` use the fixed-size
  `eig_sym_3x3` instead of LAPACK: OGDEN and the lazy isochoric stretches differ from the LAPACK
  result at the 1e-12 relative level.
- The fixed-size kernels that the perf code had copied file by file live once in
  `Simulation/Maths`: `small_dense.hpp` (3x3 products `mul3`, `mul3_AtB`, `mul3_ABt`, the 6-vector
  `dot6` and rank-one `outer6`, next to `mul6`/`mul66`), the new `small_voigt.hpp` (the Voigt index
  tables, `mat_from_voigt` / `voigt_from_mat` on the stack), `lanes.hpp` (two-lane arithmetic) and
  `rotation.hpp` (`is_identity_rotation`). `tensor.hpp` adds the typed forms `voigt_to_mat` /
  `mat_to_voigt` on `Tensor2Type`, used by `tensor2` itself and by the MODUL hardening. The local
  copies in `objective_rates.cpp`, `hyperelastic.cpp`, `contimech.cpp`, `eshelby.cpp`,
  `fastor_bridge.hpp`, `hardening.cpp` and `modular_umat.cpp` are gone; results bitwise unchanged.
- EPKCP thermomechanical: the cutting plane evaluates the power-law hardening through
  `PowerLawHardening` (C1 blend below p = 1e-6), as EPICP and its own closest point; for
  m < 1 the first plastic increment no longer meets an infinite slope.
- Faster kernels, results bitwise unchanged unless stated: reuse of the
  ModularUMAT instances, fixed-size 3x3 spectral solves and small dense LU,
  Eshelby integration and cached homogenization, closed-form
  Fischer-Burmeister for one and two mechanisms.
- The Chaboche, composite and MOORI identification examples use the least
  squares; the JSON-keys and analytical hyperelastic examples keep
  differential evolution.

### Fixed

- The copies of the phase state variables: the copy constructor of `state_variables` skipped
  `tangent_mode` and `ncp` (a phase copy fell back to the default tangent and to
  Fischer-Burmeister), `state_variables_M/T::operator=` and `copy_fields_M/T` re-listed the
  common fields instead of calling the base operations, and `construct()` fell through from
  `sv_type = 0` into the mechanical case: `sv_type` is now 1 (mechanical) or 2 (thermomechanical),
  any other value throws. The dead `isv.hpp` is removed. Test `Tstate_variables` (#149).
- A legacy kernel given too few state variables (`EPICP` with `nstatev = 1`, alone or as the
  sub-phase of a mean-field scheme) failed with Armadillo's bare `Mat::operator(): index out of
  bounds`: the UMAT dispatch now raises `std::invalid_argument` (`ValueError`) naming the law and
  its `nstatev`, as the modular adapter already did upfront.
- An invalid elastic convention ended the process with `exit(0)` (`M_iso`, `L_cubic`, `M_cubic`,
  `L_ortho`, `M_ortho`, `L_isotrans`, `M_isotrans`; from Python the interpreter quit with a success
  code) or, in `L_iso`, returned a zero stiffness: they now throw `std::invalid_argument`
  (`ValueError` in Python).
- Small code defects without effect on results: `assert(Y_crit(n))` out of bounds in
  `Fischer_Burmeister(_limits)` (a Debug build aborted), `variant::build` did not store `g`,
  the `state_variables_T` default constructor left `Wm`/`Wm_start` to the Armadillo default,
  an unreachable `HfR = 0.5*a2* -a3` branch of the thermomechanical unified_T closest point.
- Copies of `step_meca` / `step_thermomeca` lost their loading: the copy constructor set
  `BC_mecas` from `mecas`, the assignment skipped `BC_mecas`, `BC_Ts`, `control_type` and
  `tab_data`. The copy operations of `step`, `step_meca` and `step_thermomeca` are now the
  memberwise defaults (test `Tstep`). Nothing in simcoon copied a step.
- Under finite strain the solver ran UMEXT as a no-op (unchanged, zero stress) and failed on
  UMABA with "Unknown umat name": the plugins are now refused with an explicit error.
- `Eshelby_penny` set the in-plane shear S_2323 (Voigt 23) to 1 and S_1212 to 0: now
  S_1212 = S_1313 = 1/2 (1 in engineering Voigt) and S_2323 = 0, the flat oblate limit.
  No homogenization scheme called it.
- `d2lagrange_pow_1` placed `p0` wrongly (exact only for p0 = 1); its test now checks it
  against a finite difference of `dlagrange_pow_1`. No law called it.
- Documentation: .cpp comments cut to short implementation notes, contracts (props,
  statev, theory) moved to the header Doxygen, dead code and wrong `@file` banners removed;
  the EPCHA header example read the back-strains as the backstresses (statev 20-31).
- `isochoric_pstretch` ignored `J` on a single point (it recomputed it from the
  tensor), unlike its batch form.
- Race on the Eshelby quadrature tables between concurrent solves.
- MIHEN crash (null pointer) and MISCN failure on the second block of a
  loading path.
- `n_threads` only took effect for 1 thread: it is now honoured under GCD and
  OpenMP, and `Lt_convert` takes it too.

## [2.2.0] - 2026-10-14

Tensor bases (#126, #127, #128), cached rotation operators (#130), Eshelby memo
and lazy fields (#129), closest-point return mapping (#131, #133), Johnson-Cook
(#85), analytic SMA tangents (#136), logarithmic U (#137). See the
[release notes](https://github.com/3MAH/simcoon/releases/tag/2.2.0).

## [2.1.0] - 2026-10-01

See the [release notes](https://github.com/3MAH/simcoon/releases/tag/2.1.0).

## [2.0.1] - 2026-09-09

See the [release notes](https://github.com/3MAH/simcoon/releases/tag/2.0.1).

## [2.0.0] - 2026-09-07

See the [release notes](https://github.com/3MAH/simcoon/releases/tag/2.0.0).

[Unreleased]: https://github.com/3MAH/simcoon/compare/2.4.0...HEAD
[2.4.0]: https://github.com/3MAH/simcoon/compare/2.3.0...2.4.0
[2.3.0]: https://github.com/3MAH/simcoon/compare/2.2.0...2.3.0
[2.2.0]: https://github.com/3MAH/simcoon/compare/2.1.0...2.2.0
[2.1.0]: https://github.com/3MAH/simcoon/compare/2.0.1...2.1.0
[2.0.1]: https://github.com/3MAH/simcoon/compare/2.0.0...2.0.1
[2.0.0]: https://github.com/3MAH/simcoon/compare/2.0.0b1...2.0.0
