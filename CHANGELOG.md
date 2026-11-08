# Changelog

All notable changes to simcoon are listed here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and simcoon uses
[semantic versioning](https://semver.org/). Each pull request adds its entry
under **Unreleased**; the section takes the version number at the release.

## In progress

Known gaps and planned work, not yet in a release. An item moves to
**Unreleased** with the pull request that closes it.

- `ncp="min"` serves the mechanical EPICP, EPKCP and MODUL only: every
  thermomechanical law (`_T`, `sim.umat_T`, thermomechanical blocks) is refused
  by `check_ncp_serves` and keeps Fischer-Burmeister.
- `ncp="min"` is refused with Armstrong-Frederick back stresses under the
  cutting plane (path dependent); whether the closest point, which restarts
  from the start of the step, lifts this is not verified.
- `integrator="closest_point"` is provided with its algorithmic tangent only
  (`tangent_mode="algorithmic"`): its continuum operator is not implemented.
- Other complementarity functions to measure against Fischer-Burmeister and
  `min`, under both integrators: active set, scaled Fischer-Burmeister.
- The numerical options of the mean-field schemes (MIMTN/MISCN closed form,
  rule, Anderson, tolerance) are still props slots: they are per phase (nested
  schemes), so they need a named builder rather than a run option.
- fedoo does not pass `ncp` and `integrator` through yet.
- The external plugins (UMEXT, UMABA, UMANS) serve the small-strain route only: a finite-strain
  entry in the plugin API (law fed the logarithmic strain, returning tau, as the internal
  kernels) would let them run under NLGEOM.
- `tangent_assembly` takes per-mechanism `std::vector<arma::vec>` / `std::vector<arma::mat>`
  (kappa, dPhi/dsigma, dLambda/dsigma): every MODUL tangent copies the 6x6 flow Hessians to the heap, the
  largest remaining allocation source in MODUL plasticity. Fixed-size or view-based arguments
  need an API change of `assemble_*_tangent`.
- The perf kernels of 2.3 reproduce Armadillo's rounding so that the solver dumps stay
  bit-for-bit (`dot_voigt` with the two interleaved partial sums of `accu` in `contimech.cpp`;
  `trace`, `powmat(X, 2)` and the I2 block split around `gamma_11` in `hyperelastic.cpp`;
  `a*Ivol + b*X` entry by entry in `constitutive.cpp`; the `calG` summation order in
  `eshelby.cpp`). Once the perf pull requests are merged, these get their natural form: the
  invariants of `contimech.cpp` (Mises, J2, norms, eta flows) and `batch_mises` share one scalar
  kernel on the engineering-Voigt self-contraction, dispatched on `Tensor2Type` (shear weight 2
  for a stress, 1/2 for a strain, declared once next to `mandel_factors`), without building a
  `tensor2`; the `vdev`/`vdev2` pair disappears. The gate becomes `dump_results.py compare` /
  `regression_baseline.py --compare` with an absolute tolerance scaled to each table (Q and r
  of the SMA_T laws are differences of large terms), the zero-tolerance pass staying available
  to pin a defect. The Fastor to Armadillo rewrite starts
  under this relaxed policy. `mul6`/`mul66` keep their fixed order: they replace the BLAS, they
  do not imitate it.
- Migration of the laws to the modular engine (MODUL), state in
  `benchmark/README.md`. Routed through the adapter (code 201): ELISO, ELIST,
  ELORT, EPKCP, EPHIL, EPTRI, EPHAC, EPANI, EPDFA, EPCHG, EPHIN. To do:
  - adapters for the hyperelastic laws whose MODUL equivalent exists (NEOHC,
    MOORI, YEOHH, ISHAH, GETHH, SWANH, HOLZA; NEOHI to confirm);
  - EPICP, EPCHA and PRONK keep their legacy kernel while it is faster;
  - no MODUL equivalent yet: OGDEN, SNTVE, HYPOO, EPJCK, LLDM0, ZENER, ZENNK,
    the SMA family, the micromechanics schemes;
  - the thermomechanical dispatch (`umat_T`) has no adapter: every `_T` law
    runs its legacy kernel.

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

[Unreleased]: https://github.com/3MAH/simcoon/compare/2.2.0...HEAD
[2.2.0]: https://github.com/3MAH/simcoon/compare/2.1.0...2.2.0
[2.1.0]: https://github.com/3MAH/simcoon/compare/2.0.1...2.1.0
[2.0.1]: https://github.com/3MAH/simcoon/compare/2.0.0...2.0.1
[2.0.0]: https://github.com/3MAH/simcoon/compare/2.0.0b1...2.0.0
