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

Removals scheduled for the next major version: the transition code of 2.x.

- The migration message for `tangent_mode=3` / `"closest_point"`: the check
  `tangent_mode in ("closest_point", 3)` of `solve` (`simcoon/solver/core.py`)
  and its C++ counterpart in `kernel_tangent_mode` (`umat_smart.cpp`). The
  value then gets the generic error of an unknown tangent mode; the internal
  kernel code 3 of the closest point is unchanged.

## [Unreleased]

### Breaking changes

- `copy_parameters`, `apply_parameters`, `copy_constants` and `apply_constants`
  are removed. `KeyTemplate(src, params, consts).write(dst)` copies the
  templates and substitutes the keys in one call, about 5 times faster
  (templates read once, each working file written once).
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

- Faster kernels, results bitwise unchanged unless stated: reuse of the
  ModularUMAT instances, fixed-size 3x3 spectral solves and small dense LU,
  Eshelby integration and cached homogenization, closed-form
  Fischer-Burmeister for one and two mechanisms.
- The Chaboche, composite and MOORI identification examples use the least
  squares; the JSON-keys and analytical hyperelastic examples keep
  differential evolution.

### Fixed

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
