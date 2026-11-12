# Benchmarks and non-regression tools

> **Temporary folder.** These scripts serve two ongoing efforts and will be
> removed once both are done; nothing in the library or the test suite uses them.
>
> 1. **Computation time.** Each optimization step is checked for its outputs
>    (`dump_results.py`: bit-identical, or its worst differences published) and
>    for its speed (`benchmark_batch.py`, `benchmark_solver.py`, interleaved A/B
>    against the previous step).
> 2. **Moving the laws to the modular engine, one at a time.** A legacy name is
>    routed to MODUL through a props adapter (code 201 in `umat_smart.cpp`) once
>    `dump_results.py` shows the same response and `bench_legacy_vs_modular.py`
>    shows no loss of speed. The table at the end tracks where each law stands.

| Script | Measures |
|---|---|
| `benchmark_solver.py` | the 0D solver, in-memory route against the path-file route |
| `benchmark_batch.py` | the batch APIs (N points per call) and fedoo-like iteration sequences |
| `dump_results.py` | output arrays of the solver and of the batch `umat` / `umat_T`, compared between two builds |
| `bench_legacy_vs_modular.py` | the legacy kernels against their MODUL twins |

## Solver benchmark

Times identical loading cases through the two shipped solver entry points and
cross-checks that they produce the same response:

- **memory** — `sim.solver.solve(...)`: blocks built in Python
  (`_core.solver_run` + memory sink), no filesystem involved.
- **file** — `sim.solver.save_path_json` / `load_path_json` then `solve`: the
  path-file route, including the path.json write and read.

Since 2.0 the C++ engine reads nothing, so both routes reach the same Newton
engine with the same blocks: the timing gap isolates the path.json write and
its read, and the equivalence column checks the JSON round-trip.

```bash
python benchmark/benchmark_solver.py           # all cases, 5 repeats
python benchmark/benchmark_solver.py -r 20     # tighter statistics
```

Cases cover linear isotropic/anisotropic elasticity (ELISO, ELIST), isotropic
and kinematic plasticity (EPICP, EPKCP) and finite-strain hyperelasticity
(NEOHC), with props mirroring the validated `testBin` materials.

## Batch API benchmark

`benchmark_batch.py` times the batch Python APIs on N material points per call
(`sim.umat`, `sim.umat_T`, `objective_rate`, `Lt_convert`, `stress_convert`,
`Log_strain`, `Rotation`, `Tensor2`/`Tensor4`), and two sequences that mirror
one fedoo global iteration: `fedooSS_<law>` (small strain, the `umat` call as
fedoo makes it) and `fedooUL_<law>` (updated Lagrangian: `objective_rate`
log_R, frame transport, `umat` with the box tangent, two `Lt_convert`).

```bash
python benchmark/benchmark_batch.py run -l master -o master.json
python benchmark/benchmark_batch.py run -l branch -o branch.json --N 1e4,1e5
python benchmark/benchmark_batch.py summary master.json branch.json
```

Each point is the minimum wall time per call over adaptive repetitions. The
summary fits `t = a + b·N` on the serial runs: `a` is the fixed cost of a call,
`b` the cost of one point (the figure a kernel optimization should move); the
parallel runs (`n_threads=0`) are tabulated per point. For an A/B comparison,
interleave runs of the two builds and pass all the files to `summary`.

## Non-regression dumps

`dump_results.py` saves every output array of a fixed set of computations and
compares two such dumps:

- `solver`: all result fields of the `benchmark_solver.py` cases, of each family
  of `bench_legacy_vs_modular.py` (legacy name and MODUL twin), and of every
  `solve` call of five gallery examples (EPCHA, MODUL closest point, SMA_T,
  SMA_TR, thermomechanical SMA_T) — 1435 arrays;
- `batch`: `sim.umat` on 997 points over six loading/unloading increments, for
  ELISO, EPICP, EPKCP, MODUL and EPCHG, serial and parallel — 240 arrays;
- `batch_T`: `sim.umat_T` on 251 points over six increments with a temperature
  ramp, for EPICP_T, EPKCP_T and six MODUL compositions (J2, Prager, Chaboche,
  Prony + Prager, damage, entropic neo-Hookean), cutting plane and closest
  point: stress, state, works, heat and its four derivatives.

```bash
python benchmark/dump_results.py solver ref.npz      # reference build
python benchmark/dump_results.py solver new.npz      # build under test
python benchmark/dump_results.py compare ref.npz new.npz             # bit for bit
python benchmark/dump_results.py compare ref.npz new.npz --rtol 1e-12
```

`compare` exits with status 1 on any difference beyond `--rtol` (default 0) and
lists the largest relative difference per case. A performance change that keeps
the floating-point operations (caching, layout, dead code) must compare
bit-identical; one that reorders them (another linear solver, a closed form)
publishes its worst differences.

## Legacy kernels against their MODUL twins

`bench_legacy_vs_modular.py` times a dedicated (closed-form) kernel against the
equivalent MODUL configuration (interleaved A/B, minimum over repeats). Three
dedicated kernels remain, under the names `EPICP_closed`, `EPCHA_closed` and
`PRONK_closed`; every legacy name, these three included, is routed to the
modular engine by a props adapter (code 201 in `umat_smart.cpp`), so for the
other families (ELISO, ELIST, ELORT, EPKCP, EPHIL, EPHAC, EPANI, EPDFA, EPCHG,
EPHIN) the benchmark compares the modular engine with itself. By default it
runs the three true pairs; `--all` adds the adapter families, flagged as such.

## Hotspots of the February 2026 profiling: current state

The (superseded) PR #63 profiling listed five hotspots of the C++ engine. Their
state was measured on 2.2.0 (October 2026, `benchmark_solver.py` cases and
per-component timings of an ELISO increment, about 3 µs in Python), then updated
with the optimization lots that followed:

| # | Hotspot (February) | State | Measured cost (2.2.0) |
|---|---|---|---|
| 1 | UMAT dispatch maps rebuilt on every call (`select_umat_*`, `umat_smart.cpp`) | **fixed**: `static const` maps | 24 ns per lookup (~1 % of an increment); rebuilt, the map would cost ~1.9 µs |
| 2 | `inv(K)` then multiply instead of `solve` in the global Newton (`solver.cpp`) | **fixed** in the mixed-control Newton: `newton_correction`, a fixed-size LU on the stack (`small_lu<7>`), no LAPACK call; the RNL path (`solver_type` 1) still uses `inv(K)` | ~650 ns per iteration, 20 % of a mixed-control ELISO increment; `arma::solve` gains only 0–30 % on a 6×6 (LAPACK call overhead) |
| 3 | Allocations inside the Newton loop | **absent** on the solver side | 0 malloc per iteration; the allocations of an increment are in the constitutive kernels (fixed-size matrices to use there) |
| 4 | Rotations applied to isotropic materials (`global2local`/`local2global`) | products already skipped; **copies remain** | ~250 ns per UMAT call |
| 5 | State-variable packing with copies (`statev_2_phases`/`phases_2_statev`) | **outside the 0D path** (Abaqus interface and plugins only) | 0 on the measured loads |

The costs that dominated a small-strain increment next have been removed as
well: the modular engine is no longer rebuilt at every call (a per-thread cache
of the `ModularUMAT` instance; it cost ~0.9 µs of an ELISO increment), the result
sink reserves its records up front, and `Kirchoff2PKII` is skipped when F = I
exactly. What remains is in the constitutive kernels: heap temporaries and the
LAPACK inverse of the algorithmic tangent. Measure with `benchmark_solver.py` (memory route, high `ninc`) and
`benchmark_batch.py` before and after any change, and check the outputs with
`dump_results.py`.

## Migration of the laws to the modular engine

State of both dispatches (`umat_smart.cpp`): the thermomechanical one
(`umat_T`) serves `MODUL` by name and the legacy names EPICP, EPKCP and PRONK
through the `_T` adapter, their dedicated kernels answering the `_closed`
names. Measured on `sim.umat_T`, 2000 points, 60 increments of a cyclic plastic
path, 1 thread, median over the plastic increments, `MODUL` costs 1.86 us per
point against 1.46 for EPICP_closed (1.28x) and 2.30 against 2.16 for
EPKCP_closed (1.07x) under the cutting plane; under the closest point 4.9
against 5.1 and 7.5 against 8.1 (0.9x, both on the generic solver) (Nov 2026,
after the natural-form lots of the modular engine and the perf pass of the
thermomechanical faces: the cycle took `MODUL` from 3.34 / 3.94 us and the
legacy thermomechanical kernels from 1.68 / 2.35 through the shared tangent
assembly and `mul6_At`). The thermal face is not
where the time goes: it costs about 0.3 us on top of the mechanical `MODUL` run
(0.4 us for EPICP_T on top of EPICP, whose kernel is the closed radial-return
tangent), since it reads the multiplier assembly `compute_tangent` kept; the
J2 gap is the mechanical engine, 1.37 us against 1.0 for EPICP on the same
increments: the generic algorithmic assembly (`small_lu` of
`M = I + Dp L dLambda/dsigma`) that EPICP replaces by its closed form, which
MODUL keeps generic by design. Profile with `sample` on a `sim.umat_T` loop;
the min-over-increments estimator picks the elastic increments and hides this.

| Status | Laws |
|---|---|
| Routed to MODUL (adapter, code 201) | ELISO, ELIST, ELORT, EPICP, EPKCP, EPCHA, PRONK, EPHIL, EPTRI, EPHAC, EPANI, EPDFA, EPCHG, EPHIN (thermomechanical face: EPICP, EPKCP, PRONK) |
| MODUL equivalent exists, adapter to write | NEOHC, MOORI, YEOHH, ISHAH, GETHH, SWANH, HOLZA (`HYPER_INVARIANTS` block, same `HyperPotential`); NEOHI to confirm (neo-Hookean with the quadratic volumetric term, kappa = 2/D1) |
| Dedicated closed-form kernel kept under a `_closed` name (benchmark, comparison) | EPICP_closed, EPCHA_closed, PRONK_closed; EPKCP_closed on the thermomechanical face |
| MODUL equivalent by parameter map, adapter to write | ZENER (one Maxwell branch: `X_1^M = X_0^2/(X_0 + X_1)`, `eta^M = eta^K X_0^2/(X_0 + X_1)^2` per isotropic mode, see `docs/simulation/modular_python.rst`) |
| No MODUL equivalent | OGDEN (principal stretches), SNTVE, HYPOO, EPJCK, LLDM0 (anisotropic damage), ZENNK (Kelvin chains, N >= 2), SMA family, micromechanics (MIHEN, MIMTN, MISCN, MIPLN) |
