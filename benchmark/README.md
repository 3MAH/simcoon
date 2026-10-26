# Benchmarks and non-regression tools

| Script | Measures |
|---|---|
| `benchmark_solver.py` | the 0D solver, in-memory route against the path-file route |
| `benchmark_batch.py` | the batch APIs (N points per call) and fedoo-like iteration sequences |
| `dump_results.py` | output arrays of the solver and of the batch `umat`, compared between two builds |
| `../bench/bench_legacy_vs_modular.py` | the legacy kernels against their MODUL twins |

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
  of `bench/bench_legacy_vs_modular.py` (legacy name and MODUL twin), and of every
  `solve` call of five gallery examples (EPCHA, MODUL closest point, SMA_T,
  SMA_TR, thermomechanical SMA_T) — 1435 arrays;
- `batch`: `sim.umat` on 997 points over six loading/unloading increments, for
  ELISO, EPICP, EPKCP, MODUL and EPCHG, serial and parallel — 240 arrays.

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

`bench/bench_legacy_vs_modular.py` times a legacy UMAT against the equivalent
MODUL configuration (interleaved A/B, minimum over repeats). Only EPICP, EPCHA
and PRONK still run a legacy kernel: ELISO, ELIST, ELORT, EPKCP, EPHIL, EPHAC,
EPANI, EPDFA, EPCHG and EPHIN are routed to the modular engine by a props
adapter (code 201 in `umat_smart.cpp`), so for them the benchmark compares the
modular engine with itself. By default it runs the three true pairs; `--all`
adds the adapter families, flagged as such.

## Hotspots of the February 2026 profiling: measured state

The (superseded) PR #63 profiling listed five hotspots of the C++ engine. Their
state, measured on 2.2.0 (October 2026, `benchmark_solver.py` cases and
per-component timings of an ELISO increment, about 3 µs in Python):

| # | Hotspot (February) | Measured state | Measured cost |
|---|---|---|---|
| 1 | UMAT dispatch maps rebuilt on every call (`select_umat_*`, `umat_smart.cpp`) | **fixed**: `static const` maps | 24 ns per lookup (~1 % of an increment); rebuilt, the map would cost ~1.9 µs |
| 2 | `inv(K)` then multiply instead of `solve` in the global Newton (`solver.cpp`) | **present** | ~650 ns per iteration, 20 % of a mixed-control ELISO increment; `arma::solve` gains only 0–30 % on a 6×6 (LAPACK call overhead), a fixed-size stack solver is needed |
| 3 | Allocations inside the Newton loop | **absent** on the solver side | 0 malloc per iteration; the allocations of an increment are in the constitutive kernels (fixed-size matrices to use there) |
| 4 | Rotations applied to isotropic materials (`global2local`/`local2global`) | products already skipped; **copies remain** | ~250 ns per UMAT call |
| 5 | State-variable packing with copies (`statev_2_phases`/`phases_2_statev`) | **outside the 0D path** (Abaqus interface and plugins only) | 0 on the measured loads |

The dominant costs of a small-strain increment are now elsewhere: the modular
engine rebuilt at every call for the names routed through it (~0.9 µs of an
ELISO increment), the result sink, and `Kirchoff2PKII` evaluated with F = I.
Measure with `benchmark_solver.py` (memory route, high `ninc`) and
`benchmark_batch.py` before and after any change, and check the outputs with
`dump_results.py`.
