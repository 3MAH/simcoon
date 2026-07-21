# Solver benchmark

Times identical loading cases through the two shipped solver entry points and
cross-checks that they produce the same response:

- **memory** — `sim.solver.solve(...)`: the 2.0 API (`_core.solver_run` +
  memory sink), no filesystem involved.
- **file** — `sim._core.solver(...)`: the legacy path-file driver, including
  its path.txt write and results-file parse.

Both wrap the same C++ Newton engine (`solver_run`), so the timing gap
isolates the file round-trip overhead, and the equivalence column is a
file-vs-memory regression check.

```bash
python benchmark/benchmark_solver.py           # all cases, 5 repeats
python benchmark/benchmark_solver.py -r 20     # tighter statistics
```

Cases cover linear isotropic/anisotropic elasticity (ELISO, ELIST), isotropic
and kinematic plasticity (EPICP, EPKCP) and finite-strain hyperelasticity
(NEOHC), with props mirroring the validated `testBin` materials.

## Known optimization opportunities in the C++ engine

Profiling done on the (superseded) PR #63 branch identified hotspots that are
still present in the current engine — verified against today's sources. The
speedup figures are the original (Feb 2026) profiling estimates, to be
re-measured with this benchmark before acting on any of them:

1. **UMAT dispatch maps rebuilt on every call** — `std::map` literals
   `list_umat` are function-local in each `select_umat_*` dispatcher
   (`src/Continuum_mechanics/Umat/umat_smart.cpp`), reconstructed once per
   increment per iteration (est. 2–3×on small models). `static const` maps
   (or a shared registry) would build them once.
2. **`inv(K)` then multiply instead of `solve`** — the global Newton uses
   `inv(invK, K); Delta = -invK * residual` (`src/Simulation/Solver/solver.cpp`).
   `arma::solve(K, residual)` is cheaper and better conditioned (est. 2–4× on
   the 6×6 solve itself).
3. **Allocations inside the Newton loop** — `zeros(6,6)` / temporary vectors
   recreated per iteration rather than reused (est. 2–3×).
4. **Rotations applied for isotropic materials** — `global2local`/`local2global`
   run regardless of whether the phase orientation is identity (est. 2–3× on
   multiphase models).
5. **State-variable packing with copies** — `statev_2_phases`/`phases_2_statev`
   copy subvectors per phase per call rather than aliasing (est. 2×).

These numbers compound only in the per-increment hot path; measure with this
benchmark (memory route, high `ninc`) before and after any change.
