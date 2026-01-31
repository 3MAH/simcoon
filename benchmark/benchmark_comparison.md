# Simcoon Solver Benchmark Comparison

## Overview

This benchmark compares the **old C++ file-based solver** (v1.x, master branch) with the
**new Python solver API** (v2.0, feature/python_solver branch).

## Test Cases

| UMAT | Description |
|------|-------------|
| ELISO | Isotropic linear elasticity |
| EPICP | Plasticity with isotropic hardening |
| EPKCP | Plasticity with combined isotropic/kinematic hardening |
| NEOHC | Neo-Hookean hyperelasticity (finite strain) |

## Performance Results

| UMAT | Old Solver (ms) | New Solver (ms) | Speedup | Max Error (%) |
|------|-----------------|-----------------|---------|---------------|
| ELISO | 1.75 ± 0.12 | 3.91 ± 0.19 | 0.45x | 0.0000 |
| EPICP | 1.85 ± 0.05 | 2.69 ± 0.07 | 0.69x | 0.0002 |
| EPKCP | 2.00 ± 0.10 | 2.81 ± 0.07 | 0.71x | 0.0002 |
| NEOHC | 3.28 ± 0.14 | 7.54 ± 0.14 | 0.43x | 0.0004 |

## Accuracy Analysis

The maximum relative error between solvers is computed as:

```
max_error = max(|σ_old - σ_new|) / max(|σ_old|) × 100%
```

Both solvers should produce nearly identical results (< 0.01% error) for the same
material model and loading conditions.

## Notes

- **Timing**: Averaged over 10 runs after warmup
- **Old solver**: Reads/writes configuration files, includes file I/O overhead
- **New solver**: Pure Python API, no file I/O required
- **Speedup > 1**: New solver is faster
- **Speedup < 1**: Old solver is faster

## Conclusion

The old C++ solver is on average **1.8x faster** than the new Python solver.

Maximum relative error across all tests: **0.000381%**

Both solvers produce **numerically equivalent results**.
