"""
Compare benchmark results from old and new solvers.

Run after executing both:
- benchmark_old_solver.py (on master branch)
- benchmark_new_solver.py (on feature/python_solver branch)

Generates benchmark_comparison.md
"""

import numpy as np
from pathlib import Path

CASE_NAMES = ["ELISO", "EPICP", "EPKCP", "NEOHC"]


def load_results(filename):
    """Load benchmark results from npz file."""
    data = np.load(filename, allow_pickle=True)
    results = {}
    for case in CASE_NAMES:
        try:
            results[case] = {
                "strain": data[f"{case}_strain"],
                "stress": data[f"{case}_stress"],
                "times": data[f"{case}_times"],
                "n_points": int(data[f"{case}_n_points"])
            }
        except KeyError:
            results[case] = None
    return results


def compare_results(old_results, new_results):
    """Compare old and new solver results."""
    comparison = {}

    for case in CASE_NAMES:
        old = old_results.get(case)
        new = new_results.get(case)

        if old is None and new is None:
            comparison[case] = {"status": "both_failed"}
            continue
        elif old is None:
            comparison[case] = {"status": "old_failed"}
            continue
        elif new is None:
            comparison[case] = {"status": "new_failed"}
            continue

        # Timing comparison
        old_mean = np.mean(old["times"]) * 1000
        old_std = np.std(old["times"]) * 1000
        new_mean = np.mean(new["times"]) * 1000
        new_std = np.std(new["times"]) * 1000
        speedup = old_mean / new_mean if new_mean > 0 else 0

        # Accuracy comparison - interpolate to common strain values
        # Use the result with fewer points as reference
        if len(old["strain"]) <= len(new["strain"]):
            ref_strain = old["strain"]
            ref_stress = old["stress"]
            cmp_stress = np.interp(ref_strain, new["strain"], new["stress"])
        else:
            ref_strain = new["strain"]
            ref_stress = new["stress"]
            cmp_stress = np.interp(ref_strain, old["strain"], old["stress"])

        # Compute relative error
        max_stress = np.max(np.abs(ref_stress))
        if max_stress > 0:
            rel_error = np.max(np.abs(ref_stress - cmp_stress)) / max_stress * 100
        else:
            rel_error = 0.0

        comparison[case] = {
            "status": "success",
            "old_time_ms": old_mean,
            "old_time_std": old_std,
            "new_time_ms": new_mean,
            "new_time_std": new_std,
            "speedup": speedup,
            "max_rel_error_pct": rel_error,
            "old_n_points": old["n_points"],
            "new_n_points": new["n_points"],
        }

    return comparison


def generate_markdown(comparison):
    """Generate markdown report."""
    md = """# Simcoon Solver Benchmark Comparison

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
"""

    for case in CASE_NAMES:
        c = comparison[case]
        if c["status"] == "success":
            md += f"| {case} | {c['old_time_ms']:.2f} ± {c['old_time_std']:.2f} | {c['new_time_ms']:.2f} ± {c['new_time_std']:.2f} | {c['speedup']:.2f}x | {c['max_rel_error_pct']:.4f} |\n"
        elif c["status"] == "old_failed":
            md += f"| {case} | FAILED | - | - | - |\n"
        elif c["status"] == "new_failed":
            md += f"| {case} | - | FAILED | - | - |\n"
        else:
            md += f"| {case} | FAILED | FAILED | - | - |\n"

    md += """
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

"""

    # Summary statistics
    successful = [c for c in comparison.values() if c["status"] == "success"]
    if successful:
        avg_speedup = np.mean([c["speedup"] for c in successful])
        max_error = max([c["max_rel_error_pct"] for c in successful])

        if avg_speedup > 1:
            md += f"The new Python solver is on average **{avg_speedup:.1f}x faster** than the old C++ solver.\n"
        else:
            md += f"The old C++ solver is on average **{1/avg_speedup:.1f}x faster** than the new Python solver.\n"

        md += f"\nMaximum relative error across all tests: **{max_error:.6f}%**\n"

        if max_error < 0.01:
            md += "\nBoth solvers produce **numerically equivalent results**.\n"
    else:
        md += "No successful comparisons could be made.\n"

    return md


def main():
    benchmark_dir = Path(__file__).parent

    old_file = benchmark_dir / "benchmark_results_old.npz"
    new_file = benchmark_dir / "benchmark_results_new.npz"

    if not old_file.exists():
        print(f"Error: {old_file} not found. Run benchmark_old_solver.py on master branch first.")
        return
    if not new_file.exists():
        print(f"Error: {new_file} not found. Run benchmark_new_solver.py on feature branch first.")
        return

    print("Loading results...")
    old_results = load_results(old_file)
    new_results = load_results(new_file)

    print("Comparing...")
    comparison = compare_results(old_results, new_results)

    print("\nGenerating report...")
    md = generate_markdown(comparison)

    output_file = benchmark_dir / "benchmark_comparison.md"
    with open(output_file, 'w') as f:
        f.write(md)

    print(f"Report saved to {output_file}")
    print("\n" + md)


if __name__ == "__main__":
    main()
