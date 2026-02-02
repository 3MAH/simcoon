"""
Benchmark for the new Python solver (v2.0).
Run on feature/python_solver branch.

Saves results to benchmark_results_new.npz
"""

import numpy as np
import time
from pathlib import Path

# Test cases with parameters that work reliably
BENCHMARK_CASES = {
    "ELISO": {
        "props": np.array([210000.0, 0.3, 1e-5]),
        "nstatev": 1,
        "strain_max": 0.01,
        "n_increments": 100,
        "control_type": "small_strain",
        "control": ['strain', 'stress', 'stress', 'stress', 'stress', 'stress'],
        "description": "Isotropic linear elasticity"
    },
    "EPICP": {
        # E, nu, alpha, sigmaY, k, m, kx1, Dx1, kx2, Dx2, kx3, Dx3
        "props": np.array([70000.0, 0.3, 1e-5, 300.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
        "nstatev": 14,
        "strain_max": 0.02,
        "n_increments": 100,
        "control_type": "small_strain",
        # Use fully strain-controlled for comparable results
        "control": ['strain', 'strain', 'strain', 'strain', 'strain', 'strain'],
        "description": "Plasticity with isotropic hardening"
    },
    "EPKCP": {
        # E, nu, alpha, sigmaY, k, m, kx1, Dx1, kx2, Dx2, kx3, Dx3
        "props": np.array([70000.0, 0.3, 1e-5, 300.0, 500.0, 1.0, 10000.0, 100.0, 0.0, 0.0, 0.0, 0.0]),
        "nstatev": 14,
        "strain_max": 0.03,
        "n_increments": 100,
        "control_type": "small_strain",
        # Use fully strain-controlled for comparable results
        "control": ['strain', 'strain', 'strain', 'strain', 'strain', 'strain'],
        "description": "Plasticity with kinematic hardening"
    },
    "NEOHC": {
        # E, nu, alpha (same as ELISO for Neo-Hookean in simcoon)
        "props": np.array([70000.0, 0.3, 1e-5]),
        "nstatev": 1,
        "strain_max": 0.1,
        "n_increments": 200,
        "control_type": "logarithmic",
        # Fully strain-controlled for stability in finite strain
        "control": ['strain', 'strain', 'strain', 'strain', 'strain', 'strain'],
        "description": "Neo-Hookean hyperelasticity"
    },
}

N_REPEATS = 10


def run_benchmark():
    from simcoon.solver import Solver, Block, StepMeca

    results = {}

    for case_name, cfg in BENCHMARK_CASES.items():
        print(f"\nBenchmarking {case_name}: {cfg['description']}")

        step = StepMeca(
            DEtot_end=np.array([cfg["strain_max"], 0, 0, 0, 0, 0]),
            Dsigma_end=np.zeros(6),
            control=cfg["control"],
            Dn_init=cfg["n_increments"],  # Same as old solver
            Dn_mini=cfg["n_increments"],  # Fixed (no adaptive reduction)
            Dn_inc=cfg["n_increments"],   # Fixed (no sub-stepping)
            time=1.0
        )

        block = Block(
            steps=[step],
            umat_name=case_name,
            props=cfg["props"],
            nstatev=cfg["nstatev"],
            control_type=cfg["control_type"]
        )

        solver = Solver(blocks=[block], max_iter=50, tol=1e-9)

        # Warmup
        try:
            history = solver.solve()
        except Exception as e:
            print(f"  FAILED: {e}")
            continue

        # Timing runs
        times = []
        for i in range(N_REPEATS):
            # Re-create solver for fresh state
            solver = Solver(blocks=[block], max_iter=50, tol=1e-9)
            start = time.perf_counter()
            history = solver.solve()
            elapsed = time.perf_counter() - start
            times.append(elapsed)

        strain = np.array([h.Etot[0] for h in history])
        stress = np.array([h.sigma[0] for h in history])

        results[case_name] = {
            "strain": strain,
            "stress": stress,
            "times": np.array(times),
            "n_points": len(history)
        }

        print(f"  {len(history)} points, {np.mean(times)*1000:.2f} ± {np.std(times)*1000:.2f} ms")

    # Save results
    output_file = Path(__file__).parent / "benchmark_results_new.npz"
    np.savez(output_file, **{f"{k}_{v2}": v1[v2] for k, v1 in results.items() for v2 in v1})
    print(f"\nResults saved to {output_file}")

    return results


if __name__ == "__main__":
    print("=" * 60)
    print("Simcoon v2.0 - New Python Solver Benchmark")
    print("=" * 60)
    run_benchmark()
