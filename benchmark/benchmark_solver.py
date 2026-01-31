"""
Benchmark script for comparing old C++ solver vs new Python solver.

Run this script on both branches to generate benchmark results:
1. On feature/python_solver: python benchmark_solver.py --solver new
2. On master: python benchmark_solver.py --solver old

Results are saved to benchmark_results_<solver>.json
"""

import argparse
import json
import time
import numpy as np
from pathlib import Path
import sys
import os

# Benchmark configuration
BENCHMARK_CASES = {
    # Linear elasticity
    "ELISO": {
        "props": [210000.0, 0.3, 1e-5],  # E, nu, alpha
        "nstatev": 1,
        "strain_max": 0.02,
        "n_increments": 100,
        "control_type": 1,  # small_strain
        "description": "Isotropic linear elasticity"
    },
    "ELIST": {
        "props": [3000, 1000, 0.4, 0.3, 700, 1e-5, 1e-5],
        "nstatev": 1,
        "strain_max": 0.02,
        "n_increments": 100,
        "control_type": 1,
        "description": "Transversely isotropic elasticity"
    },
    # Plasticity (nonlinear - more iterations expected)
    "EPICP": {
        "props": [210000.0, 0.3, 0.0, 300.0, 1000.0, 0.3],
        "nstatev": 8,
        "strain_max": 0.05,
        "n_increments": 200,
        "control_type": 1,
        "description": "Plasticity with isotropic hardening"
    },
    "EPKCP": {
        "props": [210000.0, 0.3, 0.0, 300.0, 500.0, 0.3, 20000.0, 200.0],
        "nstatev": 10,
        "strain_max": 0.05,
        "n_increments": 200,
        "control_type": 1,
        "description": "Plasticity with kinematic hardening"
    },
    # Finite strain (geometric nonlinearity)
    "NEOHC": {
        "props": [1.0, 100.0],  # C1, K (bulk modulus)
        "nstatev": 1,
        "strain_max": 0.5,
        "n_increments": 100,
        "control_type": 3,  # logarithmic
        "description": "Neo-Hookean hyperelasticity (compressible)"
    },
}

# Number of repetitions for timing
N_REPEATS = 5


def run_new_solver(case_name, case_config):
    """Run benchmark using new Python solver."""
    from simcoon.solver import Solver, Block, StepMeca

    props = np.array(case_config["props"])
    nstatev = case_config["nstatev"]
    strain_max = case_config["strain_max"]
    n_increments = case_config["n_increments"]
    control_type = case_config["control_type"]

    # Map control_type integer to string
    control_type_map = {1: 'small_strain', 2: 'green_lagrange', 3: 'logarithmic'}
    control_type_str = control_type_map.get(control_type, 'small_strain')

    step = StepMeca(
        DEtot_end=np.array([strain_max, 0, 0, 0, 0, 0]),
        Dsigma_end=np.zeros(6),
        control=['strain', 'stress', 'stress', 'stress', 'stress', 'stress'],
        Dn_init=1,
        Dn_inc=n_increments,
        time=1.0
    )

    block = Block(
        steps=[step],
        umat_name=case_name,
        props=props,
        nstatev=nstatev,
        control_type=control_type_str
    )

    solver = Solver(blocks=[block], max_iter=50, tol=1e-9)

    # Time the solve (exclude import overhead)
    start = time.perf_counter()
    history = solver.solve()
    elapsed = time.perf_counter() - start

    # Extract results
    strain = np.array([h.Etot[0] for h in history])
    stress = np.array([h.sigma[0] for h in history])

    return {
        "strain": strain.tolist(),
        "stress": stress.tolist(),
        "n_increments_actual": len(history),
        "elapsed_time": elapsed
    }


def run_old_solver(case_name, case_config):
    """Run benchmark using old C++ solver (file-based)."""
    import simcoon as sim

    props = np.array(case_config["props"])
    nstatev = case_config["nstatev"]
    strain_max = case_config["strain_max"]
    n_increments = case_config["n_increments"]
    control_type = case_config["control_type"]

    # Create temporary data directory
    data_dir = Path("benchmark_data")
    results_dir = Path("benchmark_results_temp")
    data_dir.mkdir(exist_ok=True)
    results_dir.mkdir(exist_ok=True)

    # Write path file
    dn_inc = 1.0 / n_increments
    path_content = f"""#Initial_temperature
290
#Number_of_blocks
1

#Block
1
#Loading_type
1
#Control_type(NLGEOM)
{control_type}
#Repeat
1
#Steps
1

#Mode
1
#Dn_init 1.
#Dn_mini 1.
#Dn_inc {dn_inc}
#time
1
#prescribed_mechanical_state
E {strain_max}
S 0 S 0
S 0 S 0 S 0
#prescribed_temperature_state
T 290
"""

    path_file = data_dir / f"path_{case_name}.txt"
    with open(path_file, 'w') as f:
        f.write(path_content)

    # Run solver
    start = time.perf_counter()
    sim.solver(
        case_name,
        props,
        nstatev,
        0.0, 0.0, 0.0,  # Euler angles
        0,  # solver_type
        2,  # corate_type (logarithmic)
        str(data_dir),
        str(results_dir),
        f"path_{case_name}.txt",
        f"results_{case_name}.txt"
    )
    elapsed = time.perf_counter() - start

    # Read results
    output_file = results_dir / f"results_{case_name}_global-0.txt"
    data = np.loadtxt(output_file)

    # Columns: time, T, Wm, Wm_r, Wm_ir, Wm_d, Wm_d_meca, Wm_d_therm,
    #          e11, e22, e33, e12, e13, e23, s11, s22, s33, s12, s13, s23
    strain = data[:, 8]  # e11
    stress = data[:, 14]  # s11

    return {
        "strain": strain.tolist(),
        "stress": stress.tolist(),
        "n_increments_actual": len(strain),
        "elapsed_time": elapsed
    }


def run_benchmarks(solver_type):
    """Run all benchmark cases."""
    results = {
        "solver_type": solver_type,
        "n_repeats": N_REPEATS,
        "cases": {}
    }

    run_func = run_new_solver if solver_type == "new" else run_old_solver

    for case_name, case_config in BENCHMARK_CASES.items():
        print(f"\nBenchmarking {case_name}: {case_config['description']}")

        times = []
        result = None

        for i in range(N_REPEATS):
            try:
                result = run_func(case_name, case_config)
                times.append(result["elapsed_time"])
                print(f"  Run {i+1}/{N_REPEATS}: {result['elapsed_time']*1000:.2f} ms")
            except Exception as e:
                import traceback
                print(f"  Run {i+1}/{N_REPEATS}: FAILED - {e}")
                traceback.print_exc()
                continue

        if result and times:
            results["cases"][case_name] = {
                "description": case_config["description"],
                "config": case_config,
                "strain": result["strain"],
                "stress": result["stress"],
                "n_increments_actual": result["n_increments_actual"],
                "timing": {
                    "mean_ms": np.mean(times) * 1000,
                    "std_ms": np.std(times) * 1000,
                    "min_ms": np.min(times) * 1000,
                    "max_ms": np.max(times) * 1000,
                    "all_ms": [t * 1000 for t in times]
                }
            }
            print(f"  Average: {np.mean(times)*1000:.2f} ± {np.std(times)*1000:.2f} ms")
        else:
            results["cases"][case_name] = {"error": "All runs failed"}

    return results


def main():
    parser = argparse.ArgumentParser(description="Benchmark simcoon solver")
    parser.add_argument("--solver", choices=["new", "old"], required=True,
                       help="Which solver to benchmark (new=Python, old=C++ file-based)")
    args = parser.parse_args()

    print(f"=" * 60)
    print(f"Simcoon Solver Benchmark - {args.solver.upper()} solver")
    print(f"=" * 60)

    results = run_benchmarks(args.solver)

    # Save results
    output_file = Path(f"benchmark_results_{args.solver}.json")
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)

    print(f"\nResults saved to {output_file}")

    # Print summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    for case_name, case_result in results["cases"].items():
        if "error" in case_result:
            print(f"{case_name}: FAILED")
        else:
            timing = case_result["timing"]
            print(f"{case_name}: {timing['mean_ms']:.2f} ± {timing['std_ms']:.2f} ms")


if __name__ == "__main__":
    main()
