"""
Benchmark for the old C++ file-based solver (v1.x).
Run on master branch.

Saves results to benchmark_results_old.npz
"""

import numpy as np
import time
import os
from pathlib import Path
import shutil

# Test cases - same as new solver benchmark
BENCHMARK_CASES = {
    "ELISO": {
        "props": np.array([210000.0, 0.3, 1e-5]),
        "nstatev": 1,
        "strain_max": 0.01,
        "n_increments": 100,
        "control_type": 1,  # small_strain
        "control_str": "E {strain}\nS 0 S 0\nS 0 S 0 S 0",
        "description": "Isotropic linear elasticity"
    },
    "EPICP": {
        # E, nu, alpha, sigmaY, k, m, kx1, Dx1, kx2, Dx2, kx3, Dx3
        "props": np.array([70000.0, 0.3, 1e-5, 300.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
        "nstatev": 14,
        "strain_max": 0.02,
        "n_increments": 100,
        "control_type": 1,
        # Fully strain-controlled for comparable results
        "control_str": "E {strain}\nE 0 E 0\nE 0 E 0 E 0",
        "description": "Plasticity with isotropic hardening"
    },
    "EPKCP": {
        # E, nu, alpha, sigmaY, k, m, kx1, Dx1, kx2, Dx2, kx3, Dx3
        "props": np.array([70000.0, 0.3, 1e-5, 300.0, 500.0, 1.0, 10000.0, 100.0, 0.0, 0.0, 0.0, 0.0]),
        "nstatev": 14,
        "strain_max": 0.03,
        "n_increments": 100,
        "control_type": 1,
        # Fully strain-controlled for comparable results
        "control_str": "E {strain}\nE 0 E 0\nE 0 E 0 E 0",
        "description": "Plasticity with kinematic hardening"
    },
    "NEOHC": {
        # E, nu, alpha
        "props": np.array([70000.0, 0.3, 1e-5]),
        "nstatev": 1,
        "strain_max": 0.1,
        "n_increments": 200,
        "control_type": 3,  # logarithmic
        # Fully strain-controlled for stability
        "control_str": "E {strain}\nE 0 E 0\nE 0 E 0 E 0",
        "description": "Neo-Hookean hyperelasticity"
    },
}

N_REPEATS = 10


def create_path_file(data_dir, case_name, strain_max, n_increments, control_type, control_str):
    """Create path.txt file for old solver."""
    dn_inc = 1.0 / n_increments
    # Format control string with actual strain value
    meca_state = control_str.format(strain=strain_max)
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
{meca_state}
#prescribed_temperature_state
T 290
"""
    path_file = data_dir / f"path_{case_name}.txt"
    with open(path_file, 'w') as f:
        f.write(path_content)
    return f"path_{case_name}.txt"


def run_benchmark():
    import simcoon as sim

    # Setup directories
    benchmark_dir = Path(__file__).parent
    data_dir = benchmark_dir / "data_old"
    results_dir = benchmark_dir / "results_old"

    data_dir.mkdir(exist_ok=True)
    results_dir.mkdir(exist_ok=True)

    results = {}

    for case_name, cfg in BENCHMARK_CASES.items():
        print(f"\nBenchmarking {case_name}: {cfg['description']}")

        # Create path file
        path_file = create_path_file(
            data_dir, case_name,
            cfg["strain_max"], cfg["n_increments"], cfg["control_type"], cfg["control_str"]
        )

        # Warmup run
        try:
            sim.solver(
                case_name,
                cfg["props"],
                cfg["nstatev"],
                0.0, 0.0, 0.0,  # Euler angles
                0,  # solver_type
                2,  # corate_type
                str(data_dir),
                str(results_dir),
                path_file,
                f"results_{case_name}.txt"
            )
        except Exception as e:
            print(f"  FAILED: {e}")
            continue

        # Timing runs
        times = []
        for i in range(N_REPEATS):
            start = time.perf_counter()
            sim.solver(
                case_name,
                cfg["props"],
                cfg["nstatev"],
                0.0, 0.0, 0.0,
                0, 2,
                str(data_dir),
                str(results_dir),
                path_file,
                f"results_{case_name}.txt"
            )
            elapsed = time.perf_counter() - start
            times.append(elapsed)

        # Read results
        output_file = results_dir / f"results_{case_name}_global-0.txt"
        data = np.loadtxt(output_file)

        # Columns: time, T, Wm, Wm_r, Wm_ir, Wm_d, Wm_d_meca, Wm_d_therm,
        #          e11, e22, e33, e12, e13, e23, s11, s22, s33, s12, s13, s23
        strain = data[:, 8]   # e11
        stress = data[:, 14]  # s11

        results[case_name] = {
            "strain": strain,
            "stress": stress,
            "times": np.array(times),
            "n_points": len(strain)
        }

        print(f"  {len(strain)} points, {np.mean(times)*1000:.2f} ± {np.std(times)*1000:.2f} ms")

    # Save results
    output_file = benchmark_dir / "benchmark_results_old.npz"
    np.savez(output_file, **{f"{k}_{v2}": v1[v2] for k, v1 in results.items() for v2 in v1})
    print(f"\nResults saved to {output_file}")

    # Cleanup
    shutil.rmtree(data_dir, ignore_errors=True)
    shutil.rmtree(results_dir, ignore_errors=True)

    return results


if __name__ == "__main__":
    print("=" * 60)
    print("Simcoon v1.x - Old C++ Solver Benchmark")
    print("=" * 60)
    run_benchmark()
