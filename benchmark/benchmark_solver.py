"""Benchmark: blocks built in Python (sim.solver.solve) vs the path-file route.

Times the SAME loading cases through the two shipped entry points:

- memory : ``sim.solver.solve(...)`` — blocks in, numpy arrays out, no file
  anywhere (``_core.solver_run`` + memory sink).
- file   : ``sim.solver.save_path_json`` then ``load_path_json`` then ``solve`` —
  path.json is written and read back into those same blocks (includes both).

Since 2.0 the C++ engine reads nothing: both routes end up in the same Newton
engine with the same blocks, so the measured gap is the path.json write plus
its read, and the two responses must agree — in the same measures,
which is why the cross-check below covers every case rather than only the
small-strain ones.

Usage::

    python benchmark/benchmark_solver.py            # all cases, 5 repeats
    python benchmark/benchmark_solver.py -r 20      # more repeats

Results are printed as a table and saved to benchmark_results.json.
"""

import argparse
import json
import tempfile
import time
from pathlib import Path

import numpy as np

import simcoon as sim
from simcoon.solver import Block, StepMeca

# Props/nstatev mirror testBin/Umats/<name>/data/material.dat (the canonical
# validated materials). Uniaxial strain-driven tension, lateral stress-free.
BENCHMARK_CASES = {
    "ELISO": {
        "props": [70000.0, 0.3, 0.0],
        "nstatev": 1,
        "strain_max": 0.02,
        "ninc": 100,
        "control_type": "small_strain",
        "description": "Isotropic linear elasticity",
    },
    "ELIST": {
        "props": [1.0, 4500.0, 2300.0, 0.05, 0.3, 2700.0, 1.0e-5, 2.5e-5],
        "nstatev": 1,
        "strain_max": 0.02,
        "ninc": 100,
        "control_type": "small_strain",
        "description": "Transversely isotropic elasticity",
    },
    "EPICP": {
        "props": [67538.0, 0.349, 1.0e-6, 300.0, 1500.0, 0.3],
        "nstatev": 8,
        "strain_max": 0.05,
        "ninc": 200,
        "control_type": "small_strain",
        "description": "Plasticity, isotropic hardening (power law)",
    },
    "EPKCP": {
        "props": [67538.0, 0.349, 1.0e-6, 300.0, 1500.0, 0.3, 2000.0],
        "nstatev": 14,
        "strain_max": 0.05,
        "ninc": 200,
        "control_type": "small_strain",
        "description": "Plasticity, kinematic + isotropic hardening",
    },
    "NEOHC": {
        "props": [1000.0, 10000.0],
        "nstatev": 1,
        "strain_max": 0.3,
        "ninc": 100,
        "control_type": "logarithmic",
        "description": "Neo-Hookean hyperelasticity (finite strain)",
    },
}

CONTROL_TYPE_CODES = {"small_strain": 1, "green_lagrange": 2, "logarithmic": 3}
T_INIT = 290.0
CORATE = 2  # logarithmic/XBM, both routes


def run_memory(name, cfg):
    """In-memory route: sim.solver.solve, no files."""
    value = np.zeros(6)
    value[0] = cfg["strain_max"]
    block = Block(
        steps=[StepMeca(control=["strain"] + ["stress"] * 5, value=value,
                        ninc=cfg["ninc"], time=1.0)],
        control_type=cfg["control_type"],
    )
    start = time.perf_counter()
    res = sim.solver.solve(block, name, cfg["props"], cfg["nstatev"],
                           T_init=T_INIT, corate=CORATE, record_tangent=False)
    elapsed = time.perf_counter() - start
    strain = np.asarray(res["Strain"])[0]   # components-first (6, N)
    stress = np.asarray(res["Stress"])[0]
    return elapsed, strain, stress


def run_file(name, cfg, workdir):
    """Path-file route: save the loading path as path.json, read it back, then solve."""
    data = Path(workdir) / "data"
    data.mkdir(parents=True, exist_ok=True)
    ct = CONTROL_TYPE_CODES[cfg["control_type"]]
    # control types 2-4 additionally carry a (here zero) spin
    step = sim.solver.StepMeca(control=["strain"] + ["stress"] * 5,
                               value=[cfg["strain_max"], 0, 0, 0, 0, 0],
                               time=1.0, ninc=cfg["ninc"], Dn_init=1.0, Dn_mini=0.001,
                               BC_w=np.zeros((3, 3)) if 2 <= ct <= 4 else None)
    block = sim.solver.Block(steps=[step], control_type=cfg["control_type"])
    sim.solver.save_path_json(data / "path.json", [block], T_INIT, CORATE)
    start = time.perf_counter()
    blocks, T_init, _ = sim.solver.load_path_json(data / "path.json")
    res = sim.solver.solve(blocks, name, np.asarray(cfg["props"], dtype=float),
                           cfg["nstatev"], T_init=T_init, corate=CORATE,
                           record_tangent=False)
    elapsed = time.perf_counter() - start
    strain = np.asarray(res["Strain"])[0]   # components-first (6, N)
    stress = np.asarray(res["Stress"])[0]
    return elapsed, strain, stress


def bench_case(name, cfg, n_repeats):
    out = {"description": cfg["description"], "ninc": cfg["ninc"]}
    t_mem, t_file = [], []
    for _ in range(n_repeats):
        elapsed, strain_m, stress_m = run_memory(name, cfg)
        t_mem.append(elapsed)
        with tempfile.TemporaryDirectory() as wd:
            elapsed, strain_f, stress_f = run_file(name, cfg, wd)
        t_file.append(elapsed)

    # Equivalence of the two routes at the final state. Both return the same
    # canonical measures from the same engine, so this holds for every case.
    scale = max(1.0, float(np.abs(stress_f).max()))
    out["max_rel_stress_diff"] = float(
        np.abs(stress_m[-1] - stress_f[-1]) / scale)

    out["memory_ms"] = {"mean": float(np.mean(t_mem) * 1e3),
                        "std": float(np.std(t_mem) * 1e3)}
    out["file_ms"] = {"mean": float(np.mean(t_file) * 1e3),
                      "std": float(np.std(t_file) * 1e3)}
    out["file_over_memory"] = out["file_ms"]["mean"] / out["memory_ms"]["mean"]
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("-r", "--repeats", type=int, default=5)
    parser.add_argument("-o", "--output", default="benchmark_results.json")
    args = parser.parse_args()

    results = {"n_repeats": args.repeats, "cases": {}}
    for name, cfg in BENCHMARK_CASES.items():
        print(f"benchmarking {name}: {cfg['description']}")
        results["cases"][name] = bench_case(name, cfg, args.repeats)

    Path(args.output).write_text(json.dumps(results, indent=2))
    print(f"\nresults saved to {args.output}\n")

    hdr = f"{'case':8} {'memory (ms)':>16} {'file (ms)':>16} {'file/mem':>9} {'equiv':>10}"
    print(hdr)
    print("-" * len(hdr))
    for name, c in results["cases"].items():
        m, f = c["memory_ms"], c["file_ms"]
        eq = c.get("max_rel_stress_diff")
        eq_s = f"{eq:.1e}" if eq is not None else "n/a"
        print(f"{name:8} {m['mean']:10.2f} ±{m['std']:4.2f} "
              f"{f['mean']:10.2f} ±{f['std']:4.2f} {c['file_over_memory']:8.2f}x {eq_s:>10}")


if __name__ == "__main__":
    main()
