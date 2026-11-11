"""Non-regression dumps: run a fixed set of computations, save every output array,
and compare two dumps bit for bit (or within a tolerance).

Two dumps:

- ``solver`` — every result field (``scalar_data`` and ``field_data``) of:

  * the ``benchmark_solver.py`` cases (with ``record_tangent=True``);
  * each family of ``benchmark/bench_legacy_vs_modular.py``, under its legacy
    name and as its MODUL twin, on ``examples/data/MODUL_path.json``;
  * every ``sim.solver.solve`` call made by a few gallery examples (EPCHA,
    MODUL closest point, SMA_T, SMA_TR, thermomechanical SMA_T): the examples are
    run as scripts with plotting disabled and each solve's result is recorded.

- ``batch`` — ``sim.umat`` on (6, N) arrays (N = 997 random directions) over
  six increments, loading then unloading, for ELISO, EPICP, EPKCP, MODUL and
  EPCHG, serial (``n_threads=1``) and parallel (``n_threads=0``): sigma, statev,
  Wm and Lt after each increment.

Usage::

    python benchmark/dump_results.py solver ref_solver.npz   # with the reference build
    python benchmark/dump_results.py batch  ref_batch.npz
    ...                                                       # switch builds
    python benchmark/dump_results.py solver new_solver.npz
    python benchmark/dump_results.py compare ref_solver.npz new_solver.npz
    python benchmark/dump_results.py compare ref_batch.npz new_batch.npz --rtol 1e-12

``compare`` prints, per case, either BIT-IDENTICAL or the largest relative
difference ``max|a - b| / max|a|`` over its arrays, and exits with status 1 if an
array is missing, changes shape, or differs by more than ``--rtol`` (default 0:
bit-identical). A change that is meant to be bit-identical (memory layout,
caching, dead code) must pass with the default; a change of floating-point
operation order (other linear solver, closed forms) lists its worst differences
with ``--rtol``.

The build compared is whatever ``import simcoon`` resolves: dump with one build,
switch ``PYTHONPATH`` (or reinstall), dump again. Keep the two runs on the same
machine: LAPACK/BLAS and the compiler's FMA contraction differ across platforms.
"""

import argparse
import contextlib
import io
import os
import runpy
import sys
import warnings
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
EXAMPLES = REPO / "examples"

#: Gallery examples whose solver calls enter the ``solver`` dump (they cover the
#: legacy EPCHA kernel, the modular closest point return, and the SMA laws).
DEFAULT_EXAMPLES = (
    "mechanical/EPCHA.py",
    "mechanical/MODUL_closest_point.py",
    "mechanical/SMA_T.py",
    "mechanical/SMA_TR.py",
    "thermomechanical/SMA_T.py",
)


# --------------------------------------------------------------------------- solver


def _fields(res):
    """All numeric result fields of a SolverResults, as float arrays."""
    out = {}
    for k, v in list(res.scalar_data.items()) + list(res.field_data.items()):
        try:
            out[k] = np.asarray(v, dtype=float)
        except (TypeError, ValueError):
            pass
    return out


def _run_example(path):
    """Run one example script; return the results of its ``sim.solver.solve`` calls."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import simcoon as sim
    import simcoon.solver as ssolver
    from simcoon.solver import core as score

    results = []
    orig = score.solve

    def recording_solve(*args, **kwargs):
        res = orig(*args, **kwargs)
        results.append(res)
        return res

    saved = (ssolver.solve, score.solve, plt.show, plt.savefig, os.getcwd())
    ssolver.solve = score.solve = sim.solver.solve = recording_solve
    plt.show = plt.savefig = lambda *a, **k: None
    os.chdir(path.parent)
    try:
        with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
            warnings.simplefilter("ignore")
            runpy.run_path(str(path), run_name="__main__")
    except Exception as e:   # a failure after the solves (plotting...) is not ours to judge
        print(f"  {path.name} raised {type(e).__name__}: {str(e)[:120]}")
    finally:
        ssolver.solve = score.solve = sim.solver.solve = saved[0]
        plt.show, plt.savefig = saved[2], saved[3]
        os.chdir(saved[4])
        plt.close("all")
    return results


def solver_cases():
    """(name, thunk) pairs; each thunk returns a SolverResults."""
    import simcoon as sim
    from simcoon.solver import Block, StepMeca

    sys.path.insert(0, str(REPO / "benchmark"))
    import benchmark_solver as bs
    import bench_legacy_vs_modular as blm

    cases = []
    for name, cfg in bs.BENCHMARK_CASES.items():
        value = np.zeros(6)
        value[0] = cfg["strain_max"]
        blk = Block(steps=[StepMeca(control=["strain"] + ["stress"] * 5, value=value,
                                    ninc=cfg["ninc"], time=1.0)],
                    control_type=cfg["control_type"])
        cases.append((f"benchmark/{name}", lambda blk=blk, name=name, cfg=cfg: sim.solver.solve(
            blk, name, cfg["props"], cfg["nstatev"], T_init=290.0, corate=2,
            record_tangent=True)))

    blocks, T_init, _ = sim.solver.load_path_json(str(EXAMPLES / "data" / "MODUL_path.json"))
    for fam, (lname, lprops, lnst, _path, mat) in blm.FAMILIES.items():
        cases.append((f"lvm/{fam}/legacy", lambda lname=lname, lprops=lprops, lnst=lnst:
                      sim.solver.solve(blocks, lname, lprops, lnst, T_init=T_init,
                                       record_tangent=True)))
        cases.append((f"lvm/{fam}/modul", lambda mat=mat: sim.solver.solve(
            blocks, "MODUL", mat.props, mat.nstatev, T_init=T_init, record_tangent=True)))
    return cases


def dump_solver(out, examples):
    store = {}

    def keep(name, res):
        for k, v in _fields(res).items():
            store[f"{name}|{k}"] = v

    for name, thunk in solver_cases():
        try:
            keep(name, thunk())
        except Exception as e:
            print(f"FAIL {name}: {type(e).__name__}: {str(e)[:120]}")
            store[f"{name}|__error__"] = np.array([1.0])
    for rel in examples:
        path = (EXAMPLES / rel).resolve()
        tag = "_".join(Path(rel).with_suffix("").parts)
        for i, res in enumerate(_run_example(path)):
            keep(f"ex/{tag}/{i}", res)
    np.savez(out, **store)
    print(f"saved {len(store)} arrays to {out}")


# --------------------------------------------------------------------------- batch


def dump_batch(out, N=997, n_inc=6):
    import simcoon as sim
    sys.path.insert(0, str(REPO / "benchmark"))
    from benchmark_batch import modular_props

    F = np.asfortranarray
    laws = {
        "ELISO": ([200000.0, 0.3, 1e-5], 1),
        "EPICP": ([200000.0, 0.3, 1e-5, 300.0, 1000.0, 0.5], 8),
        "EPKCP": ([200000.0, 0.3, 0.0, 300.0, 1000.0, 0.5, 10000.0], 14),
        "MODUL": modular_props(),
        "EPCHG": ([200000.0, 0.3, 80000.0, 1e-5, 300.0, 2, 2, 0, 100.0, 10.0, 50.0, 5.0,
                   20000.0, 200.0, 5000.0, 50.0], 20),
    }
    rng = np.random.default_rng(3)
    store = {}
    for name, (props, nstatev) in laws.items():
        for nt in (1, 0):
            props_a = F(np.asarray(props, float).reshape(-1, 1))
            etot = F(np.zeros((6, N)))
            sigma = F(np.zeros((6, N)))
            statev = F(np.zeros((nstatev, N)))
            Wm = F(np.zeros((4, N)))
            DR = F(np.repeat(np.eye(3)[:, :, None], N, axis=2))
            direction = rng.standard_normal((6, N))
            for inc in range(n_inc):
                # 4 loading increments past yield, then unloading
                Detot = F(2e-3 * direction * (1 if inc < 4 else -1))
                try:
                    sigma, statev, Wm, Lt = sim.umat(
                        name, etot, Detot, np.array([]), np.array([]), sigma, DR, props_a,
                        statev, inc * 0.1, 0.1, Wm, temp=np.array([290.0]), n_threads=nt,
                        start=(inc == 0))
                except Exception as e:
                    print(f"FAIL {name} nt={nt} inc {inc}: {type(e).__name__}: {str(e)[:100]}")
                    store[f"batch/{name}/nt{nt}/inc{inc}|__error__"] = np.array([1.0])
                    break
                etot = F(etot + Detot)
                for k, v in (("sigma", sigma), ("statev", statev), ("Wm", Wm), ("Lt", Lt)):
                    store[f"batch/{name}/nt{nt}/inc{inc}|{k}"] = np.array(v)
    np.savez(out, **store)
    print(f"saved {len(store)} arrays to {out}")


# --------------------------------------------------------------------------- compare


def compare(path_a, path_b, rtol=0.0, verbose=False, field_rtol=None):
    """Compare two dumps; return True when every array matches within rtol.

    *field_rtol* maps a field name (the part after '|') to its own tolerance: the heat tables
    Q and r of the thermomechanical laws are differences of large terms.

    Prints the cases that differ (all of them with *verbose*) and a summary line.
    """
    a, b = np.load(path_a), np.load(path_b)
    field_rtol = field_rtol or {}
    keys = sorted(set(a.files) | set(b.files))
    worst = {}     # case -> (relative difference, field)
    ok = True
    n_exact = n_cmp = 0
    for k in keys:
        case, _, field = k.rpartition("|")
        if k not in a.files or k not in b.files:
            print(f"MISSING in {'b' if k in a.files else 'a'}: {k}")
            ok = False
            continue
        x, y = a[k], b[k]
        if x.shape != y.shape:
            print(f"SHAPE {k}: {x.shape} vs {y.shape}")
            worst[case] = (np.inf, field)
            ok = False
            continue
        n_cmp += 1
        if np.array_equal(x, y, equal_nan=True):
            n_exact += 1
            worst.setdefault(case, (0.0, ""))
            continue
        scale = max(float(np.nanmax(np.abs(x))) if x.size else 0.0, 1e-300)
        rel = float(np.nanmax(np.abs(x - y)) / scale)
        if np.isnan(rel):   # NaN on one side only
            rel = np.inf
        if rel > worst.get(case, (0.0, ""))[0]:
            worst[case] = (rel, field)
        if rel > field_rtol.get(field, rtol):
            ok = False
    for case, (rel, field) in sorted(worst.items()):
        if rel == 0 and not verbose:
            continue
        status = "BIT-IDENTICAL" if rel == 0 else f"max rel diff {rel:.2e} ({field})"
        print(f"  {case:55s} {status}")
    print(f"{n_exact}/{n_cmp} arrays bit-identical, {len(keys)} keys"
          f"{'' if rtol == 0 else f', tolerance {rtol:g}'}: {'PASS' if ok else 'FAIL'}")
    return ok


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("solver", help="dump the 0D solver results")
    p.add_argument("output")
    p.add_argument("--examples", default=",".join(DEFAULT_EXAMPLES),
                   help="comma-separated example scripts, relative to examples/ ('' for none)")
    p = sub.add_parser("batch", help="dump the batch sim.umat results")
    p.add_argument("output")
    p = sub.add_parser("compare", help="compare two dumps")
    p.add_argument("a")
    p.add_argument("b")
    p.add_argument("--rtol", type=float, default=0.0,
                   help="largest relative difference accepted (default 0: bit-identical)")
    p.add_argument("-v", "--verbose", action="store_true", help="also list identical cases")
    p.add_argument("--field-rtol", default="",
                   help="per-field tolerances, e.g. Q=1e-5,r=1e-5 (fields not listed use --rtol)")
    args = parser.parse_args()

    if args.command == "solver":
        dump_solver(args.output, [e for e in args.examples.split(",") if e])
    elif args.command == "batch":
        dump_batch(args.output)
    else:
        field_rtol = {}
        for item in filter(None, args.field_rtol.split(",")):
            name, sep, value = item.partition("=")
            try:
                field_rtol[name] = float(value)
            except ValueError:
                parser.error(f"--field-rtol expects NAME=TOL items, got {item!r}")
            if not sep:
                parser.error(f"--field-rtol expects NAME=TOL items, got {item!r}")
        sys.exit(0 if compare(args.a, args.b, args.rtol, args.verbose, field_rtol) else 1)


if __name__ == "__main__":
    main()
