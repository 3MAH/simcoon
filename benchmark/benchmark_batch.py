"""Benchmark: the batch (N material points per call) Python APIs of simcoon.

Times the calls a finite-element code such as fedoo makes once per global
iteration on all its integration points: ``sim.umat``/``sim.umat_T`` on
(6, N) arrays, ``objective_rate``, ``Lt_convert``, ``stress_convert``,
``Log_strain``, the rotation helpers, and two composite sequences:

- ``fedooSS_<law>`` — small strain: the ``umat`` call alone, with the
  arguments fedoo passes (``F0 = F1 = []``, ``DR = I``, ``temp``);
- ``fedooUL_<law>`` — updated Lagrangian, one global iteration:
  ``objective_rate("log_R")``, ``Rotation.from_matrix`` + ``apply_strain``/
  ``apply_stress``, ``umat`` with ``corate=3`` and the box tangent, then the two
  ``Lt_convert`` (``DsigmaDe_2_DSDE``, ``DSDE_2_Dsigma_LieDD``).

Every (case, N, n_threads) is timed as the MINIMUM wall time per call over
adaptive repetitions (least-noise estimator on a loaded machine). The
``summary`` command fits ``t = a + b N`` on the serial points (relative
least squares): ``a`` is the fixed cost of one call (argument parsing,
allocation of the outputs), ``b`` the marginal cost of one point, which is what
an optimization of a kernel should move. Parallel rows (``n_threads=0``, all
cores) are reported per point.

Usage::

    python benchmark/benchmark_batch.py run -l master -o master.json
    python benchmark/benchmark_batch.py run -l branch -o branch.json \\
        --cases fedooSS_ELISO,fedooUL_ELISO --N 1e4,1e5
    python benchmark/benchmark_batch.py summary master.json branch.json

For an A/B comparison, interleave the two builds (run A, run B, run A, run B,
...) and give all the files to ``summary``: it keeps the minimum per label.
The simcoon module timed is whatever ``import simcoon`` resolves (installed
package, or a build directory put first on ``PYTHONPATH``). Set
``OMP_NUM_THREADS``/``VECLIB_MAXIMUM_THREADS`` explicitly when the machine is
shared.
"""

import argparse
import json
import sys
import time
from collections import defaultdict

import numpy as np

import simcoon as sim

F = np.asfortranarray

ELISO = ([200000.0, 0.3, 0.0], 1)
EPICP = ([200000.0, 0.3, 0.0, 300.0, 1000.0, 0.5], 8)
EPKCP = ([200000.0, 0.3, 0.0, 300.0, 1000.0, 0.5, 10000.0], 14)


def stack_eye(N):
    return F(np.repeat(np.eye(3)[:, :, None], N, axis=2))


def rand_F(N, amp=0.05, seed=0):
    """N deformation gradients I + amp * gaussian noise, shape (3, 3, N)."""
    rng = np.random.default_rng(seed)
    return F(np.eye(3)[:, :, None] + amp * rng.standard_normal((3, 3, N)))


def modular_props():
    """props/nstatev of a MODUL J2 material (isotropic elasticity + power-law hardening)."""
    from simcoon.modular import (IsotropicElasticity, ModularMaterial, Plasticity,
                                 PowerLawHardening, VonMisesYield)
    m = ModularMaterial(
        elasticity=IsotropicElasticity(C1=200000.0, C2=0.3, alpha=0.0),
        mechanisms=[Plasticity(sigma_Y=300.0, yield_criterion=VonMisesYield(),
                               isotropic_hardening=PowerLawHardening(k=1000.0, m=0.5))])
    return np.asarray(m.props, float), int(m.nstatev)


def law(name):
    return {"ELISO": ELISO, "EPICP": EPICP, "EPKCP": EPKCP, "MODUL": modular_props()}[name]


def umat_case(name, props, nstatev, N, regime, finite=False, n_threads=1, temp=False):
    """One ``sim.umat`` call on N identical points; regime "el" stays elastic,
    "pl" is far past yield (sigma_Y = 300, E = 2e5)."""
    etot = F(np.zeros((6, N)))
    Detot = np.zeros((6, N))
    if regime == "el":
        Detot[0, :] = 1e-4
    else:
        Detot[0, :] = 5e-3
        Detot[1, :] = -2.5e-3
        Detot[2, :] = -2.5e-3
    Detot = F(Detot)
    sigma = F(np.zeros((6, N)))
    DR = stack_eye(N)
    statev = F(np.zeros((nstatev, N)))
    Wm = F(np.zeros((4, N)))
    props = F(np.asarray(props, float).reshape(-1, 1))
    if finite:
        F0 = stack_eye(N)
        F1 = F0.copy(order="F")
        F1[0, 0, :] = 1.0 + Detot[0, 0]
        F1[0, 1, :] = 0.01
    else:
        F0 = F1 = np.array([])
    kw = dict(n_threads=n_threads, start=False)
    if temp:
        kw["temp"] = np.array([290.0])

    def call():
        return sim.umat(name, etot, Detot, F0, F1, sigma, DR, props, statev, 1.0, 0.01, Wm, **kw)
    return call


def umat_T_case(name, props, nstatev, N, regime, n_threads=1):
    etot = F(np.zeros((6, N)))
    Detot = np.zeros((6, N))
    Detot[0, :] = 1e-4 if regime == "el" else 5e-3
    Detot = F(Detot)
    sigma = F(np.zeros((6, N)))
    DR = stack_eye(N)
    statev = F(np.zeros((nstatev, N)))
    Wm = F(np.zeros((4, N)))
    Wt = F(np.zeros((3, N)))
    T = np.full(N, 290.0)
    DT = np.zeros(N)
    props = F(np.asarray(props, float).reshape(-1, 1))

    def call():
        return sim.umat_T(name, etot, Detot, sigma, DR, props, statev, 1.0, 0.01, Wm, Wt, T, DT,
                          n_threads=n_threads, start=False)
    return call


def fedoo_ul_case(name, N, n_threads):
    """The simcoon calls of one fedoo global iteration in updated Lagrangian
    (corate log_R, box tangent converted twice by Lt_convert)."""
    props, nstatev = law(name)
    props = F(np.asarray(props, float).reshape(-1, 1))
    F0 = rand_F(N, 0.02, seed=1)
    F1 = F(F0 + 0.004 * np.eye(3)[:, :, None] * np.array([1.0, -0.5, -0.5])[:, None, None])
    e0 = F(np.zeros((6, N)))
    s0 = F(np.zeros((6, N)))
    statev = F(np.zeros((nstatev, N)))
    Wm = F(np.zeros((4, N)))

    def call():
        de, _D, DR, _Omega = sim.objective_rate("log_R", F0, F1, 0.01, True, n_threads)
        rot = sim.Rotation.from_matrix(DR.transpose(2, 0, 1))
        e = rot.apply_strain(e0)
        s = rot.apply_stress(s0)
        sig, _sv, _wm, Lt = sim.umat(name, e, de, F0, F1, s, DR, props, statev, 1.0, 0.01, Wm,
                                     temp=np.array([290.0]), ndi=3, n_threads=n_threads,
                                     tangent_mode=2, corate=3, work_correction=True,
                                     tangent_output="box", start=False)
        dsde = sim.Lt_convert(Lt, F1, sig, "DsigmaDe_2_DSDE")
        return sim.Lt_convert(dsde, F1, sig, "DSDE_2_Dsigma_LieDD")
    return call


def build(case, N, n_threads):
    """Return a zero-argument callable timing *case* on N points."""
    if case == "umat_ELISO":
        return umat_case("ELISO", *ELISO, N, "el", n_threads=n_threads)
    if case == "umat_ELISO_T":   # with temp, as fedoo passes it
        return umat_case("ELISO", *ELISO, N, "el", n_threads=n_threads, temp=True)
    if case in ("umat_EPICP_el", "umat_EPICP_pl"):
        return umat_case("EPICP", *EPICP, N, case[-2:], n_threads=n_threads)
    if case == "umat_EPKCP_pl":
        return umat_case("EPKCP", *EPKCP, N, "pl", n_threads=n_threads)
    if case in ("umat_MODUL_el", "umat_MODUL_pl"):
        return umat_case("MODUL", *modular_props(), N, case[-2:], n_threads=n_threads)
    if case == "umat_NEOHC":
        return umat_case("NEOHC", [1000.0, 10000.0], 1, N, "el", finite=True, n_threads=n_threads)
    if case == "umat_ELISO_finite":
        return umat_case("ELISO", *ELISO, N, "el", finite=True, n_threads=n_threads)
    if case == "umatT_ELISO":
        return umat_T_case("ELISO", [7800.0, 500.0, 200000.0, 0.3, 0.0], 1, N, "el",
                           n_threads=n_threads)
    if case == "umatT_EPICP_pl":
        return umat_T_case("EPICP", [7800.0, 500.0, 200000.0, 0.3, 0.0, 300.0, 1000.0, 0.5], 8,
                           N, "pl", n_threads=n_threads)
    if case.startswith("objrate_"):
        rate = case.split("_", 1)[1]
        F0, F1 = stack_eye(N), rand_F(N)
        return lambda: sim.objective_rate(rate, F0, F1, 0.01, True, n_threads)
    if case.startswith("Ltconv_"):   # no n_threads argument
        key = case.split("_", 1)[1]
        Lt = F(np.repeat((np.eye(6) * 1000.0)[:, :, None], N, axis=2))
        F1, s = rand_F(N), F(np.ones((6, N)))
        return lambda: sim.Lt_convert(Lt, F1, s, key)
    if case.startswith("stressconv_"):
        key = case.split("_", 1)[1]
        F1, s = rand_F(N), F(np.ones((6, N)))
        return lambda: sim.stress_convert(s, F1, key, copy=False)
    if case == "Log_strain":
        F1 = rand_F(N)
        return lambda: sim.Log_strain(F1, True, False)
    if case == "rot_from_matrix_apply_strain":   # fedoo mechanical3d, finite strain
        DR = np.ascontiguousarray(rand_F(N, 0.0).transpose(2, 0, 1))
        e = F(np.ones((6, N)))

        def call():
            r = sim.Rotation.from_matrix(DR)
            return r.apply_strain(e), r.apply_stress(e)
        return call
    if case == "rot_apply_stiffness":
        from scipy.spatial.transform import Rotation as R
        r = sim.Rotation.from_quat(R.random(max(N, 2), random_state=0).as_quat())
        L = F(np.repeat((np.eye(6) * 1000.0)[:, :, None], len(r), axis=2))
        return lambda: r.apply_stiffness(L)
    if case == "tensor2_rotate":
        from scipy.spatial.transform import Rotation as R
        t = sim.Tensor2.stress(np.ones((N, 6)))
        r = sim.Rotation.from_matrix(R.random(N, random_state=0).as_matrix())
        return lambda: t.rotate(r)
    if case == "tensor4_push_forward":
        Fm = np.ascontiguousarray(rand_F(N).transpose(2, 0, 1))
        L = sim.Tensor4.stiffness(np.repeat((np.eye(6) * 1000.0)[None], N, axis=0))
        return lambda: L.push_forward(Fm)
    if case.startswith("fedooUL_"):
        return fedoo_ul_case(case.split("_", 1)[1], N, n_threads)
    if case.startswith("fedooSS_"):
        name = case.split("_", 1)[1]
        return umat_case(name, *law(name), N, "el", n_threads=n_threads, temp=True)
    raise KeyError(f"unknown case {case!r}")


ALL_CASES = [
    "umat_ELISO", "umat_ELISO_T", "umat_EPICP_el", "umat_EPICP_pl", "umat_EPKCP_pl",
    "umat_MODUL_el", "umat_MODUL_pl", "umat_NEOHC", "umat_ELISO_finite", "umatT_ELISO",
    "umatT_EPICP_pl", "objrate_log_R", "objrate_jaumann", "objrate_log",
    "Ltconv_DsigmaDe_2_DSDE", "Ltconv_DSDE_2_Dsigma_LieDD", "stressconv_Cauchy2PKII",
    "Log_strain", "rot_from_matrix_apply_strain", "rot_apply_stiffness", "tensor2_rotate",
    "tensor4_push_forward",
    "fedooSS_ELISO", "fedooSS_EPICP", "fedooSS_MODUL",
    "fedooUL_ELISO", "fedooUL_EPICP", "fedooUL_MODUL",
]


def time_call(call, budget=0.25, min_reps=5, max_reps=2000):
    """Minimum wall time of call() in microseconds, over adaptive repetitions."""
    call()   # warm-up (first-touch allocations, lazy imports)
    best = float("inf")
    t_start = time.perf_counter()
    reps = 0
    while reps < min_reps or (time.perf_counter() - t_start < budget and reps < max_reps):
        t0 = time.perf_counter_ns()
        call()
        best = min(best, time.perf_counter_ns() - t0)
        reps += 1
    return best / 1e3, reps


def cmd_run(args):
    cases = args.cases.split(",") if args.cases else ALL_CASES
    Ns = [int(float(x)) for x in args.N.split(",")]
    threads = [int(x) for x in args.threads.split(",")]
    rows = []
    for case in cases:
        for nt in threads:
            for N in Ns:
                if nt != 1 and N <= 100:
                    continue   # the batch APIs run serially below 100 points anyway
                try:
                    us, reps = time_call(build(case, N, nt), budget=args.budget)
                except Exception as e:   # one broken case must not void the run
                    print(f"{case} N={N} nt={nt}: ERROR {type(e).__name__}: {e}", file=sys.stderr)
                    continue
                rows.append(dict(label=args.label, case=case, N=N, n_threads=nt, us=us, reps=reps))
                print(f"{args.label:8s} {case:30s} nt={nt} N={N:>7d} {us:12.1f} us "
                      f"{us / N:9.3f} us/pt", flush=True)
    with open(args.output, "w") as f:
        json.dump(rows, f, indent=1)
    print(f"\nresults saved to {args.output}")


def fit_affine(N, t):
    """t = a + b N by least squares on the relative residuals (each point weighs
    the same whatever its magnitude, so N = 1 pins a and N = 1e5 pins b)."""
    N, t = np.asarray(N, float), np.asarray(t, float)
    A = np.vstack([np.ones_like(N), N]).T / t[:, None]
    a, b = np.linalg.lstsq(A, np.ones_like(t), rcond=None)[0]
    return a, b


def cmd_summary(args):
    rows = []
    for p in args.files:
        with open(p) as f:
            rows += json.load(f)
    best = {}
    cases = []
    for r in rows:
        k = (r["label"], r["case"], r["N"], r["n_threads"])
        best[k] = min(best.get(k, float("inf")), r["us"])
        if r["case"] not in cases:
            cases.append(r["case"])
    labels = sorted({k[0] for k in best}, key=lambda l: [r["label"] for r in rows].index(l))

    print("## serial (n_threads=1): t = a + b N  [a in us/call, b in us/point]")
    print(f"{'case':30s} " + " ".join(f"{l + ' a':>12s} {l + ' b':>12s}" for l in labels))
    for c in cases:
        line = f"{c:30s} "
        for l in labels:
            pts = sorted((k[2], v) for k, v in best.items() if k[:2] == (l, c) and k[3] == 1)
            if len(pts) < 3:
                line += f"{'-':>12s} {'-':>12s} "
                continue
            a, b = fit_affine(*zip(*pts))
            line += f"{a:12.2f} {b:12.3f} "
        print(line)

    Ns = sorted({k[2] for k in best})
    print("\n## time per point (us/pt) by N; n_threads=0 is the parallel run")
    print(f"{'':30s} {'':6s} " + " ".join(f"{N:>9d}" for N in Ns))
    for c in cases:
        for nt in sorted({k[3] for k in best}, reverse=True):
            for l in labels:
                vals = [best.get((l, c, N, nt)) for N in Ns]
                if all(v is None for v in vals):
                    continue
                cells = " ".join(f"{v / N:9.3f}" if v is not None else f"{'':9s}"
                                 for v, N in zip(vals, Ns))
                print(f"{c:30s} nt={nt:<3d} {cells}  {l}")


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="command", required=True)
    run = sub.add_parser("run", help="time the cases and save the raw minima (JSON)")
    run.add_argument("-l", "--label", default="current", help="build label stored in the rows")
    run.add_argument("-o", "--output", default="benchmark_batch.json")
    run.add_argument("--cases", default="", help="comma-separated subset of: " + ",".join(ALL_CASES))
    run.add_argument("--N", default="1,10,1000,10000,100000", help="points per call")
    run.add_argument("--threads", default="1,0", help="n_threads values (0 = all cores)")
    run.add_argument("--budget", type=float, default=0.25, help="seconds of repetitions per point")
    run.set_defaults(func=cmd_run)
    summ = sub.add_parser("summary", help="fit a + b N and tabulate one or more runs")
    summ.add_argument("files", nargs="+")
    summ.set_defaults(func=cmd_summary)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
