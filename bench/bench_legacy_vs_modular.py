"""Benchmark legacy UMATs against their MODUL (modular UMAT) twin configurations.

ZENER/ZENNK are intentionally absent: Zener_fast/Zener_Nfast are generalized
KELVIN chains (branch driving force sigma - L_i EV_i), a different rheological
model from the modular Prony (generalized Maxwell) — no twin exists.

Methodology: INTERLEAVED A/B — each repetition runs the legacy kernel then the
modular twin back-to-back, so load spikes hit both arms equally; the reported
estimator is the MINIMUM over repetitions (least-noise wall time; naive medians
proved 5x sensitive to background load). A max-deviation column doubles as an
equivalence smoke check.

Usage:
    python bench/bench_legacy_vs_modular.py [--families ELISO,EPCHA,...]
                                            [--repeats 8]

Run from anywhere; paths resolve relative to the repo. Requires the simcoon
python package (editable install). Results land in a temp dir and are removed.
"""
import argparse
import os
import statistics
import tempfile
import time
from pathlib import Path

import numpy as np

import simcoon as sim
from simcoon.modular import (
    ModularMaterial, IsotropicElasticity, CubicElasticity,
    TransverseIsotropicElasticity, OrthotropicElasticity,
    Plasticity, VoceHardening, PowerLawHardening,
    HillYield, DFAYield, AnisotropicYield,
    PragerHardening, ChabocheHardening, Viscoelasticity,
)

REPO = Path(__file__).resolve().parents[1]
EXAMPLES = REPO / "examples" / "mechanical"

VE_TERMS = [(1500., 0.35, 3000., 1200.),
            (800., 0.35, 30000., 12000.),
            (400., 0.35, 300000., 120000.)]


def _mat(elasticity, mechanisms=()):
    return ModularMaterial(elasticity=elasticity, mechanisms=list(mechanisms))


# family -> (legacy_name, legacy_props, legacy_nstatev, path_file, ModularMaterial)
FAMILIES = {
    "ELISO": ("ELISO", [210000., 0.3, 1.2e-5], 1, "MODUL_path.txt",
              _mat(IsotropicElasticity(C1=210000., C2=0.3, alpha=1.2e-5))),
    "ELIST": ("ELIST", [3, 230000., 15000., 0.02, 0.4, 50000., 0., 0.], 1,
              "MODUL_path.txt",
              _mat(TransverseIsotropicElasticity(EL=230000., ET=15000.,
                                                 nuTL=0.02, nuTT=0.4,
                                                 GLT=50000., axis=3))),
    "ELORT": ("ELORT", [70000., 30000., 15000., 0.3, 0.3, 0.3,
                        8000., 6000., 5000., 1e-5, 2e-5, 3e-5], 1,
              "MODUL_path.txt",
              _mat(OrthotropicElasticity(C1=70000., C2=30000., C3=15000.,
                                         C4=0.3, C5=0.3, C6=0.3,
                                         C7=8000., C8=6000., C9=5000.,
                                         alpha1=1e-5, alpha2=2e-5, alpha3=3e-5))),
    "EPICP": ("EPICP", [210000., 0.3, 0., 300., 1000., 0.3], 8, "MODUL_path.txt",
              _mat(IsotropicElasticity(C1=210000., C2=0.3),
                   [Plasticity(sigma_Y=300.,
                               isotropic_hardening=PowerLawHardening(k=1000., m=0.3))])),
    "EPKCP": ("EPKCP", [210000., 0.3, 0., 300., 1000., 1.0, 20000.], 33,
              "MODUL_path.txt",
              _mat(IsotropicElasticity(C1=210000., C2=0.3),
                   [Plasticity(sigma_Y=300.,
                               isotropic_hardening=PowerLawHardening(k=1000., m=1.0),
                               kinematic_hardening=PragerHardening(C=1.5 * 20000.))])),
    "EPCHA": ("EPCHA", [210000., 0.3, 0., 300., 200., 20.,
                        30000., 172., 19500., 301.], 33, "MODUL_path.txt",
              _mat(IsotropicElasticity(C1=210000., C2=0.3),
                   [Plasticity(sigma_Y=300.,
                               isotropic_hardening=VoceHardening(Q=200., b=20.),
                               kinematic_hardening=ChabocheHardening(
                                   terms=((30000., 172.), (19500., 301.))))])),
    "EPHIL": ("EPHIL", [210000., 0.3, 0., 300., 5000., 1.0,
                        0.5, 0.4, 0.6, 1.5, 1.5, 1.5], 33, "MODUL_path.txt",
              _mat(IsotropicElasticity(C1=210000., C2=0.3),
                   [Plasticity(sigma_Y=300.,
                               yield_criterion=HillYield(F=0.5, G=0.4, H=0.6,
                                                         L=1.5, M=1.5, N=1.5),
                               isotropic_hardening=PowerLawHardening(k=5000., m=1.0))])),
    "EPHAC": ("EPHAC", [210000., 0.3, 85000., 0., 300., 200., 20.,
                        30000., 172., 19500., 301.,
                        0.5, 0.4, 0.6, 1.5, 1.5, 1.5], 33, "MODUL_path.txt",
              _mat(CubicElasticity(C1=210000., C2=0.3, C3=85000.),
                   [Plasticity(sigma_Y=300.,
                               yield_criterion=HillYield(F=0.5, G=0.4, H=0.6,
                                                         L=1.5, M=1.5, N=1.5),
                               isotropic_hardening=VoceHardening(Q=200., b=20.),
                               kinematic_hardening=ChabocheHardening(
                                   terms=((30000., 172.), (19500., 301.))))])),
    "EPANI": ("EPANI", [210000., 0.3, 85000., 0., 300., 200., 20.,
                        30000., 172., 19500., 301.,
                        1.2, 1.1, 1.1, -0.6, -0.6, -0.5, 1.6, 1.5, 1.4], 33,
              "MODUL_path.txt",
              _mat(CubicElasticity(C1=210000., C2=0.3, C3=85000.),
                   [Plasticity(sigma_Y=300.,
                               yield_criterion=AnisotropicYield(
                                   P11=1.2, P22=1.1, P33=1.1, P12=-0.6,
                                   P13=-0.6, P23=-0.5, P44=1.6, P55=1.5,
                                   P66=1.4),
                               isotropic_hardening=VoceHardening(Q=200., b=20.),
                               kinematic_hardening=ChabocheHardening(
                                   terms=((30000., 172.), (19500., 301.))))])),
    "EPDFA": ("EPDFA", [210000., 0.3, 85000., 0., 300., 200., 20.,
                        30000., 172., 19500., 301.,
                        0.5, 0.4, 0.6, 1.5, 1.5, 1.5, 0.1], 33, "MODUL_path.txt",
              _mat(CubicElasticity(C1=210000., C2=0.3, C3=85000.),
                   [Plasticity(sigma_Y=300.,
                               yield_criterion=DFAYield(F=0.5, G=0.4, H=0.6,
                                                        L=1.5, M=1.5, N=1.5,
                                                        K=0.1),
                               isotropic_hardening=VoceHardening(Q=200., b=20.),
                               kinematic_hardening=ChabocheHardening(
                                   terms=((30000., 172.), (19500., 301.))))])),
    "EPCHG": ("EPCHG", [210000., 0.3, 85000., 0., 300., 2, 2, 0,
                        150., 15., 50., 40., 30000., 172., 19500., 301.], 33,
              "MODUL_path.txt",
              _mat(CubicElasticity(C1=210000., C2=0.3, C3=85000.),
                   [Plasticity(sigma_Y=300.,
                               # legacy N-term iso couples through a single Hp
                               # -> ONE effective Voce (b_eff=sum b_i,
                               # Q_eff=sum(b_i Q_i)/sum(b_i))
                               isotropic_hardening=VoceHardening(
                                   Q=(15. * 150. + 40. * 50.) / 55., b=55.),
                               kinematic_hardening=ChabocheHardening(
                                   terms=((30000., 172.), (19500., 301.))))])),
    # EPHIN legacy is defective for N >= 2 (NaN even for identical/inactive
    # second surface) — bench the provable N = 1 case.
    "EPHIN": ("EPHIN", [210000., 0.3, 0., 1,
                        300., 3000., 1.0, 0.5, 0.4, 0.6, 1.5, 1.5, 1.5], 33,
              "MODUL_path.txt",
              _mat(IsotropicElasticity(C1=210000., C2=0.3),
                   [Plasticity(sigma_Y=300.,
                               yield_criterion=HillYield(F=0.5, G=0.4, H=0.6,
                                                         L=1.5, M=1.5, N=1.5),
                               isotropic_hardening=PowerLawHardening(k=3000., m=1.0))])),
    "PRONK": ("PRONK", [3000., 0.35, 0., 3] + [x for t in VE_TERMS for x in t],
              7 + 7 * 3, "PRONK_path.txt",
              _mat(IsotropicElasticity(C1=3000., C2=0.35),
                   [Viscoelasticity(terms=VE_TERMS)])),
}


def bench_family(key, results_dir, repeats):
    name, props, nstatev, path_file, mat = FAMILIES[key]
    props = np.asarray(props, dtype=float)
    t_leg, t_mod = [], []
    out_l, out_m = f"b_{key}_l.txt", f"b_{key}_m.txt"
    for _ in range(repeats):
        t0 = time.perf_counter()
        sim._core.solver(name, props, nstatev, 0., 0., 0., 0, 1,
                   "../data", results_dir, path_file, out_l)
        t_leg.append(time.perf_counter() - t0)
        t0 = time.perf_counter()
        sim._core.solver("MODUL", mat.props, mat.nstatev, 0., 0., 0., 0, 1,
                   "../data", results_dir, path_file, out_m)
        t_mod.append(time.perf_counter() - t0)
    # discard the first (warmup) pair; min = least-noise estimator
    tl, tm = min(t_leg[1:]), min(t_mod[1:])
    res = Path(EXAMPLES) / results_dir
    a = np.loadtxt(res / out_l.replace(".txt", "_global-0.txt"))[:, 14]
    b = np.loadtxt(res / out_m.replace(".txt", "_global-0.txt"))[:, 14]
    dev = np.max(np.abs(a - b)) / max(np.max(np.abs(a)), 1e-12)
    return tl * 1e3, tm * 1e3, dev


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--families", default=",".join(FAMILIES),
                    help="comma-separated subset of: " + ",".join(FAMILIES))
    ap.add_argument("--repeats", type=int, default=8)
    args = ap.parse_args()

    os.chdir(EXAMPLES)
    with tempfile.TemporaryDirectory() as tmp:
        link = EXAMPLES / f"_bench_{os.getpid()}"
        link.symlink_to(tmp, target_is_directory=True)
        try:
            print(f"| family | legacy [ms] | MODUL [ms] | ratio | max dev |")
            print(f"|---|---|---|---|---|")
            for key in args.families.split(","):
                key = key.strip()
                if key not in FAMILIES:
                    print(f"| {key} | (unknown family) | | | |")
                    continue
                tl, tm, dev = bench_family(key, link.name, args.repeats)
                print(f"| {key} | {tl:.1f} | {tm:.1f} | {tm/tl:.2f} | {dev:.2e} |")
        finally:
            link.unlink()


if __name__ == "__main__":
    main()
