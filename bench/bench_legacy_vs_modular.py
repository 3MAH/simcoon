"""Benchmark legacy UMATs against their MODUL (modular UMAT) twin configurations.

ZENER/ZENNK are intentionally absent: Zener_fast/Zener_Nfast are generalized
KELVIN chains (branch driving force sigma - L_i EV_i), a different rheological
model from the modular Prony (generalized Maxwell) — no twin exists.

Methodology: INTERLEAVED A/B — each repetition runs the legacy kernel then the
modular twin back-to-back, so load spikes hit both arms equally; the reported
estimator is the MINIMUM over repetitions (least-noise wall time; naive medians
proved 5x sensitive to background load). A max-deviation column doubles as an
equivalence smoke check.

Note that the load path is parsed once per arm OUTSIDE the timed region, so what is timed
is the constitutive kernel alone rather than the path-file parse the two arms
share.

Usage:
    python bench/bench_legacy_vs_modular.py [--families ELISO,EPCHA,...]
                                            [--repeats 8]

Run from anywhere; paths resolve relative to the repo. Requires the simcoon
python package (editable install). Results land in a temp dir and are removed.
"""

import argparse
import time
from pathlib import Path

import numpy as np

import simcoon as sim
from simcoon.modular import (
    ModularMaterial,
    IsotropicElasticity,
    CubicElasticity,
    TransverseIsotropicElasticity,
    OrthotropicElasticity,
    Plasticity,
    VoceHardening,
    PowerLawHardening,
    HillYield,
    DFAYield,
    AnisotropicYield,
    PragerHardening,
    ChabocheHardening,
    Viscoelasticity,
)

REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "examples" / "data"

VE_TERMS = [
    (1500.0, 0.35, 3000.0, 1200.0),
    (800.0, 0.35, 30000.0, 12000.0),
    (400.0, 0.35, 300000.0, 120000.0),
]


def _mat(elasticity, mechanisms=()):
    return ModularMaterial(elasticity=elasticity, mechanisms=list(mechanisms))


# family -> (legacy_name, legacy_props, legacy_nstatev, path_file, ModularMaterial)
FAMILIES = {
    "ELISO": (
        "ELISO",
        [210000.0, 0.3, 1.2e-5],
        1,
        "MODUL_path.json",
        _mat(IsotropicElasticity(C1=210000.0, C2=0.3, alpha=1.2e-5)),
    ),
    "ELIST": (
        "ELIST",
        [3, 230000.0, 15000.0, 0.02, 0.4, 50000.0, 0.0, 0.0],
        1,
        "MODUL_path.json",
        _mat(
            TransverseIsotropicElasticity(
                EL=230000.0, ET=15000.0, nuTL=0.02, nuTT=0.4, GLT=50000.0, axis=3
            )
        ),
    ),
    "ELORT": (
        "ELORT",
        [
            70000.0,
            30000.0,
            15000.0,
            0.3,
            0.3,
            0.3,
            8000.0,
            6000.0,
            5000.0,
            1e-5,
            2e-5,
            3e-5,
        ],
        1,
        "MODUL_path.json",
        _mat(
            OrthotropicElasticity(
                C1=70000.0,
                C2=30000.0,
                C3=15000.0,
                C4=0.3,
                C5=0.3,
                C6=0.3,
                C7=8000.0,
                C8=6000.0,
                C9=5000.0,
                alpha1=1e-5,
                alpha2=2e-5,
                alpha3=3e-5,
            )
        ),
    ),
    "EPICP": (
        "EPICP",
        [210000.0, 0.3, 0.0, 300.0, 1000.0, 0.3],
        8,
        "MODUL_path.json",
        _mat(
            IsotropicElasticity(C1=210000.0, C2=0.3),
            [
                Plasticity(
                    sigma_Y=300.0,
                    isotropic_hardening=PowerLawHardening(k=1000.0, m=0.3),
                )
            ],
        ),
    ),
    "EPKCP": (
        "EPKCP",
        [210000.0, 0.3, 0.0, 300.0, 1000.0, 1.0, 20000.0],
        33,
        "MODUL_path.json",
        _mat(
            IsotropicElasticity(C1=210000.0, C2=0.3),
            [
                Plasticity(
                    sigma_Y=300.0,
                    isotropic_hardening=PowerLawHardening(k=1000.0, m=1.0),
                    kinematic_hardening=PragerHardening(C=1.5 * 20000.0),
                )
            ],
        ),
    ),
    "EPCHA": (
        "EPCHA",
        [210000.0, 0.3, 0.0, 300.0, 200.0, 20.0, 30000.0, 172.0, 19500.0, 301.0],
        33,
        "MODUL_path.json",
        _mat(
            IsotropicElasticity(C1=210000.0, C2=0.3),
            [
                Plasticity(
                    sigma_Y=300.0,
                    isotropic_hardening=VoceHardening(Q=200.0, b=20.0),
                    kinematic_hardening=ChabocheHardening(
                        terms=((30000.0, 172.0), (19500.0, 301.0))
                    ),
                )
            ],
        ),
    ),
    "EPHIL": (
        "EPHIL",
        [210000.0, 0.3, 0.0, 300.0, 5000.0, 1.0, 0.5, 0.4, 0.6, 1.5, 1.5, 1.5],
        33,
        "MODUL_path.json",
        _mat(
            IsotropicElasticity(C1=210000.0, C2=0.3),
            [
                Plasticity(
                    sigma_Y=300.0,
                    yield_criterion=HillYield(F=0.5, G=0.4, H=0.6, L=1.5, M=1.5, N=1.5),
                    isotropic_hardening=PowerLawHardening(k=5000.0, m=1.0),
                )
            ],
        ),
    ),
    "EPHAC": (
        "EPHAC",
        [
            210000.0,
            0.3,
            85000.0,
            0.0,
            300.0,
            200.0,
            20.0,
            30000.0,
            172.0,
            19500.0,
            301.0,
            0.5,
            0.4,
            0.6,
            1.5,
            1.5,
            1.5,
        ],
        33,
        "MODUL_path.json",
        _mat(
            CubicElasticity(C1=210000.0, C2=0.3, C3=85000.0),
            [
                Plasticity(
                    sigma_Y=300.0,
                    yield_criterion=HillYield(F=0.5, G=0.4, H=0.6, L=1.5, M=1.5, N=1.5),
                    isotropic_hardening=VoceHardening(Q=200.0, b=20.0),
                    kinematic_hardening=ChabocheHardening(
                        terms=((30000.0, 172.0), (19500.0, 301.0))
                    ),
                )
            ],
        ),
    ),
    "EPANI": (
        "EPANI",
        [
            210000.0,
            0.3,
            85000.0,
            0.0,
            300.0,
            200.0,
            20.0,
            30000.0,
            172.0,
            19500.0,
            301.0,
            1.2,
            1.1,
            1.1,
            -0.6,
            -0.6,
            -0.5,
            1.6,
            1.5,
            1.4,
        ],
        33,
        "MODUL_path.json",
        _mat(
            CubicElasticity(C1=210000.0, C2=0.3, C3=85000.0),
            [
                Plasticity(
                    sigma_Y=300.0,
                    yield_criterion=AnisotropicYield(
                        P11=1.2,
                        P22=1.1,
                        P33=1.1,
                        P12=-0.6,
                        P13=-0.6,
                        P23=-0.5,
                        P44=1.6,
                        P55=1.5,
                        P66=1.4,
                    ),
                    isotropic_hardening=VoceHardening(Q=200.0, b=20.0),
                    kinematic_hardening=ChabocheHardening(
                        terms=((30000.0, 172.0), (19500.0, 301.0))
                    ),
                )
            ],
        ),
    ),
    "EPDFA": (
        "EPDFA",
        [
            210000.0,
            0.3,
            85000.0,
            0.0,
            300.0,
            200.0,
            20.0,
            30000.0,
            172.0,
            19500.0,
            301.0,
            0.5,
            0.4,
            0.6,
            1.5,
            1.5,
            1.5,
            0.1,
        ],
        33,
        "MODUL_path.json",
        _mat(
            CubicElasticity(C1=210000.0, C2=0.3, C3=85000.0),
            [
                Plasticity(
                    sigma_Y=300.0,
                    yield_criterion=DFAYield(
                        F=0.5, G=0.4, H=0.6, L=1.5, M=1.5, N=1.5, K=0.1
                    ),
                    isotropic_hardening=VoceHardening(Q=200.0, b=20.0),
                    kinematic_hardening=ChabocheHardening(
                        terms=((30000.0, 172.0), (19500.0, 301.0))
                    ),
                )
            ],
        ),
    ),
    "EPCHG": (
        "EPCHG",
        [
            210000.0,
            0.3,
            85000.0,
            0.0,
            300.0,
            2,
            2,
            0,
            150.0,
            15.0,
            50.0,
            40.0,
            30000.0,
            172.0,
            19500.0,
            301.0,
        ],
        33,
        "MODUL_path.json",
        _mat(
            CubicElasticity(C1=210000.0, C2=0.3, C3=85000.0),
            [
                Plasticity(
                    sigma_Y=300.0,
                    # legacy N-term iso couples through a single Hp
                    # -> ONE effective Voce (b_eff=sum b_i,
                    # Q_eff=sum(b_i Q_i)/sum(b_i))
                    isotropic_hardening=VoceHardening(
                        Q=(15.0 * 150.0 + 40.0 * 50.0) / 55.0, b=55.0
                    ),
                    kinematic_hardening=ChabocheHardening(
                        terms=((30000.0, 172.0), (19500.0, 301.0))
                    ),
                )
            ],
        ),
    ),
    # EPHIN legacy is defective for N >= 2 (NaN even for identical/inactive
    # second surface) — bench the provable N = 1 case.
    "EPHIN": (
        "EPHIN",
        [210000.0, 0.3, 0.0, 1, 300.0, 3000.0, 1.0, 0.5, 0.4, 0.6, 1.5, 1.5, 1.5],
        33,
        "MODUL_path.json",
        _mat(
            IsotropicElasticity(C1=210000.0, C2=0.3),
            [
                Plasticity(
                    sigma_Y=300.0,
                    yield_criterion=HillYield(F=0.5, G=0.4, H=0.6, L=1.5, M=1.5, N=1.5),
                    isotropic_hardening=PowerLawHardening(k=3000.0, m=1.0),
                )
            ],
        ),
    ),
    "PRONK": (
        "PRONK",
        [3000.0, 0.35, 0.0, 3] + [x for t in VE_TERMS for x in t],
        7 + 7 * 3,
        "PRONK_path.json",
        _mat(
            IsotropicElasticity(C1=3000.0, C2=0.35), [Viscoelasticity(terms=VE_TERMS)]
        ),
    ),
}


def bench_family(key, repeats):
    name, props, nstatev, path_file, mat = FAMILIES[key]
    props = np.asarray(props, dtype=float)
    t_leg, t_mod = [], []
    for _ in range(repeats):
        # re-parse per arm so neither one inherits blocks the other consumed
        blocks, T_init = sim.solver.load_path_json(str(DATA / path_file))[:2]
        t0 = time.perf_counter()
        res_l = sim.solver.solve(blocks, name, props, nstatev, T_init=T_init, corate=1)
        t_leg.append(time.perf_counter() - t0)
        blocks, T_init = sim.solver.load_path_json(str(DATA / path_file))[:2]
        t0 = time.perf_counter()
        res_m = sim.solver.solve(
            blocks, "MODUL", mat.props, mat.nstatev, T_init=T_init, corate=1
        )
        t_mod.append(time.perf_counter() - t0)
    # discard the first (warmup) pair; min = least-noise estimator
    tl, tm = min(t_leg[1:]), min(t_mod[1:])
    a = np.asarray(res_l["Stress"])[0]
    b = np.asarray(res_m["Stress"])[0]
    dev = np.max(np.abs(a - b)) / max(np.max(np.abs(a)), 1e-12)
    return tl * 1e3, tm * 1e3, dev


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--families",
        default=",".join(FAMILIES),
        help="comma-separated subset of: " + ",".join(FAMILIES),
    )
    ap.add_argument("--repeats", type=int, default=8)
    args = ap.parse_args()

    print("| family | legacy [ms] | MODUL [ms] | ratio | max dev |")
    print("|---|---|---|---|---|")
    for key in args.families.split(","):
        key = key.strip()
        if key not in FAMILIES:
            print(f"| {key} | (unknown family) | | | |")
            continue
        tl, tm, dev = bench_family(key, args.repeats)
        print(f"| {key} | {tl:.1f} | {tm:.1f} | {tm / tl:.2f} | {dev:.2e} |")


if __name__ == "__main__":
    main()
