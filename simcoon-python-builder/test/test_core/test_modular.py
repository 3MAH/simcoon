"""Regression tests for the modular UMAT Python API.

These tests exercise the end-to-end `ModularMaterial → sim.solver("MODUL", ...)`
path for the cases that would silently break if the modular C++ orchestrator
or Python props-serialization regressed.
"""

import os
from pathlib import Path

import numpy as np
import pytest

import simcoon as sim
from simcoon.modular import (
    ArmstrongFrederickHardening,
    IsotropicElasticity,
    LinearIsotropicHardening,
    ModularMaterial,
    Plasticity,
    VoceHardening,
    Viscoelasticity,
)


EXAMPLES_DIR = Path(__file__).resolve().parents[3] / "examples" / "mechanical"


@pytest.fixture
def work_in_examples(tmp_path, monkeypatch):
    """sim.solver reads/writes relative paths — chdir to the examples folder
    so `data/MODUL_path.txt` resolves, but point results into a pytest tmpdir."""
    monkeypatch.chdir(EXAMPLES_DIR)
    results = tmp_path / "results"
    results.mkdir()
    # sim.solver's 'results' folder arg is relative; create a symlink for this
    # test run so its output lands in the pytest tmp area.
    link = EXAMPLES_DIR / f"_pytest_results_{tmp_path.name}"
    link.symlink_to(results, target_is_directory=True)
    yield link.name
    link.unlink()


def _run_solver(mat: ModularMaterial, results_dir: str, outfile: str) -> np.ndarray:
    """Run the MODUL solver through the standard path file; return the
    (n_steps, 2) array of (eps_11, sigma_11)."""
    sim.solver(
        mat.umat_name, mat.props, mat.nstatev,
        0.0, 0.0, 0.0,        # psi_rve, theta_rve, phi_rve
        0, 1,                 # solver_type, corate_type
        "../data", results_dir, "MODUL_path.txt", outfile,
    )
    out = Path(EXAMPLES_DIR) / results_dir / outfile.replace(".txt", "_global-0.txt")
    return np.loadtxt(out, usecols=(8, 14))


def test_af_with_list_args_rejected():
    """Regression: AF takes scalars only; lists should fail fast with a
    helpful message pointing to ChabocheHardening."""
    with pytest.raises(TypeError, match="ChabocheHardening"):
        ArmstrongFrederickHardening(C=[30000.0, 195000.0], D=[172.0, 3012.0])


def test_viscoelasticity_legacy_tuple_rejected():
    """The old (g, tau) 2-tuple form must be rejected with a pointer to the
    new Prony_Nfast (E, nu, etaB, etaS) 4-tuple layout."""
    with pytest.raises(TypeError, match=r"\(E, nu, etaB, etaS\)"):
        Viscoelasticity(terms=[(0.3, 1.0), (0.2, 10.0)])


def test_two_plasticity_equivalent_to_one_linear(work_in_examples):
    """Two identical linear-hardening plasticity mechanisms with modulus 2H
    each must give the same stress as a single mechanism with modulus H, by
    superposition of associated flow. This exercises the cross-mechanism
    Jacobian assembly in ModularUMAT::return_mapping — without it the FB
    iteration diverges to nonphysical stress."""
    E, nu, sigma_Y, H = 210000.0, 0.3, 300.0, 10000.0

    mat_single = ModularMaterial(
        elasticity=IsotropicElasticity(E=E, nu=nu),
        mechanisms=[Plasticity(
            sigma_Y=sigma_Y,
            isotropic_hardening=LinearIsotropicHardening(H=H),
        )],
    )
    mat_two = ModularMaterial(
        elasticity=IsotropicElasticity(E=E, nu=nu),
        mechanisms=[
            Plasticity(sigma_Y=sigma_Y,
                       isotropic_hardening=LinearIsotropicHardening(H=2 * H)),
            Plasticity(sigma_Y=sigma_Y,
                       isotropic_hardening=LinearIsotropicHardening(H=2 * H)),
        ],
    )

    hist_single = _run_solver(mat_single, work_in_examples, "res_single.txt")
    hist_two    = _run_solver(mat_two,    work_in_examples, "res_two.txt")

    # Same number of increments, same loading history → stresses should match
    # within solver tolerance at every step.
    assert hist_single.shape == hist_two.shape
    # Compare peak |sigma_11| — the most load-bearing quantity for plasticity.
    peak_single = np.max(np.abs(hist_single[:, 1]))
    peak_two    = np.max(np.abs(hist_two[:, 1]))
    assert peak_single > 0, "sanity"
    # Relative tolerance covers small solver-stepping differences from the
    # doubled constraint count; coupling fix brings this from 100%+ to ~1.5%.
    rel_err = abs(peak_two - peak_single) / peak_single
    assert rel_err < 0.03, (
        f"two-plasticity peak {peak_two:.2f} differs from single {peak_single:.2f} "
        f"by {rel_err:.3%}; cross-mechanism Jacobian coupling may have regressed"
    )


def test_voce_plasticity_runs_end_to_end(work_in_examples):
    """Smoke test for the baseline MODUL path — ensures the props layout and
    statev allocation are in sync between Python and C++ after any modular
    refactor."""
    mat = ModularMaterial(
        elasticity=IsotropicElasticity(E=210000.0, nu=0.3, alpha=1.2e-5),
        mechanisms=[Plasticity(
            sigma_Y=300.0,
            isotropic_hardening=VoceHardening(Q=200.0, b=10.0),
        )],
    )
    hist = _run_solver(mat, work_in_examples, "res_voce.txt")
    peak = np.max(np.abs(hist[:, 1]))
    # With sigma_Y + Q = 500 and b=10 at 2% strain, peak sits around 460 MPa
    # after the cyclic path saturates hardening.
    assert 400.0 < peak < 500.0, f"Voce peak {peak:.1f} outside expected band"


def test_viscoelastic_matches_pronk_reference(work_in_examples):
    """MODUL with a Viscoelasticity mechanism must reproduce the reference
    PRONK (Prony_Nfast) UMAT through the solver on a relaxation path.

    Regression for the consistent-tangent bug where tangent_contribution
    rebuilt K(i,i) without the -1/DTime term (Lt went singular and the
    mixed-control solver diverged)."""
    terms = [(1500.0, 0.35, 3000.0, 1200.0),
             (800.0, 0.35, 30000.0, 12000.0),
             (400.0, 0.35, 300000.0, 120000.0)]
    E0, nu0 = 3000.0, 0.35

    pronk_props = np.array([E0, nu0, 0.0, len(terms)]
                           + [x for t in terms for x in t])
    sim.solver("PRONK", pronk_props, 7 + 7 * len(terms),
               0.0, 0.0, 0.0, 0, 1,
               "../data", work_in_examples, "PRONK_path.txt", "res_pronk.txt")
    ref = np.loadtxt(Path(EXAMPLES_DIR) / work_in_examples
                     / "res_pronk_global-0.txt", usecols=(8, 14))

    mat = ModularMaterial(
        elasticity=IsotropicElasticity(E=E0, nu=nu0),
        mechanisms=[Viscoelasticity(terms=terms)],
    )
    sim.solver(mat.umat_name, mat.props, mat.nstatev,
               0.0, 0.0, 0.0, 0, 1,
               "../data", work_in_examples, "PRONK_path.txt", "res_veq.txt")
    hist = np.loadtxt(Path(EXAMPLES_DIR) / work_in_examples
                      / "res_veq_global-0.txt", usecols=(8, 14))

    assert hist.shape == ref.shape
    peak = np.max(np.abs(ref[:, 1]))
    assert peak > 1.0, "sanity: relaxation path produced stress"
    max_diff = np.max(np.abs(hist[:, 1] - ref[:, 1]))
    assert max_diff / peak < 1e-3, (
        f"MODUL viscoelastic deviates from PRONK by {max_diff/peak:.3%} "
        "(consistent-tangent regression?)"
    )


def test_damage_softens_end_to_end(work_in_examples):
    """LINEAR damage through sim.solver: the damaged run must (a) develop
    0 < D-driven softening vs the undamaged elastic run at the same strain,
    and (b) keep stress positive. Regression for the orchestrator ignoring
    stiffness_reduction (D evolved but stress was never softened)."""
    from simcoon.modular import Damage, elastic_model

    mat_el = elastic_model(E=210000.0, nu=0.3)
    hist_el = _run_solver(mat_el, work_in_examples, "res_el.txt")

    mat_dm = ModularMaterial(
        elasticity=IsotropicElasticity(E=210000.0, nu=0.3),
        mechanisms=[Damage(Y_0=0.05, Y_c=10.0)],
    )
    hist_dm = _run_solver(mat_dm, work_in_examples, "res_dm.txt")

    assert hist_dm.shape == hist_el.shape
    peak_el = np.max(np.abs(hist_el[:, 1]))
    peak_dm = np.max(np.abs(hist_dm[:, 1]))
    assert peak_el > 0
    assert peak_dm > 0
    assert peak_dm < 0.95 * peak_el, (
        f"damage peak {peak_dm:.1f} not softened vs elastic {peak_el:.1f}"
    )


def test_tangent_mode_1_same_converged_response(work_in_examples):
    """tangent_mode=1 (Simo-Hughes algorithmic) must reproduce the mode-0
    converged response — the tangent steers the global Newton, not the
    residual. Guards the mode-1 assembly (Bhat = -B sign convention and the
    coupled sub-block extraction) at solver level."""
    mat = ModularMaterial(
        elasticity=IsotropicElasticity(E=210000.0, nu=0.3),
        mechanisms=[Plasticity(
            sigma_Y=300.0,
            isotropic_hardening=VoceHardening(Q=200.0, b=20.0),
            kinematic_hardening=ArmstrongFrederickHardening(C=30000.0, D=172.0),
        )],
    )

    outs = {}
    for mode in (0, 1):
        out = f"res_tg{mode}.txt"
        sim.solver(mat.umat_name, mat.props, mat.nstatev,
                   0.0, 0.0, 0.0, 0, 1,
                   "../data", work_in_examples, "MODUL_path.txt", out,
                   mode)
        outs[mode] = np.loadtxt(Path(EXAMPLES_DIR) / work_in_examples
                                / out.replace(".txt", "_global-0.txt"))

    assert outs[0].shape == outs[1].shape
    peak = np.max(np.abs(outs[0][:, 14]))
    diff = np.max(np.abs(outs[0][:, 14] - outs[1][:, 14]))
    assert peak > 100.0
    assert diff / peak < 1e-5, (
        f"mode-1 response deviates from mode-0 by {diff/peak:.2e} "
        "(algorithmic tangent must not change the converged solution)"
    )


def test_chaboche_matches_epcha_reference(work_in_examples):
    """MODUL (Voce + 2-term Chaboche) must reproduce the reference EPCHA UMAT
    through the solver on the cyclic path.

    Regression for the FB-Jacobian hardening-sign bug: B carried +H_total
    instead of -H_total (dPhi/dp = -H), a wrong Newton slope that made the FB
    error metric report convergence prematurely — 4.4% peak-stress error vs
    EPCHA on this path, invisible to every monotonic test."""
    from simcoon.modular import ChabocheHardening

    epcha_props = np.array([210000.0, 0.3, 0.0,
                            300.0, 200.0, 20.0,
                            30000.0, 172.0, 19500.0, 301.0])
    sim.solver("EPCHA", epcha_props, 33, 0.0, 0.0, 0.0, 0, 1,
               "../data", work_in_examples, "MODUL_path.txt", "res_epcha.txt")
    ref = np.loadtxt(Path(EXAMPLES_DIR) / work_in_examples
                     / "res_epcha_global-0.txt", usecols=(8, 14))

    mat = ModularMaterial(
        elasticity=IsotropicElasticity(E=210000.0, nu=0.3),
        mechanisms=[Plasticity(
            sigma_Y=300.0,
            isotropic_hardening=VoceHardening(Q=200.0, b=20.0),
            kinematic_hardening=ChabocheHardening(
                terms=((30000.0, 172.0), (19500.0, 301.0))),
        )],
    )
    sim.solver(mat.umat_name, mat.props, mat.nstatev, 0.0, 0.0, 0.0, 0, 1,
               "../data", work_in_examples, "MODUL_path.txt", "res_mchab.txt")
    hist = np.loadtxt(Path(EXAMPLES_DIR) / work_in_examples
                      / "res_mchab_global-0.txt", usecols=(8, 14))

    assert hist.shape == ref.shape
    peak = np.max(np.abs(ref[:, 1]))
    max_diff = np.max(np.abs(hist[:, 1] - ref[:, 1]))
    assert peak > 500.0, "sanity: cyclic path reached the hardened regime"
    assert max_diff / peak < 2e-3, (
        f"MODUL Chaboche deviates from EPCHA by {max_diff/peak:.3%} "
        "(FB Jacobian sign regression?)"
    )


def test_hill_matches_ephil_reference(work_in_examples):
    """MODUL (Hill yield + linear hardening) must reproduce the reference
    EPHIL UMAT bit-for-bit through the solver.

    End-to-end coverage for the anisotropic-criterion path (Hill params ->
    P_Hill dispatch -> flow direction). Linear hardening (power-law exponent
    m = 1) is used deliberately: the general power-law tangent
    dR/dp = m k p^(m-1) is singular at p = 0, so the two CCP loops differ by
    a transient (~1.4%) at yield onset that re-converges — with m = 1 that
    singularity is absent and the integrators agree to machine precision,
    isolating the Hill criterion itself."""
    from simcoon.modular import HillYield, PowerLawHardening

    ephil_props = np.array([210000., 0.3, 0., 300., 5000., 1.0,
                            0.5, 0.4, 0.6, 1.5, 1.5, 1.5])
    sim.solver("EPHIL", ephil_props, 33, 0.0, 0.0, 0.0, 0, 1,
               "../data", work_in_examples, "MODUL_path.txt", "res_ephil.txt")
    ref = np.loadtxt(Path(EXAMPLES_DIR) / work_in_examples
                     / "res_ephil_global-0.txt", usecols=(8, 14))

    mat = ModularMaterial(
        elasticity=IsotropicElasticity(E=210000., nu=0.3),
        mechanisms=[Plasticity(
            sigma_Y=300.,
            yield_criterion=HillYield(F=0.5, G=0.4, H=0.6, L=1.5, M=1.5, N=1.5),
            isotropic_hardening=PowerLawHardening(k=5000., m=1.0),
        )],
    )
    sim.solver(mat.umat_name, mat.props, mat.nstatev, 0.0, 0.0, 0.0, 0, 1,
               "../data", work_in_examples, "MODUL_path.txt", "res_mhill.txt")
    hist = np.loadtxt(Path(EXAMPLES_DIR) / work_in_examples
                      / "res_mhill_global-0.txt", usecols=(8, 14))

    assert hist.shape == ref.shape
    peak = np.max(np.abs(ref[:, 1]))
    assert peak > 500.0, "sanity: reached the hardened regime"
    assert np.max(np.abs(hist[:, 1] - ref[:, 1])) < 1e-6, (
        "MODUL Hill deviates from EPHIL (criterion params/dispatch regression?)"
    )
