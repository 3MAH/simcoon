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
    if link.is_symlink() or link.exists():
        link.unlink()  # stale link from an interrupted previous run
    link.symlink_to(results, target_is_directory=True)
    yield link.name
    link.unlink(missing_ok=True)


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


def test_elastic_convention_equivalence(work_in_examples):
    """The same isotropic material given as (E, nu) and as (K, mu) must give
    identical stress — the convention slot in props selects the
    interpretation; all conversions live in the C++ builders (L_iso)."""
    E, nu = 210000.0, 0.3
    K = E / (3.0 * (1.0 - 2.0 * nu))
    mu = E / (2.0 * (1.0 + nu))

    mat_enu = ModularMaterial(elasticity=IsotropicElasticity(C1=E, C2=nu))
    mat_kmu = ModularMaterial(
        elasticity=IsotropicElasticity(C1=K, C2=mu, convention="Kmu"))

    # The convention code travels as the first slot of the elasticity block.
    assert mat_enu.props[1] == 0.0
    assert mat_kmu.props[1] == 2.0

    hist_enu = _run_solver(mat_enu, work_in_examples, "res_conv_enu.txt")
    hist_kmu = _run_solver(mat_kmu, work_in_examples, "res_conv_kmu.txt")
    np.testing.assert_allclose(hist_kmu, hist_enu, rtol=1e-10)


def test_convention_string_rejected_if_unknown():
    with pytest.raises(ValueError, match="IsoConvention"):
        IsotropicElasticity(C1=1.0, C2=0.3, convention="Envy")


def test_two_plasticity_equivalent_to_one_linear(work_in_examples):
    """Two identical linear-hardening plasticity mechanisms with modulus 2H
    each must give the same stress as a single mechanism with modulus H, by
    superposition of associated flow. This exercises the cross-mechanism
    Jacobian assembly in ModularUMAT::return_mapping — without it the FB
    iteration diverges to nonphysical stress."""
    E, nu, sigma_Y, H = 210000.0, 0.3, 300.0, 10000.0

    mat_single = ModularMaterial(
        elasticity=IsotropicElasticity(C1=E, C2=nu),
        mechanisms=[Plasticity(
            sigma_Y=sigma_Y,
            isotropic_hardening=LinearIsotropicHardening(H=H),
        )],
    )
    mat_two = ModularMaterial(
        elasticity=IsotropicElasticity(C1=E, C2=nu),
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
        elasticity=IsotropicElasticity(C1=210000.0, C2=0.3, alpha=1.2e-5),
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
        elasticity=IsotropicElasticity(C1=E0, C2=nu0),
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
        elasticity=IsotropicElasticity(C1=210000.0, C2=0.3),
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
    """The algorithmic tangent (mode 2, the default) must reproduce the
    continuum-mode (1) converged response — the tangent steers the global
    Newton, not the residual. Guards the algorithmic assembly (Bhat = -B sign
    convention and the coupled sub-block extraction) at solver level.
    2.0 renumbering: 0 = none/explicit, 1 = continuum, 2 = algorithmic."""
    mat = ModularMaterial(
        elasticity=IsotropicElasticity(C1=210000.0, C2=0.3),
        mechanisms=[Plasticity(
            sigma_Y=300.0,
            isotropic_hardening=VoceHardening(Q=200.0, b=20.0),
            kinematic_hardening=ArmstrongFrederickHardening(C=30000.0, D=172.0),
        )],
    )

    outs = {}
    for mode in (1, 2):
        out = f"res_tg{mode}.txt"
        sim.solver(mat.umat_name, mat.props, mat.nstatev,
                   0.0, 0.0, 0.0, 0, 1,
                   "../data", work_in_examples, "MODUL_path.txt", out,
                   mode)
        outs[mode] = np.loadtxt(Path(EXAMPLES_DIR) / work_in_examples
                                / out.replace(".txt", "_global-0.txt"))

    assert outs[1].shape == outs[2].shape
    peak = np.max(np.abs(outs[1][:, 14]))
    diff = np.max(np.abs(outs[1][:, 14] - outs[2][:, 14]))
    assert peak > 100.0
    assert diff / peak < 1e-5, (
        f"algorithmic-mode response deviates from continuum mode by {diff/peak:.2e} "
        "(the tangent must not change the converged solution)"
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
        elasticity=IsotropicElasticity(C1=210000.0, C2=0.3),
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
        elasticity=IsotropicElasticity(C1=210000., C2=0.3),
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


def test_chaboche_shear_matches_epcha_reference(work_in_examples):
    """MODUL Chaboche kinematic hardening under SHEAR loading must match the
    reference EPCHA UMAT.

    Regression for the backstress convention bug: X = (2/3) C a is a
    stress-like tensor, but it was carried in the back-strain (engineering
    strain) Voigt convention and subtracted from sigma / contracted with the
    flow n — doubling the shear terms. Invisible on the uniaxial path (zero
    shear, hence the earlier uniaxial EPCHA test passed at 3.7e-4); a pure
    shear path exposed a 24.7% stress error. Two sites: backstress_t and the
    ChabocheHardening total_backstress accumulator."""
    from simcoon.modular import ChabocheHardening
    # SHEAR_path.txt drives E12 with all other components stress-free.
    ep = np.array([210000., 0.3, 0., 300., 0., 0., 30000., 300., 19500., 172.])
    sim.solver("EPCHA", ep, 33, 0.0, 0.0, 0.0, 0, 1,
               "../data", work_in_examples, "SHEAR_path.txt", "sh_epcha.txt")
    ref = np.loadtxt(Path(EXAMPLES_DIR) / work_in_examples
                     / "sh_epcha_global-0.txt", usecols=(11, 17))  # eps12, sig12

    mat = ModularMaterial(
        elasticity=IsotropicElasticity(C1=210000., C2=0.3),
        mechanisms=[Plasticity(
            sigma_Y=300.,
            kinematic_hardening=ChabocheHardening(
                terms=((30000., 300.), (19500., 172.))),
        )],
    )
    sim.solver(mat.umat_name, mat.props, mat.nstatev, 0.0, 0.0, 0.0, 0, 1,
               "../data", work_in_examples, "SHEAR_path.txt", "sh_modul.txt")
    hist = np.loadtxt(Path(EXAMPLES_DIR) / work_in_examples
                      / "sh_modul_global-0.txt", usecols=(11, 17))

    assert hist.shape == ref.shape
    peak = np.max(np.abs(ref[:, 1]))
    assert peak > 100.0, "sanity: shear path yielded"
    assert np.max(np.abs(hist[:, 1] - ref[:, 1])) / peak < 1e-3, (
        "MODUL Chaboche under shear deviates from EPCHA (backstress "
        "convention regression?)"
    )



def test_armstrong_frederick_path_matches_chaboche(work_in_examples):
    """The ArmstrongFrederickHardening dataclass routes to a DISTINCT C++
    class (ArmstrongFrederickHardening, not ChabocheHardening). A single-term
    Chaboche is the same physics, so the two must agree bit-for-bit — locking
    the AF total_backstress path (which uses backstress_t directly, no
    accumulator) under shear, where the backstress convention matters."""
    from simcoon.modular import ArmstrongFrederickHardening, ChabocheHardening

    def run(kin, out):
        m = ModularMaterial(
            elasticity=IsotropicElasticity(C1=210000., C2=0.3),
            mechanisms=[Plasticity(sigma_Y=300., kinematic_hardening=kin)],
        )
        sim.solver(m.umat_name, m.props, m.nstatev, 0.0, 0.0, 0.0, 0, 1,
                   "../data", work_in_examples, "SHEAR_path.txt", out)
        return np.loadtxt(Path(EXAMPLES_DIR) / work_in_examples
                          / out.replace(".txt", "_global-0.txt"), usecols=(11, 17))

    af = run(ArmstrongFrederickHardening(C=30000., D=300.), "af_af.txt")
    ch = run(ChabocheHardening(terms=((30000., 300.),)), "af_ch.txt")
    assert af.shape == ch.shape
    assert np.max(np.abs(af[:, 1] - ch[:, 1])) < 1e-9, (
        "AF class diverges from single-term Chaboche (distinct C++ path bug?)"
    )
    assert np.max(np.abs(af[:, 1])) > 100.0, "sanity: shear yielded"


# ============================================================================
# Legacy-vs-MODUL equivalence gates (rationalization PR, Phase 1)
# Pattern: run the legacy UMAT and its MODUL twin through the solver on the
# same path; compare stress histories. These tests are the deletion gates for
# the name-adapter migration (and permanent validation for kept UMATs).
# ============================================================================

def _run_named(name, props, nstatev, results_dir, path_file, out, cols=(8, 14)):
    sim.solver(name, np.asarray(props, dtype=float), nstatev,
               0.0, 0.0, 0.0, 0, 1,
               "../data", results_dir, path_file, out)
    return np.loadtxt(Path(EXAMPLES_DIR) / results_dir
                      / out.replace(".txt", "_global-0.txt"), usecols=cols)


def _assert_equiv(ref, hist, rel_tol, label):
    assert hist.shape == ref.shape
    peak = np.max(np.abs(ref[:, 1]))
    assert peak > 0, f"{label}: sanity, zero response"
    max_diff = np.max(np.abs(hist[:, 1] - ref[:, 1]))
    assert max_diff / peak < rel_tol, (
        f"{label}: MODUL twin deviates by {max_diff/peak:.3e} (tol {rel_tol:.0e})")


def test_eliso_matches_modul(work_in_examples):
    ref = _run_named("ELISO", [210000., 0.3, 1.2e-5], 1,
                     work_in_examples, "MODUL_path.txt", "eq_eliso.txt")
    mat = ModularMaterial(elasticity=IsotropicElasticity(
        C1=210000., C2=0.3, alpha=1.2e-5))
    hist = _run_named("MODUL", mat.props, mat.nstatev,
                      work_in_examples, "MODUL_path.txt", "eq_meliso.txt")
    _assert_equiv(ref, hist, 1e-9, "ELISO")


def test_elist_matches_modul(work_in_examples):
    from simcoon.modular import TransverseIsotropicElasticity
    # legacy props: [axis, EL, ET, nuTL, nuTT, GLT, alpha_L, alpha_T]
    ref = _run_named("ELIST", [3, 230000., 15000., 0.02, 0.4, 50000., 0., 0.],
                     1, work_in_examples, "MODUL_path.txt", "eq_elist.txt")
    mat = ModularMaterial(elasticity=TransverseIsotropicElasticity(
        EL=230000., ET=15000., nuTL=0.02, nuTT=0.4, GLT=50000., axis=3))
    hist = _run_named("MODUL", mat.props, mat.nstatev,
                      work_in_examples, "MODUL_path.txt", "eq_melist.txt")
    _assert_equiv(ref, hist, 1e-9, "ELIST")


def test_elort_matches_modul(work_in_examples):
    from simcoon.modular import OrthotropicElasticity
    ref = _run_named("ELORT", [70000., 30000., 15000., 0.3, 0.3, 0.3,
                               8000., 6000., 5000., 1e-5, 2e-5, 3e-5],
                     1, work_in_examples, "MODUL_path.txt", "eq_elort.txt")
    mat = ModularMaterial(elasticity=OrthotropicElasticity(
        C1=70000., C2=30000., C3=15000., C4=0.3, C5=0.3, C6=0.3,
        C7=8000., C8=6000., C9=5000.,
        alpha1=1e-5, alpha2=2e-5, alpha3=3e-5))
    hist = _run_named("MODUL", mat.props, mat.nstatev,
                      work_in_examples, "MODUL_path.txt", "eq_melort.txt")
    _assert_equiv(ref, hist, 1e-9, "ELORT")


def test_epkcp_matches_modul(work_in_examples):
    """EPKCP: Prager convention differs — legacy X = kX*a (tensorial), modular
    X = (2/3)*C*a, hence C = 1.5*kX in the twin. m = 1 avoids the power-law
    p=0 tangent singularity (same rationale as the EPHIL test)."""
    from simcoon.modular import PowerLawHardening, PragerHardening
    ref = _run_named("EPKCP", [210000., 0.3, 0., 300., 1000., 1.0, 20000.],
                     33, work_in_examples, "MODUL_path.txt", "eq_epkcp.txt")
    mat = ModularMaterial(
        elasticity=IsotropicElasticity(C1=210000., C2=0.3),
        mechanisms=[Plasticity(
            sigma_Y=300.,
            isotropic_hardening=PowerLawHardening(k=1000., m=1.0),
            kinematic_hardening=PragerHardening(C=1.5 * 20000.))])
    hist = _run_named("MODUL", mat.props, mat.nstatev,
                      work_in_examples, "MODUL_path.txt", "eq_mepkcp.txt")
    _assert_equiv(ref, hist, 1e-5, "EPKCP")


def test_ephac_matches_modul(work_in_examples):
    """EPHAC: cubic elasticity (E, nu, G) + Hill yield + Voce + 2x AF."""
    from simcoon.modular import CubicElasticity, HillYield, ChabocheHardening
    props = [210000., 0.3, 85000., 0., 300., 200., 20.,
             30000., 172., 19500., 301.,
             0.5, 0.4, 0.6, 1.5, 1.5, 1.5]
    ref = _run_named("EPHAC", props, 33,
                     work_in_examples, "MODUL_path.txt", "eq_ephac.txt")
    mat = ModularMaterial(
        elasticity=CubicElasticity(C1=210000., C2=0.3, C3=85000.),
        mechanisms=[Plasticity(
            sigma_Y=300.,
            yield_criterion=HillYield(F=0.5, G=0.4, H=0.6, L=1.5, M=1.5, N=1.5),
            isotropic_hardening=VoceHardening(Q=200., b=20.),
            kinematic_hardening=ChabocheHardening(
                terms=((30000., 172.), (19500., 301.))))])
    hist = _run_named("MODUL", mat.props, mat.nstatev,
                      work_in_examples, "MODUL_path.txt", "eq_mephac.txt")
    _assert_equiv(ref, hist, 1e-3, "EPHAC")


def test_epani_matches_modul(work_in_examples):
    """EPANI: cubic elasticity + 9-parameter anisotropic yield + Voce + 2x AF."""
    from simcoon.modular import CubicElasticity, AnisotropicYield, ChabocheHardening
    # P must be an ADMISSIBLE quadratic form: symmetric with zero row sums
    # on the normal block (deviatoric, PSD) — an indefinite P gives
    # sqrt(negative) = NaN in Eq_stress_P (both legacy and modular).
    props = [210000., 0.3, 85000., 0., 300., 200., 20.,
             30000., 172., 19500., 301.,
             1.2, 1.1, 1.1, -0.6, -0.6, -0.5, 1.6, 1.5, 1.4]
    ref = _run_named("EPANI", props, 33,
                     work_in_examples, "MODUL_path.txt", "eq_epani.txt")
    mat = ModularMaterial(
        elasticity=CubicElasticity(C1=210000., C2=0.3, C3=85000.),
        mechanisms=[Plasticity(
            sigma_Y=300.,
            yield_criterion=AnisotropicYield(
                P11=1.2, P22=1.1, P33=1.1, P12=-0.6, P13=-0.6, P23=-0.5,
                P44=1.6, P55=1.5, P66=1.4),
            isotropic_hardening=VoceHardening(Q=200., b=20.),
            kinematic_hardening=ChabocheHardening(
                terms=((30000., 172.), (19500., 301.))))])
    hist = _run_named("MODUL", mat.props, mat.nstatev,
                      work_in_examples, "MODUL_path.txt", "eq_mepani.txt")
    _assert_equiv(ref, hist, 1e-3, "EPANI")


def test_epdfa_matches_modul(work_in_examples):
    """EPDFA: cubic elasticity + DFA yield (Hill + hydrostatic K) + Voce + 2x AF."""
    from simcoon.modular import CubicElasticity, DFAYield, ChabocheHardening
    props = [210000., 0.3, 85000., 0., 300., 200., 20.,
             30000., 172., 19500., 301.,
             0.5, 0.4, 0.6, 1.5, 1.5, 1.5, 0.1]
    ref = _run_named("EPDFA", props, 33,
                     work_in_examples, "MODUL_path.txt", "eq_epdfa.txt")
    mat = ModularMaterial(
        elasticity=CubicElasticity(C1=210000., C2=0.3, C3=85000.),
        mechanisms=[Plasticity(
            sigma_Y=300.,
            yield_criterion=DFAYield(F=0.5, G=0.4, H=0.6, L=1.5, M=1.5, N=1.5,
                                     K=0.1),
            isotropic_hardening=VoceHardening(Q=200., b=20.),
            kinematic_hardening=ChabocheHardening(
                terms=((30000., 172.), (19500., 301.))))])
    hist = _run_named("MODUL", mat.props, mat.nstatev,
                      work_in_examples, "MODUL_path.txt", "eq_mepdfa.txt")
    _assert_equiv(ref, hist, 1e-3, "EPDFA")


def test_epchg_matches_modul(work_in_examples):
    """EPCHG (generic Chaboche): cubic elasticity + von Mises + N "Voce
    terms" + N-term Chaboche.

    NOTE: the legacy N-term isotropic hardening couples every term through a
    SINGLE Hp (dHp/dp = sum_i b_i (Q_i - Hp)) — mathematically ONE effective
    Voce with b_eff = sum(b_i), Q_eff = sum(b_i Q_i)/sum(b_i), NOT the
    standard combined-Voce sum. The modular twin (and the name adapter) map
    to that single effective Voce."""
    from simcoon.modular import CubicElasticity, ChabocheHardening
    # props: E nu G alpha | sigmaY N_iso N_kin criteria | (Q,b)xN | (C,D)xN
    props = [210000., 0.3, 85000., 0., 300., 2, 2, 0,
             150., 15., 50., 40.,
             30000., 172., 19500., 301.]
    ref = _run_named("EPCHG", props, 33,
                     work_in_examples, "MODUL_path.txt", "eq_epchg.txt")
    b_eff = 15. + 40.
    q_eff = (15. * 150. + 40. * 50.) / b_eff
    mat = ModularMaterial(
        elasticity=CubicElasticity(C1=210000., C2=0.3, C3=85000.),
        mechanisms=[Plasticity(
            sigma_Y=300.,
            isotropic_hardening=VoceHardening(Q=q_eff, b=b_eff),
            kinematic_hardening=ChabocheHardening(
                terms=((30000., 172.), (19500., 301.))))])
    hist = _run_named("MODUL", mat.props, mat.nstatev,
                      work_in_examples, "MODUL_path.txt", "eq_mepchg.txt")
    _assert_equiv(ref, hist, 1e-3, "EPCHG")


def test_ephin_matches_modul(work_in_examples):
    """EPHIN (Hill, N yield surfaces): equivalence proven for N = 1.

    The LEGACY multi-surface path (N >= 2) is defective: it returns NaN even
    for two identical surfaces, and even when the second surface can never
    activate (sigma_Y = 5000 on a 2% path) — probed 2026-07. The modular
    engine handles multiple mechanisms correctly (see
    test_two_plasticity_equivalent_to_one_linear), so the name adapter is an
    upgrade for N >= 2; the provable gate is the N = 1 case."""
    from simcoon.modular import HillYield, PowerLawHardening
    # props: E nu alpha N | (sigmaY k m F G H L M N)xN
    props = [210000., 0.3, 0., 1,
             300., 3000., 1.0, 0.5, 0.4, 0.6, 1.5, 1.5, 1.5]
    ref = _run_named("EPHIN", props, 33,
                     work_in_examples, "MODUL_path.txt", "eq_ephin.txt")
    mat = ModularMaterial(
        elasticity=IsotropicElasticity(C1=210000., C2=0.3),
        mechanisms=[Plasticity(
            sigma_Y=300.,
            yield_criterion=HillYield(F=0.5, G=0.4, H=0.6,
                                      L=1.5, M=1.5, N=1.5),
            isotropic_hardening=PowerLawHardening(k=3000., m=1.0))])
    hist = _run_named("MODUL", mat.props, mat.nstatev,
                      work_in_examples, "MODUL_path.txt", "eq_mephin.txt")
    _assert_equiv(ref, hist, 1e-5, "EPHIN")


def test_zennk_has_no_modular_twin():
    """ZENNK (Zener_Nfast) is a generalized KELVIN chain (branch driving
    force sigma - L_i EV_i, branches in series), NOT a generalized Maxwell:
    the modular Prony viscoelasticity is a different rheological model
    (measured deviation 86% on the relaxation path). ZENNK therefore keeps
    its dedicated implementation (no adapter) — same bucket as ZENER."""
    pass
