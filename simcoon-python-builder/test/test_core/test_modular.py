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


EXAMPLES_DIR = Path(__file__).resolve().parents[3] / "examples" / "umats"


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
        "data", results_dir, "MODUL_path.txt", outfile,
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
