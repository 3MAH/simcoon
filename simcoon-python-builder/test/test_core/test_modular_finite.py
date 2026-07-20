"""MODUL under finite strain (NLGEOM): Hencky hyperelastic composition.

Registered in the finite dispatcher on the log-strain measures with a hard
corate_type == 3 (log_R) requirement: only there does the accumulated
corotational strain equal ln V exactly, making sigma = L:(ln V - sum eps_inel)
a genuine stored-energy law. These tests run the file solver on self-contained
temp path files (control_type 3: log-strain / Kirchhoff conjugate pair).
"""

import numpy as np
import pytest

from simcoon.modular import (
    ModularMaterial,
    IsotropicElasticity,
    Plasticity,
    VonMisesYield,
    VoceHardening,
)
from solver_harness import (
    C_TIME,
    S_STRAIN,
    S_STRESS,
    S_WM,
    path_file,
    run_path,
)


def _run_finite(base_dir, umat_name, props, nstatev, corate, targets,
                control_type=3):
    return run_path(base_dir, umat_name, props, nstatev, corate,
                    path_file(targets, control_type))


def _elastic_modul():
    return ModularMaterial(
        elasticity=IsotropicElasticity(C1=100000.0, C2=0.3, alpha=0.0,
                                       convention="Enu"),
    )


def test_modul_finite_matches_legacy_eliso(tmp_path):
    """Elastic-only MODUL and legacy ELISO (modular adapter) are the same
    log-strain Kirchhoff box; under NLGEOM ct3 + corate 3 their histories must
    coincide, and the controlled Kirchhoff stress must map to e11 = tau11/E
    exactly (pins the Kirchhoff-route output: a spurious Cauchy->Kirchhoff
    conversion would scale the response by J)."""
    mat = _elastic_modul()
    h_mod = _run_finite(tmp_path / "mod", mat.umat_name, mat.props,
                        mat.nstatev, 3, [("S", 15000.0)])
    h_leg = _run_finite(tmp_path / "leg", "ELISO", [100000.0, 0.3, 0.0], 1, 3,
                        [("S", 15000.0)])

    assert abs(h_mod[-1, C_TIME] - 1.0) < 1e-6
    assert abs(h_leg[-1, C_TIME] - 1.0) < 1e-6
    np.testing.assert_allclose(h_mod[-1, S_STRESS], h_leg[-1, S_STRESS],
                               atol=1e-6 * 15000.0)
    np.testing.assert_allclose(h_mod[-1, S_STRAIN], h_leg[-1, S_STRAIN],
                               atol=1e-8)
    # Hencky elasticity: log strain e11 = tau11 / E, exactly, at finite stretch
    assert abs(h_mod[-1, S_STRAIN][0] - 0.15) < 1e-6


def test_modul_finite_hyperelastic_closed_cycle(tmp_path):
    """Load to 15% log strain and back to zero stress: a hyperelastic law
    leaves no residual strain and no residual work. This is the hyper/hypo
    consistency property corate 3 buys — a Jaumann-style rate would not
    return to zero."""
    mat = _elastic_modul()
    hist = _run_finite(tmp_path, mat.umat_name, mat.props, mat.nstatev, 3,
                       [("S", 15000.0), ("S", 0.0)])

    assert abs(hist[-1, C_TIME] - 2.0) < 1e-6
    peak_e11 = np.max(np.abs(hist[:, S_STRAIN][:, 0]))
    assert peak_e11 > 0.14  # genuinely finite strain
    final = hist[-1]
    assert np.max(np.abs(final[S_STRAIN])) < 1e-8
    assert np.max(np.abs(final[S_STRESS])) < 1e-4
    assert np.max(np.abs(final[S_WM])) < 1e-6 * 0.5 * 15000.0 * 0.15


def test_modul_finite_rejects_non_log_corate(tmp_path):
    """Any corate other than 3 (log_R) degrades the composition to a
    non-integrable hypoelastic rate: rejected up front."""
    mat = _elastic_modul()
    for corate in (0, 1, 2, 5):
        with pytest.raises(RuntimeError, match="corate_type = 3"):
            _run_finite(tmp_path / f"co{corate}", mat.umat_name, mat.props,
                        mat.nstatev, corate, [("S", 15000.0)])


def test_modul_finite_plasticity_dissipates(tmp_path):
    """Voce elasto-plasticity under NLGEOM: a strain-controlled cycle to 5%
    log strain and back yields on loading, re-yields in compression on the
    way back, and leaves strictly positive dissipation.

    Strain-controlled on purpose: this test pins the plasticity physics
    deterministically; the stress-controlled unloading path (which needs the
    solver's Newton divergence guards) is covered by
    test_solver_robustness.py."""
    mat = ModularMaterial(
        elasticity=IsotropicElasticity(C1=210000.0, C2=0.3, alpha=0.0,
                                       convention="Enu"),
        mechanisms=[
            Plasticity(
                sigma_Y=300.0,
                yield_criterion=VonMisesYield(),
                isotropic_hardening=VoceHardening(Q=200.0, b=10.0),
            ),
        ],
    )
    hist = _run_finite(tmp_path, mat.umat_name, mat.props, mat.nstatev, 3,
                       [("E", 0.05), ("E", 0.0)])

    assert abs(hist[-1, C_TIME] - 2.0) < 1e-6
    peak_s11 = np.max(hist[:, S_STRESS][:, 0])
    assert peak_s11 > 300.0                    # yielded on loading
    final = hist[-1]
    assert abs(final[S_STRAIN][0]) < 1e-8      # strain driven back to zero
    assert final[S_STRESS][0] < -300.0         # reverse yield in compression
    assert final[S_WM][3] > 0.0                # Wm_d > 0
    # No kinematic hardening -> no backstress storage -> Wm_ir stays 0
    assert abs(final[S_WM][2]) < 1e-12
