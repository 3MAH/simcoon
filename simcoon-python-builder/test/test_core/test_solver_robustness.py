"""Solver Newton divergence guards.

Stress-controlled UNLOADING through the plastic->elastic branch flip used to
diverge geometrically when the hardening tangent is soft and SATURATING: the
predictor (elastoplastic tangent) overshoots past the reverse-yield span, the
trial re-yields in compression, and the undamped Newton oscillates with
exploding corrections (traced DE11: 8e-3 -> 0.83 -> 433 -> 4e13). The run died
on a singular inverse (tangent assembly / RU_decomposition).

Guards under test:
- ModularUMAT rejects a pathological plastic multiplier (dp > 1 per increment)
  with tnew_dt = 0.5 instead of committing a hardening-saturated garbage state,
  so the solver bisects the increment away (modular_umat.cpp);
- tangent assembly falls back to the elastic operator on a degenerate inverse
  instead of aborting (tangent_assembly.cpp);
- singular global Jacobian K -> step cut, then the existing inforce path at the
  minimal fraction (never a hard throw); exception_inv / exception_det from the
  kinematics -> step cut, rethrown only at Dn_mini (solver.cpp).

Control type 1 on purpose: reproduces without the finite-strain MODUL
registration. The fix is in the shared solver Newton loop, so it is not
modular-specific (see the note below the test on why the EPICP illustration
was dropped).
"""

from simcoon.modular import (
    ModularMaterial,
    IsotropicElasticity,
    Plasticity,
    VonMisesYield,
    VoceHardening,
)
import numpy as np
import pytest

from simcoon.solver import Block, StepMeca, solve
from solver_harness import C_TIME, S_STRESS, S_WM, path_file, run_path


def _run_stress_cycle(base_dir, umat_name, props, nstatev, targets):
    """Small-strain (ct1), fully stress-controlled uniaxial cycle."""
    return run_path(base_dir, umat_name, props, nstatev, 1,
                    path_file([("S", t) for t in targets], 1))


def test_modul_voce_stress_unload_cycle(tmp_path):
    """MODUL + saturating Voce, load past yield then stress-unload to zero:
    used to die on 'inv(): matrix is singular'. Must now complete the cycle
    with the stress driven back to zero and positive plastic dissipation."""
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
    hist = _run_stress_cycle(tmp_path, mat.umat_name, mat.props, mat.nstatev,
                             [400.0, 0.0])
    final = hist[-1]
    assert abs(final[C_TIME] - 2.0) < 1e-6
    assert abs(final[S_STRESS][0]) < 1e-3     # unloaded to zero stress
    assert final[S_WM][3] > 1.0               # plastic dissipation happened

# NOTE: an EPICP (legacy CCP) soft-hardening variant of this test was dropped.
# The crash fix lives in the shared solver Newton loop (not the modular engine),
# so it is not modular-specific — but at the EPICP soft-hardening extreme
# (k=200, m=0.3) the stress-unload branch flip is genuinely non-convergent and
# resolves differently across LAPACK backends (completes on macOS, does not on
# Linux/Windows). Asserting convergence there tests the platform, not the fix;
# the modular Voce case above is the portable regression guard.


# ---------------------------------------------------------------------------
# a prescribed state the material cannot reach
# ---------------------------------------------------------------------------

E, NU, SIGMA_Y = 70000.0, 0.3, 300.0


def _cycle(target, ninc):
    """One stress-driven cycle to +/- target, two steps of unit duration."""
    return [
        Block(
            steps=[
                StepMeca(control=["stress"] * 6, value=[target, 0, 0, 0, 0, 0], ninc=ninc),
                StepMeca(control=["stress"] * 6, value=[-target, 0, 0, 0, 0, 0], ninc=ninc),
            ],
            ncycle=1,
        )
    ]


@pytest.mark.parametrize("ninc", [25, 100])
def test_unreachable_stress_target_aborts_instead_of_reporting_success(ninc):
    """A perfectly plastic law cannot be driven above its plateau.

    ``inforce`` is a *correction*: it closes an increment that did not quite converge
    and carries the residual into the next one, which absorbs it. When the target is
    physically out of reach the residual is never absorbed, so that path used to close
    every remaining increment and return status 0 on a state the run never reached --
    with a time axis inflated by ``1/Dn_mini`` on top, because ``DTime`` was only
    refreshed inside the branch that calls the UMAT.
    """
    props = [E, NU, 1.0e-5, SIGMA_Y, 0.0, 0.3]          # k = 0: no hardening
    blocks = _cycle(2.0 * SIGMA_Y, ninc)                # far above the plateau

    res = solve(blocks, "EPICP", props, 8, raise_on_abort=False)
    assert res.status != 0, "an unreachable target must not be reported as converged"
    # and the time axis stays inside the two unit-duration steps
    assert float(np.asarray(res["Time"])[-1]) <= 2.0 + 1e-9
    # the default contract turns that status into an exception
    with pytest.raises(RuntimeError):
        solve(blocks, "EPICP", props, 8)


@pytest.mark.parametrize("ninc", [25, 100])
def test_reachable_stress_target_is_unaffected(ninc):
    """The guard must not disturb a target the material can reach.

    Same cycle, once inside the plateau of the perfectly plastic law and once on a
    hardening law that reaches it: both run to completion with an exact time axis.
    """
    for props, target in (
        ([E, NU, 1.0e-5, SIGMA_Y, 0.0, 0.3], 0.8 * SIGMA_Y),      # below the plateau
        ([E, NU, 1.0e-5, SIGMA_Y, 1000.0, 0.3], 2.0 * SIGMA_Y),   # hardening reaches it
    ):
        res = solve(_cycle(target, ninc), "EPICP", props, 8)
        assert res.status == 0
        assert float(np.asarray(res["Time"])[-1]) == pytest.approx(2.0, abs=1e-9)
        assert res["Stress"][0, -1] == pytest.approx(-target, rel=1e-6)
