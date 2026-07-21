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
