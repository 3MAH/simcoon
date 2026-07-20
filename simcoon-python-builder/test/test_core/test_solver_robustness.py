"""Solver Newton divergence guards.

Stress-controlled UNLOADING through the plastic->elastic branch flip used to
diverge geometrically when the hardening tangent is soft and SATURATING: the
predictor (elastoplastic tangent) overshoots past the reverse-yield span, the
trial re-yields in compression, and the undamped Newton oscillates with
exploding corrections (traced DE11: 8e-3 -> 0.83 -> 433 -> 4e13) while the
residual stays bounded (the yield surface caps the stress). The run died on a
singular inverse (tangent assembly / RU_decomposition), uncaught by the
step-cut machinery.

Guards under test (solver.cpp, tangent_assembly.cpp, modular_umat.cpp):
- Newton correction-norm cap -> step cut before the box integrates garbage;
- singular global Jacobian K -> step cut, loud at the minimal fraction;
- exception_inv / exception_det from kinematics or tangent assembly -> step
  cut (rethrown only at Dn_mini);
- ModularUMAT rejects pathological plastic multipliers (dp > 1 per increment,
  or committed negative dp) with tnew_dt = 0.5 instead of committing a
  hardening-saturated garbage state.

Control type 1 on purpose: reproduces without the finite-strain MODUL
registration, and EPICP shows the failure was never modular-specific.
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


def test_epicp_soft_hardening_stress_unload_cycle(tmp_path):
    """EPICP with a soft power law fails the same way (the divergence was
    never modular-specific) — the solver guards must rescue it too."""
    hist = _run_stress_cycle(tmp_path, "EPICP",
                             [210000.0, 0.3, 0.0, 300.0, 200.0, 0.3], 8,
                             [400.0, 0.0])
    final = hist[-1]
    assert abs(final[C_TIME] - 2.0) < 1e-6
    assert abs(final[S_STRESS][0]) < 1e-3
    assert final[S_WM][3] > 1.0
