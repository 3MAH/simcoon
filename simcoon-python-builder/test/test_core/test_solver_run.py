"""Tests for the in-memory solver API (simcoon.solver / _core.solver_run).

These used to compare every run against the file-driven solver `sim._core.solver`,
which ran the same C++ engine through solver_file_sink. That binding left with the
2.0 JSON-only migration, and an equivalence gate between two implementations is
vacuous once only one remains. The cases were kept — control types 2/3/4 crossed
with corates, spin, sinusoidal and tabular modes, thermomechanical blocks, legacy
path parsing — and their assertions rewritten as invariants that hold on their own:
the prescribed state is reached, the response is finite and non-trivial, and what
the physics imposes (free thermal expansion, cycle bookkeeping, plastic flow,
the rotation history) is checked directly.
"""

import json
import subprocess
import sys
import textwrap
from pathlib import Path

import numpy as np
import pytest

import simcoon as sim
from simcoon.modular import elastic_model
from simcoon.solver import Block, StepMeca, StepThermomeca, solve

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def assert_ran_and_responded(res, n_expected=None):
    """The run completed and produced a finite, non-trivial mechanical response."""
    assert res.status == 0
    if n_expected is not None:
        assert len(res) == n_expected
    assert np.isfinite(res["Stress"]).all()
    assert np.isfinite(res["Strain"]).all()
    # the PEAK of the history, not its last increment: several cases unload back to
    # zero stress by construction.
    assert np.abs(res["Stress"]).max() > 1.0



#: path-file component order (11, 12, 22, 13, 23, 33) -> Voigt index
_FILE_ORDER = [0, 3, 1, 4, 5, 2]


ELISO_PROPS = [70000.0, 0.3, 1.0e-5]
EPICP_PROPS = [70000.0, 0.3, 1.0e-5, 300.0, 1000.0, 0.3]  # E nu alpha sigmaY k m
EPICP_NSTATEV = 8
#: a purely elastic MODUL: what translate_ELISO makes of ELISO_PROPS
MODUL_ELASTIC_PROPS = elastic_model(*ELISO_PROPS).props.tolist()
SNTVE_PROPS = [70000.0, 0.3, 1.0e-5]

_UNIAXIAL = ["strain"] + ["stress"] * 5


# ---------------------------------------------------------------------------
# small strain
# ---------------------------------------------------------------------------


def test_eliso_uniaxial_analytic():
    E, nu, _ = ELISO_PROPS
    step = StepMeca(control=_UNIAXIAL, value=[0.01, 0, 0, 0, 0, 0], ninc=50)
    res = solve(step, "ELISO", ELISO_PROPS, 1, T_init=290.0)
    assert res.status == 0
    np.testing.assert_allclose(res["Stress"][0, -1], E * 0.01, rtol=1e-10)
    np.testing.assert_allclose(res["Strain"][1, -1], -nu * 0.01, atol=1e-12)
    np.testing.assert_allclose(res["Stress"][1:, -1], 0.0, atol=1e-8)
    # tangent is the isotropic stiffness
    L = sim.L_iso([E, nu], "Enu")
    np.testing.assert_allclose(res["TangentMatrix"][:, :, -1], L, rtol=1e-8)


def test_eliso_load_unload():
    E = ELISO_PROPS[0]
    s1 = StepMeca(control=_UNIAXIAL, value=[0.02, 0, 0, 0, 0, 0], ninc=100)
    s2 = StepMeca(control=_UNIAXIAL, value=[0.0, 0, 0, 0, 0, 0], ninc=100)
    res = solve(Block(steps=[s1, s2]), "ELISO", ELISO_PROPS, 1, T_init=290.0)

    assert_ran_and_responded(res, n_expected=200)
    # elastic load then unload: the peak follows Hooke and nothing is left behind
    np.testing.assert_allclose(res["Stress"][0].max(), E * 0.02, rtol=1e-8)
    np.testing.assert_allclose(res["Strain"][0, -1], 0.0, atol=1e-12)
    np.testing.assert_allclose(res["Stress"][0, -1], 0.0, atol=1e-6)


def test_epicp_mixed_cyclic():
    # stress-controlled uniaxial cycling into the plastic range
    s1 = StepMeca(control="stress", value=[400.0, 0, 0, 0, 0, 0], ninc=100)
    s2 = StepMeca(control="stress", value=[0.0, 0, 0, 0, 0, 0], ninc=100)
    res = solve(Block(steps=[s1, s2], ncycle=2), "EPICP", EPICP_PROPS,
                EPICP_NSTATEV, T_init=290.0)

    assert_ran_and_responded(res, n_expected=2 * 2 * 100)
    # cycle index bookkeeping
    assert res["Cycle"].max() == 1
    # the prescribed stress is reached at each peak, and plasticity happened
    np.testing.assert_allclose(res["Stress"][0].max(), 400.0, atol=1e-4)
    assert res["Statev"][0].max() > 1.0e-4


# ---------------------------------------------------------------------------
# finite strain: control types 2/3/4 x corates
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("umat,props,nstatev", [
    ("ELISO", ELISO_PROPS, 1),
    ("EPICP", EPICP_PROPS, EPICP_NSTATEV),
    ("SNTVE", SNTVE_PROPS, 1),
])
@pytest.mark.parametrize("control_type", [2, 3, 4])
@pytest.mark.parametrize("corate", [0, 2, 5])
def test_finite_strain_controls(umat, props, nstatev, control_type, corate):
    if umat == "EPICP" and control_type == 4 and corate == 5:
        # pre-existing engine limitation (identical for the file solver): plasticity
        # under strain-controlled Biot loading with the log_F rate fails in
        # logarithmic_F (singular DF inversion)
        pytest.skip("EPICP + Biot control + logarithmic_F: known engine limitation")
    # Biot control drives the STRETCH U11 itself, not a strain: asking for 0.08 there
    # would prescribe a 92 % compression (and the engine would deliver it).
    if control_type == 4:
        target = 1.05 if umat == "SNTVE" else 1.08
    elif umat == "SNTVE":
        target = 0.05  # keep the hyperelastic case well-conditioned
    else:
        target = 0.08
    step = StepMeca(control=_UNIAXIAL, value=[target, 0, 0, 0, 0, 0], ninc=50,
                    BC_w=np.zeros((3, 3)))
    ct = {2: "green_lagrange", 3: "logarithmic", 4: "biot"}[control_type]
    res = solve(Block(steps=[step], control_type=ct), umat, props, nstatev,
                T_init=290.0, corate=corate)

    assert_ran_and_responded(res, n_expected=50)
    # The target is prescribed in the measure of THIS control type — Green-Lagrange,
    # logarithmic or Biot — so the invariant is put on the stretch itself, which every
    # measure agrees on: the bar elongates, monotonically.
    F11 = res["F"][0, 0]
    assert F11[-1] > 1.0
    assert np.all(np.diff(F11) > -1e-12)
    # lateral faces stay stress-free
    assert np.abs(res["Stress"][1:, -1]).max() < 1e-3 * max(abs(res["Stress"][0, -1]), 1.0)


def test_ct3_spin():
    # logarithmic control with a superimposed rotation rate (BC_w)
    BC_w = np.array([[0.0, 0.2, 0.0], [-0.2, 0.0, 0.0], [0.0, 0.0, 0.0]])
    step = StepMeca(control=_UNIAXIAL, value=[0.05, 0, 0, 0, 0, 0], ninc=50, BC_w=BC_w)
    res = solve(Block(steps=[step], control_type="logarithmic"), "ELISO",
                ELISO_PROPS, 1, T_init=290.0, corate="logarithmic")

    assert_ran_and_responded(res, n_expected=50)
    # The target is a LOG strain, and it is prescribed in the COROTATIONAL frame: with a
    # superimposed spin its 11 component in the lab frame drifts (0.0491 here, 1.8 % off,
    # of order theta^2/2). The invariant is the principal strain, which no rotation moves.
    # ascontiguousarray: a column of a (6, N) history is not contiguous, and carma
    # refuses to borrow such an array.
    e_end = sim.v2t_strain(np.ascontiguousarray(res["LogStrain"][:, -1]))
    # Under the spin the state is not exactly uniaxial in the lab frame — the two
    # lateral principal strains differ (-0.015001 vs -0.014784) — so no simple
    # quantity equals the target exactly: the principal strain measures 0.049786.
    # The invariant is that the driven magnitude is delivered to ~0.5 %; a broken
    # rotation path moves it by far more.
    np.testing.assert_allclose(np.linalg.eigvalsh(e_end).max(), 0.05, rtol=5e-3)
    # the rotation history is captured, and R stays a rotation
    R_end = res["R"][:, :, -1]
    assert np.abs(R_end - np.eye(3)).max() > 1e-3
    np.testing.assert_allclose(R_end @ R_end.T, np.eye(3), atol=1e-10)


# ct5 is the only fully kinematic path (nK == 0) and the only caller of
# step_meca::generate_kin, whose incremental F is built from arma::logmat /
# arma::expmat -- code no other control type touches. Cover more than one
# target shape and corate through it: a diagonal target makes log(F_target)
# symmetric, which takes a different branch of expmat than the shear one.
# Ordered least-exotic first on purpose: a native crash ends the whole pytest
# process, so the cases that exercise the plainest code path must run before the
# ones that add a branch, otherwise a failure downstream hides everything.
# "stretch" keeps log(F_target) symmetric and its norm under the inverse
# scaling-and-squaring cutoff; "big_stretch" is still symmetric but crosses the
# cutoff, so the Denman-Beavers square-root loop runs; "shear" is the
# non-symmetric case.
@pytest.mark.parametrize("corate", ["jaumann", "logarithmic"])
@pytest.mark.parametrize("shape", ["stretch", "big_stretch", "shear"])
def test_ct5_F_variants(shape, corate):
    if shape == "shear":
        F_target = np.eye(3)
        F_target[0, 1] = 0.2
    elif shape == "big_stretch":
        F_target = np.diag([1.6, 1.0, 1.0])
    else:
        F_target = np.diag([1.05, 1.0, 1.0])

    step = StepMeca(control="F", value=F_target.ravel(), ninc=10)
    res = solve(Block(steps=[step], control_type="F"), "ELISO", ELISO_PROPS, 1,
                T_init=290.0, corate=corate)

    assert res.status == 0
    # the prescribed gradient is reached exactly at the end of the step
    np.testing.assert_allclose(res["F"][:, :, -1], F_target, atol=1e-10)
    stress = res["Stress"][:, -1]
    assert np.isfinite(stress).all()
    assert np.abs(stress).max() > 1.0


def test_strain_keys_at_finite_strain():
    """'Strain' is the logarithmic strain (alias 'LogStrain'); 'GreenLagrange' is E."""
    F_target = np.diag([1.6, 1.0, 1.0])
    step = StepMeca(control="F", value=F_target.ravel(), ninc=10)
    res = solve(Block(steps=[step], control_type="F"), "ELISO", ELISO_PROPS, 1,
                T_init=290.0, corate="logarithmic_R")
    assert res.status == 0
    assert res["Strain"] is res["LogStrain"]
    F = res["F"][:, :, -1]
    # the log strain is integrated increment by increment (corotational rate): the
    # endpoint carries the integration error of 10 increments, not a measure mix-up
    np.testing.assert_allclose(res["Strain"][0, -1], np.log(1.6), rtol=1e-3)
    np.testing.assert_allclose(res["GreenLagrange"][0, -1], 0.5 * (1.6 ** 2 - 1.0), atol=1e-10)
    E = 0.5 * (F.T @ F - np.eye(3))
    np.testing.assert_allclose(res["GreenLagrange"][:3, -1], np.diag(E), atol=1e-10)
    cols = res.to_dataframe().columns
    assert "Strain_11" in cols and "GreenLagrange_11" in cols and "LogStrain_11" not in cols


def test_ct5_F_control():
    # simple shear driven by the full deformation gradient
    F_target = np.eye(3)
    F_target[0, 1] = 0.2
    step = StepMeca(control="F", value=F_target.ravel(), ninc=50)
    res = solve(Block(steps=[step], control_type="F"), "SNTVE", SNTVE_PROPS, 1,
                T_init=290.0, corate="logarithmic")
    assert res.status == 0
    np.testing.assert_allclose(res["F"][0, 1, -1], 0.2, atol=1e-10)
    assert abs(res["Stress"][3, -1]) > 1.0  # shear stress developed


# ---------------------------------------------------------------------------
# loading modes
# ---------------------------------------------------------------------------

def test_sinusoidal_mode():
    step = StepMeca(control=_UNIAXIAL, value=[0.01, 0, 0, 0, 0, 0], ninc=100,
                    mode="sinusoidal")
    res = solve(step, "ELISO", ELISO_PROPS, 1, T_init=290.0)

    assert_ran_and_responded(res, n_expected=100)
    np.testing.assert_allclose(res["Strain"][0, -1], 0.01, rtol=1e-8)
    # sinusoidal profile: increments are not uniform
    de = np.diff(res["Strain"][0])
    assert de.max() / de.min() > 1.5


def test_tabular_memory():
    t = np.linspace(0.02, 1.0, 50)
    e11 = 0.015 * np.sin(np.pi * t)
    table = np.column_stack([t, e11])

    step = StepMeca(control=["strain"] + ["zero"] * 5, mode="tabular", tabular=table)
    res = solve(step, "ELISO", ELISO_PROPS, 1, T_init=290.0)

    assert_ran_and_responded(res, n_expected=len(t))
    # the table is followed exactly, in its own absolute time
    np.testing.assert_allclose(res["Time"], t, atol=1e-12)
    np.testing.assert_allclose(res["Strain"][0], e11, atol=1e-10)


# ---------------------------------------------------------------------------
# thermomechanical blocks
# ---------------------------------------------------------------------------

ELISO_T_PROPS = [1.0e-9, 1.0, 70000.0, 0.3, 1.0e-5]  # rho c_p E nu alpha


def test_thermomeca_temperature_ramp():
    step = StepThermomeca(control="stress", value=[0.0] * 6, ninc=50, T_final=340.0)
    res = solve(step, "ELISO", ELISO_T_PROPS, 1, T_init=290.0)

    assert res.status == 0
    assert len(res) == 50
    # free thermal expansion, and the ramp lands on its target
    np.testing.assert_allclose(res["Temp"][-1], 340.0, rtol=1e-10)
    np.testing.assert_allclose(res["Strain"][0, -1], 1.0e-5 * 50.0, rtol=1e-8)
    np.testing.assert_allclose(res["Stress"][:, -1], 0.0, atol=1e-6)
    assert "Q" in res and "Wt" in res and "dSdE" in res


def test_thermomeca_heat_flux_and_convection():
    # prescribed heat flux drives the temperature
    step = StepThermomeca(control="stress", value=[0.0] * 6, ninc=50,
                          thermal_control="heat_flux", Q=1.0e-6)
    res = solve(step, "ELISO", ELISO_T_PROPS, 1, T_init=290.0)
    assert res.status == 0
    assert res["Temp"][-1] != pytest.approx(290.0)

    # 0D convection relaxes towards T_init
    step = StepThermomeca(control="stress", value=[0.0] * 6, ninc=50,
                          thermal_control="convection", q_conv=1.0e-6)
    res = solve(step, "ELISO", ELISO_T_PROPS, 1, T_init=290.0)
    assert res.status == 0


def test_thermomeca_epicp():
    props = [1.0e-9, 1.0] + EPICP_PROPS  # rho c_p then mechanical props
    step = StepThermomeca(control=_UNIAXIAL, value=[0.01, 0, 0, 0, 0, 0], ninc=50,
                          T_final=290.0)
    res = solve(step, "EPICP", props, EPICP_NSTATEV, T_init=290.0)

    assert_ran_and_responded(res, n_expected=50)
    np.testing.assert_allclose(res["Strain"][0, -1], 0.01, rtol=1e-8)
    # guards the thermomechanical stress output fix (sigma, not the unset tau route)
    assert res["Stress"][0, -1] > 300.0


# ---------------------------------------------------------------------------
# modular UMAT / adapters / tangent modes
# ---------------------------------------------------------------------------

def test_modul_and_adapter_match():
    from simcoon.modular import (ModularMaterial, IsotropicElasticity, Plasticity,
                                 VonMisesYield, PowerLawHardening)

    mat = ModularMaterial(
        elasticity=IsotropicElasticity(C1=70000.0, C2=0.3, alpha=1.0e-5),
        mechanisms=[Plasticity(sigma_Y=300.0,
                               yield_criterion=VonMisesYield(),
                               isotropic_hardening=PowerLawHardening(k=1000.0, m=0.3))],
    )
    step = StepMeca(control=_UNIAXIAL, value=[0.02, 0, 0, 0, 0, 0], ninc=100)
    res_mod = solve(step, mat.umat_name, mat.props, mat.nstatev, T_init=290.0)
    # EPICP is a dedicated kernel; the modular engine integrates the same model
    # with a different local scheme -> close but not identical in the transient
    res_leg = solve(step, "EPICP", EPICP_PROPS, EPICP_NSTATEV, T_init=290.0)
    np.testing.assert_allclose(res_mod["Stress"][0], res_leg["Stress"][0], rtol=3e-2)
    # the adapter-served legacy name ELISO and a pure-elastic MODUL are the SAME
    # engine and must match exactly
    from simcoon.modular import ModularMaterial as MM, IsotropicElasticity as IE
    mat_el = MM(elasticity=IE(C1=70000.0, C2=0.3, alpha=1.0e-5))
    step_el = StepMeca(control=_UNIAXIAL, value=[0.01, 0, 0, 0, 0, 0], ninc=20)
    res_el_mod = solve(step_el, mat_el.umat_name, mat_el.props, mat_el.nstatev, T_init=290.0)
    res_el_leg = solve(step_el, "ELISO", ELISO_PROPS, 1, T_init=290.0)
    np.testing.assert_allclose(res_el_mod["Stress"], res_el_leg["Stress"], atol=1e-9)


def test_solver_is_safe_as_the_first_call_of_a_process():
    """The engine runs with the GIL released, and Armadillo allocates through numpy's
    allocator (in _core everywhere, in libsimcoon too on Windows). numpy's C-API table
    used to be imported lazily by the first allocation of each translation unit, a
    Python call made without the GIL: an access violation on Windows whenever the solver
    was the first binding of the process to touch a law (Sep 2026, feature/micro CI). It
    is now imported once at `import simcoon`. A fresh interpreter makes the solver the
    first caller, on every platform, for the modular engine (MODUL, and ELISO through
    its adapter) and a dedicated kernel (EPICP)."""
    cases = [("MODUL", MODUL_ELASTIC_PROPS, 1), ("ELISO", ELISO_PROPS, 1),
             ("EPICP", EPICP_PROPS, EPICP_NSTATEV)]
    code = textwrap.dedent(f"""
        import numpy as np
        from simcoon.solver import StepMeca, solve
        step = StepMeca(control=["strain"] + ["stress"] * 5,
                        value=[0.002, 0, 0, 0, 0, 0], ninc=2)
        for name, props, nstatev in {cases!r}:
            res = solve(step, name, np.asarray(props), nstatev, T_init=290.0)
            assert res.status == 0, name
        print("ok")
    """)
    proc = subprocess.run([sys.executable, "-X", "faulthandler", "-c", code],
                          capture_output=True, text=True, timeout=600)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "ok" in proc.stdout


@pytest.mark.parametrize("name, props", [("MODUL", MODUL_ELASTIC_PROPS), ("ELISO", ELISO_PROPS)])
def test_modular_rejects_too_few_statev(name, props):
    """ModularUMAT::initialize counts the state variables its mechanisms register and
    refuses a shorter statev, on the direct route and through the legacy adapter."""
    step = StepMeca(control=_UNIAXIAL, value=[0.002, 0, 0, 0, 0, 0], ninc=1)
    with pytest.raises(Exception, match="nstatev"):
        solve(step, name, np.asarray(props), 0, T_init=290.0)


def test_orientation_is_given_in_degrees():
    """solve() takes the RVE Euler angles in degrees, as material.dat and the docs do,
    and converts them for the C++ side. A transversely isotropic material is invariant
    under a 180 degree rotation about the third axis and not under 90 degrees; 180 rad
    would be neither."""
    elist = [1.0, 4500.0, 2300.0, 0.05, 0.3, 2700.0, 0.0, 0.0]   # axis E_L E_T nu_TL nu_TT G_LT alphas (examples/mechanical/ELIST.py)
    step = StepMeca(control=_UNIAXIAL, value=[0.002, 0, 0, 0, 0, 0], ninc=2)
    run = lambda psi: solve(step, "ELIST", elist, 1, T_init=290.0,
                            orientation=(psi, 0.0, 0.0))["Stress"][0, -1]
    assert run(180.0) == pytest.approx(run(0.0), rel=1e-9)
    assert abs(run(90.0) - run(0.0)) > 0.1 * abs(run(0.0))


def test_record_tangent_false_omits_tangent_history():
    step = StepMeca(control=_UNIAXIAL, value=[0.002, 0, 0, 0, 0, 0], ninc=2)
    res = solve(step, "ELISO", ELISO_PROPS, 1, T_init=290.0, record_tangent=False)
    assert res.status == 0
    assert "TangentMatrix" not in res
    assert "TangentMatrix" in solve(step, "ELISO", ELISO_PROPS, 1, T_init=290.0)


@pytest.mark.parametrize("mode", ["none", "continuum", "algorithmic"])
def test_tangent_modes_converge_to_same_state(mode):
    step = StepMeca(control="stress", value=[400.0, 0, 0, 0, 0, 0], ninc=50)
    res = solve(step, "EPICP", EPICP_PROPS, EPICP_NSTATEV, T_init=290.0,
                tangent_mode=mode)
    assert res.status == 0
    np.testing.assert_allclose(res["Stress"][0, -1], 400.0, atol=1e-4)


# ---------------------------------------------------------------------------
# error paths
# ---------------------------------------------------------------------------

def test_thermomechanical_rejects_finite_strain_control():
    step = StepThermomeca(control="stress", value=[0.0] * 6, T_final=300.0)
    with pytest.raises(ValueError):
        solve(Block(steps=[step], control_type="logarithmic"), "ELISO",
              ELISO_T_PROPS, 1)


def test_bad_control_flag():
    with pytest.raises(ValueError):
        StepMeca(control=["bad_flag"] * 6, value=[0.0] * 6).to_dict(1, 290.0)


def test_tabular_requires_table():
    step = StepMeca(control="strain", mode="tabular")
    with pytest.raises(ValueError):
        solve(step, "ELISO", ELISO_PROPS, 1)


def test_stress_control_rejected_for_F_control():
    step = StepMeca(control=["stress"] * 9, value=np.eye(3).ravel())
    with pytest.raises(ValueError):
        solve(Block(steps=[step], control_type="F"), "ELISO", ELISO_PROPS, 1)


def test_mixed_step_types_rejected():
    b = Block(steps=[StepMeca(value=[0.0] * 6),
                     StepThermomeca(value=[0.0] * 6)])
    with pytest.raises(ValueError):
        b.to_dict(290.0)


def test_mixed_block_types_rejected():
    # a mechanical block followed by a thermomechanical one must fail up front
    # (the engine constructs the state variables once, for the first block)
    b1 = Block(steps=[StepMeca(control=_UNIAXIAL, value=[0.01, 0, 0, 0, 0, 0], ninc=5)])
    b2 = Block(steps=[StepThermomeca(control="stress", value=[0.0] * 6,
                                     T_final=340.0, ninc=5)])
    with pytest.raises(ValueError):
        solve([b1, b2], "ELISO", ELISO_T_PROPS, 1, T_init=290.0)


def test_tabular_heat_flux():
    # thermomechanical tabular step with a prescribed heat-flux column (Q = 0:
    # adiabatic); strain driven so the mechanical response is deterministic
    t = np.linspace(0.05, 1.0, 20)
    e11 = 0.002 * t
    table = np.column_stack([t, np.zeros_like(t), e11])
    step = StepThermomeca(control=["strain"] + ["zero"] * 5, mode="tabular",
                          tabular=table, thermal_control="heat_flux")
    res = solve(step, "ELISO", ELISO_T_PROPS, 1, T_init=290.0)
    assert res.status == 0
    np.testing.assert_allclose(res["Strain"][0], e11, atol=1e-9)
    assert "Q" in res


def test_tabular_cycling_rejected():
    # the time column of a table is absolute: cycling a tabular step is ill-defined
    table = np.column_stack([np.linspace(0.1, 1.0, 10), np.linspace(0, 0.01, 10)])
    step = StepMeca(control=["strain"] + ["zero"] * 5, mode="tabular", tabular=table)
    with pytest.raises(ValueError):
        solve(Block(steps=[step], ncycle=2), "ELISO", ELISO_PROPS, 1)


def test_tabular_restarting_time_rejected():
    # a second-step table whose ABSOLUTE time column restarts at 0 would feed
    # negative DTime to the UMAT
    s1 = StepMeca(control=_UNIAXIAL, value=[0.005, 0, 0, 0, 0, 0], ninc=10, time=1.0)
    table = np.column_stack([np.linspace(0.1, 1.0, 10), np.linspace(0.005, 0.01, 10)])
    s2 = StepMeca(control=["strain"] + ["zero"] * 5, mode="tabular", tabular=table)
    with pytest.raises(RuntimeError):
        solve(Block(steps=[s1, s2]), "ELISO", ELISO_PROPS, 1)


# ---------------------------------------------------------------------------
# thermomechanical tabular + parameter spellings
# ---------------------------------------------------------------------------

def test_thermomechanical_tabular_constant_T():
    # default (constant temperature) thermomechanical tabular step: guards the
    # cBC_T >= 2 column-count fix in step_thermomeca::generate
    t = np.linspace(0.05, 1.0, 20)
    e11 = 0.005 * t
    step = StepThermomeca(control=["strain"] + ["zero"] * 5, mode="tabular",
                          tabular=np.column_stack([t, e11]))
    res = solve(step, "ELISO", ELISO_T_PROPS, 1, T_init=290.0)
    assert res.status == 0
    np.testing.assert_allclose(res["Strain"][0], e11, atol=1e-10)
    np.testing.assert_allclose(res["Temp"], 290.0, atol=1e-10)


def test_tabular_T_chains_to_next_step():
    # a tabular temperature ramp must chain its FINAL temperature into the next
    # step's hold value (no spurious ramp back to T_init)
    t = np.linspace(0.05, 1.0, 20)
    T_col = np.linspace(290.0 + 3.0, 350.0, 20)
    s1 = StepMeca(control=["zero"] * 6, mode="tabular", tabular_T=True,
                  tabular=np.column_stack([t, T_col]))
    s2 = StepMeca(control="stress", value=[0.0] * 6, ninc=10, time=1.0)  # T_final=None
    res = solve(Block(steps=[s1, s2]), "ELISO", ELISO_PROPS, 1, T_init=290.0)
    np.testing.assert_allclose(res["Temp"][-10:], 350.0, atol=1e-9)


def test_lambda_solver_param():
    step = StepMeca(control=_UNIAXIAL, value=[0.01, 0, 0, 0, 0, 0], ninc=10)
    res = solve(step, "ELISO", ELISO_PROPS, 1, T_init=290.0, lambda_solver=5.0e4)
    assert res.status == 0
    np.testing.assert_allclose(res["Stress"][0, -1], 700.0, rtol=1e-8)


