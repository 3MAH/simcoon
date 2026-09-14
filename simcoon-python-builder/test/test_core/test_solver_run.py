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

import numpy as np
import pytest

import simcoon as sim
from simcoon.solver import Block, StepMeca, StepThermomeca, from_file, solve

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


# The legacy path-file grammar. Nothing in simcoon reads it any more — but
# simcoon.solver.from_file parses it in Python, and the tests of that parser need
# a file to parse, so the writers stay here as fixtures.

#: path-file component order (11, 12, 22, 13, 23, 33) -> Voigt index
_FILE_ORDER = [0, 3, 1, 4, 5, 2]


def _meca_state_lines(flags, values):
    """Emit the #prescribed_mechanical_state lines in path-file order."""
    toks = []
    for k in _FILE_ORDER:
        toks.append(f"{flags[k]} {values[k]}")
    return "{}\n{} {}\n{} {} {}".format(*toks)


def _step_text(flags, values, time=1.0, ninc=100, mode=1, T=290.0, BC_w=None):
    txt = f"""#Mode
{mode}
#Dn_init 1.
#Dn_mini 1.
#Dn_inc {1.0/ninc}
#time
{time}
#prescribed_mechanical_state
{_meca_state_lines(flags, values)}
"""
    if BC_w is not None:
        BC_w = np.asarray(BC_w)
        txt += "#Rotation\n" + "\n".join(
            " ".join(str(BC_w[i, j]) for j in range(3)) for i in range(3)
        ) + "\n"
    txt += f"#prescribed_temperature_state\nT {T}\n"
    return txt


def _thermo_step_text(flags, values, thermal_token, time=1.0, ninc=50):
    return f"""#Mode
1
#Dn_init 1.
#Dn_mini 1.
#Dn_inc {1.0/ninc}
#time
{time}
#prescribed_mechanical_state
{_meca_state_lines(flags, values)}
#prescribed_thermal_state
{thermal_token}
"""


def _path_text(steps_text, control_type=1, loading_type=1, ncycle=1, T_init=290.0):
    return f"""#Initial_temperature
{T_init}
#Number_of_blocks
1

#Block
1
#Loading_type
{loading_type}
#Control_type(NLGEOM)
{control_type}
#Repeat
{ncycle}
#Steps
{len(steps_text)}

""" + "\n".join(steps_text)


def write_path_file(tmp_path, path_text, extra_files=None):
    """Drop a legacy path.txt (and its tabular files) in tmp_path/data."""
    data = tmp_path / "data"
    data.mkdir(exist_ok=True)
    (data / "path.txt").write_text(path_text)
    for name, content in (extra_files or {}).items():
        (data / name).write_text(content)
    return str(data)


ELISO_PROPS = [70000.0, 0.3, 1.0e-5]
EPICP_PROPS = [70000.0, 0.3, 1.0e-5, 300.0, 1000.0, 0.3]  # E nu alpha sigmaY k m
EPICP_NSTATEV = 8
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
    # logarithmic or Biot — which is not the canonical strain the results carry, so the
    # invariant is put on the stretch itself: the bar elongates, monotonically.
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


# ---------------------------------------------------------------------------
# legacy file parsing (from_file / material_from_file)
# ---------------------------------------------------------------------------

def test_from_file_mechanical(tmp_path):
    # parse a legacy path.txt into Blocks and run the programme it describes
    flags = ["E"] + ["S"] * 5
    steps = [_step_text(flags, [0.02, 0, 0, 0, 0, 0], ninc=50),
             _step_text(flags, [0.0, 0, 0, 0, 0, 0], ninc=50)]
    data = write_path_file(tmp_path, _path_text(steps, ncycle=2))

    blocks, T_init = from_file(data, "path.txt")
    assert T_init == 290.0
    assert len(blocks) == 1 and blocks[0].ncycle == 2 and len(blocks[0].steps) == 2

    res = solve(blocks, "ELISO", ELISO_PROPS, 1, T_init=T_init)
    assert_ran_and_responded(res, n_expected=2 * 2 * 50)
    # the parsed programme loads to 0.02 and unloads to 0, twice
    np.testing.assert_allclose(res["Strain"][0].max(), 0.02, rtol=1e-8)
    np.testing.assert_allclose(res["Strain"][0, -1], 0.0, atol=1e-12)


def test_from_file_thermomechanical(tmp_path):
    flags = ["S"] * 6
    steps = [_thermo_step_text(flags, [0.0] * 6, "T 340")]
    data = write_path_file(tmp_path, _path_text(steps, loading_type=2))

    blocks, T_init = from_file(data, "path.txt")
    s = blocks[0].steps[0]
    assert isinstance(s, StepThermomeca) and s.thermal_control == "temperature"
    assert s.T_final == 340.0

    res = solve(blocks, "ELISO", ELISO_T_PROPS, 1, T_init=T_init)
    assert res.status == 0
    # the parsed ramp lands on its target, stress-free, with free thermal expansion
    np.testing.assert_allclose(res["Temp"][-1], 340.0, rtol=1e-10)
    np.testing.assert_allclose(res["Strain"][0, -1], 1.0e-5 * 50.0, rtol=1e-8)
    np.testing.assert_allclose(res["Stress"][:, -1], 0.0, atol=1e-6)


def test_from_file_tabular(tmp_path):
    t = np.linspace(0.02, 1.0, 50)
    e11 = 0.015 * np.sin(np.pi * t)
    tab_lines = "\n".join(f"{i+1} {t[i]:.16g} {e11[i]:.16g}" for i in range(len(t)))
    step3 = """#Mode
3
#File
tab.txt
#Dn_init 1.
#Dn_mini 1.
#prescribed_mechanical_state
E
0 0
0 0 0
#T_is_set
0
"""
    data = write_path_file(tmp_path, _path_text([step3]),
                           extra_files={"tab.txt": tab_lines})

    blocks, T_init = from_file(data, "path.txt")
    s = blocks[0].steps[0]
    assert s.mode == 3 and not s.tabular_T
    assert s.control[0] == "strain" and s.control[1:] == ["zero"] * 5
    np.testing.assert_allclose(s.tabular[:, 0], t, atol=1e-14)

    res = solve(blocks, "ELISO", ELISO_PROPS, 1, T_init=T_init)
    assert_ran_and_responded(res, n_expected=len(t))
    # the table read from disk drives the run, in its own absolute time
    np.testing.assert_allclose(res["Time"], t, atol=1e-12)
    np.testing.assert_allclose(res["Strain"][0], e11, atol=1e-10)


def test_material_from_file(tmp_path):
    from simcoon.solver import material_from_file
    (tmp_path / "material.dat").write_text("""Material
Name\tELISO
Number_of_material_parameters\t3
Number_of_internal_variables\t1

#Orientation
psi\t0.1
theta\t0.2
phi\t0.3

#Mechanical
E 70000.
nu 0.3
alpha 1.E-5
""")
    kw = material_from_file(str(tmp_path), "material.dat")
    assert kw["umat_name"] == "ELISO"
    assert kw["nstatev"] == 1
    np.testing.assert_allclose(kw["props"], [70000.0, 0.3, 1.0e-5])
    assert kw["orientation"] == (0.1, 0.2, 0.3)
    res = solve(StepMeca(control=_UNIAXIAL, value=[0.01, 0, 0, 0, 0, 0], ninc=10),
                T_init=290.0, **kw)
    assert res.status == 0


# The file-driven binding this module used to cross-check against (sim._core.solver)
# is gone, so from_file can no longer be compared run-for-run against the C++ reader.
# The grammar it reproduces is still the documented one and still implemented in C++
# (test/support/file_readers.cpp), so pin the parse itself: every field below is a
# place where a silent drift would produce a different loading programme that the
# "status == 0 and the stress is finite" assertions would not notice.
PATH_FIXTURE = """#Initial_temperature
300.0
#Number_of_blocks
2

#Block
1
#Loading_type
1
#Control_type(NLGEOM)
1
#Repeat
3
#Steps
2

#Mode
1
#Dn_init 1.
#Dn_mini 0.1
#Dn_inc 0.02
#time
7.5
#Consigne
E 0.11
S 0.12 E 0.22
S 0.13 S 0.23 E 0.33
#Consigne_T
T 305.

#Mode
2
#Dn_init 0.5
#Dn_mini 0.01
#Dn_inc 0.25
#time
2.
#Consigne
S 0. E 0. S 0. E 0. S 0. E 0.
#Consigne_T
T 310.

#Block
2
#Loading_type
2
#Control_type(NLGEOM)
1
#Repeat
1
#Steps
1

#Mode
1
#Dn_init 1.
#Dn_mini 0.05
#Dn_inc 0.1
#time
4.
#Consigne
E 0.01
S 0. S 0.
S 0. S 0. S 0.
#Consigne_T
Q 1500.
"""


def test_from_file_parses_the_documented_grammar(tmp_path):
    """Field-by-field pin of the legacy path-file grammar as from_file reads it."""
    (tmp_path / "path.txt").write_text(PATH_FIXTURE)
    blocks, T_init = sim.solver.from_file(str(tmp_path), "path.txt")

    assert T_init == 300.0
    assert len(blocks) == 2

    mech = blocks[0]
    assert mech.ncycle == 3 and mech.control_type == 1
    assert len(mech.steps) == 2

    first = mech.steps[0]
    assert first.mode in (1, "linear")
    assert first.time == 7.5
    assert first.ninc == 50                       # round(1 / Dn_inc)
    assert first.Dn_init == 1.0 and first.Dn_mini == 0.1
    assert first.T_final == 305.0
    # the file lists 11, 12, 22, 13, 23, 33 (lower triangle, row-wise); Voigt is
    # 11, 22, 33, 12, 13, 23. Swapping the two would scramble the shear components.
    assert list(first.control) == ["strain", "strain", "strain",
                                   "stress", "stress", "stress"]
    np.testing.assert_allclose(np.asarray(first.value, dtype=float),
                               [0.11, 0.22, 0.33, 0.12, 0.13, 0.23])

    second = mech.steps[1]
    assert second.mode in (2, "sinusoidal")
    assert second.ninc == 4 and second.time == 2.0
    # S E S E S E read as 11, 12, 22, 13, 23, 33 lands in Voigt as
    # 11=S, 22=S, 33=E, 12=E, 13=E, 23=S — the interleaving is the point
    assert list(second.control) == ["stress", "stress", "strain",
                                    "strain", "strain", "stress"]

    thermo = blocks[1]
    assert thermo.ncycle == 1
    heat = thermo.steps[0]
    assert heat.ninc == 10 and heat.time == 4.0
    # a 'Q' thermal condition is a heat flux, not a temperature target
    assert getattr(heat, "thermal_control", None) in ("heat_flux", 1)
    assert float(getattr(heat, "Q", 0.0)) == 1500.0
