"""Tests for the in-memory solver API (simcoon.solver / _core.solver_run).

The reference is the file-driven solver `sim._core.solver` (same C++ engine
through solver_file_sink): for every control type / corate / loading mode
exercised here, the in-memory run must reproduce the file run to numerical
identity.
"""

import numpy as np
import pytest

import simcoon as sim
from simcoon.solver import Block, StepMeca, StepThermomeca, from_file, solve

# ---------------------------------------------------------------------------
# helpers: file-driven reference runs
# ---------------------------------------------------------------------------

#: path-file component order (11, 12, 22, 13, 23, 33) -> Voigt index
_FILE_ORDER = [0, 3, 1, 4, 5, 2]

# default file output columns (no output.dat):
# 0 block | 1 cycle | 2 step | 3 inc | 4 time | 5 T | 6 Q | 7 r |
# 8:14 strain (Green-Lagrange) | 14:20 stress (Cauchy) | 20:24 Wm
_C_TIME, _C_T = 4, 5
_S_STRAIN, _S_STRESS, _S_WM = slice(8, 14), slice(14, 20), slice(20, 24)


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


def run_file_solver(tmp_path, path_text, umat, props, nstatev, corate=2,
                    solver_type=0, extra_files=None):
    data = tmp_path / "data"
    data.mkdir(exist_ok=True)
    results = tmp_path / "results"
    results.mkdir(exist_ok=True)
    (data / "path.txt").write_text(path_text)
    for name, content in (extra_files or {}).items():
        (data / name).write_text(content)
    sim._core.solver(umat, np.asarray(props, dtype=float), nstatev, 0.0, 0.0, 0.0,
               solver_type, corate, str(data), str(results), "path.txt", "res.txt")
    return np.loadtxt(results / "res_global-0.txt")


def _close_to_text(mem, txt):
    """Compare in-memory values against the same values parsed from the result
    text file, whose writer prints ~6 significant digits."""
    scale = max(np.abs(txt).max(), 1e-30)
    np.testing.assert_allclose(mem, txt, rtol=2e-5, atol=2e-6 * scale)


def assert_matches_file(res, out):
    """Compare an in-memory SolverResults against a file-solver output table."""
    assert len(res) == out.shape[0]
    _close_to_text(res["Time"], out[:, _C_TIME])
    _close_to_text(res["Temp"], out[:, _C_T])
    _close_to_text(res["Strain"].T, out[:, _S_STRAIN])
    _close_to_text(res["Stress"].T, out[:, _S_STRESS])
    _close_to_text(res["Wm"].T, out[:, _S_WM])


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


def test_eliso_load_unload_vs_file(tmp_path):
    flags = ["E"] + ["S"] * 5
    steps = [_step_text(flags, [0.02, 0, 0, 0, 0, 0], ninc=100),
             _step_text(flags, [0.0, 0, 0, 0, 0, 0], ninc=100)]
    out = run_file_solver(tmp_path, _path_text(steps), "ELISO", ELISO_PROPS, 1)

    s1 = StepMeca(control=_UNIAXIAL, value=[0.02, 0, 0, 0, 0, 0], ninc=100)
    s2 = StepMeca(control=_UNIAXIAL, value=[0.0, 0, 0, 0, 0, 0], ninc=100)
    res = solve(Block(steps=[s1, s2]), "ELISO", ELISO_PROPS, 1, T_init=290.0)
    assert_matches_file(res, out)


def test_epicp_mixed_cyclic_vs_file(tmp_path):
    # stress-controlled uniaxial cycling into the plastic range
    flags = ["S"] * 6
    steps = [_step_text(flags, [400.0, 0, 0, 0, 0, 0], ninc=100),
             _step_text(flags, [0.0, 0, 0, 0, 0, 0], ninc=100)]
    out = run_file_solver(tmp_path, _path_text(steps, ncycle=2), "EPICP",
                          EPICP_PROPS, EPICP_NSTATEV)

    s1 = StepMeca(control="stress", value=[400.0, 0, 0, 0, 0, 0], ninc=100)
    s2 = StepMeca(control="stress", value=[0.0, 0, 0, 0, 0, 0], ninc=100)
    res = solve(Block(steps=[s1, s2], ncycle=2), "EPICP", EPICP_PROPS,
                EPICP_NSTATEV, T_init=290.0)
    assert_matches_file(res, out)
    # cycle index bookkeeping
    assert res["Cycle"].max() == 1
    assert len(res) == 2 * 2 * 100
    # plasticity happened
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
def test_finite_strain_vs_file(tmp_path, umat, props, nstatev, control_type, corate):
    if umat == "EPICP" and control_type == 4 and corate == 5:
        # pre-existing engine limitation (identical for the file solver): plasticity
        # under strain-controlled Biot loading with the log_F rate fails in
        # logarithmic_F (singular DF inversion)
        pytest.skip("EPICP + Biot control + logarithmic_F: known engine limitation")
    if umat == "SNTVE" and control_type == 4:
        target = 0.05  # keep the Biot-controlled hyperelastic case well-conditioned
    else:
        target = 0.08
    flags = ["E"] + ["S"] * 5
    steps = [_step_text(flags, [target, 0, 0, 0, 0, 0], ninc=50,
                        BC_w=np.zeros((3, 3)))]
    out = run_file_solver(tmp_path, _path_text(steps, control_type=control_type),
                          umat, props, nstatev, corate=corate)

    step = StepMeca(control=_UNIAXIAL, value=[target, 0, 0, 0, 0, 0], ninc=50,
                    BC_w=np.zeros((3, 3)))
    ct = {2: "green_lagrange", 3: "logarithmic", 4: "biot"}[control_type]
    res = solve(Block(steps=[step], control_type=ct), umat, props, nstatev,
                T_init=290.0, corate=corate)
    assert_matches_file(res, out)


def test_ct3_spin_vs_file(tmp_path):
    # logarithmic control with a superimposed rotation rate (BC_w)
    BC_w = np.array([[0.0, 0.2, 0.0], [-0.2, 0.0, 0.0], [0.0, 0.0, 0.0]])
    flags = ["E"] + ["S"] * 5
    steps = [_step_text(flags, [0.05, 0, 0, 0, 0, 0], ninc=50, BC_w=BC_w)]
    out = run_file_solver(tmp_path, _path_text(steps, control_type=3), "ELISO",
                          ELISO_PROPS, 1, corate=2)

    step = StepMeca(control=_UNIAXIAL, value=[0.05, 0, 0, 0, 0, 0], ninc=50, BC_w=BC_w)
    res = solve(Block(steps=[step], control_type="logarithmic"), "ELISO",
                ELISO_PROPS, 1, T_init=290.0, corate="logarithmic")
    assert_matches_file(res, out)
    # the rotation history is captured
    R_end = res["R"][:, :, -1]
    assert np.abs(R_end - np.eye(3)).max() > 1e-3


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

def test_sinusoidal_vs_file(tmp_path):
    flags = ["E"] + ["S"] * 5
    steps = [_step_text(flags, [0.01, 0, 0, 0, 0, 0], ninc=100, mode=2)]
    out = run_file_solver(tmp_path, _path_text(steps), "ELISO", ELISO_PROPS, 1)

    step = StepMeca(control=_UNIAXIAL, value=[0.01, 0, 0, 0, 0, 0], ninc=100,
                    mode="sinusoidal")
    res = solve(step, "ELISO", ELISO_PROPS, 1, T_init=290.0)
    assert_matches_file(res, out)
    # sinusoidal profile: increments are not uniform
    de = np.diff(res["Strain"][0])
    assert de.max() / de.min() > 1.5


def test_tabular_memory_vs_file(tmp_path):
    # same table driven from disk (mode-3 file) and from memory (tab_data)
    t = np.linspace(0.02, 1.0, 50)
    e11 = 0.015 * np.sin(np.pi * t)
    table = np.column_stack([t, e11])
    tab_lines = "\n".join(
        f"{i+1} {t[i]:.16g} {e11[i]:.16g}" for i in range(len(t))
    )
    # mode-3 file step: flags only (no values), thermal token '0' = constant T
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
    out = run_file_solver(tmp_path, _path_text([step3]), "ELISO", ELISO_PROPS, 1,
                          extra_files={"tab.txt": tab_lines})

    step = StepMeca(control=["strain"] + ["zero"] * 5, mode="tabular", tabular=table)
    res = solve(step, "ELISO", ELISO_PROPS, 1, T_init=290.0)
    assert_matches_file(res, out)
    np.testing.assert_allclose(res["Time"], t, atol=1e-12)


# ---------------------------------------------------------------------------
# thermomechanical blocks
# ---------------------------------------------------------------------------

ELISO_T_PROPS = [1.0e-9, 1.0, 70000.0, 0.3, 1.0e-5]  # rho c_p E nu alpha


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


def test_thermomeca_temperature_ramp_vs_file(tmp_path):
    flags = ["S"] * 6
    steps = [_thermo_step_text(flags, [0.0] * 6, "T 340")]
    out = run_file_solver(tmp_path, _path_text(steps, loading_type=2), "ELISO",
                          ELISO_T_PROPS, 1)

    step = StepThermomeca(control="stress", value=[0.0] * 6, ninc=50, T_final=340.0)
    res = solve(step, "ELISO", ELISO_T_PROPS, 1, T_init=290.0)
    assert_matches_file(res, out)
    # free thermal expansion
    np.testing.assert_allclose(res["Strain"][0, -1], 1.0e-5 * 50.0, rtol=1e-8)
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


def test_thermomeca_epicp_vs_file(tmp_path):
    props = [1.0e-9, 1.0] + EPICP_PROPS  # rho c_p then mechanical props
    flags = ["E"] + ["S"] * 5
    steps = [_thermo_step_text(flags, [0.01, 0, 0, 0, 0, 0], "T 290")]
    out = run_file_solver(tmp_path, _path_text(steps, loading_type=2), "EPICP",
                          props, EPICP_NSTATEV)

    step = StepThermomeca(control=_UNIAXIAL, value=[0.01, 0, 0, 0, 0, 0], ninc=50,
                          T_final=290.0)
    res = solve(step, "EPICP", props, EPICP_NSTATEV, T_init=290.0)
    assert_matches_file(res, out)
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

def test_from_file_mechanical_vs_file(tmp_path):
    # parse a legacy path.txt into Blocks and solve in memory: must reproduce
    # the file-driven run of the very same file
    flags = ["E"] + ["S"] * 5
    steps = [_step_text(flags, [0.02, 0, 0, 0, 0, 0], ninc=50),
             _step_text(flags, [0.0, 0, 0, 0, 0, 0], ninc=50)]
    path_text = _path_text(steps, ncycle=2)
    out = run_file_solver(tmp_path, path_text, "ELISO", ELISO_PROPS, 1)

    blocks, T_init = from_file(str(tmp_path / "data"), "path.txt")
    assert T_init == 290.0
    assert len(blocks) == 1 and blocks[0].ncycle == 2 and len(blocks[0].steps) == 2
    res = solve(blocks, "ELISO", ELISO_PROPS, 1, T_init=T_init)
    assert_matches_file(res, out)


def test_from_file_thermomechanical_vs_file(tmp_path):
    flags = ["S"] * 6
    steps = [_thermo_step_text(flags, [0.0] * 6, "T 340")]
    path_text = _path_text(steps, loading_type=2)
    out = run_file_solver(tmp_path, path_text, "ELISO", ELISO_T_PROPS, 1)

    blocks, T_init = from_file(str(tmp_path / "data"), "path.txt")
    s = blocks[0].steps[0]
    assert isinstance(s, StepThermomeca) and s.thermal_control == "temperature"
    assert s.T_final == 340.0
    res = solve(blocks, "ELISO", ELISO_T_PROPS, 1, T_init=T_init)
    assert_matches_file(res, out)


def test_from_file_tabular_vs_file(tmp_path):
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
    path_text = _path_text([step3])
    out = run_file_solver(tmp_path, path_text, "ELISO", ELISO_PROPS, 1,
                          extra_files={"tab.txt": tab_lines})

    blocks, T_init = from_file(str(tmp_path / "data"), "path.txt")
    s = blocks[0].steps[0]
    assert s.mode == 3 and not s.tabular_T
    assert s.control[0] == "strain" and s.control[1:] == ["zero"] * 5
    np.testing.assert_allclose(s.tabular[:, 0], t, atol=1e-14)
    res = solve(blocks, "ELISO", ELISO_PROPS, 1, T_init=T_init)
    assert_matches_file(res, out)


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
