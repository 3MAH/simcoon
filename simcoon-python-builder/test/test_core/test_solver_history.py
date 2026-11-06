"""The recorded history of _core.solver_run, field by field.

The binding records the solver state into its own sink (solver_run.cpp, history_sink) rather
than through simcoon::solver_memory_sink, the sink of the C++ API. Both capture the same
converged state at every output point; these tests pin that contract on the arrays the binding
returns: one record per increment, the Cauchy stress as tau / det(F1), matrices flattened row
by row, the fields of each state-variable type, and the (N, n) Fortran layout SolverResults
takes views of.
"""

import numpy as np
import pytest

import simcoon as sim
from simcoon import _core
from simcoon.solver import Block, StepMeca, StepThermomeca, SolverResults

MIXED = ["strain"] + ["stress"] * 5
EPICP_PROPS = [70000.0, 0.3, 1.0e-5, 300.0, 1000.0, 0.3]
EPICP_T_PROPS = [7800.0, 500.0, 70000.0, 0.3, 1.0e-5, 300.0, 1000.0, 0.3]
VECTORS_M = {"Etot": 6, "etot": 6, "PKII": 6, "tau": 6, "sigma": 6, "Wm": 4}
MATRICES = {"F1": 9, "R": 9, "DR": 9}


def run(blocks, umat, props, nstatev, record_tangent=True, corate=3, **params):
    blocks_py = [b.to_dict(290.0) for b in blocks]
    return _core.solver_run(blocks_py, 290.0, umat, np.asarray(props, dtype=float), nstatev, None, 0,
                            corate, dict(tangent_mode=2, **params), record_tangent, None)


def check_common(raw, n, nstatev):
    for key in ("block", "cycle", "step", "inc", "time", "T"):
        assert raw[key].shape == (n,)
    for key, ncomp in {**VECTORS_M, **MATRICES}.items():
        a = raw[key]
        assert a.shape == (n, ncomp), key
        assert a.flags.f_contiguous, key
    assert raw["statev"].shape == (n, nstatev)
    # the Cauchy stress of every record is tau / det(F1)
    F = raw["F1"].reshape(n, 3, 3)
    np.testing.assert_allclose(raw["sigma"], raw["tau"] / np.linalg.det(F)[:, None], rtol=1e-14, atol=0.0)


CASES = {
    "small_strain": (lambda: [Block(steps=[StepMeca(control=MIXED, value=[0.01, 0, 0, 0, 0, 0], ninc=12),
                                           StepMeca(control=MIXED, value=[-0.005, 0, 0, 0, 0, 0], ninc=8)],
                                    ncycle=2)], 40),
    "logarithmic": (lambda: [Block(steps=[StepMeca(control=MIXED, value=[0.2, 0, 0, 0, 0, 0], ninc=15)],
                                   control_type="logarithmic")], 15),
}


@pytest.mark.parametrize("case", sorted(CASES))
@pytest.mark.parametrize("record_tangent", [True, False])
def test_mechanical_history(case, record_tangent):
    make, n = CASES[case]
    raw = run(make(), "EPICP", EPICP_PROPS, 8, record_tangent)
    assert raw["status"] == 0 and raw["sv_type"] == 1
    check_common(raw, n, 8)
    assert "dSdE" not in raw and "Q" not in raw
    if record_tangent:
        assert raw["Lt"].shape == (n, 36) and raw["Lt"].flags.f_contiguous
    else:
        assert "Lt" not in raw
    res = SolverResults(raw, [case != "small_strain"], 3)
    assert res["Stress"].shape == (6, n) and res["Stress"].flags.c_contiguous
    assert np.shares_memory(res["Stress"], raw["sigma"])


def test_matrices_are_flattened_row_by_row():
    # a non-symmetric deformation gradient prescribed in full: the last F1 is the target,
    # written row by row in BC_meca and in the record alike
    target = np.array([[1.2, 0.1, 0.0], [0.05, 0.9, 0.02], [0.0, 0.03, 0.95]])
    raw = run([Block(steps=[StepMeca(control="F", value=target.ravel(), ninc=10)], control_type="F")],
              "NEOHC", [1000.0, 10000.0], 1)
    assert raw["status"] == 0
    check_common(raw, 10, 1)
    np.testing.assert_allclose(raw["F1"][-1].reshape(3, 3), target, rtol=0.0, atol=1e-12)
    res = SolverResults(raw, [True], 3)
    np.testing.assert_array_equal(res["F"][:, :, -1], raw["F1"][-1].reshape(3, 3))
    np.testing.assert_array_equal(res["TangentMatrix"][:, :, 3], raw["Lt"][3].reshape(6, 6))


@pytest.mark.parametrize("record_tangent", [True, False])
def test_thermomechanical_history(record_tangent):
    n = 10
    step = StepThermomeca(control=MIXED, value=[0.01, 0, 0, 0, 0, 0], ninc=n, T_final=320.0)
    raw = run([Block(steps=[step])], "EPICP", EPICP_T_PROPS, 8, record_tangent)
    assert raw["status"] == 0 and raw["sv_type"] == 2
    check_common(raw, n, 8)
    assert raw["Q"].shape == (n,) and raw["r"].shape == (n,) and raw["Wt"].shape == (n, 3)
    assert "Lt" not in raw
    if record_tangent:
        assert raw["dSdE"].shape == (n, 36)
        for key in ("dSdT", "drdE"):
            assert raw[key].shape == (n, 6)
        assert raw["drdT"].shape == (n, 1)
    else:
        assert "dSdE" not in raw


def test_aborted_run_keeps_its_records():
    # a stress beyond the saturation of a perfectly plastic law, without inforce: the run stops
    # early and every history holds the same records
    step = StepMeca(control=["stress"] * 6, value=[900.0, 0, 0, 0, 0, 0], ninc=50, Dn_mini=0.25)
    raw = run([Block(steps=[step])], "EPICP", [70000.0, 0.3, 1.0e-5, 300.0, 0.0, 1.0], 8, inforce=0)
    assert raw["status"] == 1
    n = raw["time"].shape[0]
    assert 0 < n < 50
    check_common(raw, n, 8)
    assert raw["Lt"].shape == (n, 36)
