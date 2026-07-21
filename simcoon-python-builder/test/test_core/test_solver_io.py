"""Tests for the solver JSON configuration and results persistence."""

import numpy as np
import pytest

from simcoon.solver import (Block, SolverResults, StepMeca, StepThermomeca,
                            load_material_json, load_path_json,
                            load_simulation_json, save_material_json,
                            save_path_json, solve)

ELISO_PROPS = [70000.0, 0.3, 1.0e-5]


def _uniaxial_block(ninc=20):
    step = StepMeca(control=["strain"] + ["stress"] * 5,
                    value=[0.01, 0, 0, 0, 0, 0], ninc=ninc)
    return Block(steps=[step])


def test_material_roundtrip(tmp_path):
    f = tmp_path / "material.json"
    save_material_json(str(f), "ELISO", ELISO_PROPS, 1, orientation=(0.1, 0.2, 0.3))
    kw = load_material_json(str(f))
    assert kw["umat_name"] == "ELISO"
    np.testing.assert_allclose(kw["props"], ELISO_PROPS)
    assert kw["nstatev"] == 1
    assert kw["orientation"] == (0.1, 0.2, 0.3)


def test_path_roundtrip_mechanical(tmp_path):
    f = tmp_path / "path.json"
    b = _uniaxial_block()
    save_path_json(str(f), [b], T_init=300.0, corate="green_naghdi")
    blocks, T_init, corate = load_path_json(str(f))
    assert T_init == 300.0
    assert corate == "green_naghdi"
    assert len(blocks) == 1 and len(blocks[0].steps) == 1
    s = blocks[0].steps[0]
    assert isinstance(s, StepMeca) and not s._thermomechanical
    np.testing.assert_allclose(s.value, [0.01, 0, 0, 0, 0, 0])
    assert s.ninc == 20


def test_path_roundtrip_thermomechanical(tmp_path):
    f = tmp_path / "path.json"
    step = StepThermomeca(control="stress", value=[0.0] * 6, ninc=10,
                          thermal_control="heat_flux", Q=2.5)
    save_path_json(str(f), [Block(steps=[step])], T_init=290.0)
    blocks, _, _ = load_path_json(str(f))
    s = blocks[0].steps[0]
    assert isinstance(s, StepThermomeca)
    assert s.thermal_control == "heat_flux"
    assert s.Q == 2.5


def test_path_roundtrip_tabular(tmp_path):
    f = tmp_path / "path.json"
    table = np.column_stack([np.linspace(0.1, 1.0, 10), np.linspace(0, 0.01, 10)])
    step = StepMeca(control=["strain"] + ["zero"] * 5, mode="tabular", tabular=table)
    save_path_json(str(f), [Block(steps=[step])], T_init=290.0)
    blocks, _, _ = load_path_json(str(f))
    np.testing.assert_allclose(blocks[0].steps[0].tabular, table)


def test_simulation_json_solve_equivalence(tmp_path):
    mf, pf = tmp_path / "m.json", tmp_path / "p.json"
    save_material_json(str(mf), "ELISO", ELISO_PROPS, 1)
    save_path_json(str(pf), [_uniaxial_block()], T_init=290.0)

    res_direct = solve(_uniaxial_block(), "ELISO", ELISO_PROPS, 1, T_init=290.0)
    res_json = solve(**load_simulation_json(str(mf), str(pf)))
    np.testing.assert_allclose(res_json["Stress"], res_direct["Stress"], atol=1e-12)
    np.testing.assert_allclose(res_json["Strain"], res_direct["Strain"], atol=1e-12)


def test_results_npz_roundtrip(tmp_path):
    res = solve(_uniaxial_block(), "ELISO", ELISO_PROPS, 1, T_init=290.0)
    f = tmp_path / "results.npz"
    res.save(str(f))
    res2 = SolverResults.load(str(f))
    assert res2.status == res.status
    assert sorted(res2.keys()) == sorted(res.keys())
    for k in res.keys():
        np.testing.assert_allclose(res2[k], res[k])


def test_results_dataframe():
    pd = pytest.importorskip("pandas")
    res = solve(_uniaxial_block(ninc=5), "ELISO", ELISO_PROPS, 1, T_init=290.0)
    df = res.to_dataframe()
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 5
    assert "Stress_11" in df.columns and "Time" in df.columns


def test_results_unknown_key():
    res = solve(_uniaxial_block(ninc=5), "ELISO", ELISO_PROPS, 1, T_init=290.0)
    with pytest.raises(KeyError):
        res["NotAField"]
