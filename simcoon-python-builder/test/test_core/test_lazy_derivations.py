"""Derived data computed once and cached: the Eshelby memo of the homogenization schemes,
the rotation handles of a Basis, the lazy 'LogStrain' of the solver results.

None of these changes a number: the memo returns the tensor a fresh integration would
(bitwise key), the Basis cache holds the same rotation, LogStrain is the same formula on
the same F -- they only avoid recomputing.
"""

import pickle

import numpy as np
import pytest

import simcoon as sim
from simcoon.solver import Block, StepMeca, solve
from simcoon.solver.micromechanics import Ellipsoid

UNI = ["strain"] + ["stress"] * 5
MATRIX = dict(umat_name="ELISO", save=1, nstatev=1, props=np.array([2250.0, 0.19, 8.8e-5]))
FIBRE = dict(umat_name="ELISO", save=1, nstatev=1, props=np.array([73000.0, 0.19, 0.5e-6]))


def _phases(n_inclusions, fraction=0.3, oriented=True):
    phases = [Ellipsoid(number=0, concentration=1 - fraction, **MATRIX)]
    for i in range(1, n_inclusions + 1):
        rot = sim.Rotation.from_euler("zxz", [10.0 * i, 5.0 * i, 0.0], degrees=True) if oriented else sim.Rotation.identity()
        phases.append(Ellipsoid(number=i, concentration=fraction / n_inclusions, material_orientation=rot, **FIBRE))
    return phases


@pytest.mark.parametrize("scheme", ["MIMTN", "MISCN"])
def test_identical_inclusions_split_or_merged_give_the_same_effective_stiffness(scheme):
    # the memo shares one Eshelby integration between identical spheres in the same medium;
    # splitting a phase into identical sub-phases must not change the result
    one = sim.L_eff(scheme, [50, 50, 0], 0, orientation=(0.0, 0.0, 0.0), phases=_phases(1, oriented=False))
    ten = sim.L_eff(scheme, [50, 50, 0], 0, orientation=(0.0, 0.0, 0.0), phases=_phases(10, oriented=False))
    np.testing.assert_allclose(ten, one, rtol=1e-12, atol=1e-8)


def test_oriented_inclusions_share_the_integration_only_in_an_isotropic_medium():
    # Mori-Tanaka: the matrix is isotropic, so every orientation sees the same medium in its
    # own frame and the memo serves all of them; the result must equal the per-phase one.
    phases = _phases(10)
    L = sim.L_eff("MIMTN", [50, 50, 0], 0, orientation=(0.0, 0.0, 0.0), phases=phases)
    # an isotropic sphere is orientation-blind: same stiffness as the unoriented ten
    L_ref = sim.L_eff("MIMTN", [50, 50, 0], 0, orientation=(0.0, 0.0, 0.0), phases=_phases(10, oriented=False))
    np.testing.assert_allclose(L, L_ref, rtol=1e-10, atol=1e-6)
    # distinct geometries must not share: a prolate inclusion changes the answer
    phases[1] = Ellipsoid(number=1, concentration=0.03, a1=3.0, **FIBRE)
    L_prolate = sim.L_eff("MIMTN", [50, 50, 0], 0, orientation=(0.0, 0.0, 0.0), phases=phases)
    assert not np.allclose(L_prolate, L_ref, rtol=1e-6)


def test_solver_with_an_elastic_matrix_reuses_the_eshelby_tensor():
    # the per-inclusion memo makes a 200-increment Mori-Tanaka run cost a few Eshelby
    # integrations instead of thousands; the history must be the one of a fresh run
    phases = _phases(3)
    blk = Block(steps=[StepMeca(control=UNI, value=[0.01, 0, 0, 0, 0, 0], ninc=20)])
    r1 = solve(blk, "MIMTN", [50, 50, 0], 0, T_init=290.0, phases=phases)
    r2 = solve(blk, "MIMTN", [50, 50, 0], 0, T_init=290.0, phases=phases)
    np.testing.assert_array_equal(r1["Stress"], r2["Stress"])
    assert r1.status == 0 and np.all(np.isfinite(r1["Stress"]))


def test_basis_rotation_handles_are_cached_and_droppable():
    R = sim.Rotation.from_euler("zxz", [25.0, 40.0, -15.0], degrees=True)
    b = sim.Basis(rotation=R, name="material")
    s = sim.Tensor2.stress(np.array([100.0, -40.0, 25.0, 30.0, -12.0, 8.0]))
    np.testing.assert_allclose(s.to_basis(b).voigt, s.rotate(R, active=False).voigt, rtol=1e-13, atol=1e-10)
    assert b._inverse_rotation() is b._inverse_rotation()               # computed once
    assert b._rotation_cpp(True) is b._rotation_cpp(True)
    np.testing.assert_allclose(s.to_basis(b).to_basis(None).voigt, s.voigt, rtol=1e-13, atol=1e-10)
    clone = pickle.loads(pickle.dumps(s.to_basis(b)))                   # cache dropped, basis kept
    assert clone.basis.equals(b) and clone.basis.name == "material"
    np.testing.assert_allclose(clone.to_basis(None).voigt, s.voigt, rtol=1e-13, atol=1e-10)


def test_log_strain_is_derived_on_first_access():
    blk = Block(steps=[StepMeca(control=UNI, value=[0.05, 0, 0, 0, 0, 0], ninc=10)], control_type=3)
    res = solve(blk, "ELISO", [70000.0, 0.3, 0.0], 1, T_init=290.0)
    assert "LogStrain" in res._pending and "LogStrain" in res and "LogStrain" in res.keys()
    log = res["LogStrain"]
    assert "LogStrain" not in res._pending and res["LogStrain"] is log      # computed once
    F = res["F"]
    V2 = np.einsum("ijn,kjn->ikn", F, F)                                      # b = F F^T, ln V = 1/2 ln b
    for n in range(len(res)):
        w, Q = np.linalg.eigh(V2[:, :, n])
        lnV = Q @ np.diag(0.5 * np.log(w)) @ Q.T
        np.testing.assert_allclose(sim.Tensor2.strain(log[:, n]).mat, lnV, rtol=1e-10, atol=1e-12)
    # save/load materialise it
    import tempfile, os
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "r.npz")
        res.save(path)
        back = sim.solver.results.SolverResults.load(path)
        np.testing.assert_array_equal(back["LogStrain"], log)


def test_to_dataframe_materialises_the_pending_log_strain():
    pytest.importorskip("pandas")
    blk = Block(steps=[StepMeca(control=UNI, value=[0.05, 0, 0, 0, 0, 0], ninc=10)], control_type=3)
    res = solve(blk, "ELISO", [70000.0, 0.3, 0.0], 1, T_init=290.0)
    assert "LogStrain" in res._pending
    df = res.to_dataframe()
    assert "LogStrain_11" in df.columns
    np.testing.assert_array_equal(df["LogStrain_11"].to_numpy(), res["LogStrain"][0])
