"""The reference cases of testBin/Umats, driven from their JSON inputs against the
committed results of the historical C++ tests.

Each case used to be a gtest reading path.txt / material.dat / N<kind>.dat and comparing
the solver's result file with comparison/*.txt to 1e-6. The inputs are JSON now
(simcoon reads no text format), the engine is the same, and the references carry 6
significant digits, hence the relative tolerance.
"""

from pathlib import Path

import numpy as np
import pytest

from simcoon import solver as slv
from simcoon.solver.micromechanics import load_ellipsoids_json, load_layers_json

UMATS = Path(__file__).resolve().parents[3] / "testBin" / "Umats"

#: the reference result files: columns 8:14 strain, 14:20 stress (Voigt 11 22 33 12 13 23)
S_STRAIN, S_STRESS = slice(8, 14), slice(14, 20)


def _run(case, material="material.json", **extra):
    data = UMATS / case / "data"
    kwargs = slv.load_simulation_json(data / material, data / "path.json")
    kwargs.update(extra)
    return slv.solve(**kwargs)


def _reference(case, name="results_job_global-0.txt"):
    return np.loadtxt(UMATS / case / "comparison" / name)


def _assert_history(res, ref, rows=None):
    strain, stress = res["Strain"].T, res["Stress"].T
    if rows is not None:
        strain, stress = strain[rows], stress[rows]
    assert strain.shape[0] == ref.shape[0], (strain.shape, ref.shape)
    np.testing.assert_allclose(strain, ref[:, S_STRAIN], rtol=1e-5, atol=1e-9)
    np.testing.assert_allclose(stress, ref[:, S_STRESS], rtol=1e-5, atol=1e-6)


@pytest.mark.parametrize("case, extra", [
    ("ELISO", {}), ("ELIST", {}), ("ELORT", {}),
    ("EPICP", {"precision": 1e-5}), ("EPKCP", {"precision": 1e-5}),
])
def test_homogeneous_case(case, extra):
    _assert_history(_run(case, **extra), _reference(case))


@pytest.mark.parametrize("case", ["EPCHA", "EPCHG", "EPHAC"])
def test_cyclic_case_replays_its_table(case):
    """Three blocks (pre-cycle, alignment, tabular replay of the experiment); the
    reference holds the last rows of the history (the third block only for EPCHA and
    EPHAC, the whole run for EPCHG)."""
    res = _run(case)
    ref = _reference(case, "simul_1.txt")
    _assert_history(res, ref, rows=slice(-ref.shape[0], None))


def test_mori_tanaka_nested():
    """A two-level composite: the inner MIMTN phase carries its own ellipsoids."""
    data = UMATS / "MIMTN" / "data"
    phases = load_ellipsoids_json(data / "ellipsoids0.json")
    phases[0].phases = load_ellipsoids_json(data / "ellipsoids1.json")
    _assert_history(_run("MIMTN", phases=phases), _reference("MIMTN"))


def test_mori_tanaka_plastic_matrix():
    """MIMTN with a modular elasto-plastic matrix: the tangent family matters here."""
    data = UMATS / "MIMTP" / "data"
    _assert_history(_run("MIMTP", phases=load_ellipsoids_json(data / "ellipsoids0.json")),
                    _reference("MIMTP"))


def test_laminate():
    data = UMATS / "MIPLN" / "data"
    _assert_history(_run("MIPLN", phases=load_layers_json(data / "layers0.json")),
                    _reference("MIPLN"))


@pytest.mark.parametrize("tag", ["NH", "MR", "IS", "GT", "OG"])
def test_hyperelastic_models(tag):
    """Uniaxial stress to 3 MPa under finite strain: the reference carries the three
    isochoric principal stretches (columns 8:11, ascending) and the Cauchy stress (11:17)."""
    res = _run("HYPER", material=f"material_{tag}.json", precision=1e-5)
    ref = _reference("HYPER", f"results_{tag}.dat")
    F = np.moveaxis(res["F"], -1, 0)
    stretches = np.sqrt(np.linalg.eigvalsh(F @ F.transpose(0, 2, 1)))
    stretches /= np.cbrt(np.linalg.det(F))[:, None]
    assert stretches.shape[0] == ref.shape[0]
    np.testing.assert_allclose(stretches, ref[:, 8:11], rtol=1e-5, atol=1e-9)
    np.testing.assert_allclose(res["Stress"].T, ref[:, 11:17], rtol=1e-5, atol=1e-6)


def test_logarithmic_strain_under_F_control():
    """Control type 5 drives F to a prescribed tensor; the logarithmic strain of the last
    increment must be that of the target."""
    res = _run("LOG_int")
    F_target = np.array([[1.0, 0.5, 0.2], [0.0, 0.8, 0.0], [0.0, 0.0, 1.0]])
    F_end = np.ascontiguousarray(res["F"][:, :, -1])
    np.testing.assert_allclose(F_end, F_target, atol=1e-9)
    # the strain the solver integrated along the path (logarithmic rate) against ln V of the
    # target, computed here from F: the two are independent
    w, n = np.linalg.eigh(F_target @ F_target.T)
    lnV = 0.5 * (n * np.log(w)) @ n.T
    lnV_voigt = np.array([lnV[0, 0], lnV[1, 1], lnV[2, 2], 2 * lnV[0, 1], 2 * lnV[0, 2], 2 * lnV[1, 2]])
    np.testing.assert_allclose(res["Strain"][:, -1], lnV_voigt, atol=1e-3)
