"""Hyperelastic elasticity blocks of the modular UMAT.

The contract: MODUL composed with a potential and no mechanism IS the
standalone kernel of that potential. The bridge is b_el = exp(2 eps_el), so
the reference kernel is driven with F = exp(eps_el).
"""

import numpy as np
import pytest
from scipy.linalg import expm

import simcoon as sim
from simcoon.modular import (
    GentThomasElasticity,
    IsiharaElasticity,
    ModularMaterial,
    MooneyRivlinElasticity,
    NeoHookeanElasticity,
    SwansonElasticity,
    YeohElasticity,
)

# (modular block, standalone UMAT name, that UMAT's props)
MODELS = [
    (NeoHookeanElasticity(mu=0.5673, kappa=1000.0), "NEOHC", [0.5673, 1000.0]),
    (MooneyRivlinElasticity(C10=0.2588, C01=-0.0449, kappa=10000.0), "MOORI",
     [0.2588, -0.0449, 10000.0]),
    (YeohElasticity(C10=0.30, C20=-0.010, C30=0.0005, kappa=1000.0), "YEOHH",
     [0.30, -0.010, 0.0005, 1000.0]),
    (IsiharaElasticity(C10=0.1161, C20=0.0136, C01=0.0114, kappa=4000.0), "ISHAH",
     [0.1161, 0.0136, 0.0114, 4000.0]),
    (GentThomasElasticity(c1=0.2837, c2=2.81e-11, kappa=4000.0), "GETHH",
     [0.2837, 2.81e-11, 4000.0]),
    (SwansonElasticity(terms=((0.5, 0.1, 0.9, 0.6), (0.2, 0.05, 1.1, 0.8)),
                       kappa=4000.0), "SWANH",
     [2.0, 4000.0, 0.5, 0.1, 0.9, 0.6, 0.2, 0.05, 1.1, 0.8]),
]

STATES = {
    "generic": [0.28, -0.11, -0.09, 0.06, -0.03, 0.04],
    "two_equal": [0.15, -0.07, -0.07, 0.0, 0.0, 0.0],
    "dilatation": [0.10, 0.10, 0.10, 0.0, 0.0, 0.0],
    "ground": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
}


def _v2m(e):
    return np.array([[e[0], e[3] / 2.0, e[4] / 2.0],
                     [e[3] / 2.0, e[1], e[5] / 2.0],
                     [e[4] / 2.0, e[5] / 2.0, e[2]]])


def _umat(name, props, etot, F1, nstatev):
    n = 1
    z6 = lambda: np.zeros((6, n), order="F")
    eye = np.eye(3).reshape(3, 3, 1).copy(order="F")
    stress, sv, wm, Lt = sim.umat(
        name,
        np.asfortranarray(np.asarray(etot, dtype=float).reshape(6, 1)),
        z6(),
        eye if F1 is not None else np.array([]),
        np.asarray(F1).reshape(3, 3, 1).copy(order="F") if F1 is not None else np.array([]),
        z6(), eye,
        np.asfortranarray(np.asarray(props, dtype=float).reshape(-1, 1)),
        np.zeros((nstatev, n), order="F"), 0.0, 1.0,
        np.zeros((4, n), order="F"), n_threads=1)
    return stress[:, 0], Lt[:, :, 0]


@pytest.mark.parametrize("state", list(STATES))
@pytest.mark.parametrize("block,umat_name,umat_props",
                         MODELS, ids=[m[1] for m in MODELS])
def test_modular_hyper_matches_standalone_kernel(block, umat_name, umat_props, state):
    eps = np.array(STATES[state])
    mat = ModularMaterial(elasticity=block)

    sigma_mod, Lt_mod = _umat("MODUL", mat.props, eps, None, mat.nstatev)

    # the standalone kernels output Cauchy; MODUL is a kirchhoff_box model
    sigma_ref, Lt_ref = _umat(umat_name, umat_props, np.zeros(6),
                              expm(_v2m(eps)), 1)
    tau_ref = np.exp(eps[:3].sum()) * sigma_ref

    scale = max(1.0, np.abs(tau_ref).max())
    assert np.abs(sigma_mod - tau_ref).max() < 1e-10 * scale
    assert np.abs(Lt_mod - Lt_ref).max() < 1e-9 * np.abs(Lt_ref).max()


def test_yeoh_ground_state_is_L_iso():
    """L0 (the mechanisms' reference) is the closed-form ground state."""
    mat = ModularMaterial(
        elasticity=YeohElasticity(C10=0.30, C20=-0.010, C30=0.0005, kappa=1000.0))
    _, Lt = _umat("MODUL", mat.props, np.zeros(6), None, mat.nstatev)
    L_ref = np.asarray(sim.L_iso([1000.0, 2.0 * 0.30], "Kmu"))
    assert np.abs(Lt - L_ref).max() / np.abs(L_ref).max() < 1e-9


def test_props_roundtrip_layout():
    """[potential, n_params, params..., alpha] preceded by the elasticity type."""
    block = YeohElasticity(C10=1.0, C20=2.0, C30=3.0, kappa=4.0, alpha=5.0)
    assert block.to_props() == [2.0, 4.0, 1.0, 2.0, 3.0, 4.0, 5.0]
    assert block.nprops == 7
    mat = ModularMaterial(elasticity=block)
    # elasticity_type, then the block, then the mechanism count
    assert list(np.asarray(mat.props)) == [4.0] + block.to_props() + [0.0]


def test_swanson_rejects_malformed_terms():
    with pytest.raises(TypeError):
        SwansonElasticity(terms=((0.5, 0.1, 0.9),), kappa=1.0)
