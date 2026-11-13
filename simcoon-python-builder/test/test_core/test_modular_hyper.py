"""Hyperelastic elasticity blocks of the modular UMAT.

The contract: MODUL composed with a potential and no mechanism IS the
standalone kernel of that potential. The bridge is b_el = exp(2 eps_el), so
the reference kernel is driven with F = exp(eps_el).
"""

import dataclasses

import numpy as np
import pytest
from scipy.linalg import expm

import simcoon as sim
from simcoon.modular import (
    GentThomasElasticity,
    HolzapfelElasticity,
    IsiharaElasticity,
    ModularMaterial,
    MooneyRivlinElasticity,
    NeoHookeanElasticity,
    SwansonElasticity,
    YeohElasticity,
)

# (modular block, UMAT name, that UMAT's props): NEOHC_closed and NEOHI_closed are the kernels built on
# F (an independent oracle); the other names run the same block through their props adapter, so the
# comparison checks the translation
MODELS = [
    (NeoHookeanElasticity(mu=0.5673, kappa=1000.0), "NEOHC_closed", [0.5673, 1000.0]),
    (NeoHookeanElasticity(mu=2000.0 / (2 * 1.40), kappa=2000.0 / (3 * 0.20), volumetric="quadratic"),
     "NEOHI_closed", [2000.0, 0.40, 0.0]),
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
    (HolzapfelElasticity(C10=0.0354, k1=0.0107, k2=7.48, kappa_d=0.1,
                         fibres=sim.Rotation.from_euler(
                             "zxz", [[0.0, 0.0, 40.0], [0.0, 0.0, -40.0]], degrees=True),
                         kappa=1000.0), "HOLZA",
     [0.0354, 0.0107, 7.48, 0.1, 2.0,
      np.cos(np.deg2rad(40.0)), np.sin(np.deg2rad(40.0)), 0.0,
      np.cos(np.deg2rad(40.0)), -np.sin(np.deg2rad(40.0)), 0.0, 1000.0]),
]

STATES = {
    "generic": [0.28, -0.11, -0.09, 0.06, -0.03, 0.04],
    "two_equal": [0.15, -0.07, -0.07, 0.0, 0.0, 0.0],
    "dilatation": [0.10, 0.10, 0.10, 0.0, 0.0, 0.0],
    "ground": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
}


def _umat(name, props, etot, F1, nstatev, Detot=None):
    n = 1
    z6 = lambda: np.zeros((6, n), order="F")
    eye = np.eye(3).reshape(3, 3, 1).copy(order="F")
    stress, sv, wm, Lt = sim.umat(
        name,
        np.asfortranarray(np.asarray(etot, dtype=float).reshape(6, 1)),
        np.asfortranarray(np.asarray(Detot, dtype=float).reshape(6, 1)) if Detot is not None else z6(),
        eye if F1 is not None else np.array([]),
        np.asarray(F1).reshape(3, 3, 1).copy(order="F") if F1 is not None else np.array([]),
        z6(), eye,
        np.asfortranarray(np.asarray(props, dtype=float).reshape(-1, 1)),
        np.zeros((nstatev, n), order="F"), 0.0, 1.0,
        np.zeros((4, n), order="F"), n_threads=1)
    return stress[:, 0], Lt[:, :, 0]


@pytest.mark.parametrize("volumetric", ["log", "quadratic"])
@pytest.mark.parametrize("state", list(STATES))
@pytest.mark.parametrize("block,umat_name,umat_props",
                         MODELS, ids=[m[1] for m in MODELS])
def test_modular_hyper_matches_standalone_kernel(block, umat_name, umat_props, state, volumetric):
    eps = np.array(STATES[state])
    if umat_name == "NEOHI_closed":
        if volumetric == "log":
            pytest.skip("NEOHI is the quadratic volumetric form")
    else:
        block = dataclasses.replace(block, volumetric=volumetric)
        # the law selects U(J) by an optional trailing prop (absent = log)
        umat_props = list(umat_props) + ([1.0] if volumetric == "quadratic" else [])
    mat = ModularMaterial(elasticity=block)

    sigma_mod, Lt_mod = _umat("MODUL", mat.props, eps, None, mat.nstatev)

    # through sim.umat with F the laws output Cauchy; MODUL on the small-strain arguments returns tau.
    # A box law reads Detot (the kernels built on F ignore it).
    sigma_ref, Lt_ref = _umat(umat_name, umat_props, np.zeros(6),
                              expm(sim.v2t_strain(eps)), 1, Detot=eps)
    tau_ref = np.exp(eps[:3].sum()) * sigma_ref

    scale = max(1.0, np.abs(tau_ref).max())
    assert np.abs(sigma_mod - tau_ref).max() < 1e-10 * scale
    assert np.abs(Lt_mod - Lt_ref).max() < 1e-9 * np.abs(Lt_ref).max()


@pytest.mark.parametrize("umat_name,props,kappa", [
    ("NEOHC", [0.5673, 1000.0], 1000.0), ("YEOHH", [0.30, -0.010, 0.0005, 1000.0], 1000.0),
    ("OGDEN", [2.0, 1000.0, 0.4, 1.3, 0.1, 5.0], 1000.0),
])
@pytest.mark.parametrize("volumetric,dUdJ,dU2dJ2", [
    ("log", lambda k, J: k * np.log(J), lambda k, J: k / J),
    ("quadratic", lambda k, J: k * (J - 1.0), lambda k, J: k),
])
def test_volumetric_potential_under_pure_dilatation(umat_name, props, kappa, volumetric, dUdJ, dU2dJ2):
    """Pure dilatation isolates U(J): sigma = U'(J) I, and the hydrostatic tangent is
    3 J (U' + J U'') on the Kirchhoff box; a finite difference confirms it."""
    props = list(props) + ([1.0] if volumetric == "quadratic" else [])
    J = 1.2
    F = np.cbrt(J) * np.eye(3)
    ln_V = lambda Jx: np.log(np.cbrt(Jx)) * np.array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0])   # the box reads Detot
    sigma, Lt = _umat(umat_name, props, np.zeros(6), F, 1, Detot=ln_V(J))
    np.testing.assert_allclose(sigma[:3], dUdJ(kappa, J), rtol=1e-10, atol=1e-12)
    np.testing.assert_allclose(sigma[3:], 0.0, atol=1e-12)
    np.testing.assert_allclose(Lt[0, :3].sum(), 3.0 * J * (dUdJ(kappa, J) + J * dU2dJ2(kappa, J)),
                               rtol=1e-9)
    # finite difference of tau_11 along a hydrostatic log-strain increment
    d = 1e-6
    tau = lambda Jx: Jx * _umat(umat_name, props, np.zeros(6), np.cbrt(Jx) * np.eye(3), 1, Detot=ln_V(Jx))[0][0]
    fd = (tau(J * np.exp(3 * d)) - tau(J * np.exp(-3 * d))) / (2 * d)
    np.testing.assert_allclose(Lt[0, :3].sum(), fd, rtol=1e-6)


def test_volumetric_selector_is_validated():
    with pytest.raises(ValueError):
        NeoHookeanElasticity(mu=0.5, kappa=1000.0, volumetric="cubic").to_props()
    with pytest.raises(ValueError, match="volumetric potential"):
        _umat("NEOHC", [0.5673, 1000.0, 2.0], np.zeros(6), np.eye(3), 1)


def test_yeoh_ground_state_is_L_iso():
    """L0 (the mechanisms' reference) is the closed-form ground state."""
    mat = ModularMaterial(
        elasticity=YeohElasticity(C10=0.30, C20=-0.010, C30=0.0005, kappa=1000.0))
    _, Lt = _umat("MODUL", mat.props, np.zeros(6), None, mat.nstatev)
    L_ref = np.asarray(sim.L_iso([1000.0, 2.0 * 0.30], "Kmu"))
    assert np.abs(Lt - L_ref).max() / np.abs(L_ref).max() < 1e-9


def test_props_roundtrip_layout():
    """[potential, n_params, params..., alpha, omega] preceded by the elasticity type."""
    block = YeohElasticity(C10=1.0, C20=2.0, C30=3.0, kappa=4.0, alpha=5.0)
    assert block.to_props() == [2.0, 5.0, 1.0, 2.0, 3.0, 4.0, 0.0, 5.0, 0.0]
    assert block.nprops == 9
    quad = YeohElasticity(C10=1.0, C20=2.0, C30=3.0, kappa=4.0, volumetric="quadratic")
    assert quad.to_props() == [2.0, 5.0, 1.0, 2.0, 3.0, 4.0, 1.0, 0.0, 0.0]
    mat = ModularMaterial(elasticity=block)
    # elasticity_type, then the block, then the mechanism count
    assert list(np.asarray(mat.props)) == [4.0] + block.to_props() + [0.0]


def test_swanson_rejects_malformed_terms():
    with pytest.raises(TypeError):
        SwansonElasticity(terms=((0.5, 0.1, 0.9),), kappa=1.0)


def test_potential_parameters_are_required_and_alpha_keyword_only():
    """Positional arguments fill the potential's parameters, never alpha."""
    block = NeoHookeanElasticity(0.5673, 1000.0)
    assert (block.mu, block.kappa, block.alpha) == (0.5673, 1000.0, 0.0)
    with pytest.raises(TypeError):
        NeoHookeanElasticity(0.5673, 1000.0, 1.2e-5)
    with pytest.raises(TypeError):
        YeohElasticity(C10=0.30)


def test_hyper_blocks_are_exported():
    import simcoon.modular as md

    for name in ("HyperPotential", "VolumetricPotential", "NeoHookeanElasticity", "MooneyRivlinElasticity",
                 "YeohElasticity", "IsiharaElasticity", "GentThomasElasticity",
                 "SwansonElasticity"):
        assert name in md.__all__


# ------------------------------------------------------------------ the adapters

ADAPTERS = ["NEOHC", "NEOHI", "MOORI", "YEOHH", "ISHAH", "GETHH", "SWANH", "HOLZA"]


def test_hyperelastic_names_run_the_modular_engine_under_log_R_only():
    """The plain names run the HYPER_INVARIANTS block, a potential of the logarithmic strain: the
    solver serves them under corate 3 and the total form; another corate is refused by name, the
    message pointing at the kernel built on F when one exists (NEOHC_closed, NEOHI_closed)."""
    from simcoon.solver import Block, StepMeca
    st = StepMeca(control="F", value=[1.05, 0.02, 0., 0., 0.98, 0., 0., 0., 1.01], time=1., ninc=2)
    props = {"NEOHC": [1000., 10000.], "NEOHI": [2000., 0.40, 0.], "YEOHH": [0.30, -0.010, 0.0005, 1000.0]}
    for name, p in props.items():
        for form in ({"corate": 3}, {"kinematic_form": "total"}):
            sim.solver.solve(Block(steps=[st], control_type="F"), name, p, 1, T_init=290., **form)
        with pytest.raises(Exception, match="corate_type = 3"):
            sim.solver.solve(Block(steps=[st], control_type="F"), name, p, 1, T_init=290., corate=0)
    with pytest.raises(Exception, match="NEOHC_closed"):
        sim.solver.solve(Block(steps=[st], control_type="F"), "NEOHC", props["NEOHC"], 1, T_init=290., corate=4)
    for name in ("NEOHC_closed", "NEOHI_closed"):   # built on F: every corate
        for corate in (0, 1, 2, 3, 4, 5):
            sim.solver.solve(Block(steps=[st], control_type="F"), name, props[name[:5]], 1, T_init=290., corate=corate)


def test_neohi_adapter_is_the_quadratic_neo_hookean():
    """NEOHI (E, nu, alpha): C10 = mu/2 and 1/D1 = kappa/2 with mu = E/2(1+nu), kappa = E/3(1-2nu),
    i.e. the NEOHC potential with the quadratic U(J). The adapter and the kernel built on F agree
    at round-off under corate 3 (stress, box tangent, Wm)."""
    eps = np.array(STATES["generic"])
    F1 = expm(sim.v2t_strain(eps))
    s_a, Lt_a = _umat("NEOHI", [2000.0, 0.40, 0.0], np.zeros(6), F1, 1, Detot=eps)
    s_c, Lt_c = _umat("NEOHI_closed", [2000.0, 0.40, 0.0], np.zeros(6), F1, 1)
    np.testing.assert_allclose(s_a, s_c, rtol=0, atol=1e-12 * np.abs(s_c).max())
    np.testing.assert_allclose(Lt_a, Lt_c, rtol=0, atol=1e-11 * np.abs(Lt_c).max())
    s_n, _ = _umat("NEOHC", [2000.0 / 2.8, 2000.0 / 0.6, 1.0], np.zeros(6), F1, 1, Detot=eps)
    np.testing.assert_allclose(s_a, s_n, rtol=1e-14)   # the props differ by the rounding of kappa
