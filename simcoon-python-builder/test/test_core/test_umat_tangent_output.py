"""sim.umat(tangent_output=...) returns the tangent a finite-element coupler assembles.

'material' and 'spatial' must equal the box tangent converted by Lt_convert afterwards, which is
what couplers did before (one pass for dS/dE, two for the Lie tangent), with the conversion that
matches the corate the law ran with.
"""

import numpy as np
import pytest

import simcoon as sim

N = 150  # past the parallel cutoff (100)

# box -> dS/dE key of each corate in Lt_convert (corate 5 has no key: it is checked through 'spatial')
BOX_TO_DSDE = {0: "DsigmaDe_JaumannDD_2_DSDE", 1: "DsigmaDe_GreenNaghdiDD_2_DSDE",
               2: "DsigmaDe_2_DSDE", 3: "DsigmaDe_2_DSDE", 4: "Dsigma_LieDD_2_DSDE"}

LAWS = {  # name: (props, nstatev, strain increment)
    "EPICP": ([200000.0, 0.3, 0.0, 300.0, 1000.0, 0.5], 8, [4e-3, -1e-3, -1e-3, 2e-3, 0.0, 1e-3]),
    "NEOHC": ([1000.0, 5000.0], 1, [4e-2, -1e-2, -1e-2, 2e-2, 0.0, 1e-2]),
}


def _col(a, n=N):
    return np.asfortranarray(np.tile(np.asarray(a, dtype=float).reshape(-1, 1), (1, n)))


def _call(name, corate, **kwargs):
    props, nstatev, de = LAWS[name]
    rng = np.random.default_rng(1)
    F0 = np.asfortranarray(np.eye(3)[:, :, None] + 0.05 * rng.standard_normal((3, 3, N)))
    F1 = np.asfortranarray(F0 + 0.01 * rng.standard_normal((3, 3, N)))
    DR = np.asfortranarray(np.tile(np.eye(3)[:, :, None], (1, 1, N)))
    statev = np.zeros((nstatev, N), order="F")
    statev[0] = 290.0
    out = sim.umat(name, _col(np.zeros(6)), _col(de), F0, F1, _col([150.0, 20.0, 0, 10.0, 0, 0]), DR,
                   _col(props), statev, 0.5, 1.0, _col(np.zeros(4)), corate=corate, **kwargs)
    return out, F1


@pytest.mark.parametrize("name", sorted(LAWS))
@pytest.mark.parametrize("corate", range(6))
def test_fused_tangent_matches_lt_convert(name, corate):
    (sigma, _, _, box), F1 = _call(name, corate)
    (s_m, _, _, material), _ = _call(name, corate, tangent_output="material")
    (s_s, _, _, spatial), _ = _call(name, corate, tangent_output="spatial")
    np.testing.assert_array_equal(s_m, sigma)
    np.testing.assert_array_equal(s_s, sigma)
    scale = np.abs(box).max()
    if corate in BOX_TO_DSDE:
        np.testing.assert_allclose(material, sim.Lt_convert(box, F1, sigma, BOX_TO_DSDE[corate]),
                                   rtol=0, atol=1e-10 * scale)
    np.testing.assert_allclose(spatial, sim.Lt_convert(material, F1, sigma, "DSDE_2_Dsigma_LieDD"),
                               rtol=0, atol=1e-10 * scale)


def test_default_is_the_box():
    (_, _, _, box), _ = _call("EPICP", 3)
    (_, _, _, explicit), _ = _call("EPICP", 3, tangent_output="box")
    np.testing.assert_array_equal(explicit, box)


def test_unknown_tangent_output_raises():
    with pytest.raises(ValueError, match="tangent_output"):
        _call("EPICP", 3, tangent_output="lie")


def test_converted_tangent_needs_one_F_per_point():
    props, nstatev, de = LAWS["EPICP"]
    eye = np.asfortranarray(np.tile(np.eye(3)[:, :, None], (1, 1, N)))
    with pytest.raises(ValueError, match="F1"):
        sim.umat("EPICP", _col(np.zeros(6)), _col(de), eye, eye[:, :, :1].copy(order="F"),
                 _col(np.zeros(6)), eye, _col(props), _col([290.0] + [0.0] * 7), 0.5, 1.0,
                 _col(np.zeros(4)), tangent_output="spatial")
