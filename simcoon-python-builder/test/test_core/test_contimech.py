"""sph / dev at the python boundary: every accepted input shape and layout.

Regression for the 3x3 branches, which converted the matrix with the column converter
and raised for every 3x3 input, and for the Voigt sph path (libsimcoon has no Voigt sph;
the wrapper takes the exact complement of dev).
"""

import numpy as np
import pytest

import simcoon as sim

M = np.array([[1.0, 2.0, 3.0], [2.0, 4.0, 5.0], [3.0, 5.0, 6.0]])
V = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])


def test_sph_dev_3x3():
    tr = np.trace(M) / 3.0
    assert np.allclose(sim.sph(M), tr * np.eye(3))
    assert np.allclose(sim.dev(M), M - tr * np.eye(3))
    assert np.allclose(sim.sph(M) + sim.dev(M), M)


@pytest.mark.parametrize("shape", [(6,), (6, 1)])
def test_sph_dev_voigt(shape):
    v = V.reshape(shape)
    sph = np.asarray(sim.sph(v)).ravel()
    dev = np.asarray(sim.dev(v)).ravel()
    assert np.allclose(sph, [2.0, 2.0, 2.0, 0.0, 0.0, 0.0])
    assert np.allclose(dev, [-1.0, 0.0, 1.0, 4.0, 5.0, 6.0])
    assert np.allclose(sph + dev, V)


@pytest.mark.parametrize(
    "make",
    [
        lambda: np.ascontiguousarray(M),
        lambda: np.asfortranarray(M),
        lambda: np.arange(18.0).reshape(3, 6)[:, ::2],  # strided view
        lambda: np.asfortranarray(M)[::-1, :],  # reversed
    ],
    ids=["C", "F", "strided", "reversed"],
)
def test_dev_any_layout(make):
    m = make()
    ref = np.ascontiguousarray(m)
    assert np.allclose(sim.dev(m), ref - np.trace(ref) / 3.0 * np.eye(3))


def test_dev_readonly_input_untouched():
    m = np.asfortranarray(M)
    m.flags.writeable = False
    before = m.copy()
    out = sim.dev(m)
    assert np.allclose(out, M - np.trace(M) / 3.0 * np.eye(3))
    assert np.array_equal(m, before)
    assert not m.flags.writeable


@pytest.mark.parametrize("bad", [np.zeros((2, 2)), np.zeros(5), np.zeros((3, 3, 3))])
def test_sph_dev_reject_bad_shapes(bad):
    with pytest.raises(ValueError):
        sim.sph(bad)
    with pytest.raises(ValueError):
        sim.dev(bad)
