"""sim.umat / sim.umat_T(start=...): who decides that a material point is fresh.

start=True re-initialises the point (statev[0] = T, stress, internal variables, Wm).
By default it is inferred from time <= 1e-9, so every call at time 0 re-initialises --
a coupler's Newton corrections of the first increment included. start=False must keep
the point's state whatever the time.
"""

import numpy as np

import simcoon as sim

N = 3
EPICP = [200000.0, 0.3, 0.0, 300.0, 1000.0, 0.5]


def _col(a, n=N):
    return np.asfortranarray(np.tile(np.asarray(a, dtype=float).reshape(-1, 1), (1, n)))


def _eye():
    return np.asfortranarray(np.repeat(np.eye(3)[:, :, None], N, axis=2))


DE = [4e-3, -1e-3, -1e-3, 0, 0, 0]


def _epicp(sigma, statev, wm, etot, de, temp, time, **kwargs):
    return sim.umat("EPICP", _col(etot), _col(de), np.empty(0), np.empty(0), sigma, _eye(),
                    _col(EPICP), statev, time, 1.0, wm, temp=np.full(N, temp), n_threads=1, **kwargs)


def _plastic_state():
    # past the yield stress: non-zero stress, plastic strain and work
    sigma, statev, wm, _ = _epicp(_col(np.zeros(6)), _col([0.0] * 8), _col(np.zeros(4)),
                                  np.zeros(6), DE, 293.15, 0.0, start=True)
    assert statev[1].min() > 0.0  # accumulated plastic strain p
    np.testing.assert_array_equal(statev[0], 293.15)
    return sigma, statev, wm


def test_start_false_keeps_the_point_at_time_zero():
    sigma, statev, wm = _plastic_state()
    s1, v1, w1, _ = _epicp(sigma.copy(order="F"), statev.copy(order="F"), wm.copy(order="F"),
                           DE, np.zeros(6), 350.0, 0.0, start=False)
    np.testing.assert_array_equal(v1[0], 293.15)  # the reference temperature survives
    np.testing.assert_allclose(v1, statev, rtol=1e-12, atol=1e-15)  # so does the plastic state
    np.testing.assert_allclose(s1, sigma, rtol=1e-12, atol=1e-9)
    np.testing.assert_allclose(w1, wm, rtol=1e-12, atol=1e-12)


def test_default_still_infers_start_from_time():
    sigma, statev, wm = _plastic_state()
    args = lambda: (sigma.copy(order="F"), statev.copy(order="F"), wm.copy(order="F"), DE, np.zeros(6), 350.0)
    s0, v0, w0, _ = _epicp(*args(), 0.0)
    np.testing.assert_array_equal(v0[0], 350.0)  # time 0: re-initialised (legacy behaviour)
    s_true, v_true, w_true, _ = _epicp(*args(), 0.0, start=True)
    for got, want in [(s0, s_true), (v0, v_true), (w0, w_true)]:
        np.testing.assert_array_equal(got, want)
    _, v2, _, _ = _epicp(*args(), 1.0)
    np.testing.assert_array_equal(v2[0], 293.15)  # time > 0: not re-initialised


def test_umat_T_start():
    props = _col([7800.0, 460.0, 200000.0, 0.3, 1e-5])

    def call(t_init, T, **kwargs):
        return sim.umat_T("ELISO", _col(np.zeros(6)), _col(np.zeros(6)), _col(np.zeros(6)), _eye(), props,
                          _col([t_init]), 0.0, 1.0, _col(np.zeros(4)), _col(np.zeros(3)),
                          np.full(N, T), np.zeros(N), n_threads=1, **kwargs)

    np.testing.assert_array_equal(call(293.15, 350.0, start=False)[1][0], 293.15)
    np.testing.assert_array_equal(call(293.15, 350.0)[1][0], 350.0)
    np.testing.assert_array_equal(call(293.15, 350.0, start=True)[1][0], 350.0)
