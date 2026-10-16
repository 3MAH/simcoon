"""sim.umat(temp=..., DT=...): the temperature at the start of the increment and its increment."""

import numpy as np
import pytest

import simcoon as sim

E, NU, ALPHA = 70000.0, 0.3, 2.0e-5
T0, DT = 293.15, 40.0


def _col(a):
    return np.asfortranarray(np.asarray(a, dtype=float).reshape(-1, 1))


def _eliso(de, n=1, **kwargs):
    DR = np.asfortranarray(np.repeat(np.eye(3)[:, :, None], n, axis=2))
    sigma, statev, wm, Lt = sim.umat("ELISO", np.zeros((6, n), order="F"), np.asfortranarray(np.tile(_col(de), n)),
                                     np.empty(0), np.empty(0), np.zeros((6, n), order="F"), DR, _col([E, NU, ALPHA]),
                                     np.zeros((1, n), order="F"), 0.0, 1.0, np.zeros((4, n), order="F"),
                                     n_threads=1, start=True, **kwargs)
    return sigma


def test_free_thermal_expansion_is_stress_free():
    sigma = _eliso(ALPHA * DT * np.array([1, 1, 1, 0, 0, 0]), temp=np.array([T0]), DT=np.array([DT]))
    np.testing.assert_allclose(sigma, 0.0, atol=1e-9)


def test_blocked_heating_gives_the_thermal_stress():
    sigma = _eliso(np.zeros(6), temp=np.array([T0]), DT=np.array([DT]))
    np.testing.assert_allclose(sigma[:3, 0], -E / (1 - 2 * NU) * ALPHA * DT, rtol=1e-12)


def test_default_DT_is_zero():
    np.testing.assert_array_equal(_eliso(np.zeros(6), temp=np.array([T0])), 0.0)


def test_temp_and_DT_shared_or_per_point():
    shared = _eliso(np.zeros(6), n=3, temp=np.array([T0]), DT=np.array([DT]))
    per_point = _eliso(np.zeros(6), n=3, temp=np.full(3, T0), DT=np.array([0.0, DT, 2 * DT]))
    np.testing.assert_allclose(shared, np.tile(shared[:, :1], 3), rtol=1e-14)
    np.testing.assert_allclose(per_point[:, 2], 2 * per_point[:, 1], rtol=1e-12)
    np.testing.assert_array_equal(per_point[:, 0], 0.0)


@pytest.mark.parametrize("kwargs", [dict(DT=np.array([DT])),                               # DT without temp
                                    dict(temp=np.full(2, T0)),                             # wrong length
                                    dict(temp=np.full(3, T0), DT=np.full(2, DT))])
def test_invalid_temperature_input(kwargs):
    with pytest.raises(ValueError):
        _eliso(np.zeros(6), n=3, **kwargs)


def test_sma_closest_point_tangent_with_phase_dependent_expansion():
    """SMADI with alpha_A != alpha_M, heated while transforming: the thermal strain of the
    transformation enters the flows as Dalpha (T + DT - T_init); the closest-point tangent is then
    the exact derivative of the returned stress (it read Dalpha (T + DT) before)."""
    props = [0, 67538., 67538., 0.349, 0.349, 1e-5, 3e-5, 0., 0.0418, 0.021, 0., 10., 10., 300., 290.,
             295., 305., 0.2, 0.2, 0.2, 0.2, 300., 0.2, 2., 1e-6, 1e-3, 1., 1e8]
    T_a, T_b, n = 320., 350., 40
    DR = np.asfortranarray(np.eye(3)[:, :, None])
    sigma, statev, wm, etot = _col(np.zeros(6)), _col(np.zeros(17)), _col(np.zeros(4)), np.zeros(6)
    de = np.array([0.02, -0.01, -0.01, 0., 0., 0.]) / n
    for i in range(n):
        sigma, statev, wm, _ = sim.umat("SMADI", _col(etot), _col(de), np.empty(0), np.empty(0), sigma, DR,
                                        _col(props), statev, float(i), 0.01, wm,
                                        temp=np.array([T_a + (T_b - T_a) * i / n]), DT=np.array([(T_b - T_a) / n]),
                                        n_threads=1, start=(i == 0), tangent_mode=3)
        etot = etot + de
    assert statev[1, 0] > 0.1                                   # transforming

    def call(De):
        s, _, _, Lt = sim.umat("SMADI", _col(etot), _col(De), np.empty(0), np.empty(0), sigma.copy(order="F"),
                               DR, _col(props), statev.copy(order="F"), 1., 0.01, wm.copy(order="F"),
                               temp=np.array([T_b]), DT=np.array([0.5]), n_threads=1, start=False, tangent_mode=3)
        return s[:, 0], Lt[:, :, 0]

    De0, h = np.array([2e-4, -6e-5, -6e-5, 1e-4, 0., 0.]), 1e-7
    Lt = call(De0)[1]
    fd = np.column_stack([(call(De0 + h * e)[0] - call(De0 - h * e)[0]) / (2 * h) for e in np.eye(6)])
    assert np.abs(Lt - fd).max() < 1e-8 * np.abs(fd).max()
