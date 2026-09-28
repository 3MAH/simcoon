"""HYPOO against ELORT: the rate and the total form of the same orthotropic stiffness.

Both are Kirchhoff-native, so any difference between them is the integration form, never
the stress measure. Along a path without rotation increments the two are identical; under
rotation the rate form differs by a first-order time-discretisation gap (present even for an
isotropic stiffness) and, for an anisotropic one, by the non-commutation of the stress
transport with L -- the effect the pair exists to show.
"""

import numpy as np
import pytest

import simcoon as sim
from simcoon.solver import Block, StepMeca

ISO = [70000.] * 3 + [0.3] * 3 + [70000. / 2.6] * 3 + [0.] * 3
ORT = [70000., 30000., 15000., 0.3, 0.3, 0.3, 8000., 6000., 5000., 0., 0., 0.]


def _uniaxial(name, props, corate):
    step = StepMeca(control=["strain"] + ["stress"] * 5, value=[0.3, 0, 0, 0, 0, 0],
                    time=1.0, ninc=50, Dn_init=1.0, Dn_mini=1e-4)
    return sim.solver.solve(Block(steps=[step], control_type="logarithmic"), name,
                            np.asarray(props, float), 1, T_init=290.0, corate=corate)


def _shear(name, props, ninc, gamma=1.0):
    step = StepMeca(control="F", value=[1., gamma, 0., 0., 1., 0., 0., 0., 1.],
                    time=1.0, ninc=ninc, Dn_init=1.0, Dn_mini=1e-4)
    res = sim.solver.solve(Block(steps=[step], control_type="F"), name,
                           np.asarray(props, float), 1, T_init=290.0, corate=3)
    return res["Kirchhoff"][:, -1]


def _gap(props, ninc):
    tau_e, tau_h = _shear("ELORT", props, ninc), _shear("HYPOO", props, ninc)
    return np.abs(tau_h - tau_e).max() / np.abs(tau_e).max()


@pytest.mark.parametrize("corate", [0, 3])
def test_hypoo_equals_elort_without_rotation(corate):
    """Same L, same transported strain, no rotation: tau, Cauchy output and Wm all agree."""
    e, h = _uniaxial("ELORT", ORT, corate), _uniaxial("HYPOO", ORT, corate)
    for key in ("Kirchhoff", "Stress", "Wm"):
        np.testing.assert_allclose(h[key], e[key], rtol=1e-12,
                                   atol=1e-12 * np.abs(e[key]).max(), err_msg=key)


def test_hypoo_shear_gap_is_discretisation_for_iso_and_anisotropy_for_ortho():
    """Isotropic L: a first-order gap that halves with the increment. Orthotropic L: a gap
    of order one, far above the discretisation part."""
    iso_100, iso_200 = _gap(ISO, 100), _gap(ISO, 200)
    assert iso_100 < 0.02
    assert iso_200 < 0.6 * iso_100, "the isotropic gap must shrink at first order"

    ort_100 = _gap(ORT, 100)
    assert ort_100 > 20.0 * iso_100, "the anisotropic gap must dominate the discretisation"
