"""HYPOO against ELORT: the rate and the total form of the same orthotropic stiffness.

Both are Kirchhoff-native, so any difference between them is the integration form, never
the stress measure. For an isotropic L they are the same law on any path, rotation included
(the transported strain and stress satisfy the same recursion). For an anisotropic L they
differ under rotation, by an amount of order one: HYPOO is kept as that reference.
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


def _gap(props, ninc=100):
    tau_e, tau_h = _shear("ELORT", props, ninc), _shear("HYPOO", props, ninc)
    return np.abs(tau_h - tau_e).max() / np.abs(tau_e).max()


@pytest.mark.parametrize("corate", [0, 3])
def test_hypoo_equals_elort_without_rotation(corate):
    """Same L, same transported strain, no rotation: tau, Cauchy output and Wm all agree."""
    e, h = _uniaxial("ELORT", ORT, corate), _uniaxial("HYPOO", ORT, corate)
    for key in ("Kirchhoff", "Stress", "Wm"):
        np.testing.assert_allclose(h[key], e[key], rtol=1e-12,
                                   atol=1e-12 * np.abs(e[key]).max(), err_msg=key)


def test_hypoo_equals_elort_in_shear_for_isotropic_L():
    """Isotropic L: the rate and the total form coincide exactly, even in simple shear.

    This used to show a ~1 % first-order gap, which was the transport defect of the finite
    route (the total form received its start strain untransported), not a property of the
    rate form.
    """
    assert _gap(ISO) < 1e-12


def test_hypoo_departs_from_elort_in_shear_for_orthotropic_L():
    """Orthotropic L: the two are different laws under rotation (reference effect)."""
    assert _gap(ORT) > 0.1
