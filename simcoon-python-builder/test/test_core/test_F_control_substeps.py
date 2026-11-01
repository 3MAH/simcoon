"""Kinematic control (F, grad u): an increment cut into sub-increments advances F along the
increment's geodesic, F1 = exp(Dtinc dL) F0 with dL the increment's generator, as the other control
types advance by Dtinc times their increment; the last sub-increment lands on the tabulated F.
"""

import numpy as np
import pytest
from scipy.linalg import expm, logm

from simcoon.solver import Block, StepMeca, StepThermomeca, solve

T0 = 293.15
# Johnson-Cook: rate-dependent yield, so the sub-increments' time fractions matter
EPJCK = np.array([200000.0, 0.33, 1.0e-5, 792.0, 510.0, 0.26, 0.014, 1.0, 1.03, 293.0, 1793.0])
ZENER_T = np.array([4.4, 0.656, 3000., 0.35, 1e-4, 1500., 0.35, 3000., 1200.])
F_END = np.array([[1.15, 0.3, 0.], [0., 0.95, 0.], [0., 0., 1.02]])
FIXED = dict(mul_tnew_dt=1.0)   # sub-increments keep their size (no growth on fast convergence)


def _last(res, key="Kirchhoff"):
    v = np.asarray(res[key])
    return v[..., -1]


def test_cut_increment_is_the_finer_path_mechanical():
    """One increment run as 4 sub-increments (Dn_init = 0.25) is the 4-increment step: same
    geodesic, same time fractions, so the same viscoplastic response."""
    cut = StepMeca(control="F", value=F_END.ravel(), ninc=1, time=0.01, Dn_init=0.25)
    fine = StepMeca(control="F", value=F_END.ravel(), ninc=4, time=0.01)
    r_cut = solve(Block(steps=[cut], control_type="F"), "EPJCK", EPJCK, 9, T_init=T0, **FIXED)
    r_fine = solve(Block(steps=[fine], control_type="F"), "EPJCK", EPJCK, 9, T_init=T0, **FIXED)
    r_whole = solve(Block(steps=[StepMeca(control="F", value=F_END.ravel(), ninc=1, time=0.01)], control_type="F"),
                    "EPJCK", EPJCK, 9, T_init=T0)
    tau = _last(r_fine)
    assert np.asarray(r_fine["Statev"])[1, -1] > 0.05   # p: the path is plastic
    np.testing.assert_allclose(_last(r_cut), tau, rtol=0, atol=1e-9 * np.abs(tau).max())
    # the sub-increments matter: the single increment is another (coarser) integration
    assert np.abs(_last(r_whole) - tau).max() > 1e-4 * np.abs(tau).max()


@pytest.mark.parametrize("form", ["rate", "total"])
def test_cut_increment_is_the_finer_path_thermomechanical(form):
    """Same for a thermomechanical viscoelastic law with a heat flux (adiabatic): the temperature
    and the stress of the cut increment are those of the finer step."""
    kw = dict(control="F", value=F_END.ravel(), time=0.5, thermal_control="heat_flux", Q=0.)
    cut = StepThermomeca(ninc=1, Dn_init=0.25, **kw)
    fine = StepThermomeca(ninc=4, **kw)
    r_cut = solve(Block(steps=[cut], control_type="F"), "ZENER", ZENER_T, 8, T_init=T0, kinematic_form=form, **FIXED)
    r_fine = solve(Block(steps=[fine], control_type="F"), "ZENER", ZENER_T, 8, T_init=T0, kinematic_form=form, **FIXED)
    tau = _last(r_fine)
    np.testing.assert_allclose(_last(r_cut), tau, rtol=0, atol=1e-9 * np.abs(tau).max())
    np.testing.assert_allclose(_last(r_cut, "Temp"), _last(r_fine, "Temp"), rtol=0, atol=1e-10)


def test_grad_u_control_is_F_minus_identity():
    """Control type 'gradU' prescribes grad u = F - I."""
    G = F_END - np.eye(3)
    r_G = solve(Block(steps=[StepMeca(control="F", value=G.ravel(), ninc=10)], control_type="gradU"), "EPICP",
                np.array([210000., 0.3, 1.2e-5, 300., 500., 0.4]), 8, T_init=T0)
    r_F = solve(Block(steps=[StepMeca(control="F", value=F_END.ravel(), ninc=10)], control_type="F"), "EPICP",
                np.array([210000., 0.3, 1.2e-5, 300., 500., 0.4]), 8, T_init=T0)
    np.testing.assert_allclose(_last(r_G), _last(r_F), rtol=0, atol=1e-10 * np.abs(_last(r_F)).max())


def test_tabular_F_on_the_geodesic_is_the_linear_step():
    """A table whose rows sit on the geodesic of the linear step gives the linear step: the
    increments' generators are read back from the rows, and a cut increment follows them."""
    n, t_end = 5, 0.01
    L = logm(F_END).real
    table = np.array([np.concatenate([[t_end * (i + 1) / n], expm(L * (i + 1) / n).ravel()]) for i in range(n)])
    tab = StepMeca(control="F", mode="tabular", tabular=table, Dn_init=0.5)
    lin = StepMeca(control="F", value=F_END.ravel(), ninc=2 * n, time=t_end)
    r_tab = solve(Block(steps=[tab], control_type="F"), "EPJCK", EPJCK, 9, T_init=T0, **FIXED)
    r_lin = solve(Block(steps=[lin], control_type="F"), "EPJCK", EPJCK, 9, T_init=T0, **FIXED)
    tau = _last(r_lin)
    np.testing.assert_allclose(_last(r_tab), tau, rtol=0, atol=1e-8 * np.abs(tau).max())
