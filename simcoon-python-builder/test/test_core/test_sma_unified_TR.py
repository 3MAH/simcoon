"""Shape-memory path (cool at zero stress below Mf, load the martensite, unload, heat above Af)
for SMRDI / SMRAI against SMADI. Guards the zero-stress active set of umat_sma_unified_TR and
the magnitude normalisation of the convergence measure of the SMA kernels (see the notes in
unified_TR.hpp and unified_T.hpp).
"""

import time

import numpy as np
import pytest

import simcoon as sim
from simcoon.solver import Block, StepMeca

# SMADI block (28 props): flagT, E_A, E_M, nu_A, nu_M, alpha_A, alpha_M, Hmin, Hmax, k1, sigmacrit,
# C_A, C_M, Ms0, Mf0, As0, Af0, n1..n4, sigmacaliber, prager_b, prager_n, c/p0/n/alpha_lambda
SMADI = [0, 70000.0, 70000.0, 0.3, 0.3, 1e-6, 1e-6, 0.0, 0.05, 0.021, 0.0, 6.0, 5.0,
         293.15, 273.15, 313.15, 333.15, 0.2, 0.2, 0.2, 0.2, 300.0, 1.4, 2.0,
         1e-6, 1e-3, 1.0, 1e8]
DFA_ISO = [0.5, 0.5, 0.5, 1.5, 1.5, 1.5, 0.0]  # F, G, H, L, M, N, K: von Mises operator
TR_VARIANTS = [("SMRDI", []), ("SMRAI", DFA_ISO)]
T_START, T_COLD, T_HOT = 353.15, 253.15, 373.15  # above Af, below Mf, above Af
EPS_MAX = 0.04
I_XI = 1


def reo(y_reo, h_reo=5000.0):
    return [y_reo, h_reo, 0.05, 1e-6, 1e-3, 1.0, 1e8]


def shape_memory_cycle(ninc=100):
    free = ["stress"] * 6
    uni = ["strain"] + ["stress"] * 5
    zero = [0.0] * 6
    return [Block(steps=[
        StepMeca(control=free, value=zero, time=1.0, ninc=ninc, T_final=T_COLD, Dn_mini=0.01),
        StepMeca(control=uni, value=[EPS_MAX, 0, 0, 0, 0, 0], time=1.0, ninc=ninc, Dn_mini=0.01),
        StepMeca(control=free, value=zero, time=1.0, ninc=ninc, Dn_mini=0.01),
        StepMeca(control=free, value=zero, time=1.0, ninc=2 * ninc, T_final=T_HOT, Dn_mini=0.01),
    ])]


def solve(name, props, nstatev):
    t0 = time.perf_counter()
    res = sim.solver.solve(shape_memory_cycle(), name, np.asarray(props, dtype=float), nstatev,
                           T_init=T_START, corate=3)
    elapsed = time.perf_counter() - t0
    assert res.status == 0
    assert np.isfinite(res["Stress"]).all()
    return res, elapsed


@pytest.fixture(scope="module")
def smadi():
    """Transformation-only reference, and the healthy wall-clock of this path on this machine."""
    return solve("SMADI", SMADI, 17)


def assert_no_crawl(elapsed, reference):
    # the regression was a ~100x slowdown (step cuts at the minimal sub-increment)
    assert elapsed < max(5.0 * reference, 2.0), "the solver crawled"


@pytest.mark.parametrize("name, extra", TR_VARIANTS)
def test_zero_stress_thermal_steps(name, extra, smadi):
    """Reorientation active: cooling at exactly zero stress, then the full shape-memory effect."""
    res, elapsed = solve(name, SMADI + extra + reo(200.0), 30)
    assert_no_crawl(elapsed, smadi[1])
    t, xi, e11 = res["Time"], res["Statev"][I_XI], res["Strain"][0]
    cooled = np.argmax(t >= 1.0) - 1
    assert xi[cooled] > 0.99 and abs(e11[cooled]) < 2e-4          # self-accommodated martensite
    assert e11[np.argmax(t >= 3.0) - 1] > 0.5 * EPS_MAX            # strain kept after unloading
    assert abs(e11[-1]) < 1e-4 and xi[-1] < 1e-3                  # recovered on heating


@pytest.mark.parametrize("name, extra", TR_VARIANTS)
def test_reorientation_off_loads_martensite_like_smadi(name, extra, smadi):
    """Y_Reo out of reach: the cooled martensite loads elastically to 2.8 GPa, as SMADI does,
    through the band where Y0t + D Hcur Mises(sigma) changes sign."""
    ref, t_ref = smadi
    res, elapsed = solve(name, SMADI + extra + reo(1e10), 30)
    assert_no_crawl(elapsed, t_ref)
    assert res["Stress"][0].max() > 2500.0
    np.testing.assert_allclose(res["Stress"], ref["Stress"], rtol=0, atol=0.5)
    np.testing.assert_allclose(res["Statev"][I_XI], ref["Statev"][I_XI], rtol=0, atol=1e-4)


@pytest.mark.parametrize("name, extra", TR_VARIANTS)
@pytest.mark.parametrize("y_reo", [100.0, 70.0, 30.0])
def test_superelastic_unloading_with_low_reorientation_limit(name, extra, y_reo):
    """Above Af with a low Y_Reo, the reverse reorientation surface is reached while the reverse
    transformation is running: both mechanisms are active together and their local Jacobian is
    not symmetric. The stress-controlled unloading must follow, and the superelastic loop close
    (regression: transposed multiplier elimination in assemble_continuum_tangent)."""
    free = ["stress"] * 6
    uni = ["strain"] + ["stress"] * 5
    blocks = [Block(steps=[
        StepMeca(control=uni, value=[EPS_MAX, 0, 0, 0, 0, 0], time=1.0, ninc=200, Dn_mini=0.01),
        StepMeca(control=free, value=[0.0] * 6, time=1.0, ninc=200, Dn_mini=0.01),
    ])]
    res = sim.solver.solve(blocks, name, np.asarray(SMADI + extra + reo(y_reo), dtype=float), 30,
                           T_init=T_START, corate=3)
    assert res.status == 0
    assert res["Time"][-1] == pytest.approx(2.0)
    assert res["Statev"][I_XI].max() > 0.5       # the transformation did take place
    assert res["Statev"][I_XI][-1] < 1e-3        # and reversed completely
    assert abs(res["Strain"][0][-1]) < 1e-4      # the loop closes
