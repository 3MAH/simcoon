"""Concurrent mean-field solves with different Eshelby quadratures.

The Gauss tables of the Eshelby integrals used to be process-wide statics, rewritten at the
start of every ellipsoidal solve and every L_eff while solver_run releases the GIL: two
threads solving with different (mp, np) read each other's tables (IndexError from the
Armadillo bounds check, or silently wrong results). Each (mp, np) now has its own immutable
tables, read by every pass of the scheme.
"""

from concurrent.futures import ThreadPoolExecutor

import numpy as np

import simcoon as sim
from simcoon import solver as slv
from simcoon.solver.micromechanics import Ellipsoid


def _phases(elastic=False):
    matrix = Ellipsoid(number=0, umat_name="EPICP", save=1, concentration=0.7, nstatev=8,
                       props=np.array([2250.0, 0.19, 8.8e-5, 30.0, 300.0, 0.4]))
    if elastic:
        matrix = Ellipsoid(number=0, umat_name="ELISO", save=1, concentration=0.7, nstatev=1,
                           props=np.array([2250.0, 0.19, 8.8e-5]))
    fibre = Ellipsoid(number=1, umat_name="ELISO", save=1, concentration=0.3, nstatev=1,
                      props=np.array([73000.0, 0.22, 0.5e-6]), a1=2.0)
    return [matrix, fibre]


_STEP = slv.StepMeca(control=["strain"] + ["stress"] * 5, value=[0.03, 0, 0, 0, 0, 0],
                     time=1.0, ninc=50)


def _solve(mp):
    return slv.solve(_STEP, "MIMTN", np.array([mp, mp, 0.0]), 1, phases=_phases(),
                     tangent_mode="continuum")["Stress"][0, -1]


def _l_eff(mp):
    return np.asarray(sim.L_eff("MISCN", np.array([mp, mp, 0.0]), 1, phases=_phases(elastic=True)))


def test_threaded_solves_match_the_sequential_ones():
    ref = {mp: _solve(mp) for mp in (10, 40)}
    mps = [10, 40] * 8
    with ThreadPoolExecutor(8) as ex:
        out = list(ex.map(_solve, mps))
    assert [o == ref[m] for o, m in zip(out, mps)] == [True] * len(mps)


def test_threaded_l_eff_and_solves_do_not_interfere():
    ref_solve = {mp: _solve(mp) for mp in (10, 40)}
    ref_leff = {mp: _l_eff(mp) for mp in (12, 30)}
    jobs = [("solve", 10), ("leff", 12), ("solve", 40), ("leff", 30)] * 4
    def run(job):
        kind, mp = job
        return _solve(mp) if kind == "solve" else _l_eff(mp)
    with ThreadPoolExecutor(8) as ex:
        out = list(ex.map(run, jobs))
    for (kind, mp), o in zip(jobs, out):
        expected = ref_solve[mp] if kind == "solve" else ref_leff[mp]
        np.testing.assert_array_equal(o, expected)
