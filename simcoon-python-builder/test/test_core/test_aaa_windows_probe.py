"""TEMPORARY Windows narrowing probe (v2) — delete once the access violation is pinned.

v1 result: probe 1 already died — ELISO, blocks built in Python, 5 increments. So the
crash is NOT MODUL-specific: the first ``solver_run`` call dies, whatever it is. The
earlier crashes landed on test_modular only because it was the first caller collected.

v1 had a design flaw: all five cases passed ``corate=1``, so it could not separate the
corate from the route. That matters, because master's Windows job runs ``solver_run``
happily 62+ times — and its corate parametrization is ``[0, 2, 5]`` while every other
solver_run test omits ``corate`` (default ``logarithmic_R`` = 3). Green-Naghdi (1) is
never routed through ``solver_run`` on master, yet every crashing call on this branch
uses it (test_modular's ``_run_case`` and all of v1).

Known objection, which this probe also tests: master DID run corate 1 with control
type 1 through the FILE binding (test_modular::_run_solver), reaching the same engine
and the same ``set_start(1)``, and it passed on Windows. So corate 1 alone cannot be
the whole story. Probes 2-5 vary ONLY the corate, on an otherwise master-identical call.

The numbering starts at 2: a probe for ``sim.umat`` alone was dropped as redundant,
because ``run_test.py::test_umat_ogden_*`` already PASSED at 1-2 % in the very run that
crashed — the module loads and the constitutive call works on Windows.

Reading (a native crash ends the process, so the last PASSED line names the culprit):
  all pass ................ neither the route nor the corate; look at the fixture
  dies at 2 ............... the memory route itself, on the exact call master passes
  2-4 pass, dies at 5 ..... corate 1 (Green-Naghdi) through solver_run
  5 passes, dies at 6 ..... the blocks parsed from the legacy path file
"""

from pathlib import Path

import numpy as np

import simcoon as sim
from simcoon.solver import Block, StepMeca

DATA_DIR = Path(__file__).resolve().parents[3] / "examples" / "data"

ELISO_PROPS = np.array([210000.0, 0.3, 1.2e-5])
UNIAXIAL = ["strain"] + ["stress"] * 5


def _block():
    """One short small-strain uniaxial step, built in Python (no file)."""
    return [Block(steps=[StepMeca(control=UNIAXIAL, value=[0.002, 0, 0, 0, 0, 0],
                                  ninc=5, time=1.0)], ncycle=1)]


def test_probe_2_solve_default_corate():
    """Master-identical call: solve() with no corate argument (default = 3)."""
    res = sim.solver.solve(_block(), "ELISO", ELISO_PROPS, 1, T_init=290.0)
    assert len(res) > 0


def test_probe_3_solve_corate_0_jaumann():
    """corate 0 — exercised by master's parametrization [0, 2, 5]."""
    res = sim.solver.solve(_block(), "ELISO", ELISO_PROPS, 1, T_init=290.0, corate=0)
    assert len(res) > 0


def test_probe_4_solve_corate_2_logarithmic():
    """corate 2 — also exercised by master."""
    res = sim.solver.solve(_block(), "ELISO", ELISO_PROPS, 1, T_init=290.0, corate=2)
    assert len(res) > 0


def test_probe_5_solve_corate_1_green_naghdi():
    """corate 1 — NEVER routed through solver_run on master. The suspect."""
    res = sim.solver.solve(_block(), "ELISO", ELISO_PROPS, 1, T_init=290.0, corate=1)
    assert len(res) > 0


def test_probe_6_solve_corate_1_blocks_from_file():
    """Same corate, blocks parsed from the legacy path file: isolates the parsing."""
    blocks, T_init = sim.solver.from_file(str(DATA_DIR), "MODUL_path.txt")
    res = sim.solver.solve(blocks, "ELISO", ELISO_PROPS, 1, T_init=T_init, corate=1)
    assert len(res) > 0
