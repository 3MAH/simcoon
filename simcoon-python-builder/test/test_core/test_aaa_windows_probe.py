"""TEMPORARY Windows narrowing probe — delete once the access violation is pinned.

On this branch the first ``_core.solver_run`` call dies on windows-latest (conda netlib),
while master runs it 62+ times green on the same job and Linux/macOS stay green throughout.
v1 ruled out MODUL, v2 ruled out the corate.

This one dissects the inside of ``solver_run`` rather than varying its inputs: probes 2 and 3
must raise a CLEAN Python exception (``phases=[]`` on a homogeneous model reaches
``make_sub_phases``'s throw; a bad ``cBC_meca`` flag reaches the block-parsing throw). If the
exception surfaces, that layer is sound; if the process dies there instead, it is the culprit.

A native crash ends pytest, so each case is its own test, ordered plainest-first, in a module
named to collect before the rest. The last PASSED line names the culprit:

  1 → pure Python parsing   2 → argument marshalling   3 → block parsing
  4 → the engine            5 → what solve() adds on top of _core.solver_run
"""

from pathlib import Path

import numpy as np
import pytest

import simcoon as sim
from simcoon.solver import Block, StepMeca

DATA_DIR = Path(__file__).resolve().parents[3] / "examples" / "data"

#the constants master's passing tests use
ELISO_PROPS = [70000.0, 0.3, 1.0e-5]
UNIAXIAL = ["strain"] + ["stress"] * 5
T_INIT = 290.0


def _blocks_py(ninc=1):
    """The dicts solve() hands to the binding, built through the same code."""
    step = StepMeca(control=UNIAXIAL, value=[0.002, 0, 0, 0, 0, 0], ninc=ninc, time=1.0)
    return [Block(steps=[step], ncycle=1).to_dict(T_INIT)]


def test_probe_1_from_file_only_no_cpp():
    """Parse a legacy path file; call nothing in _core."""
    blocks, T_init = sim.solver.from_file(str(DATA_DIR), "MODUL_path.txt")
    assert len(blocks) > 0 and T_init > 0.0


def test_probe_2_make_sub_phases_throw_is_reached():
    """`phases=[]` on a homogeneous model must raise, not crash."""
    with pytest.raises(Exception):
        sim._core.solver_run(_blocks_py(), T_INIT, "ELISO",
                             np.asarray(ELISO_PROPS, dtype=float), 1,
                             phases=[])


def test_probe_3_block_parsing_throw_is_reached():
    """An invalid cBC_meca flag must raise from the block parsing."""
    bad = _blocks_py()
    bad[0]["steps"][0]["cBC_meca"] = np.array([9, 1, 1, 1, 1, 1])
    with pytest.raises(Exception):
        sim._core.solver_run(bad, T_INIT, "ELISO",
                             np.asarray(ELISO_PROPS, dtype=float), 1)


def test_probe_4_core_solver_run_minimal():
    """The engine: one increment, no tangent history, straight to _core."""
    res = sim._core.solver_run(_blocks_py(ninc=1), T_INIT, "ELISO",
                               np.asarray(ELISO_PROPS, dtype=float), 1,
                               record_tangent=False)
    assert res["status"] == 0


def test_probe_5_solve_master_identical():
    """Master's own call shape: a bare StepMeca, ninc=50, props as a plain list."""
    step = StepMeca(control=UNIAXIAL, value=[0.01, 0, 0, 0, 0, 0], ninc=50)
    res = sim.solver.solve(step, "ELISO", ELISO_PROPS, 1, T_init=T_INIT)
    assert len(res) > 0
