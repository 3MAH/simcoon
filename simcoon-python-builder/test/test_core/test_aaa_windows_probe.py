"""TEMPORARY Windows narrowing probe — delete once the access violation is pinned.

On this branch the first ``_core.solver_run`` call dies on windows-latest (conda netlib),
while master runs it 62+ times green on the same job and Linux/macOS stay green throughout.

Ruled out so far, each by a probe that PASSED on Windows: MODUL (v1), the corate (v2), and
the whole binding prologue (v3 — argument marshalling, block parsing, ``solver_output`` and
``check_path_output`` all run, since ``make_sub_phases``'s throw sits after them and was
reached cleanly). v3's probe 4 also used ``record_tangent=False``, so the tangent history is
out. What remains is ``simcoon::solver_run`` itself, under ``gil_scoped_release``, or the
result assembly after it.

This one walks a ladder of throws that sit deeper and deeper inside the engine. Each must
raise a CLEAN Python exception; the one that kills the process instead marks how far the
engine got. All three were measured locally first — a probe is only an instrument if its
expected behaviour is known.

  1 tangent_mode=99 → throws in solver_run's first lines ......... engine entered
  2 umat_name unknown → throws at the first UMAT call ............ inside the increment loop
  3 nstatev=0 → throws building the constitutive kernel .......... deeper still
  4 the minimal valid run ....................................... completes and assembles

Reading: the last PASSED line names how far it got. If 1-3 pass and only 4 dies, the engine
runs as far as the kernel and the fault is in the rest of the increment loop or in the result
assembly (rows_to_arr / scalars_to_arr), which no probe can isolate directly.
"""

import numpy as np
import pytest

import simcoon as sim
from simcoon.solver import Block, StepMeca

#the constants master's passing tests use
ELISO_PROPS = np.array([70000.0, 0.3, 1.0e-5])
UNIAXIAL = ["strain"] + ["stress"] * 5
T_INIT = 290.0


def _blocks_py(ninc=1):
    """The dicts solve() hands to the binding, built through the same code."""
    step = StepMeca(control=UNIAXIAL, value=[0.002, 0, 0, 0, 0, 0], ninc=ninc, time=1.0)
    return [Block(steps=[step], ncycle=1).to_dict(T_INIT)]


def test_probe_1_engine_entered():
    """solver_run validates tangent_mode in its first lines."""
    with pytest.raises(Exception):
        sim._core.solver_run(_blocks_py(), T_INIT, "ELISO", ELISO_PROPS, 1,
                             params={"tangent_mode": 99})


def test_probe_2_first_umat_call_reached():
    """An unknown law throws from the small-strain dispatch, inside the increment loop."""
    with pytest.raises(Exception):
        sim._core.solver_run(_blocks_py(), T_INIT, "XXXXX", ELISO_PROPS, 1,
                             record_tangent=False)


def test_probe_3_kernel_construction_reached():
    """nstatev = 0 throws while the constitutive kernel is being built."""
    with pytest.raises(Exception):
        sim._core.solver_run(_blocks_py(), T_INIT, "ELISO", ELISO_PROPS, 0,
                             record_tangent=False)


def test_probe_4_minimal_valid_run():
    """The case that dies on Windows: one increment, no tangent history."""
    res = sim._core.solver_run(_blocks_py(ninc=1), T_INIT, "ELISO", ELISO_PROPS, 1,
                               record_tangent=False)
    assert res["status"] == 0
