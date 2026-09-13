"""TEMPORARY Windows narrowing probe — delete once the access violation is pinned.

``windows-latest`` dies with an access violation on the FIRST MODUL run through
``solver_run`` (``test_modular::test_elastic_convention_equivalence``), three times
identically, while Linux and macOS pass the whole suite and the C++ gtests pass on
Windows through the file driver. Everything that can be checked without Windows has
been eliminated: the executed C++ is identical to master for this case, ``so`` is
sized and filled correctly by the binding, ``ncycle > 1`` is already exercised
elsewhere, and the carma/allocator invariant is now satisfied by all 21 TUs.

The one discriminator left standing is **MODUL driven through solver_run**: on master
every MODUL run went through the file route, and the solver_run tests that pass on
Windows use ELISO/EPICP/SNTVE/NEOHC.

Method (the one that pinned the two earlier Windows crashes in this repo): a native
crash ends the pytest process, so each case is its OWN test and they are ordered
plainest-first. The last ``PASSED`` line in the CI log names the culprit. The module
is named ``test_aaa_*`` so it is collected before ``test_modular.py``, i.e. before the
process dies.

Expected reading:
  all four PASS ......... the crash needs something else in test_modular's fixture
  dies at probe 3 ....... MODUL itself, through solver_run
  dies at probe 4 ....... MODUL only with blocks parsed from the legacy path file
  dies at probe 2 ....... the parsed blocks, independently of MODUL
"""

from pathlib import Path

import numpy as np

import simcoon as sim
from simcoon.solver import Block, StepMeca
from simcoon.modular import ModularMaterial, IsotropicElasticity

DATA_DIR = Path(__file__).resolve().parents[3] / "examples" / "data"

ELISO_PROPS = np.array([210000.0, 0.3, 1.2e-5])
UNIAXIAL = ["strain"] + ["stress"] * 5

# The same isotropic material as test_elastic_convention_equivalence, elasticity only.
_MAT = ModularMaterial(elasticity=IsotropicElasticity(C1=210000.0, C2=0.3))


def _in_memory_block():
    """One short small-strain uniaxial step, built in Python (no file)."""
    return [Block(steps=[StepMeca(control=UNIAXIAL, value=[0.002, 0, 0, 0, 0, 0],
                                  ninc=5, time=1.0)], ncycle=1)]


def _parsed_blocks():
    """The very blocks test_modular uses, parsed from the legacy path file."""
    return sim.solver.from_file(str(DATA_DIR), "MODUL_path.txt")


def test_probe_1_eliso_blocks_in_memory():
    """Baseline: a legacy kernel through solver_run. Master does this and passes."""
    res = sim.solver.solve(_in_memory_block(), "ELISO", ELISO_PROPS, 1, corate=1)
    assert len(res) > 0


def test_probe_2_eliso_blocks_from_file():
    """Same kernel, but the blocks now come from MODUL_path.txt: isolates parsing."""
    blocks, T_init = _parsed_blocks()
    res = sim.solver.solve(blocks, "ELISO", ELISO_PROPS, 1, T_init=T_init, corate=1)
    assert len(res) > 0


def test_probe_3_modul_blocks_in_memory():
    """MODUL through solver_run, simplest possible path: isolates the kernel."""
    res = sim.solver.solve(_in_memory_block(), "MODUL", _MAT.props, _MAT.nstatev,
                           corate=1)
    assert len(res) > 0


def test_probe_4_modul_blocks_from_file():
    """The exact combination test_modular crashes on."""
    blocks, T_init = _parsed_blocks()
    res = sim.solver.solve(blocks, "MODUL", _MAT.props, _MAT.nstatev,
                           T_init=T_init, corate=1)
    assert len(res) > 0


def test_probe_5_modul_from_file_without_tangent():
    """Same as 4 with the tangent left unrecorded: isolates the tangent history."""
    blocks, T_init = _parsed_blocks()
    res = sim.solver.solve(blocks, "MODUL", _MAT.props, _MAT.nstatev,
                           T_init=T_init, corate=1, record_tangent=False)
    assert len(res) > 0
