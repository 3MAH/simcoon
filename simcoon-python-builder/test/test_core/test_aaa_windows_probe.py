"""TEMPORARY Windows narrowing probe — delete once the access violation is pinned.

``_core.solver_run`` dies on windows-latest (conda netlib) on this branch; master runs it 62+
times green on the same job, and Linux/macOS stay green throughout.

Ruled out so far, each by a probe that PASSED on Windows: MODUL as a law (v1), the corate (v2),
the whole binding prologue (v3), the engine entry and first UMAT call (v4). v5 was decisive:
EPICP, a DEDICATED kernel, passes BOTH its nstatev guard and a full increment on Windows, while
ELISO — routed to the ``umat_legacy_modular`` adapter (``select_umat_M`` maps it to 201) — dies.
EPICP's full run also exonerates the increment loop, the memory sink and the result assembly,
by a case that passes rather than by elimination.

So the fault sits in the adapter route. This probe splits that route in two, and the split is
clean because ``translate_ELISO`` returns exactly ``{0, 0, E, nu, alpha, 0}`` — the very props
a purely elastic MODUL material carries. ``umat_modular`` therefore receives the SAME input
either way; only the road differs:

    MODUL  -> umat_modular directly (no registry, no translator)
    ELISO  -> registry() lookup -> translate_ELISO -> umat_modular

Both reach the same throw (``ModularUMAT::initialize``, modular_umat.cpp:170). Markers measured
locally first:

    MODUL nstatev=0 -> RuntimeError    MODUL valid -> dict, status 0
    ELISO nstatev=0 -> RuntimeError    ELISO valid -> dict, status 0

MODUL comes FIRST: a native crash ends the process, and the case known to kill Windows must not
take the informative ones with it.

Reading:
  MODUL passes, ELISO dies .... the adapter road only: registry()'s function-local static or
                               the translator — both trivial, which would point at static
                               initialisation inside the DLL
  MODUL dies too ............. ModularUMAT construction itself, adapter exonerated
  everything passes .......... the fault needs the fuller call solve() makes

Standing tension, not to be papered over: this code is identical to master (umat_smart.cpp and
tangent_assembly.cpp unchanged; only viscoelastic_mechanism.cpp moved 6 lines), and master
drives ELISO through this exact path green on Windows. Expect a latent fault made fatal by the
new binary layout, as in the 2.0.1 incident where the OOB was heap-layout dependent.
"""

import numpy as np
import pytest

import simcoon as sim
from simcoon.solver import Block, StepMeca

#ELISO's legacy props, and what translate_ELISO turns them into
ELISO_PROPS = np.array([70000.0, 0.3, 1.0e-5])
MODUL_PROPS = np.array([0.0, 0.0, 70000.0, 0.3, 1.0e-5, 0.0])
UNIAXIAL = ["strain"] + ["stress"] * 5
T_INIT = 290.0


def _blocks_py(ninc=1):
    """The dicts solve() hands to the binding, built through the same code."""
    step = StepMeca(control=UNIAXIAL, value=[0.002, 0, 0, 0, 0, 0], ninc=ninc, time=1.0)
    return [Block(steps=[step], ncycle=1).to_dict(T_INIT)]


def test_probe_1_modul_direct_rejects_nstatev():
    """MODUL reaches ModularUMAT::initialize without registry or translator."""
    with pytest.raises(Exception):
        sim._core.solver_run(_blocks_py(), T_INIT, "MODUL", MODUL_PROPS, 0,
                             record_tangent=False)


def test_probe_2_modul_direct_runs():
    """A full increment through ModularUMAT, reached directly."""
    res = sim._core.solver_run(_blocks_py(ninc=1), T_INIT, "MODUL", MODUL_PROPS, 1,
                               record_tangent=False)
    assert res["status"] == 0


def test_probe_3_adapter_rejects_nstatev():
    """Same throw, reached through registry() + translate_ELISO."""
    with pytest.raises(Exception):
        sim._core.solver_run(_blocks_py(), T_INIT, "ELISO", ELISO_PROPS, 0,
                             record_tangent=False)


def test_probe_4_adapter_runs():
    """Same increment, reached through the adapter: the original crash."""
    res = sim._core.solver_run(_blocks_py(ninc=1), T_INIT, "ELISO", ELISO_PROPS, 1,
                               record_tangent=False)
    assert res["status"] == 0
