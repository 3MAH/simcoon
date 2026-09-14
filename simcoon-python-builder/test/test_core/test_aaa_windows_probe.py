"""TEMPORARY Windows narrowing probe — delete once the access violation is pinned.

On this branch ``_core.solver_run`` dies on windows-latest (conda netlib) while master runs it
62+ times green on the same job, and Linux/macOS stay green throughout.

Ruled out so far, each by a probe that PASSED on Windows: MODUL (v1), the corate (v2), the
whole binding prologue (v3), and — from the v4 ladder — the engine entry and the first UMAT
call inside the increment loop. v4 localised the fault: ``nstatev=0``, which on macOS raises
``RuntimeError: 'ELISO' (modular adapter): ModularUMAT: nstatev (0) < required (1)`` from
``legacy_adapters.cpp:224``, CRASHES on Windows instead of raising.

``umat_legacy_modular`` builds a ``ModularUMAT`` on every call. ELISO is routed to it
(``select_umat_M`` maps ELISO -> 201, verified in umat_smart.cpp:274), while EPICP is a
dedicated kernel (-> 6) that never touches the adapter. This probe opposes the two on the
same call, with markers measured locally first:

    EPICP nstatev=0 -> IndexError      EPICP valid run -> dict, status 0
    ELISO nstatev=0 -> RuntimeError    ELISO valid run -> dict, status 0

Both EPICP cases come FIRST: a native crash ends the process, and the case already known to
kill Windows must not take the informative ones with it.

Reading:
  EPICP cases pass, ELISO cases die .... the modular adapter path (ModularUMAT construction)
  EPICP dies too ....................... the increment loop itself, whatever the kernel
  everything passes .................... the fault needs the fuller call solve() makes

Caveat worth keeping in view: the branch changed almost nothing in this subtree
(umat_smart.cpp and tangent_assembly.cpp are identical to master, only
viscoelastic_mechanism.cpp moved by 6 lines), and ELISO -> adapter -> ModularUMAT exists
unchanged on master, where it passes on Windows. So "the adapter is broken" cannot be the
whole story; a layout-sensitive latent fault would fit the evidence better.
"""

import numpy as np
import pytest

import simcoon as sim
from simcoon.solver import Block, StepMeca

#the constants master's passing tests use
ELISO_PROPS = np.array([70000.0, 0.3, 1.0e-5])
EPICP_PROPS = np.array([70000.0, 0.3, 1.0e-5, 300.0, 1000.0, 0.3])
UNIAXIAL = ["strain"] + ["stress"] * 5
T_INIT = 290.0


def _blocks_py(ninc=1):
    """The dicts solve() hands to the binding, built through the same code."""
    step = StepMeca(control=UNIAXIAL, value=[0.002, 0, 0, 0, 0, 0], ninc=ninc, time=1.0)
    return [Block(steps=[step], ncycle=1).to_dict(T_INIT)]


def test_probe_1_dedicated_kernel_rejects_nstatev():
    """EPICP is a dedicated kernel: its nstatev guard must raise, not crash."""
    with pytest.raises(Exception):
        sim._core.solver_run(_blocks_py(), T_INIT, "EPICP", EPICP_PROPS, 0,
                             record_tangent=False)


def test_probe_2_dedicated_kernel_runs():
    """A full increment through a dedicated kernel, no adapter involved."""
    res = sim._core.solver_run(_blocks_py(ninc=1), T_INIT, "EPICP", EPICP_PROPS, 8,
                               record_tangent=False)
    assert res["status"] == 0


def test_probe_3_adapter_rejects_nstatev():
    """ELISO goes through umat_legacy_modular: this is the call that dies on Windows."""
    with pytest.raises(Exception):
        sim._core.solver_run(_blocks_py(), T_INIT, "ELISO", ELISO_PROPS, 0,
                             record_tangent=False)


def test_probe_4_adapter_runs():
    """A full increment through the adapter: the original crash."""
    res = sim._core.solver_run(_blocks_py(ninc=1), T_INIT, "ELISO", ELISO_PROPS, 1,
                               record_tangent=False)
    assert res["status"] == 0
