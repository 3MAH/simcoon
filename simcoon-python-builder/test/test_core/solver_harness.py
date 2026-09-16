"""Shared solver harness for solver-level tests.

One copy of the uniaxial loading programme and of the history layout used by
test_modular_finite.py and test_solver_robustness.py (not collected by pytest:
no test_ prefix).

The cases used to be written out as a path.txt and run through the file-driven
binding. That binding left with the 2.0 JSON-only migration, so the steps are
built as objects and run in memory; the returned array keeps the column layout
of the old res_global-0.txt, so the constants below still address the same
quantities and the assertions of the callers are untouched.
"""

import numpy as np

import simcoon as sim
from simcoon.solver import Block, StepMeca

# Historical column layout of res_global-0.txt:
# 3 inc | 4 time | 8:14 log strain | 14:20 Kirchhoff stress | 20:29 R |
# 29:38 F | 38:42 Wm, Wm_r, Wm_ir, Wm_d
C_TIME = 4
S_STRAIN = slice(8, 14)
S_STRESS = slice(14, 20)
S_WM = slice(38, 42)

_N_COLUMNS = 42
_T_HOLD = 290.0

# The one-letter codes of the old path grammar.
_CONTROL = {"E": "strain", "S": "stress"}


def run_path(umat_name, props, nstatev, corate, targets, control_type=1):
    """Run a uniaxial case and return its history, laid out like res_global-0.txt.

    ``targets`` is a sequence of (control, value) pairs, one step each: ("S", v)
    drives the 11 stress to v, ("E", v) the 11 strain to v; the lateral components
    stay stress-free.
    """
    steps = [
        StepMeca(
            control=[_CONTROL[ctrl]] + ["stress"] * 5,
            value=[target, 0.0, 0.0, 0.0, 0.0, 0.0],
            time=1.0,
            ninc=50,          # the old "#Dn_inc 0.02"
            Dn_init=1.0,
            Dn_mini=0.001,
        )
        for ctrl, target in targets
    ]

    res = sim.solver.solve(
        Block(steps=steps, control_type=control_type),
        umat_name,
        np.asarray(props, dtype=float),
        nstatev,
        T_init=_T_HOLD,
        corate=corate,
    )

    # The old output.dat asked for strain_type 3 / stress_type 3, i.e. the LOGARITHMIC
    # strain and the KIRCHHOFF stress: "LogStrain" (= "Strain") and "Kirchhoff", not
    # the Cauchy stress "Stress" carries. (At a log strain of 0.15 the Green-Lagrange
    # measure, "GreenLagrange", would read 0.174929.)
    hist = np.zeros((len(res), _N_COLUMNS))
    hist[:, C_TIME] = res["Time"]
    hist[:, S_STRAIN] = res["LogStrain"].T
    hist[:, S_STRESS] = res["Kirchhoff"].T
    hist[:, S_WM] = res["Wm"].T
    return hist


def call_pyumat_batch(etot, Detot, sigma, DR, Wm, nstatev=0, time=0.5, dtime=1.0, **kw):
    """``sim.umat("PYEXT", ...)`` on Fortran-ordered ``(., n)`` batches with empty F0/F1 and
    props (the registered Python law holds its own parameters)."""
    n = etot.shape[1]
    return sim.umat("PYEXT", etot, Detot, np.array([]), np.array([]), sigma, DR,
                    np.zeros((0, 1), order="F"), np.zeros((nstatev, n), order="F"),
                    time, dtime, Wm, **kw)
