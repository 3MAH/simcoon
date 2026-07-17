"""solve(): drive the C++ material-point solver in memory."""

from __future__ import annotations

from typing import List, Optional, Sequence, Union

import numpy as np

import simcoon._core as _core

from .blocks import Block, StepMeca
from .maps import CORATE_TYPES, TANGENT_MODES, as_code, tangent_default
from .results import SolverResults


def solve(
    blocks: Union[Block, StepMeca, Sequence[Union[Block, StepMeca]]],
    umat_name: str,
    props: Sequence[float],
    nstatev: int,
    T_init: float = 293.15,
    corate: Union[str, int] = "logarithmic",
    tangent_mode: Union[str, int] = tangent_default,
    solver_type: int = 0,
    orientation: Sequence[float] = (0.0, 0.0, 0.0),
    record_tangent: bool = True,
    raise_on_abort: bool = True,
    **params,
) -> SolverResults:
    """Solve a homogeneous loading path with the C++ simcoon solver, in memory.

    Parameters
    ----------
    blocks : Block, StepMeca or sequence of them
        The loading path. Bare steps are wrapped in a small-strain Block.
    umat_name : str
        Constitutive model name (5 characters, e.g. 'ELISO', 'EPICP', 'MODUL').
    props : array-like
        Material properties.
    nstatev : int
        Number of internal state variables.
    T_init : float
        Initial temperature.
    corate : str or int
        Objective rate for the finite-strain control types (see CORATE_TYPES).
    tangent_mode : str or int
        Tangent operator mode: 'none', 'continuum', 'algorithmic' (default)
        or 'closest_point' (reserved).
    solver_type : int
        0 = classic Newton-Raphson (default), 1 = RNL (control_type 1 only).
    orientation : sequence of 3 floats
        Euler angles (psi, theta, phi) of the material orientation (rad).
    record_tangent : bool
        Capture the tangent operator history ('TangentMatrix' or the coupled
        thermomechanical tangents).
    raise_on_abort : bool
        Raise a RuntimeError when the solver aborts early (status != 0)
        instead of returning the partial history.
    **params
        Numeric solver controls forwarded to the C++ loop: div_tnew_dt,
        mul_tnew_dt, miniter, maxiter, inforce, precision, lambda_solver
        (penalty stiffness of the strain-driven components).

    Returns
    -------
    SolverResults
        History of the converged increments (fedoo-style data layout).
    """
    if isinstance(blocks, (Block, StepMeca)):
        blocks = [blocks]
    blocks = [b if isinstance(b, Block) else Block(steps=[b]) for b in blocks]

    corate_code = as_code(corate, CORATE_TYPES, "corate")
    run_params = dict(params)
    run_params["tangent_mode"] = as_code(tangent_mode, TANGENT_MODES, "tangent mode")

    blocks_py = []
    T_run = float(T_init)
    for b in blocks:
        blocks_py.append(b.to_dict(T_run))
        T_run = b.T_end(T_run)

    psi, theta, phi = (float(x) for x in orientation)
    raw = _core.solver_run(
        blocks_py,
        float(T_init),
        umat_name,
        np.asarray(props, dtype=float).ravel(),
        int(nstatev),
        psi, theta, phi,
        int(solver_type),
        corate_code,
        run_params,
        bool(record_tangent),
    )
    res = SolverResults(raw)
    if raise_on_abort and res.status != 0:
        raise RuntimeError(
            f"the solver aborted early after {len(res)} recorded increments "
            "(status=1); pass raise_on_abort=False to inspect the partial history"
        )
    return res
