"""solve(): drive the C++ material-point solver in memory."""

from __future__ import annotations

import contextlib
from typing import Any, Callable, Optional, Sequence, Union

import numpy as np

import simcoon._core as _core
from simcoon.pyumat import UMAT_NAME, registered

from .blocks import Block, StepMeca, StepThermomeca
from .maps import TANGENT_MODES, as_code, kinematic_code, tangent_default
from .micromechanics import to_phase_dicts, euler_angles
from .results import SolverResults


def solve(
    blocks: Union[Block, StepMeca, Sequence[Union[Block, StepMeca]]],
    umat_name: Union[str, Any, Callable],
    props: Optional[Sequence[float]] = None,
    nstatev: Optional[int] = None,
    T_init: float = 293.15,
    corate: Union[str, int, None] = None,
    tangent_mode: Union[str, int] = tangent_default,
    solver_type: int = 0,
    orientation: Sequence[float] = (0.0, 0.0, 0.0),
    phases: Optional[Sequence[Any]] = None,
    record_tangent: bool = True,
    raise_on_abort: bool = True,
    kinematic_form: Union[str, int, None] = None,
    **params,
) -> SolverResults:
    """Solve a homogeneous loading path with the C++ simcoon solver, in memory.

    Parameters
    ----------
    blocks : Block, StepMeca or sequence of them
        The loading path. Bare steps are wrapped in a small-strain Block.
    umat_name : str, PythonUMAT or callable
        Constitutive model: either the name of a built-in model (5 characters,
        e.g. 'ELISO', 'EPICP', 'MODUL') or a constitutive law written in Python
        (a :class:`simcoon.PythonUMAT` instance, or any callable with its
        ``integrate`` keyword signature). A Python law is registered under the
        ``PYEXT`` name for the duration of the call and integrated by the C++
        solver exactly like a built-in kernel.
    props : array-like, optional
        Material properties. Required for a built-in model; defaults to the
        ``props`` attribute of a Python law.
    nstatev : int, optional
        Number of internal state variables. Required for a built-in model;
        defaults to the ``nstatev`` attribute of a Python law.
    T_init : float
        Initial temperature.
    corate : str, int or None
        Former name of ``kinematic_form``: the same value under either name.
    kinematic_form : str, int or None
        Kinematic form of the finite-strain constitutive update, one value:

        - ``"total"``: the strain-space framework of Miehe, Apel and Lambrecht
          (2002); the law works on the material logarithmic strain ln U, a
          state function of F, and on the stress T dual to it. No rate is
          involved; log-strain box laws and hyperelastic laws built on F; no
          ``"logarithmic"`` nor ``"small_strain"`` control.
        - an objective rate (see CORATE_TYPES, name or code): the rate form
          with that rate; the law works on ln V (logarithmic rates) and the
          Kirchhoff stress tau.
        - ``"rate"``, or None (default): the rate form with 'logarithmic_R' —
          the exact polar rotation, whose frame increment DR = R1 R0^T makes
          the tangent transport exact even with rotated internal-variable
          history (the XBM 'logarithmic' rate keeps a small tangent residual
          there). Under the rate form MODUL requires log_R (the modular Hencky
          composition).
    tangent_mode : str or int
        Tangent operator mode: 'none', 'continuum', 'algorithmic' (default)
        or 'closest_point' (closest-point projection integrator with its exact
        operator where a law has one; the algorithmic operator elsewhere).
    solver_type : int
        0 = classic Newton-Raphson (default), 1 = RNL (control_type 1 only).
    orientation : simcoon.Rotation, dict or sequence of 3 floats
        Orientation of the material frame: a ``Rotation``, or its Euler angles
        ``(psi, theta, phi)`` in degrees (see
        :func:`~simcoon.solver.micromechanics.as_rotation`); applied actively,
        material frame to global frame.
    phases : sequence, optional
        Sub-phases of a mean-field model (MIMTN, MISCN, MIHEN, MIPLN): the
        Ellipsoid / Layer objects of :mod:`simcoon.solver.micromechanics`, or the
        dicts they convert to. Leave it None for every single-phase model.
    record_tangent : bool
        Capture the tangent operator history ('TangentMatrix' or the coupled
        thermomechanical tangents).
    raise_on_abort : bool
        Raise a RuntimeError when the solver aborts early (status != 0)
        instead of returning the partial history. The solver aborts when the
        Newton loop does not converge at the minimal increment, or when the
        increment falls below ``Dn_mini`` with ``inforce=0``.
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

    if isinstance(umat_name, str):
        if props is None or nstatev is None:
            raise TypeError(
                "props and nstatev are required when umat_name is a built-in model name"
            )
        law_ctx = contextlib.nullcontext()
    else:
        # Python law: served under PYEXT for the duration of the solve. The thermomechanical
        # dispatch (select_umat_T) has no PYEXT entry, so reject those blocks here with a clear
        # message instead of a C++ "Unknown umat name" from deep inside the solve.
        for b in blocks:
            if any(isinstance(st, StepThermomeca) for st in b.steps):
                raise TypeError(
                    "a constitutive law written in Python (PYEXT) cannot serve a thermomechanical "
                    "block: only the mechanical dispatch supports it. Use a built-in "
                    "thermomechanical UMAT name, or drive the coupling from Python."
                )
        if props is None:
            props = getattr(umat_name, "props", np.zeros(0))
        if nstatev is None:
            nstatev = getattr(umat_name, "nstatev", 0)
        law_ctx = registered(umat_name)
        umat_name = UMAT_NAME

    corate_type = kinematic_code(kinematic_form, corate)
    run_params = dict(params)
    run_params["tangent_mode"] = as_code(tangent_mode, TANGENT_MODES, "tangent mode")

    blocks_py = []
    T_run = float(T_init)
    for b in blocks:
        blocks_py.append(b.to_dict(T_run))
        T_run = b.T_end(T_run)

    # Sub-phases of a mean-field model (MIMTN, MISCN, MIHEN, MIPLN): dataclasses are accepted
    # as readily as the dicts the binding reads.
    phases_py = None if phases is None else to_phase_dicts(phases)

    with law_ctx:
        raw = _core.solver_run(
            blocks_py,
            float(T_init),
            umat_name,
            np.asarray(props, dtype=float).ravel(),
            int(nstatev),
            euler_angles(orientation),
            int(solver_type),
            corate_type,
            run_params,
            bool(record_tangent),
            phases_py,
        )
    res = SolverResults(raw, [bd["control_type"] > 1 for bd in blocks_py], corate_type)
    if raise_on_abort and res.status != 0:
        raise RuntimeError(
            f"the solver aborted early after {len(res)} recorded increments "
            "(status=1); pass raise_on_abort=False to inspect the partial history"
        )
    return res
