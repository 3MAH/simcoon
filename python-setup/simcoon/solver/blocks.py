"""Loading definition (blocks and steps) for the in-memory simcoon solver.

Components of strain/stress vectors follow the simcoon Voigt convention
[11, 22, 33, 12, 13, 23] (engineering shear strains). For the kinematic
control types ('F', 'gradU'), the 9 components are the row-major entries of
the deformation gradient (resp. displacement gradient) tensor.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Sequence, Union

import numpy as np

from .maps import BLOCK_TYPES, CONTROL_TYPES, STEP_MODES, THERMAL_CONTROL, as_code

#: Aliases accepted in the per-component `control` list
_MECA_FLAGS = {
    "strain": 0, "E": 0, "e": 0,
    "stress": 1, "S": 1, "s": 1,
    "zero": 2, "0": 2,
    "F": 0, "L": 0,  # kinematic components (control types 'F' / 'gradU')
}


def _control_flags(control, size_meca):
    """Convert the per-component control specification to cBC_meca integer flags."""
    if isinstance(control, str):
        control = [control] * size_meca
    flags = [as_code(c, _MECA_FLAGS, "control flag") for c in control]
    if len(flags) != size_meca:
        raise ValueError(
            f"control must have {size_meca} components; got {len(flags)}"
        )
    return flags


@dataclass
class StepMeca:
    """A mechanical loading step.

    Parameters
    ----------
    control : str or sequence of str
        Per-component control: 'strain' (kinematic) or 'stress' (static);
        'zero' holds a component at zero for tabular steps. A single string
        applies to all components.
    value : array-like, optional
        Absolute target values at the end of the step, per component, in
        Voigt order (6 components; 9 row-major for control types 'F'/'gradU').
        Not used for tabular steps.
    time : float
        Duration of the step.
    ninc : int
        Number of increments (linear and sinusoidal modes).
    mode : str or int
        'linear', 'sinusoidal' or 'tabular'.
    Dn_init, Dn_mini : float, optional
        Initial and minimal sub-increment fraction of one increment
        (adaptive stepping). Defaults: 1.0 and 1e-3.
    BC_w : array-like, optional
        3x3 spin matrix (rotation rate) applied during the step, for the
        mixed finite-strain control types ('green_lagrange', 'logarithmic',
        'biot').
    T_final : float, optional
        Target temperature at the end of the step (temperature ramp applied
        to the mechanical UMAT). None (default) holds the temperature.
    tabular : ndarray, optional
        Mode-3 table, one row per increment. Columns: [time, (thermal column,
        see tabular_T and StepThermomeca.thermal_control), controlled
        components in Voigt order]. Required when mode='tabular'. The time
        column is ABSOLUTE simulation time and must continue from the previous
        step's end time (a table restarting at 0 after an earlier step is
        rejected: it would produce negative time increments). Tabular steps
        cannot be cycled (Block.ncycle must be 1).
    tabular_T : bool
        Whether the tabular table contains a temperature column (after the
        time column). Default False (constant temperature). For
        thermomechanical heat-flux steps the thermal column is the flux Q and
        is always required, regardless of this flag.
    """

    control: Union[str, Sequence[Union[str, int]]] = "strain"
    value: Optional[Sequence[float]] = None
    time: float = 1.0
    ninc: int = 100
    mode: Union[str, int] = "linear"
    Dn_init: float = 1.0
    Dn_mini: float = 1.0e-3
    BC_w: Optional[Sequence[Sequence[float]]] = None
    T_final: Optional[float] = None
    tabular: Optional[np.ndarray] = None
    tabular_T: bool = False

    _thermomechanical = False

    def to_dict(self, control_type: int, T_hold: float) -> dict:
        """Marshal to the dict consumed by simcoon._core.solver_run.

        Parameters
        ----------
        control_type : int
            The control type of the enclosing block (drives the 6/9 sizing).
        T_hold : float
            Temperature to hold when T_final is None (running block value).
        """
        size_meca = 9 if control_type >= 5 else 6
        mode = as_code(self.mode, STEP_MODES, "step mode")
        d = {
            "mode": mode,
            "Dn_init": float(self.Dn_init),
            "Dn_mini": float(self.Dn_mini),
            "cBC_meca": _control_flags(self.control, size_meca),
        }
        if mode < 3:
            if self.ninc < 1:
                raise ValueError("ninc must be >= 1")
            d["Dn_inc"] = 1.0 / int(self.ninc)
            d["time"] = float(self.time)
            if self.value is None:
                raise ValueError("value is required for linear/sinusoidal steps")
            value = np.asarray(self.value, dtype=float).ravel()
            if value.size != size_meca:
                raise ValueError(
                    f"value must have {size_meca} components for this control type; got {value.size}"
                )
            d["BC_meca"] = value
        else:
            if self.tabular is None:
                raise ValueError("tabular steps require the `tabular` table")
            d["tab_data"] = np.ascontiguousarray(self.tabular, dtype=float)
        if self.BC_w is not None:
            d["BC_w"] = np.asarray(self.BC_w, dtype=float).reshape(3, 3)
        self._thermal_to_dict(d, mode, T_hold)
        return d

    def _thermal_to_dict(self, d: dict, mode: int, T_hold: float) -> None:
        if mode == 3:
            d["cBC_T"] = 0 if self.tabular_T else 2
            d["BC_T"] = 0.0
        else:
            d["cBC_T"] = 0
            d["BC_T"] = float(self.T_final) if self.T_final is not None else float(T_hold)

    def T_end(self, T_hold: float) -> float:
        """Temperature at the end of the step (for chaining T_final=None steps)."""
        if self.tabular_T and self.tabular is not None:
            return float(np.asarray(self.tabular)[-1, 1])  # T is the column after time
        return float(self.T_final) if self.T_final is not None else float(T_hold)


@dataclass
class StepThermomeca(StepMeca):
    """A thermomechanical loading step (coupled heat equation).

    In addition to the mechanical control of StepMeca:

    Parameters
    ----------
    thermal_control : str
        'temperature' (ramp to T_final), 'heat_flux' (prescribed flux Q) or
        'convection' (0D convection Q = -q_conv (T - T_init)).
    Q : float
        Prescribed heat flux (thermal_control='heat_flux').
    q_conv : float
        Convection coefficient rho*c_p/tau (thermal_control='convection').
    """

    thermal_control: str = "temperature"
    Q: float = 0.0
    q_conv: float = 0.0

    _thermomechanical = True

    def _thermal_to_dict(self, d: dict, mode: int, T_hold: float) -> None:
        cBC_T = as_code(self.thermal_control, THERMAL_CONTROL, "thermal control")
        if mode == 3:
            if cBC_T == 0:  # temperature: T column when tabular_T, else constant
                d["cBC_T"] = 0 if self.tabular_T else 2
                d["BC_T"] = 0.0
            elif cBC_T == 1:  # heat flux: the table carries the Q column after time
                d["cBC_T"] = 1
                d["BC_T"] = 0.0
            else:  # convection: no thermal column, q_conv as coefficient
                d["cBC_T"] = 3
                d["BC_T"] = float(self.q_conv)
            return
        d["cBC_T"] = cBC_T
        if cBC_T == 0:
            d["BC_T"] = float(self.T_final) if self.T_final is not None else float(T_hold)
        elif cBC_T == 1:
            d["BC_T"] = float(self.Q)
        else:  # convection
            d["BC_T"] = float(self.q_conv)

    def T_end(self, T_hold: float) -> float:
        """Temperature at the end of the step; flux/convection steps leave the
        chained hold temperature unchanged (the reached T is solution-dependent)."""
        if self.thermal_control != "temperature":
            return T_hold
        return super().T_end(T_hold)


@dataclass
class Block:
    """A loading block: a sequence of steps repeated ncycle times.

    Parameters
    ----------
    steps : list of StepMeca / StepThermomeca
        The loading steps. If any step is a StepThermomeca, the block is
        thermomechanical (coupled heat equation) and all steps must be.
    control_type : str or int
        Loading control (see CONTROL_TYPES). Thermomechanical blocks only
        support 'small_strain'.
    ncycle : int
        Number of repetitions of the step sequence.
    """

    steps: List[StepMeca] = field(default_factory=list)
    control_type: Union[str, int] = "small_strain"
    ncycle: int = 1

    def add_step(self, step: StepMeca) -> "Block":
        self.steps.append(step)
        return self

    @property
    def type(self) -> int:
        thermo = [s._thermomechanical for s in self.steps]
        if any(thermo) and not all(thermo):
            raise ValueError(
                "a block cannot mix StepMeca and StepThermomeca steps"
            )
        return BLOCK_TYPES["thermomechanical"] if any(thermo) else BLOCK_TYPES["mechanical"]

    def to_dict(self, T_hold: float) -> dict:
        """Marshal to the dict consumed by simcoon._core.solver_run."""
        if not self.steps:
            raise ValueError("a block requires at least one step")
        control_type = as_code(self.control_type, CONTROL_TYPES, "control type")
        btype = self.type
        if btype == 2 and control_type != 1:
            raise ValueError(
                "thermomechanical blocks only support control_type 'small_strain'"
            )
        steps = []
        T_run = T_hold
        for s in self.steps:
            step_dict = s.to_dict(control_type, T_run)
            if self.ncycle > 1 and step_dict["mode"] == 3:
                # the time column of a table is absolute, so repeating the step is
                # ill-defined (the engine rejects it too); unroll cycles explicitly
                raise ValueError(
                    "tabular steps cannot be cycled (ncycle > 1): the time column is "
                    "absolute; unroll the cycles into explicit steps"
                )
            steps.append(step_dict)
            T_run = s.T_end(T_run)
        return {
            "type": btype,
            "control_type": control_type,
            "ncycle": int(self.ncycle),
            "steps": steps,
        }

    def T_end(self, T_hold: float) -> float:
        T_run = T_hold
        for s in self.steps:
            T_run = s.T_end(T_run)
        return T_run
