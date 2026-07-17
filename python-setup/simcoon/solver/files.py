"""Parsers for the legacy simcoon text formats (path.txt, material.dat).

`from_file` turns a legacy loading path file into Block/Step objects for
:func:`~simcoon.solver.solve`; `material_from_file` reads a legacy material
definition into solve() keyword arguments. Both consume the files token by
token, mirroring the historical C++ readers (read_path / read_matprops in
src/Simulation/Solver/read.cpp), so any file the file-driven solver accepted
parses identically.
"""

from __future__ import annotations

import os
from typing import List, Tuple

import numpy as np

from .blocks import Block, StepMeca, StepThermomeca

#: path-file component order (11, 12, 22, 13, 23, 33) -> Voigt index [11,22,33,12,13,23]
_FILE_TO_VOIGT = [0, 3, 1, 4, 5, 2]

_MECA_FLAG_NAMES = {"E": "strain", "S": "stress", "0": "zero", "L": "strain", "U": "strain", "F": "strain"}


class _Tokens:
    """Sequential token consumer replicating the C++ `stream >> value` semantics."""

    def __init__(self, text: str):
        self._toks = text.split()
        self._i = 0

    def s(self) -> str:
        tok = self._toks[self._i]
        self._i += 1
        return tok

    def f(self) -> float:
        return float(self.s())

    def i(self) -> int:
        return int(self.s())

    def skip(self, n: int = 1) -> None:
        self._i += n


def _read_meca_targets(tk: _Tokens, control_type: int):
    """Read the per-component flags+targets of a linear/sinusoidal step."""
    if control_type <= 4:
        control = [None] * 6
        value = np.zeros(6)
        for k in range(6):
            flag = tk.s()
            control[_FILE_TO_VOIGT[k]] = _MECA_FLAG_NAMES[flag]
            value[_FILE_TO_VOIGT[k]] = tk.f()
        return control, value
    # control types 5/6: 9 raw kinematic components, no flags
    return ["strain"] * 9, np.array([tk.f() for _ in range(9)])


def _read_meca_flags(tk: _Tokens, control_type: int):
    """Read the per-component flags of a tabular step (no targets)."""
    if control_type <= 4:
        control = [None] * 6
        for k in range(6):
            control[_FILE_TO_VOIGT[k]] = _MECA_FLAG_NAMES[tk.s()]
        return control
    return [_MECA_FLAG_NAMES[tk.s()] for _ in range(9)]


def _read_rotation(tk: _Tokens):
    tk.skip()  # rotation label
    return np.array([[tk.f() for _ in range(3)] for _ in range(3)])


def _load_tab_file(path_data: str, filename: str) -> np.ndarray:
    """Load a mode-3 increment file: one row per increment, leading label stripped."""
    rows = []
    with open(os.path.join(path_data, filename)) as f:
        for line in f:
            toks = line.split()
            if toks:
                rows.append([float(x) for x in toks[1:]])
    return np.array(rows)


def _parse_meca_step(tk: _Tokens, control_type: int, path_data: str) -> StepMeca:
    tk.skip()  # '#Mode'
    mode = tk.i()
    if mode in (1, 2):
        tk.skip(); Dn_init = tk.f()
        tk.skip(); Dn_mini = tk.f()
        tk.skip(); Dn_inc = tk.f()
        tk.skip(); time = tk.f()
        tk.skip()  # mechanical-state label
        control, value = _read_meca_targets(tk, control_type)
        BC_w = _read_rotation(tk) if 2 <= control_type <= 4 else None
        tk.skip()  # thermal label
        if tk.s() != "T":
            raise ValueError("mechanical steps only accept a temperature (T) condition")
        T_final = tk.f()
        return StepMeca(control=control, value=value, time=time,
                        ninc=round(1.0 / Dn_inc), mode=mode,
                        Dn_init=Dn_init, Dn_mini=Dn_mini, BC_w=BC_w, T_final=T_final)
    if mode == 3:
        tk.skip(); tabfile = tk.s()
        tk.skip(); Dn_init = tk.f()
        tk.skip(); Dn_mini = tk.f()
        tk.skip()  # mechanical-state label
        control = _read_meca_flags(tk, control_type)
        tk.skip()  # thermal label
        thermal = tk.s()  # 'T' = temperature column in the table, '0' = constant
        return StepMeca(control=control, mode=3, Dn_init=Dn_init, Dn_mini=Dn_mini,
                        tabular=_load_tab_file(path_data, tabfile),
                        tabular_T=(thermal == "T"))
    raise ValueError(f"unknown step mode {mode} (1: linear, 2: sinusoidal, 3: tabular)")


def _parse_thermomeca_step(tk: _Tokens, control_type: int, path_data: str) -> StepThermomeca:
    tk.skip()  # '#Mode'
    mode = tk.i()
    if mode in (1, 2):
        tk.skip(); Dn_init = tk.f()
        tk.skip(); Dn_mini = tk.f()
        tk.skip(); Dn_inc = tk.f()
        tk.skip(); time = tk.f()
        tk.skip()  # mechanical-state label
        control, value = _read_meca_targets(tk, control_type)
        BC_w = _read_rotation(tk) if 2 <= control_type <= 3 else None
        tk.skip()  # thermal label
        thermal = tk.s()
        thermal_value = tk.f()
        kwargs = dict(control=control, value=value, time=time,
                      ninc=round(1.0 / Dn_inc), mode=mode,
                      Dn_init=Dn_init, Dn_mini=Dn_mini, BC_w=BC_w)
        if thermal == "T":
            return StepThermomeca(T_final=thermal_value, **kwargs)
        if thermal == "Q":
            return StepThermomeca(thermal_control="heat_flux", Q=thermal_value, **kwargs)
        if thermal == "C":
            return StepThermomeca(thermal_control="convection", q_conv=thermal_value, **kwargs)
        raise ValueError(f"unknown thermal condition '{thermal}' (T, Q or C)")
    if mode == 3:
        tk.skip(); tabfile = tk.s()
        tk.skip(); Dn_init = tk.f()
        tk.skip(); Dn_mini = tk.f()
        tk.skip()  # mechanical-state label
        control = _read_meca_flags(tk, control_type)
        tk.skip()  # thermal label
        thermal = tk.s()  # 'T'/'Q' column in the table, '0' constant, 'C' convection
        kwargs = dict(control=control, mode=3, Dn_init=Dn_init, Dn_mini=Dn_mini,
                      tabular=_load_tab_file(path_data, tabfile))
        if thermal == "T":
            return StepThermomeca(tabular_T=True, **kwargs)
        if thermal == "0":
            return StepThermomeca(**kwargs)
        if thermal == "Q":
            return StepThermomeca(thermal_control="heat_flux", **kwargs)
        if thermal == "C":
            return StepThermomeca(thermal_control="convection", q_conv=tk.f(), **kwargs)
        raise ValueError(f"unknown tabular thermal condition '{thermal}' (T, Q, C or 0)")
    raise ValueError(f"unknown step mode {mode} (1: linear, 2: sinusoidal, 3: tabular)")


def from_file(path_data: str = "data", pathfile: str = "path.txt") -> Tuple[List[Block], float]:
    """Parse a legacy loading path file into Block objects.

    Parameters
    ----------
    path_data : str
        Folder containing the path file (and any mode-3 increment files it
        references).
    pathfile : str
        Name of the loading path file.

    Returns
    -------
    blocks : list of Block
        The loading path, ready for :func:`~simcoon.solver.solve`.
    T_init : float
        The initial temperature declared in the file.
    """
    with open(os.path.join(path_data, pathfile)) as f:
        tk = _Tokens(f.read())

    tk.skip(); T_init = tk.f()
    tk.skip(); nblocks = tk.i()

    blocks = []
    for _ in range(nblocks):
        tk.skip(); number = tk.i()
        tk.skip(); btype = tk.i()
        tk.skip(); control_type = tk.i()
        tk.skip(); ncycle = tk.i()
        tk.skip(); nstep = tk.i()
        if btype not in (1, 2):
            raise ValueError(f"block {number}: unknown loading type {btype} (1: mechanical, 2: thermomechanical)")
        parse_step = _parse_meca_step if btype == 1 else _parse_thermomeca_step
        steps = [parse_step(tk, control_type, path_data) for _ in range(nstep)]
        blocks.append(Block(steps=steps, control_type=control_type, ncycle=ncycle))
    return blocks, T_init


def material_from_file(path_data: str = "data", materialfile: str = "material.dat") -> dict:
    """Parse a legacy material definition file into solve() keyword arguments.

    Returns
    -------
    dict
        {'umat_name', 'props', 'nstatev', 'orientation'} — mergeable into
        :func:`~simcoon.solver.solve` like the JSON loader.
    """
    with open(os.path.join(path_data, materialfile)) as f:
        tk = _Tokens(f.read())

    tk.skip(2); umat_name = tk.s()
    tk.skip(); nprops = tk.i()
    tk.skip(); nstatev = tk.i()
    tk.skip(2); psi = tk.f()
    tk.skip(); theta = tk.f()
    tk.skip(); phi = tk.f()
    tk.skip()  # section label of the properties
    props = np.zeros(nprops)
    for i in range(nprops):
        tk.skip()
        props[i] = tk.f()
    return {
        "umat_name": umat_name,
        "props": props,
        "nstatev": int(nstatev),
        "orientation": (psi, theta, phi),
    }
