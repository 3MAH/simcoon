#!/usr/bin/env python
"""Convert a pre-2.0 simcoon ``data`` directory to the JSON files simcoon 2.0 reads.

simcoon itself reads no legacy text format any more: ``path.txt``, ``material.dat``,
``tab_file_<n>.txt``, ``N<kind><n>.dat`` and ``solver_essentials.inp`` are parsed HERE,
once, and written out as ``path.json`` (+ one ``<stem>_tab<k>.csv`` per mode-3 table),
``material.json`` and ``<kind><n>.json`` (see
``simcoon.solver.load_simulation_json`` and ``simcoon.solver.micromechanics``). This script
is a migration tool kept outside the package on purpose: the token-based parsers below are
as fragile as the formats they read.

Usage::

    python scripts/legacy_to_json.py DATA_DIR [DATA_DIR ...] [--out OUT_DIR]

Requires an installed simcoon (the JSON writers and the loading objects come from it).
"""
from __future__ import annotations

import glob
import os
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import numpy as np

from simcoon.solver.blocks import Block, StepMeca, StepThermomeca
from simcoon.solver.io import save_material_json, save_path_json
from simcoon.solver.maps import CORATE_TYPES
from simcoon.solver.micromechanics import (Cylinder, Ellipsoid,
                                           Layer, Phase, Section,
                                           save_cylinders_json, save_ellipsoids_json,
                                           save_layers_json, save_phases_json,
                                           save_sections_json)

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
    _reanchor_tabular_times(blocks)
    return blocks, T_init


def _reanchor_tabular_times(blocks) -> None:
    """Shift restarting mode-3 time axes so from_file matches the C++ file path.

    The legacy FILE convention lets each mode-3 table restart its time axis at
    0 (first row = anchor at the current state); the C++ file reader shifts
    such an axis to continue from the running simulation time. In-memory
    ``tabular`` arrays are bound to the strict absolute-time contract instead,
    so a parsed legacy table must be re-anchored here — otherwise a file set
    that runs fine through the file-driven solver is rejected by solve().
    Tables whose time column already continues absolutely are left untouched
    (shift only when the first row is EARLIER than the running time, the same
    guard as the C++ reader).
    """
    t_run = 0.0
    for b in blocks:
        for _ in range(b.ncycle):
            for s in b.steps:
                if s.mode in ("tabular", 3) and s.tabular is not None:
                    t0 = float(s.tabular[0, 0])
                    if t0 < t_run - 1e-12:
                        s.tabular = s.tabular.copy()
                        s.tabular[:, 0] += t_run - t0
                    t_run = float(s.tabular[-1, 0])
                else:
                    t_run += float(s.time)


def material_from_file(path_data: str = "data", materialfile: str = "material.dat") -> dict:
    """Parse a legacy material definition file into solve() keyword arguments.

    Returns
    -------
    dict
        {'umat_name', 'props', 'nstatev', 'orientation'} — mergeable into
        :func:`~simcoon.solver.solve` like the JSON loader. The orientation is the
        file's Euler angles, in degrees, as solve() takes them.
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


# ---------------------------------------------------------------------------
# One-way conversion of a legacy data directory to JSON
# ---------------------------------------------------------------------------

def _corate_of_essentials(path_data: str, essentials: str):
    """The corate named by a solver_essentials.inp, or None when there is no such file.

    Nothing else in that file survives the migration: the solver type is an argument of
    solve(), and solver_control.inp only ever held the defaults of solve(params=...).
    """
    filename = os.path.join(path_data, essentials)
    if not os.path.isfile(filename):
        return None
    with open(filename) as f:
        tk = _Tokens(f.read())
    tk.skip(); tk.i()            # solver type
    tk.skip(); code = tk.i()     # Rate_type
    for name, value in CORATE_TYPES.items():
        if value == code:
            return name
    raise ValueError(f"{filename}: unknown Rate_type {code}")


def _is_path_file(filename: str) -> bool:
    """A loading path file opens with #Initial_temperature; tables and data do not."""
    try:
        with open(filename) as f:
            for line in f:
                if line.strip():
                    return line.strip().startswith("#Initial_temperature")
    except (OSError, UnicodeDecodeError):
        return False
    return False


def convert_to_json(path_data: str = "data", out_dir: Optional[str] = None,
                    materialfile: str = "material.dat",
                    essentials: str = "solver_essentials.inp") -> List[str]:
    """Convert a legacy ``data`` directory to the JSON files simcoon 2.0 reads.

    Every loading path file of the directory (a text file opening with
    ``#Initial_temperature``, whatever its name) becomes ``<stem>.json`` through
    :func:`save_path_json`, the tables its mode-3 steps referenced rewritten as
    ``<stem>_tab<k>.csv`` next to it (time column re-anchored to absolute time) and
    the corate of ``solver_essentials.inp`` recorded; ``material.dat`` becomes
    ``material.json``; every ``N<kind><n>.dat`` sub-phase file becomes ``<kind><n>.json``
    (see :func:`~simcoon.solver.micromechanics.convert_dat_to_json`). ``output.dat``
    and ``solver_control.inp`` have no JSON counterpart: the results come back in
    memory and the control parameters are the defaults of :func:`solve`.

    The legacy files are left in place. Returns the paths written.
    """
    out_dir = path_data if out_dir is None else out_dir
    os.makedirs(out_dir, exist_ok=True)
    corate = _corate_of_essentials(path_data, essentials) or "logarithmic_R"
    written = []

    for name in sorted(os.listdir(path_data)):
        src = os.path.join(path_data, name)
        stem, ext = os.path.splitext(name)
        if ext == ".txt" and _is_path_file(src):
            blocks, T_init = from_file(path_data, name)
            dst = os.path.join(out_dir, stem + ".json")
            save_path_json(dst, blocks, T_init, corate)
            written.append(dst)
            written += sorted(glob.glob(os.path.join(out_dir, stem + "_tab*.csv")))
        elif ext == ".dat" and name == materialfile:
            material = material_from_file(path_data, name)
            dst = os.path.join(out_dir, "material.json")
            save_material_json(dst, material["umat_name"], material["props"],
                               material["nstatev"], material["orientation"])
            written.append(dst)
        elif ext == ".dat":
            try:
                kind_from_dat_name(name)
            except ValueError:
                continue     # not a sub-phase file (raw data, an identification template)
            json_name = (stem[1:] if stem[:1].lower() == "n" else stem).lower() + ".json"
            try:
                written.append(str(convert_dat_to_json(src, os.path.join(out_dir, json_name))))
            except ValueError as exc:
                #an identification template or an Abaqus deck under a sub-phase file name
                warnings.warn(f"{src}: not converted ({exc})", UserWarning, stacklevel=2)
    return written
# =============================================================================
# Legacy .dat input (read-only)
# =============================================================================
#
# The historical tab-separated files (Nphases0.dat, Nlayers0.dat,
# Nellipsoids0.dat, Ncylinders0.dat, Nsections0.dat) used to be parsed in C++ by
# src/Simulation/Phase/read.cpp. They are read here instead, the way path.txt and
# material.dat are read by solver/files.py, so the C++ side never touches the
# filesystem. Reading only: JSON is the format written from now on, and
# convert_dat_to_json() is the one-way door.
#
# Every row holds the fixed columns of its kind, then nprops, nstatev, then the
# nprops property values. Parsing is token-based because the files are aligned
# with ragged tabs, so column positions cannot be trusted.

# Number of fixed columns per kind, up to and including nstatev.
_DAT_LAYOUTS = {
    'phases': 9,
    'layers': 12,
    'ellipsoids': 16,
    'cylinders': 15,
    'sections': 8,
}


def _dat_rows(filepath: Union[str, Path], kind: str) -> List[List[str]]:
    """Tokenise a legacy .dat file, dropping its header line and blank lines.

    Raises
    ------
    ValueError
        If a row does not hold exactly the columns its own nprops announces.
    """
    n_fixed = _DAT_LAYOUTS[kind]
    with open(filepath) as f:
        lines = f.readlines()
    if not lines:
        raise ValueError(f"{filepath}: empty file, a header line was expected")

    rows = []
    for lineno, line in enumerate(lines[1:], start=2):
        tokens = line.split()
        if not tokens:
            continue
        if tokens[0].startswith('*'):
            raise ValueError(
                f"{filepath}:{lineno}: this is an Abaqus deck (*Material / *Solid Section), "
                f"not a tabular {kind} file. The C++ side never read those either"
            )
        if len(tokens) < n_fixed:
            raise ValueError(
                f"{filepath}:{lineno}: {len(tokens)} columns in a {kind} row, "
                f"at least {n_fixed} expected"
            )
        try:
            nprops = int(tokens[n_fixed - 2])
        except ValueError:
            raise ValueError(
                f"{filepath}:{lineno}: nprops column is {tokens[n_fixed - 2]!r}, not an integer"
            ) from None
        if len(tokens) != n_fixed + nprops:
            raise ValueError(
                f"{filepath}:{lineno}: nprops={nprops} announces {n_fixed + nprops} "
                f"columns, {len(tokens)} found"
            )
        rows.append(tokens)
    return rows


def _dat_props(tokens: List[str], n_fixed: int) -> np.ndarray:
    """The property values of a row, which follow the fixed columns."""
    return np.array([float(t) for t in tokens[n_fixed:]], dtype=float)


def load_phases_dat(filepath: Union[str, Path]) -> List[Phase]:
    """Load phases from a legacy ``Nphases<N>.dat`` file."""
    n = _DAT_LAYOUTS['phases']
    return [
        Phase(
            number=int(t[0]),
            umat_name=t[1],
            save=int(t[2]),
            concentration=float(t[3]),
            material_orientation=(float(t[4]), float(t[5]), float(t[6])),
            nstatev=int(t[8]),
            props=_dat_props(t, n),
        )
        for t in _dat_rows(filepath, 'phases')
    ]


def load_layers_dat(filepath: Union[str, Path]) -> List[Layer]:
    """Load layers from a legacy ``Nlayers<N>.dat`` file."""
    n = _DAT_LAYOUTS['layers']
    return [
        Layer(
            number=int(t[0]),
            umat_name=t[1],
            save=int(t[2]),
            concentration=float(t[3]),
            material_orientation=(float(t[4]), float(t[5]), float(t[6])),
            geometry_orientation=(float(t[7]), float(t[8]), float(t[9])),
            nstatev=int(t[11]),
            props=_dat_props(t, n),
        )
        for t in _dat_rows(filepath, 'layers')
    ]


def load_ellipsoids_dat(filepath: Union[str, Path]) -> List[Ellipsoid]:
    """Load ellipsoidal inclusions from a legacy ``Nellipsoids<N>.dat`` file."""
    n = _DAT_LAYOUTS['ellipsoids']
    return [
        Ellipsoid(
            number=int(t[0]),
            coatingof=int(t[1]),
            umat_name=t[2],
            save=int(t[3]),
            concentration=float(t[4]),
            material_orientation=(float(t[5]), float(t[6]), float(t[7])),
            a1=float(t[8]),
            a2=float(t[9]),
            a3=float(t[10]),
            geometry_orientation=(float(t[11]), float(t[12]), float(t[13])),
            nstatev=int(t[15]),
            props=_dat_props(t, n),
        )
        for t in _dat_rows(filepath, 'ellipsoids')
    ]


def load_cylinders_dat(filepath: Union[str, Path]) -> List[Cylinder]:
    """Load cylindrical inclusions from a legacy ``Ncylinders<N>.dat`` file."""
    n = _DAT_LAYOUTS['cylinders']
    return [
        Cylinder(
            number=int(t[0]),
            coatingof=int(t[1]),
            umat_name=t[2],
            save=int(t[3]),
            concentration=float(t[4]),
            material_orientation=(float(t[5]), float(t[6]), float(t[7])),
            L=float(t[8]),
            R=float(t[9]),
            geometry_orientation=(float(t[10]), float(t[11]), float(t[12])),
            nstatev=int(t[14]),
            props=_dat_props(t, n),
        )
        for t in _dat_rows(filepath, 'cylinders')
    ]


def load_sections_dat(filepath: Union[str, Path]) -> List[Section]:
    """Load textile sections from a legacy ``Nsections<N>.dat`` file.

    Only the tabular form is read. ``Nsections1.dat``-style Abaqus decks
    (``*Material`` / ``*Solid Section``) were never read by the C++ side either.
    """
    n = _DAT_LAYOUTS['sections']
    return [
        Section(
            number=int(t[0]),
            name=t[1],
            umat_name=t[2],
            material_orientation=(float(t[3]), float(t[4]), float(t[5])),
            nstatev=int(t[7]),
            props=_dat_props(t, n),
        )
        for t in _dat_rows(filepath, 'sections')
    ]


_DAT_CONVERTERS = {
    'phases': (load_phases_dat, save_phases_json),
    'layers': (load_layers_dat, save_layers_json),
    'ellipsoids': (load_ellipsoids_dat, save_ellipsoids_json),
    'cylinders': (load_cylinders_dat, save_cylinders_json),
    'sections': (load_sections_dat, save_sections_json),
}


def kind_from_dat_name(name: str) -> str:
    """The kind ('phases', 'ellipsoids', ...) a legacy file name announces."""
    stem = Path(name).stem.lower()
    for kind in _DAT_LAYOUTS:
        if stem.startswith('n' + kind):
            return kind
    raise ValueError(
        f"{name!r}: cannot tell which kind this is, expected a name starting with "
        + ", ".join('N' + k for k in _DAT_LAYOUTS)
    )


def convert_dat_to_json(filepath: Union[str, Path],
                        json_path: Union[str, Path] = None,
                        kind: str = None,
                        prop_names: List[str] = None) -> Path:
    """Convert one legacy .dat file to its JSON equivalent.

    Parameters
    ----------
    filepath : str or Path
        The ``N<kind><N>.dat`` file to read.
    json_path : str or Path, optional
        Where to write. Defaults to the same directory, with the leading ``N``
        dropped and lowercased: ``Nellipsoids0.dat`` -> ``ellipsoids0.json``.
    kind : str, optional
        Override the kind instead of inferring it from the file name.
    prop_names : list of str, optional
        Names for the property columns, stored in the JSON as a dict.

    Returns
    -------
    Path
        The JSON file written.
    """
    filepath = Path(filepath)
    kind = kind or kind_from_dat_name(filepath.name)
    if kind not in _DAT_CONVERTERS:
        raise ValueError(f"unknown kind {kind!r}, expected one of {', '.join(_DAT_CONVERTERS)}")
    if json_path is None:
        stem = filepath.stem
        json_path = filepath.with_name((stem[1:] if stem[:1].lower() == 'n' else stem).lower() + '.json')

    load, save = _DAT_CONVERTERS[kind]
    save(json_path, load(filepath), prop_names=prop_names)
    return Path(json_path)




# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------

def main(argv=None) -> int:
    import argparse
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("data_dirs", nargs="+", help="legacy data directories to convert")
    parser.add_argument("--out", default=None, help="write the JSON files there instead of next to the sources")
    args = parser.parse_args(argv)
    for d in args.data_dirs:
        for written in convert_to_json(d, args.out):
            print(written)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
