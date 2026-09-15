"""JSON round-trip for materials and loading paths.

Schema (after the feature/python_solver PR #63 design)::

    material.json: {"name": "ELISO", "props": [...], "nstatev": 1,
                    "orientation": {"psi": 0, "theta": 0, "phi": 0}}

    path.json: {"initial_temperature": 293.15,
                "corate": "logarithmic_R",
                "blocks": [{"type": "mechanical", "control_type": "small_strain",
                            "ncycle": 1,
                            "steps": [{"mode": "linear", "control": [...],
                                       "value": [...], "time": 1.0, "ninc": 100,
                                       ...}]}]}

Mode-3 tables are the one input that is not JSON: a step's ``"tabular"`` entry
names a CSV file next to the path JSON (``<stem>_tab<k>.csv``, one row per
increment, a ``#`` header naming the columns), written by :func:`save_path_json`
and read back by :func:`load_path_json`. Whitespace-separated ``.txt`` tables are
read the same way; a nested list in place of the filename is still accepted.
"""

from __future__ import annotations

import json
import os
from typing import List, Tuple

import numpy as np

from .blocks import Block, StepMeca, StepThermomeca

_STEP_SCALARS = ("time", "ninc", "mode", "Dn_init", "Dn_mini", "T_final",
                 "thermal_control", "Q", "q_conv", "tabular_T")


def save_material_json(filename: str, umat_name: str, props, nstatev: int,
                       orientation=(0.0, 0.0, 0.0)) -> None:
    """Write a material definition to JSON."""
    psi, theta, phi = (float(x) for x in orientation)
    payload = {
        "name": umat_name,
        "props": np.asarray(props, dtype=float).ravel().tolist(),
        "nstatev": int(nstatev),
        "orientation": {"psi": psi, "theta": theta, "phi": phi},
    }
    with open(filename, "w") as f:
        json.dump(payload, f, indent=2)


def load_material_json(filename: str) -> dict:
    """Read a material definition; returns kwargs for solve()."""
    with open(filename) as f:
        payload = json.load(f)
    ori = payload.get("orientation", {})
    return {
        "umat_name": payload["name"],
        "props": np.asarray(payload["props"], dtype=float),
        "nstatev": int(payload["nstatev"]),
        "orientation": (ori.get("psi", 0.0), ori.get("theta", 0.0), ori.get("phi", 0.0)),
    }


_VOIGT_NAMES = ("11", "22", "33", "12", "13", "23")
_FULL_NAMES = ("11", "12", "13", "21", "22", "23", "31", "32", "33")
_CONTROL_LETTER = {"strain": "E", "stress": "S"}


def _table_columns(step: StepMeca, ncols: int) -> List[str]:
    """Header names of a mode-3 table: time, thermal column, controlled components."""
    names = ["time"]
    if getattr(step, "thermal_control", "temperature") == "heat_flux":
        names.append("Q")
    elif step.tabular_T:
        names.append("T")
    control = step.control
    if isinstance(control, str):
        control = [control] * (ncols - len(names))
    comp = _FULL_NAMES if len(control) == 9 else _VOIGT_NAMES
    for c, ij in zip(control, comp):
        letter = _CONTROL_LETTER.get(str(c).lower()) if isinstance(c, str) else None
        if letter is not None:
            names.append(letter + ij)
    if len(names) != ncols:  # unknown layout: keep the file readable anyway
        names = ["time"] + [f"c{k}" for k in range(1, ncols)]
    return names


def write_table(filename: str, step: StepMeca) -> None:
    """Write the mode-3 table of `step` as CSV (``#`` header, one row per increment)."""
    table = np.asarray(step.tabular, dtype=float)
    if table.ndim != 2:
        raise ValueError("a tabular table must be 2-D: [time, (T/Q), components]")
    header = "# " + ", ".join(_table_columns(step, table.shape[1]))
    with open(filename, "w") as f:
        f.write(header + "\n")
        for row in table:
            f.write(", ".join(repr(float(x)) for x in row) + "\n")


def read_table(filename: str) -> np.ndarray:
    """Read a mode-3 table: comma- or whitespace-separated, ``#`` lines ignored."""
    rows = []
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            rows.append([float(x) for x in (line.split(",") if "," in line else line.split())])
    if not rows:
        raise ValueError(f"{filename}: empty table")
    return np.array(rows, dtype=float)


def _step_to_json(step: StepMeca, table_name=None) -> dict:
    d = {"thermomechanical": step._thermomechanical}
    for k in _STEP_SCALARS:
        if hasattr(step, k):
            d[k] = getattr(step, k)
    d["control"] = list(step.control) if not isinstance(step.control, str) else step.control
    if step.value is not None:
        d["value"] = np.asarray(step.value, dtype=float).ravel().tolist()
    if step.BC_w is not None:
        d["BC_w"] = np.asarray(step.BC_w, dtype=float).reshape(3, 3).tolist()
    if step.tabular is not None:
        d["tabular"] = table_name
    return d


def _step_from_json(d: dict, base_dir: str) -> StepMeca:
    cls = StepThermomeca if d.get("thermomechanical", False) else StepMeca
    kwargs = {k: d[k] for k in _STEP_SCALARS if k in d and hasattr(cls, k)}
    kwargs["control"] = d["control"]
    if "value" in d:
        kwargs["value"] = np.asarray(d["value"], dtype=float)
    if "BC_w" in d:
        kwargs["BC_w"] = np.asarray(d["BC_w"], dtype=float)
    if "tabular" in d:
        tab = d["tabular"]
        if isinstance(tab, str):
            tab = read_table(tab if os.path.isabs(tab) else os.path.join(base_dir, tab))
        kwargs["tabular"] = np.asarray(tab, dtype=float)
    return cls(**kwargs)


def save_path_json(filename: str, blocks: List[Block], T_init: float = 293.15,
                   corate="logarithmic_R") -> None:
    """Write a loading path (list of Blocks) to JSON.

    The table of every tabular step goes to its own CSV next to the JSON,
    ``<stem>_tab<k>.csv`` (k counting the tabular steps of the path from 1), and
    the JSON references it by that name.
    """
    if isinstance(blocks, Block):
        blocks = [blocks]
    filename = str(filename)
    base_dir = os.path.dirname(filename)
    stem = os.path.splitext(os.path.basename(filename))[0]
    ntab = 0
    payload = {"initial_temperature": float(T_init), "corate": corate, "blocks": []}
    for b in blocks:
        steps = []
        for s in b.steps:
            table_name = None
            if s.tabular is not None:
                ntab += 1
                table_name = f"{stem}_tab{ntab}.csv"
                write_table(os.path.join(base_dir, table_name), s)
            steps.append(_step_to_json(s, table_name))
        payload["blocks"].append({"control_type": b.control_type,
                                  "ncycle": int(b.ncycle), "steps": steps})
    with open(filename, "w") as f:
        json.dump(payload, f, indent=2)


def load_path_json(filename: str) -> Tuple[List[Block], float, object]:
    """Read a loading path; returns (blocks, T_init, corate).

    A step's ``"tabular"`` entry is a table filename relative to the JSON (or a
    nested list, read as the table itself).
    """
    filename = str(filename)
    base_dir = os.path.dirname(filename)
    with open(filename) as f:
        payload = json.load(f)
    blocks = [
        Block(
            steps=[_step_from_json(s, base_dir) for s in bd["steps"]],
            control_type=bd.get("control_type", "small_strain"),
            ncycle=bd.get("ncycle", 1),
        )
        for bd in payload["blocks"]
    ]
    return blocks, float(payload.get("initial_temperature", 293.15)), payload.get("corate", "logarithmic_R")


def load_simulation_json(material_file: str, path_file: str) -> dict:
    """Merge a material and a path JSON into kwargs for solve()."""
    kwargs = load_material_json(material_file)
    blocks, T_init, corate = load_path_json(path_file)
    kwargs.update({"blocks": blocks, "T_init": T_init, "corate": corate})
    return kwargs
