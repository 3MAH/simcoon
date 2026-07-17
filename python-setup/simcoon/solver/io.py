"""JSON round-trip for materials and loading paths.

Schema (after the feature/python_solver PR #63 design)::

    material.json: {"name": "ELISO", "props": [...], "nstatev": 1,
                    "orientation": {"psi": 0, "theta": 0, "phi": 0}}

    path.json: {"initial_temperature": 293.15,
                "corate": "logarithmic",
                "blocks": [{"type": "mechanical", "control_type": "small_strain",
                            "ncycle": 1,
                            "steps": [{"mode": "linear", "control": [...],
                                       "value": [...], "time": 1.0, "ninc": 100,
                                       ...}]}]}

Tabular step tables are stored as nested lists.
"""

from __future__ import annotations

import json
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


def _step_to_json(step: StepMeca) -> dict:
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
        d["tabular"] = np.asarray(step.tabular, dtype=float).tolist()
    return d


def _step_from_json(d: dict) -> StepMeca:
    cls = StepThermomeca if d.get("thermomechanical", False) else StepMeca
    kwargs = {k: d[k] for k in _STEP_SCALARS if k in d and hasattr(cls, k)}
    kwargs["control"] = d["control"]
    if "value" in d:
        kwargs["value"] = np.asarray(d["value"], dtype=float)
    if "BC_w" in d:
        kwargs["BC_w"] = np.asarray(d["BC_w"], dtype=float)
    if "tabular" in d:
        kwargs["tabular"] = np.asarray(d["tabular"], dtype=float)
    return cls(**kwargs)


def save_path_json(filename: str, blocks: List[Block], T_init: float = 293.15,
                   corate="logarithmic") -> None:
    """Write a loading path (list of Blocks) to JSON."""
    if isinstance(blocks, Block):
        blocks = [blocks]
    payload = {
        "initial_temperature": float(T_init),
        "corate": corate,
        "blocks": [
            {
                "control_type": b.control_type,
                "ncycle": int(b.ncycle),
                "steps": [_step_to_json(s) for s in b.steps],
            }
            for b in blocks
        ],
    }
    with open(filename, "w") as f:
        json.dump(payload, f, indent=2)


def load_path_json(filename: str) -> Tuple[List[Block], float, object]:
    """Read a loading path; returns (blocks, T_init, corate)."""
    with open(filename) as f:
        payload = json.load(f)
    blocks = [
        Block(
            steps=[_step_from_json(s) for s in bd["steps"]],
            control_type=bd.get("control_type", "small_strain"),
            ncycle=bd.get("ncycle", 1),
        )
        for bd in payload["blocks"]
    ]
    return blocks, float(payload.get("initial_temperature", 293.15)), payload.get("corate", "logarithmic")


def load_simulation_json(material_file: str, path_file: str) -> dict:
    """Merge a material and a path JSON into kwargs for solve()."""
    kwargs = load_material_json(material_file)
    blocks, T_init, corate = load_path_json(path_file)
    kwargs.update({"blocks": blocks, "T_init": T_init, "corate": corate})
    return kwargs
