"""Results container for the in-memory simcoon solver.

The layout follows the fedoo DataSet conventions: data is accessed by name
(``res["Stress"]``), arrays are components-first (a stress history is a
``(6, N)`` array in Voigt order [11, 22, 33, 12, 13, 23], a deformation
gradient history is ``(3, 3, N)``), and scalar histories are ``(N,)``.
``fedoo.util.voigt_tensors.StressTensorList(res["Stress"])`` therefore works
directly on the returned arrays.
"""

from __future__ import annotations

from typing import Dict

import numpy as np


class SolverResults:
    """Dict-like access to the solver history.

    Attributes
    ----------
    scalar_data : dict
        Scalar histories, shape (N,): 'Time', 'Temp', 'Block', 'Cycle',
        'Step', 'Inc'; thermomechanical runs add 'Q' (heat flux) and 'r'
        (heat source).
    field_data : dict
        Tensor histories, components-first: 'Stress' (Cauchy, (6, N)),
        'Kirchhoff', 'PKII', 'Strain' (Green-Lagrange), 'LogStrain',
        'Statev' ((nstatev, N)), 'Wm' ((4, N)), 'F', 'R', 'DR' ((3, 3, N));
        'TangentMatrix' ((6, 6, N)) for mechanical runs; thermomechanical
        runs add 'Wt' ((3, N)) and the coupled tangents 'dSdE' ((6, 6, N)),
        'dSdT' ((6, N)), 'drdE' ((6, N)), 'drdT' ((N,)).
    status : int
        0 if the simulation ran to completion, 1 on early abort (the
        recorded history is then partial).
    """

    def __init__(self, raw: Dict[str, np.ndarray]):
        self.status = int(raw.get("status", 0))
        self.sv_type = int(raw.get("sv_type", 1))

        n = raw["time"].shape[0]
        self.scalar_data = {
            "Time": raw["time"],
            "Temp": raw["T"],
            "Block": raw["block"],
            "Cycle": raw["cycle"],
            "Step": raw["step"],
            "Inc": raw["inc"],
        }
        self.field_data = {
            "Stress": raw["sigma"].T,
            "Kirchhoff": raw["tau"].T,
            "PKII": raw["PKII"].T,
            "Strain": raw["Etot"].T,
            "LogStrain": raw["etot"].T,
            "Statev": raw["statev"].T,
            "Wm": raw["Wm"].T,
            "F": raw["F1"].reshape(n, 3, 3).transpose(1, 2, 0),
            "R": raw["R"].reshape(n, 3, 3).transpose(1, 2, 0),
            "DR": raw["DR"].reshape(n, 3, 3).transpose(1, 2, 0),
        }
        if "Lt" in raw:
            self.field_data["TangentMatrix"] = raw["Lt"].reshape(n, 6, 6).transpose(1, 2, 0)
        if self.sv_type == 2:
            self.scalar_data["Q"] = raw["Q"]
            self.scalar_data["r"] = raw["r"]
            self.field_data["Wt"] = raw["Wt"].T
            if "dSdE" in raw:
                self.field_data["dSdE"] = raw["dSdE"].reshape(n, 6, 6).transpose(1, 2, 0)
                self.field_data["dSdT"] = raw["dSdT"].T
                self.field_data["drdE"] = raw["drdE"].T
                self.scalar_data["drdT"] = raw["drdT"].ravel()

    # -- dict-like interface -------------------------------------------------
    def __getitem__(self, key: str) -> np.ndarray:
        if key in self.field_data:
            return self.field_data[key]
        if key in self.scalar_data:
            return self.scalar_data[key]
        raise KeyError(
            f"'{key}' not in results; available: {sorted(self.keys())}"
        )

    def __contains__(self, key: str) -> bool:
        return key in self.field_data or key in self.scalar_data

    def keys(self):
        return list(self.scalar_data) + list(self.field_data)

    def get_data(self, key: str) -> np.ndarray:
        """fedoo-style accessor (alias of __getitem__)."""
        return self[key]

    def __len__(self) -> int:
        return self.scalar_data["Time"].shape[0]

    def __repr__(self) -> str:
        kind = "thermomechanical" if self.sv_type == 2 else "mechanical"
        return (
            f"SolverResults({kind}, {len(self)} increments, "
            f"status={self.status}, fields={sorted(self.keys())})"
        )

    # -- persistence ----------------------------------------------------------
    def save(self, filename: str) -> None:
        """Save all histories to a compressed npz archive."""
        payload = {"status": np.array(self.status), "sv_type": np.array(self.sv_type)}
        for k, v in self.scalar_data.items():
            payload[f"scalar__{k}"] = v
        for k, v in self.field_data.items():
            payload[f"field__{k}"] = v
        np.savez_compressed(filename, **payload)

    @classmethod
    def load(cls, filename: str) -> "SolverResults":
        """Load a SolverResults previously written by save()."""
        data = np.load(filename)
        obj = cls.__new__(cls)
        obj.status = int(data["status"])
        obj.sv_type = int(data["sv_type"])
        obj.scalar_data = {}
        obj.field_data = {}
        for k in data.files:
            if k.startswith("scalar__"):
                obj.scalar_data[k[len("scalar__"):]] = data[k]
            elif k.startswith("field__"):
                obj.field_data[k[len("field__"):]] = data[k]
        return obj

    def to_dataframe(self):
        """Flatten scalar and 6-component histories to a pandas DataFrame."""
        import pandas as pd

        cols = {}
        for k, v in self.scalar_data.items():
            cols[k] = v
        comp = ["11", "22", "33", "12", "13", "23"]
        for k, v in self.field_data.items():
            if v.ndim == 2 and v.shape[0] == 6:
                for c in range(6):
                    cols[f"{k}_{comp[c]}"] = v[c]
            elif v.ndim == 2:
                for c in range(v.shape[0]):
                    cols[f"{k}_{c}"] = v[c]
        return pd.DataFrame(cols)
