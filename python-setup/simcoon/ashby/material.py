# -*- coding: utf-8 -*-

"""Material data model for the simcoon Ashby module."""

from __future__ import annotations

import enum
from dataclasses import dataclass, field
from typing import Any

import numpy as np


class SymmetryType(enum.Enum):
    """Crystal/material symmetry types mapped to simcoon umat names."""

    ISOTROPIC = "ELISO"
    CUBIC = "ELCUB"
    TRANSVERSE_ISOTROPIC = "ELIST"
    ORTHOTROPIC = "ELORT"
    ANISOTROPIC = "ELANI"
    UNKNOWN = "UNKNOWN"


@dataclass
class Material:
    """A single engineering material with mechanical/physical properties.

    Elastic moduli are stored in GPa, density in kg/m^3, strengths in MPa,
    CTE in 1/K.
    """

    # Identity
    name: str = ""
    formula: str = ""
    source_id: str = ""
    source: str = ""

    # Elastic (GPa)
    E: float | None = None
    nu: float | None = None
    G: float | None = None
    K: float | None = None
    elastic_tensor: np.ndarray | None = None  # optional 6x6

    # Physical
    density: float | None = None  # kg/m^3

    # Thermal
    CTE: float | None = None  # 1/K

    # Strength (MPa)
    yield_strength: float | None = None
    tensile_strength: float | None = None

    # Hardening
    hardening_type: str | None = None  # "linear", "power_law", "johnson_cook"
    hardening_params: dict[str, float] = field(default_factory=dict)

    # Classification
    symmetry: SymmetryType = SymmetryType.ISOTROPIC
    category: str = ""
    tags: list[str] = field(default_factory=list)

    def has_elastic_props(self) -> bool:
        """Return True if minimal elastic properties (E, nu) are available."""
        return self.E is not None and self.nu is not None

    def get_property(self, name: str) -> Any:
        """Return the value of an arbitrary property by attribute name."""
        return getattr(self, name, None)

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> Material:
        """Create a Material from a plain dictionary.

        Handles ``symmetry`` as a string (looked up in SymmetryType) and
        ``elastic_tensor`` as a nested list (converted to ndarray).
        """
        kw = dict(d)

        # Convert symmetry string -> enum
        sym = kw.pop("symmetry", None)
        if isinstance(sym, str):
            try:
                kw["symmetry"] = SymmetryType[sym.upper()]
            except KeyError:
                # Try matching by value (e.g. "ELISO")
                for member in SymmetryType:
                    if member.value == sym:
                        kw["symmetry"] = member
                        break
                else:
                    kw["symmetry"] = SymmetryType.UNKNOWN
        elif isinstance(sym, SymmetryType):
            kw["symmetry"] = sym

        # Convert elastic_tensor list -> ndarray
        et = kw.get("elastic_tensor")
        if et is not None and not isinstance(et, np.ndarray):
            kw["elastic_tensor"] = np.array(et)

        # Only pass recognised fields
        valid_fields = {f.name for f in cls.__dataclass_fields__.values()}
        filtered = {k: v for k, v in kw.items() if k in valid_fields}
        return cls(**filtered)


class MaterialCollection:
    """A collection of :class:`Material` objects with filtering and grouping.

    Supports iteration, indexing, slicing, ``len()``, and ``append()``.

    Parameters
    ----------
    materials : list of Material, optional
        Initial list of materials.  Defaults to an empty collection.
    """

    def __init__(self, materials: list[Material] | None = None):
        self._materials: list[Material] = list(materials) if materials else []

    # ------------------------------------------------------------------
    # Sequence-like interface
    # ------------------------------------------------------------------
    def __len__(self) -> int:
        return len(self._materials)

    def __iter__(self):
        return iter(self._materials)

    def __getitem__(self, idx):
        if isinstance(idx, slice):
            return MaterialCollection(self._materials[idx])
        return self._materials[idx]

    def append(self, mat: Material) -> None:
        self._materials.append(mat)

    def extend(self, mats) -> None:
        self._materials.extend(mats)

    def __repr__(self) -> str:
        return f"MaterialCollection({len(self._materials)} materials)"

    # ------------------------------------------------------------------
    # Filtering / grouping
    # ------------------------------------------------------------------
    def filter(self, **kwargs) -> MaterialCollection:
        """Return a new collection of materials matching all given criteria.

        Parameters
        ----------
        **kwargs
            Field name / value pairs.  A material is included only if
            every specified field equals the given value.

        Returns
        -------
        MaterialCollection

        Examples
        --------
        >>> metals = collection.filter(category="Metal")
        >>> iso_metals = collection.filter(category="Metal",
        ...                                symmetry=SymmetryType.ISOTROPIC)
        """
        def _match(mat: Material) -> bool:
            for key, val in kwargs.items():
                attr = getattr(mat, key, None)
                if attr != val:
                    return False
            return True

        return MaterialCollection([m for m in self._materials if _match(m)])

    def group_by(self, field_name: str) -> dict[str, MaterialCollection]:
        """Group materials by a field value.

        Parameters
        ----------
        field_name : str
            Attribute name of :class:`Material` to group by (e.g.
            ``"category"``).

        Returns
        -------
        dict[str, MaterialCollection]
            Mapping from each distinct field value to a sub-collection.
        """
        groups: dict[str, MaterialCollection] = {}
        for mat in self._materials:
            key = str(getattr(mat, field_name, ""))
            if key not in groups:
                groups[key] = MaterialCollection()
            groups[key].append(mat)
        return groups

    def get_property_arrays(
        self, x_prop: str, y_prop: str
    ) -> tuple[np.ndarray, np.ndarray, list[str]]:
        """Extract aligned numpy arrays for two properties.

        Materials where either property is ``None`` are skipped.

        Parameters
        ----------
        x_prop, y_prop : str
            Attribute names of :class:`Material`.

        Returns
        -------
        x : numpy.ndarray
            Values of *x_prop*.
        y : numpy.ndarray
            Values of *y_prop*.
        names : list of str
            Corresponding material names.
        """
        xs, ys, names = [], [], []
        for mat in self._materials:
            xv = getattr(mat, x_prop, None)
            yv = getattr(mat, y_prop, None)
            if xv is not None and yv is not None:
                xs.append(float(xv))
                ys.append(float(yv))
                names.append(mat.name)
        return np.array(xs), np.array(ys), names
