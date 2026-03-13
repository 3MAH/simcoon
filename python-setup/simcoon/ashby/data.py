# -*- coding: utf-8 -*-

"""Built-in dataset loaders for the simcoon Ashby module."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from simcoon.ashby.material import Material, MaterialCollection


def load_builtin() -> MaterialCollection:
    """Load the bundled curated material dataset.

    Returns a :class:`MaterialCollection` of ~60 common engineering
    materials.  Requires only numpy (no optional dependencies).
    """
    import importlib.resources as _res

    ref = _res.files("simcoon.ashby").joinpath("materials.json")
    text = ref.read_text(encoding="utf-8")
    entries: list[dict[str, Any]] = json.loads(text)
    return MaterialCollection([Material.from_dict(e) for e in entries])


def load_csv(
    path: str | Path,
    column_mapping: dict[str, str] | None = None,
) -> MaterialCollection:
    """Load materials from a user-supplied CSV file.

    Parameters
    ----------
    path : str or Path
        Path to the CSV file.  The first row must be a header.
    column_mapping : dict, optional
        Map CSV column names → :class:`Material` field names.  Any CSV
        column whose header (after mapping) matches a Material field is
        used; other columns are ignored.

    Returns
    -------
    MaterialCollection
    """
    import csv

    path = Path(path)
    mapping = column_mapping or {}
    valid_fields = {f.name for f in Material.__dataclass_fields__.values()}

    materials: list[Material] = []
    with path.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            d: dict[str, Any] = {}
            for csv_col, value in row.items():
                field_name = mapping.get(csv_col, csv_col)
                if field_name not in valid_fields or value == "":
                    continue
                # Attempt numeric conversion
                try:
                    value = float(value)  # type: ignore[assignment]
                except (ValueError, TypeError):
                    pass
                d[field_name] = value
            materials.append(Material.from_dict(d))
    return MaterialCollection(materials)
