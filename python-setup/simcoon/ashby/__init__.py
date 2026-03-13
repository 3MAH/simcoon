# -*- coding: utf-8 -*-

"""simcoon.ashby — Material property charts, databases, and simcoon bridge.

Quick start::

    from simcoon.ashby import Material, load_builtin, ashby_plot

    mats = load_builtin()               # ~60 curated engineering materials
    ashby_plot(mats, "density", "E")     # classic E-vs-density Ashby chart
"""

from __future__ import annotations

# ------------------------------------------------------------------
# Core (no optional deps)
# ------------------------------------------------------------------
from simcoon.ashby.material import Material, MaterialCollection, SymmetryType
from simcoon.ashby.data import load_builtin, load_csv
from simcoon.ashby.bridge import to_stiffness, to_compliance, to_solver_props

# ------------------------------------------------------------------
# Plotting (requires matplotlib; import errors deferred to call time)
# ------------------------------------------------------------------
from simcoon.ashby.plotting import (
    ashby_plot,
    # Legacy re-exports (backward compatibility)
    unit_vector,
    length,
    uv_2,
    poly_convexHull,
    poly_enclose,
    ellip_enclose,
)

# ------------------------------------------------------------------
# Database access (mp-api lazy — only imported when called)
# ------------------------------------------------------------------
from simcoon.ashby.database import fetch_materials

__all__ = [
    # Data model
    "Material",
    "MaterialCollection",
    "SymmetryType",
    # Dataset loaders
    "load_builtin",
    "load_csv",
    # Plotting
    "ashby_plot",
    # Bridge to simcoon
    "to_stiffness",
    "to_compliance",
    "to_solver_props",
    # Database
    "fetch_materials",
    # Legacy backward-compat
    "unit_vector",
    "length",
    "uv_2",
    "poly_convexHull",
    "poly_enclose",
    "ellip_enclose",
]
