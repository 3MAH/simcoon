# -*- coding: utf-8 -*-

"""Fetch materials from online databases (Materials Project)."""

from __future__ import annotations

from typing import Any

import numpy as np

from simcoon.ashby.material import Material, MaterialCollection, SymmetryType


def fetch_materials(
    *,
    elements: list[str] | None = None,
    num_elements: tuple[int, int] | None = None,
    has_elastic: bool = True,
    limit: int = 100,
    api_key: str | None = None,
    source: str = "materials_project",
) -> MaterialCollection:
    """Fetch materials from an online database.

    Parameters
    ----------
    elements : list of str, optional
        Chemical elements to include (e.g. ``["Al", "O"]``).
    num_elements : tuple of int, optional
        ``(min, max)`` number of elements in the formula.
    has_elastic : bool
        If ``True`` (default), only return entries with elastic data.
    limit : int
        Maximum number of results.
    api_key : str, optional
        API key.  For Materials Project this can also be set via the
        ``MP_API_KEY`` environment variable.
    source : str
        Database source.  Currently only ``"materials_project"`` is
        supported.

    Returns
    -------
    MaterialCollection
    """
    if source == "materials_project":
        return _fetch_from_mp(
            elements=elements,
            num_elements=num_elements,
            has_elastic=has_elastic,
            limit=limit,
            api_key=api_key,
        )
    raise ValueError(f"Unknown source: {source!r}.  Use 'materials_project'.")


def _fetch_from_mp(
    *,
    elements: list[str] | None,
    num_elements: tuple[int, int] | None,
    has_elastic: bool,
    limit: int,
    api_key: str | None,
) -> MaterialCollection:
    """Internal: query Materials Project via ``mp-api``."""
    try:
        from mp_api.client import MPRester
    except ImportError:
        raise ImportError(
            "The 'mp-api' package is required to fetch data from Materials "
            "Project.  Install it with:\n\n"
            "    pip install 'simcoon[ashby]'\n\n"
            "or directly:\n\n"
            "    pip install mp-api\n"
        ) from None

    kwargs: dict[str, Any] = {}
    if elements is not None:
        kwargs["elements"] = elements
    if num_elements is not None:
        kwargs["num_elements"] = num_elements

    fields = [
        "material_id",
        "formula_pretty",
        "density",
        "symmetry",
    ]
    if has_elastic:
        fields.extend([
            "bulk_modulus",
            "shear_modulus",
            "elastic_tensor",
        ])
        kwargs["has_props"] = ["elasticity"]

    with MPRester(api_key) as mpr:
        docs = mpr.materials.summary.search(
            **kwargs,
            fields=fields,
            num_chunks=1,
        )

    materials: list[Material] = []
    for doc in docs[:limit]:
        mat = _mp_doc_to_material(doc, has_elastic=has_elastic)
        if mat is not None:
            materials.append(mat)

    return MaterialCollection(materials)


def _mp_doc_to_material(doc: Any, has_elastic: bool) -> Material | None:
    """Convert a single MP summary document to a :class:`Material`."""
    kw: dict[str, Any] = {
        "name": getattr(doc, "formula_pretty", ""),
        "formula": getattr(doc, "formula_pretty", ""),
        "source_id": str(getattr(doc, "material_id", "")),
        "source": "materials_project",
        "category": "Crystal",
    }

    density = getattr(doc, "density", None)
    if density is not None:
        kw["density"] = float(density) * 1000  # g/cm³ → kg/m³

    if has_elastic:
        K = None
        G = None

        bulk = getattr(doc, "bulk_modulus", None)
        if bulk is not None:
            K = getattr(bulk, "vrh", None)
        shear = getattr(doc, "shear_modulus", None)
        if shear is not None:
            G = getattr(shear, "vrh", None)

        if K is not None and G is not None and G > 0:
            kw["K"] = float(K)  # GPa
            kw["G"] = float(G)  # GPa
            kw["E"] = 9.0 * K * G / (3.0 * K + G)
            kw["nu"] = (3.0 * K - 2.0 * G) / (2.0 * (3.0 * K + G))
        else:
            return None  # Skip entries without usable elastic data

        et = getattr(doc, "elastic_tensor", None)
        if et is not None:
            original = getattr(et, "ieee_format", None)
            if original is not None:
                kw["elastic_tensor"] = np.array(original)

    # Detect symmetry from elastic tensor if available
    kw["symmetry"] = _detect_symmetry(kw.get("elastic_tensor"))

    return Material.from_dict(kw)


def _detect_symmetry(elastic_tensor: np.ndarray | None) -> SymmetryType:
    """Detect symmetry type using ``sim.check_symetries`` when possible."""
    if elastic_tensor is None:
        return SymmetryType.ISOTROPIC  # default assumption

    try:
        import simcoon as sim

        result = sim.check_symetries(elastic_tensor, 1.0e-2)
        umat = result.get("umat_type", "")
        for member in SymmetryType:
            if member.value == umat:
                return member
    except Exception:
        pass

    return SymmetryType.UNKNOWN
