# -*- coding: utf-8 -*-

"""Bridge between Ashby material data and simcoon stiffness/compliance tensors."""

from __future__ import annotations

from typing import Any

import numpy as np

from simcoon.ashby.material import Material, SymmetryType


def to_stiffness(
    material: Material,
    convention: str = "Enu",
    unit_factor: float = 1e3,
) -> np.ndarray:
    """Convert a :class:`Material` to a 6x6 stiffness tensor.

    Parameters
    ----------
    material : Material
        The material to convert.
    convention : str
        Convention string passed to simcoon's ``L_iso`` / ``L_cubic`` etc.
        Default ``"Enu"`` (Young's modulus + Poisson's ratio).
    unit_factor : float
        Multiplier applied to moduli.  Default ``1e3`` converts GPa
        (database convention) to MPa (simcoon convention).

    Returns
    -------
    numpy.ndarray
        6x6 stiffness matrix (Voigt notation).
    """
    import simcoon as sim

    # If a full elastic tensor is provided, just scale and return it
    if material.elastic_tensor is not None:
        return np.array(material.elastic_tensor) * unit_factor

    sym = material.symmetry

    if sym == SymmetryType.ISOTROPIC:
        if material.E is None or material.nu is None:
            raise ValueError(
                f"Material '{material.name}' lacks E and nu for isotropic stiffness."
            )
        props = [material.E * unit_factor, material.nu]
        return np.array(sim.L_iso(props, convention))

    if sym == SymmetryType.CUBIC:
        if material.E is None or material.nu is None or material.G is None:
            raise ValueError(
                f"Material '{material.name}' lacks E, nu, G for cubic stiffness."
            )
        props = [material.E * unit_factor, material.nu, material.G * unit_factor]
        return np.array(sim.L_cubic(props, "EnuG"))

    if sym == SymmetryType.TRANSVERSE_ISOTROPIC:
        raise NotImplementedError(
            "Transversely isotropic materials require directional properties "
            "(EL, ET, nuTL, nuTT, GLT) which are not stored in the current "
            "Material dataclass.  Provide a full elastic_tensor instead."
        )

    if sym == SymmetryType.ORTHOTROPIC:
        raise NotImplementedError(
            "Orthotropic materials require 9 independent constants which are "
            "not stored in the current Material dataclass.  Provide a full "
            "elastic_tensor instead."
        )

    raise ValueError(f"Unsupported symmetry type: {sym}")


def to_compliance(
    material: Material,
    convention: str = "Enu",
    unit_factor: float = 1e3,
) -> np.ndarray:
    """Convert a :class:`Material` to a 6x6 compliance tensor.

    Parameters follow the same conventions as :func:`to_stiffness`.

    Returns
    -------
    numpy.ndarray
        6x6 compliance matrix (Voigt notation).
    """
    import simcoon as sim

    if material.elastic_tensor is not None:
        return np.linalg.inv(np.array(material.elastic_tensor) * unit_factor)

    sym = material.symmetry

    if sym == SymmetryType.ISOTROPIC:
        if material.E is None or material.nu is None:
            raise ValueError(
                f"Material '{material.name}' lacks E and nu for isotropic compliance."
            )
        props = [material.E * unit_factor, material.nu]
        return np.array(sim.M_iso(props, convention))

    if sym == SymmetryType.CUBIC:
        if material.E is None or material.nu is None or material.G is None:
            raise ValueError(
                f"Material '{material.name}' lacks E, nu, G for cubic compliance."
            )
        props = [material.E * unit_factor, material.nu, material.G * unit_factor]
        return np.array(sim.M_cubic(props, "EnuG"))

    if sym == SymmetryType.TRANSVERSE_ISOTROPIC:
        raise NotImplementedError(
            "Transversely isotropic materials require directional properties."
        )

    if sym == SymmetryType.ORTHOTROPIC:
        raise NotImplementedError(
            "Orthotropic materials require 9 independent constants."
        )

    raise ValueError(f"Unsupported symmetry type: {sym}")


def to_solver_props(material: Material) -> dict[str, Any]:
    """Build a dict with ``umat_name`` and ``props`` for ``sim.solver()``.

    The returned dict contains:

    - ``umat_name`` — simcoon constitutive-model name (e.g. ``"ELISO"``)
    - ``props`` — numpy array of material properties in simcoon order
    - ``nstatev`` — suggested number of state variables (1 for elastic)

    Properties are in MPa (converted from GPa via ×1000).
    """
    sym = material.symmetry
    uf = 1e3  # GPa → MPa

    if sym == SymmetryType.ISOTROPIC:
        if material.E is None or material.nu is None:
            raise ValueError(
                f"Material '{material.name}' lacks E and nu."
            )
        return {
            "umat_name": "ELISO",
            "props": np.array([material.E * uf, material.nu]),
            "nstatev": 1,
        }

    if sym == SymmetryType.CUBIC:
        if material.E is None or material.nu is None or material.G is None:
            raise ValueError(
                f"Material '{material.name}' lacks E, nu, G."
            )
        return {
            "umat_name": "ELCUB",
            "props": np.array([
                material.E * uf, material.nu, material.G * uf
            ]),
            "nstatev": 1,
        }

    raise NotImplementedError(
        f"to_solver_props not implemented for symmetry {sym}.  "
        "Provide the solver inputs manually."
    )
