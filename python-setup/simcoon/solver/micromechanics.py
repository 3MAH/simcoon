"""
Micromechanics data classes and JSON I/O for Simcoon.

This module provides standalone dataclasses and I/O functions for micromechanics
homogenization without requiring simcoon._core or solver.py. This allows users to
work with micromechanics configurations (phases, layers, ellipsoids, etc.) without
building the C++ extension module.

Classes
-------
MaterialOrientation
    Material orientation via Euler angles
GeometryOrientation
    Geometry/phase orientation via Euler angles
Phase
    Generic phase for micromechanics homogenization
Layer
    Layer phase for laminate homogenization
Ellipsoid
    Ellipsoidal inclusion for Eshelby-based homogenization
Cylinder
    Cylindrical inclusion for micromechanics
Section
    Section/yarn for textile composite homogenization

Functions
---------
load_phases_json, save_phases_json
    JSON I/O for generic phases
load_layers_json, save_layers_json
    JSON I/O for layers (laminates)
load_ellipsoids_json, save_ellipsoids_json
    JSON I/O for ellipsoidal inclusions
load_cylinders_json, save_cylinders_json
    JSON I/O for cylindrical inclusions
load_sections_json, save_sections_json
    JSON I/O for textile sections

Example
-------
>>> from simcoon.solver.micromechanics import Ellipsoid, save_ellipsoids_json
>>> import numpy as np
>>>
>>> # Create ellipsoidal phases
>>> matrix = Ellipsoid(number=0, concentration=0.7, props=np.array([3000, 0.4]))
>>> fiber = Ellipsoid(number=1, concentration=0.3, a1=50, props=np.array([70000, 0.3]))
>>>
>>> # Save to JSON
>>> save_ellipsoids_json('phases.json', [matrix, fiber])
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Union, Optional

import numpy as np


# =============================================================================
# Data Classes
# =============================================================================

@dataclass
class MaterialOrientation:
    """Material orientation via Euler angles (degrees)."""
    psi: float = 0.0    # First Euler angle (deg)
    theta: float = 0.0  # Second Euler angle (deg)
    phi: float = 0.0    # Third Euler angle (deg)


@dataclass
class GeometryOrientation:
    """Geometry/phase orientation via Euler angles (degrees)."""
    psi: float = 0.0    # First Euler angle (deg)
    theta: float = 0.0  # Second Euler angle (deg)
    phi: float = 0.0    # Third Euler angle (deg)


@dataclass
class Phase:
    """
    Generic phase for micromechanics homogenization.

    Corresponds to Nphases.dat format and C++ phase_characteristics class.

    Attributes
    ----------
    number : int
        Phase identification number
    umat_name : str
        Constitutive model name (e.g., 'ELISO', 'ELIST')
    save : int
        Save flag (1=save, 0=don't)
    concentration : float
        Volume fraction (0 to 1)
    material_orientation : MaterialOrientation
        Material orientation via Euler angles
    nstatev : int
        Number of state variables
    props : np.ndarray
        Material properties array
    """
    number: int = 0
    umat_name: str = "ELISO"
    save: int = 1
    concentration: float = 1.0
    material_orientation: MaterialOrientation = field(default_factory=MaterialOrientation)
    nstatev: int = 1
    props: np.ndarray = field(default_factory=lambda: np.array([]))

    def __post_init__(self):
        if isinstance(self.props, list):
            self.props = np.array(self.props, dtype=float)
        if isinstance(self.material_orientation, dict):
            self.material_orientation = MaterialOrientation(**self.material_orientation)


@dataclass
class Layer(Phase):
    """
    Layer phase for laminate homogenization.

    Corresponds to Nlayers.dat format and C++ layer class.
    Layers are oriented using geometry orientation angles.

    Additional Attributes
    ---------------------
    geometry_orientation : GeometryOrientation
        Geometry orientation via Euler angles
    layerup : int
        Index of layer above (-1 if none)
    layerdown : int
        Index of layer below (-1 if none)
    """
    geometry_orientation: GeometryOrientation = field(default_factory=GeometryOrientation)
    layerup: int = -1
    layerdown: int = -1

    def __post_init__(self):
        super().__post_init__()
        if isinstance(self.geometry_orientation, dict):
            self.geometry_orientation = GeometryOrientation(**self.geometry_orientation)


@dataclass
class Ellipsoid(Phase):
    """
    Ellipsoidal inclusion for Eshelby-based homogenization.

    Corresponds to Nellipsoids.dat format and C++ ellipsoid class.

    Shape types based on semi-axis ratios:
    - Sphere: a1 = a2 = a3
    - Prolate spheroid (needle): a1 > a2 = a3
    - Oblate spheroid (disc): a1 = a2 > a3
    - General ellipsoid: a1 != a2 != a3

    Additional Attributes
    ---------------------
    coatingof : int
        Index of phase this ellipsoid coats (0 if none)
    a1 : float
        First semi-axis (relative)
    a2 : float
        Second semi-axis (relative)
    a3 : float
        Third semi-axis (relative)
    geometry_orientation : GeometryOrientation
        Geometry orientation via Euler angles
    """
    coatingof: int = 0
    a1: float = 1.0
    a2: float = 1.0
    a3: float = 1.0
    geometry_orientation: GeometryOrientation = field(default_factory=GeometryOrientation)

    def __post_init__(self):
        super().__post_init__()
        if isinstance(self.geometry_orientation, dict):
            self.geometry_orientation = GeometryOrientation(**self.geometry_orientation)

    @property
    def shape_type(self) -> str:
        """Determine shape type from semi-axes."""
        tol = 1e-6
        if abs(self.a1 - self.a2) < tol and abs(self.a2 - self.a3) < tol:
            return "sphere"
        elif abs(self.a2 - self.a3) < tol and self.a1 > self.a2:
            return "prolate_spheroid"
        elif abs(self.a1 - self.a2) < tol and self.a1 > self.a3:
            return "oblate_spheroid"
        else:
            return "general_ellipsoid"


@dataclass
class Cylinder(Phase):
    """
    Cylindrical inclusion for micromechanics.

    Corresponds to Ncylinders.dat format and C++ cylinder class.

    Additional Attributes
    ---------------------
    coatingof : int
        Index of phase this cylinder coats (0 if none)
    L : float
        Length parameter
    R : float
        Radius parameter
    geometry_orientation : GeometryOrientation
        Geometry orientation via Euler angles
    """
    coatingof: int = 0
    L: float = 1.0
    R: float = 1.0
    geometry_orientation: GeometryOrientation = field(default_factory=GeometryOrientation)

    def __post_init__(self):
        super().__post_init__()
        if isinstance(self.geometry_orientation, dict):
            self.geometry_orientation = GeometryOrientation(**self.geometry_orientation)

    @property
    def aspect_ratio(self) -> float:
        """Length to radius ratio."""
        return self.L / self.R if self.R > 0 else float('inf')


@dataclass
class Section:
    """
    Section/yarn for textile composite homogenization.

    Corresponds to Nsections.dat format.

    Attributes
    ----------
    number : int
        Section identification number
    name : str
        Section name
    umat_name : str
        Constitutive model name
    material_orientation : MaterialOrientation
        Material orientation via Euler angles
    nstatev : int
        Number of state variables
    props : np.ndarray
        Material properties array
    """
    number: int = 0
    name: str = "Section"
    umat_name: str = "ELISO"
    material_orientation: MaterialOrientation = field(default_factory=MaterialOrientation)
    nstatev: int = 1
    props: np.ndarray = field(default_factory=lambda: np.array([]))

    def __post_init__(self):
        if isinstance(self.props, list):
            self.props = np.array(self.props, dtype=float)
        if isinstance(self.material_orientation, dict):
            self.material_orientation = MaterialOrientation(**self.material_orientation)


# =============================================================================
# Helper Functions
# =============================================================================

def _props_to_dict(props: np.ndarray, prop_names: List[str] = None) -> Dict[str, float]:
    """Convert props array to dict with named keys."""
    if prop_names and len(prop_names) == len(props):
        return {name: float(val) for name, val in zip(prop_names, props)}
    else:
        return {f'prop_{i}': float(val) for i, val in enumerate(props)}


# =============================================================================
# JSON I/O - Phases
# =============================================================================

def load_phases_json(filepath: Union[str, Path]) -> List[Phase]:
    """
    Load phases from a JSON file.

    Parameters
    ----------
    filepath : str or Path
        Path to the JSON phases file

    Returns
    -------
    list of Phase
        List of Phase objects

    Example JSON format
    -------------------
    ```json
    {
      "phases": [
        {
          "number": 0,
          "umat_name": "ELISO",
          "save": 1,
          "concentration": 0.8,
          "material_orientation": {"psi": 0, "theta": 0, "phi": 0},
          "nstatev": 1,
          "props": {"E": 70000, "nu": 0.3, "alpha": 1e-5}
        }
      ]
    }
    ```
    """
    with open(filepath, 'r') as f:
        data = json.load(f)

    phases = []
    for p in data.get('phases', []):
        props = p.get('props', [])
        if isinstance(props, dict):
            props = np.array(list(props.values()), dtype=float)
        else:
            props = np.array(props, dtype=float)

        phase = Phase(
            number=p.get('number', 0),
            umat_name=p.get('umat_name', 'ELISO'),
            save=p.get('save', 1),
            concentration=p.get('concentration', 1.0),
            material_orientation=MaterialOrientation(**p.get('material_orientation', {})),
            nstatev=p.get('nstatev', 1),
            props=props
        )
        phases.append(phase)

    return phases


def save_phases_json(filepath: Union[str, Path], phases: List[Phase],
                     prop_names: List[str] = None):
    """
    Save phases to a JSON file.

    Parameters
    ----------
    filepath : str or Path
        Path to save the JSON file
    phases : list of Phase
        List of Phase objects
    prop_names : list of str, optional
        Names for the properties array
    """
    phases_data = []
    for p in phases:
        props_data = _props_to_dict(p.props, prop_names)
        phase_dict = {
            'number': p.number,
            'umat_name': p.umat_name,
            'save': p.save,
            'concentration': p.concentration,
            'material_orientation': {
                'psi': p.material_orientation.psi,
                'theta': p.material_orientation.theta,
                'phi': p.material_orientation.phi
            },
            'nstatev': p.nstatev,
            'props': props_data
        }
        phases_data.append(phase_dict)

    with open(filepath, 'w') as f:
        json.dump({'phases': phases_data}, f, indent=2)


# =============================================================================
# JSON I/O - Layers
# =============================================================================

def load_layers_json(filepath: Union[str, Path]) -> List[Layer]:
    """
    Load layers from a JSON file for laminate homogenization.

    Parameters
    ----------
    filepath : str or Path
        Path to the JSON layers file

    Returns
    -------
    list of Layer
        List of Layer objects

    Example JSON format
    -------------------
    ```json
    {
      "layers": [
        {
          "number": 0,
          "umat_name": "ELISO",
          "save": 1,
          "concentration": 0.8,
          "material_orientation": {"psi": 0, "theta": 0, "phi": 0},
          "geometry_orientation": {"psi": 0, "theta": 90, "phi": -90},
          "nstatev": 1,
          "props": {"E": 70000, "nu": 0.3, "alpha": 1e-5}
        }
      ]
    }
    ```
    """
    with open(filepath, 'r') as f:
        data = json.load(f)

    layers = []
    for lyr in data.get('layers', []):
        props = lyr.get('props', [])
        if isinstance(props, dict):
            props = np.array(list(props.values()), dtype=float)
        else:
            props = np.array(props, dtype=float)

        layer = Layer(
            number=lyr.get('number', 0),
            umat_name=lyr.get('umat_name', 'ELISO'),
            save=lyr.get('save', 1),
            concentration=lyr.get('concentration', 1.0),
            material_orientation=MaterialOrientation(**lyr.get('material_orientation', {})),
            geometry_orientation=GeometryOrientation(**lyr.get('geometry_orientation', {})),
            nstatev=lyr.get('nstatev', 1),
            props=props,
            layerup=lyr.get('layerup', -1),
            layerdown=lyr.get('layerdown', -1)
        )
        layers.append(layer)

    return layers


def save_layers_json(filepath: Union[str, Path], layers: List[Layer],
                     prop_names: List[str] = None):
    """
    Save layers to a JSON file.

    Parameters
    ----------
    filepath : str or Path
        Path to save the JSON file
    layers : list of Layer
        List of Layer objects
    prop_names : list of str, optional
        Names for the properties array
    """
    layers_data = []
    for lyr in layers:
        props_data = _props_to_dict(lyr.props, prop_names)
        layer_dict = {
            'number': lyr.number,
            'umat_name': lyr.umat_name,
            'save': lyr.save,
            'concentration': lyr.concentration,
            'material_orientation': {
                'psi': lyr.material_orientation.psi,
                'theta': lyr.material_orientation.theta,
                'phi': lyr.material_orientation.phi
            },
            'geometry_orientation': {
                'psi': lyr.geometry_orientation.psi,
                'theta': lyr.geometry_orientation.theta,
                'phi': lyr.geometry_orientation.phi
            },
            'nstatev': lyr.nstatev,
            'props': props_data,
            'layerup': lyr.layerup,
            'layerdown': lyr.layerdown
        }
        layers_data.append(layer_dict)

    with open(filepath, 'w') as f:
        json.dump({'layers': layers_data}, f, indent=2)


# =============================================================================
# JSON I/O - Ellipsoids
# =============================================================================

def load_ellipsoids_json(filepath: Union[str, Path]) -> List[Ellipsoid]:
    """
    Load ellipsoids from a JSON file for Eshelby-based homogenization.

    Parameters
    ----------
    filepath : str or Path
        Path to the JSON ellipsoids file

    Returns
    -------
    list of Ellipsoid
        List of Ellipsoid objects

    Example JSON format
    -------------------
    ```json
    {
      "ellipsoids": [
        {
          "number": 0,
          "coatingof": 0,
          "umat_name": "ELISO",
          "save": 1,
          "concentration": 0.2,
          "material_orientation": {"psi": 0, "theta": 0, "phi": 0},
          "semi_axes": {"a1": 50, "a2": 1, "a3": 1},
          "geometry_orientation": {"psi": 0, "theta": 0, "phi": 0},
          "nstatev": 1,
          "props": {"E": 70000, "nu": 0.3, "alpha": 1e-5}
        }
      ]
    }
    ```
    """
    with open(filepath, 'r') as f:
        data = json.load(f)

    ellipsoids = []
    for ell in data.get('ellipsoids', []):
        props = ell.get('props', [])
        if isinstance(props, dict):
            props = np.array(list(props.values()), dtype=float)
        else:
            props = np.array(props, dtype=float)

        semi_axes = ell.get('semi_axes', {})

        ellipsoid = Ellipsoid(
            number=ell.get('number', 0),
            coatingof=ell.get('coatingof', 0),
            umat_name=ell.get('umat_name', 'ELISO'),
            save=ell.get('save', 1),
            concentration=ell.get('concentration', 1.0),
            material_orientation=MaterialOrientation(**ell.get('material_orientation', {})),
            a1=semi_axes.get('a1', ell.get('a1', 1.0)),
            a2=semi_axes.get('a2', ell.get('a2', 1.0)),
            a3=semi_axes.get('a3', ell.get('a3', 1.0)),
            geometry_orientation=GeometryOrientation(**ell.get('geometry_orientation', {})),
            nstatev=ell.get('nstatev', 1),
            props=props
        )
        ellipsoids.append(ellipsoid)

    return ellipsoids


def save_ellipsoids_json(filepath: Union[str, Path], ellipsoids: List[Ellipsoid],
                         prop_names: List[str] = None):
    """
    Save ellipsoids to a JSON file.

    Parameters
    ----------
    filepath : str or Path
        Path to save the JSON file
    ellipsoids : list of Ellipsoid
        List of Ellipsoid objects
    prop_names : list of str, optional
        Names for the properties array
    """
    ellipsoids_data = []
    for ell in ellipsoids:
        props_data = _props_to_dict(ell.props, prop_names)
        ell_dict = {
            'number': ell.number,
            'coatingof': ell.coatingof,
            'umat_name': ell.umat_name,
            'save': ell.save,
            'concentration': ell.concentration,
            'material_orientation': {
                'psi': ell.material_orientation.psi,
                'theta': ell.material_orientation.theta,
                'phi': ell.material_orientation.phi
            },
            'semi_axes': {
                'a1': ell.a1,
                'a2': ell.a2,
                'a3': ell.a3
            },
            'geometry_orientation': {
                'psi': ell.geometry_orientation.psi,
                'theta': ell.geometry_orientation.theta,
                'phi': ell.geometry_orientation.phi
            },
            'nstatev': ell.nstatev,
            'props': props_data
        }
        ellipsoids_data.append(ell_dict)

    with open(filepath, 'w') as f:
        json.dump({'ellipsoids': ellipsoids_data}, f, indent=2)


# =============================================================================
# JSON I/O - Cylinders
# =============================================================================

def load_cylinders_json(filepath: Union[str, Path]) -> List[Cylinder]:
    """
    Load cylinders from a JSON file.

    Parameters
    ----------
    filepath : str or Path
        Path to the JSON cylinders file

    Returns
    -------
    list of Cylinder
        List of Cylinder objects

    Example JSON format
    -------------------
    ```json
    {
      "cylinders": [
        {
          "number": 0,
          "coatingof": 0,
          "umat_name": "ELISO",
          "save": 1,
          "concentration": 0.2,
          "material_orientation": {"psi": 0, "theta": 0, "phi": 0},
          "geometry": {"L": 50, "R": 1},
          "geometry_orientation": {"psi": 0, "theta": 0, "phi": 0},
          "nstatev": 1,
          "props": {"E": 70000, "nu": 0.3, "alpha": 1e-5}
        }
      ]
    }
    ```
    """
    with open(filepath, 'r') as f:
        data = json.load(f)

    cylinders = []
    for cyl in data.get('cylinders', []):
        props = cyl.get('props', [])
        if isinstance(props, dict):
            props = np.array(list(props.values()), dtype=float)
        else:
            props = np.array(props, dtype=float)

        geom = cyl.get('geometry', {})

        cylinder = Cylinder(
            number=cyl.get('number', 0),
            coatingof=cyl.get('coatingof', 0),
            umat_name=cyl.get('umat_name', 'ELISO'),
            save=cyl.get('save', 1),
            concentration=cyl.get('concentration', 1.0),
            material_orientation=MaterialOrientation(**cyl.get('material_orientation', {})),
            L=geom.get('L', cyl.get('L', 1.0)),
            R=geom.get('R', cyl.get('R', 1.0)),
            geometry_orientation=GeometryOrientation(**cyl.get('geometry_orientation', {})),
            nstatev=cyl.get('nstatev', 1),
            props=props
        )
        cylinders.append(cylinder)

    return cylinders


def save_cylinders_json(filepath: Union[str, Path], cylinders: List[Cylinder],
                        prop_names: List[str] = None):
    """
    Save cylinders to a JSON file.

    Parameters
    ----------
    filepath : str or Path
        Path to save the JSON file
    cylinders : list of Cylinder
        List of Cylinder objects
    prop_names : list of str, optional
        Names for the properties array
    """
    cylinders_data = []
    for cyl in cylinders:
        props_data = _props_to_dict(cyl.props, prop_names)
        cyl_dict = {
            'number': cyl.number,
            'coatingof': cyl.coatingof,
            'umat_name': cyl.umat_name,
            'save': cyl.save,
            'concentration': cyl.concentration,
            'material_orientation': {
                'psi': cyl.material_orientation.psi,
                'theta': cyl.material_orientation.theta,
                'phi': cyl.material_orientation.phi
            },
            'geometry': {
                'L': cyl.L,
                'R': cyl.R
            },
            'geometry_orientation': {
                'psi': cyl.geometry_orientation.psi,
                'theta': cyl.geometry_orientation.theta,
                'phi': cyl.geometry_orientation.phi
            },
            'nstatev': cyl.nstatev,
            'props': props_data
        }
        cylinders_data.append(cyl_dict)

    with open(filepath, 'w') as f:
        json.dump({'cylinders': cylinders_data}, f, indent=2)


# =============================================================================
# JSON I/O - Sections
# =============================================================================

def load_sections_json(filepath: Union[str, Path]) -> List[Section]:
    """
    Load sections from a JSON file for textile composites.

    Parameters
    ----------
    filepath : str or Path
        Path to the JSON sections file

    Returns
    -------
    list of Section
        List of Section objects

    Example JSON format
    -------------------
    ```json
    {
      "sections": [
        {
          "number": 0,
          "name": "Warp_yarn",
          "umat_name": "ELISO",
          "material_orientation": {"psi": 0, "theta": 0, "phi": 0},
          "nstatev": 1,
          "props": {"E": 70000, "nu": 0.3, "alpha": 1e-5}
        }
      ]
    }
    ```
    """
    with open(filepath, 'r') as f:
        data = json.load(f)

    sections = []
    for sec in data.get('sections', []):
        props = sec.get('props', [])
        if isinstance(props, dict):
            props = np.array(list(props.values()), dtype=float)
        else:
            props = np.array(props, dtype=float)

        section = Section(
            number=sec.get('number', 0),
            name=sec.get('name', 'Section'),
            umat_name=sec.get('umat_name', 'ELISO'),
            material_orientation=MaterialOrientation(**sec.get('material_orientation', {})),
            nstatev=sec.get('nstatev', 1),
            props=props
        )
        sections.append(section)

    return sections


def save_sections_json(filepath: Union[str, Path], sections: List[Section],
                       prop_names: List[str] = None):
    """
    Save sections to a JSON file.

    Parameters
    ----------
    filepath : str or Path
        Path to save the JSON file
    sections : list of Section
        List of Section objects
    prop_names : list of str, optional
        Names for the properties array
    """
    sections_data = []
    for sec in sections:
        props_data = _props_to_dict(sec.props, prop_names)
        sec_dict = {
            'number': sec.number,
            'name': sec.name,
            'umat_name': sec.umat_name,
            'material_orientation': {
                'psi': sec.material_orientation.psi,
                'theta': sec.material_orientation.theta,
                'phi': sec.material_orientation.phi
            },
            'nstatev': sec.nstatev,
            'props': props_data
        }
        sections_data.append(sec_dict)

    with open(filepath, 'w') as f:
        json.dump({'sections': sections_data}, f, indent=2)


# =============================================================================
# Exports
# =============================================================================

__all__ = [
    # Data classes
    'MaterialOrientation',
    'GeometryOrientation',
    'Phase',
    'Layer',
    'Ellipsoid',
    'Cylinder',
    'Section',
    # JSON I/O
    'load_phases_json',
    'save_phases_json',
    'load_layers_json',
    'save_layers_json',
    'load_ellipsoids_json',
    'save_ellipsoids_json',
    'load_cylinders_json',
    'save_cylinders_json',
    'load_sections_json',
    'save_sections_json',
]
