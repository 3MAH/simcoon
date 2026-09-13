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
            material_orientation=MaterialOrientation(float(t[4]), float(t[5]), float(t[6])),
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
            material_orientation=MaterialOrientation(float(t[4]), float(t[5]), float(t[6])),
            geometry_orientation=GeometryOrientation(float(t[7]), float(t[8]), float(t[9])),
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
            material_orientation=MaterialOrientation(float(t[5]), float(t[6]), float(t[7])),
            a1=float(t[8]),
            a2=float(t[9]),
            a3=float(t[10]),
            geometry_orientation=GeometryOrientation(float(t[11]), float(t[12]), float(t[13])),
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
            material_orientation=MaterialOrientation(float(t[5]), float(t[6]), float(t[7])),
            L=float(t[8]),
            R=float(t[9]),
            geometry_orientation=GeometryOrientation(float(t[10]), float(t[11]), float(t[12])),
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
            material_orientation=MaterialOrientation(float(t[3]), float(t[4]), float(t[5])),
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


# =============================================================================
# Handing phases to the C++ side, in memory
# =============================================================================

def to_phase_dict(phase: Union[Phase, Layer, Ellipsoid, Cylinder],
                  number: Optional[int] = None) -> Dict:
    """The in-memory form of one phase, as the C++ side expects it (see ``sim.L_eff``).

    The keys are those of the JSON files, with one deliberate difference: ``props`` stays
    a plain sequence of numbers instead of the ``name -> value`` mapping ``save_*_json``
    writes, because the C++ side reads properties positionally.

    Angles are left in degrees, as in the files; the binding converts them to radians.
    """
    out = {
        'number': phase.number if number is None else number,
        'umat_name': phase.umat_name,
        'save': phase.save,
        'concentration': phase.concentration,
        'material_orientation': {
            'psi': phase.material_orientation.psi,
            'theta': phase.material_orientation.theta,
            'phi': phase.material_orientation.phi,
        },
        'nstatev': phase.nstatev,
        'props': np.asarray(phase.props, dtype=float),
    }
    geometry = getattr(phase, 'geometry_orientation', None)
    if geometry is not None:
        out['geometry_orientation'] = {
            'psi': geometry.psi, 'theta': geometry.theta, 'phi': geometry.phi,
        }
    if isinstance(phase, Ellipsoid):
        out['coatingof'] = phase.coatingof
        out['semi_axes'] = {'a1': phase.a1, 'a2': phase.a2, 'a3': phase.a3}
    elif isinstance(phase, Cylinder):
        out['coatingof'] = phase.coatingof
        out['L'] = phase.L
        out['R'] = phase.R
    elif isinstance(phase, Layer):
        out['layerup'] = phase.layerup
        out['layerdown'] = phase.layerdown
    return out


def to_phase_dicts(phases: List[Union[Phase, Layer, Ellipsoid, Cylinder]]) -> List[Dict]:
    """The whole sub-phase list, ready to be passed as ``phases=`` to ``sim.L_eff``.

    The phases are numbered by their **position** in the list, whatever each object
    carries. The mean-field schemes select the matrix with ``number == n_matrix`` and
    then index ``sub_phases[n_matrix]``, so the number has to be the position; the
    dataclass default is 0 for every phase, which would otherwise make every
    concentration tensor the identity and turn L_eff into a Voigt average.
    """
    return [to_phase_dict(p, number=i) for i, p in enumerate(phases)]


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
    # Legacy .dat input and its one-way conversion
    'load_phases_dat',
    'load_layers_dat',
    'load_ellipsoids_dat',
    'load_cylinders_dat',
    'load_sections_dat',
    'kind_from_dat_name',
    'convert_dat_to_json',
    # Handing phases to the C++ side, in memory
    'to_phase_dict',
    'to_phase_dicts',
]
