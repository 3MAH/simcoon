"""
Micromechanics data classes and JSON I/O for Simcoon.

Dataclasses and I/O functions describing the sub-phases of a mean-field model (phases,
layers, ellipsoids, ...), and their in-memory form for ``sim.L_eff`` and
``sim.solver.solve``. ``L_eff`` calls the C++ extension; the rest is pure Python.

Orientations (material frame of a phase, geometry of an inclusion or a layer) are
``simcoon.Rotation`` objects; ``as_rotation`` also takes their Euler angles in degrees.

Classes
-------
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
as_rotation, euler_angles
    An orientation as a simcoon.Rotation, and back to its (psi, theta, phi) in degrees
get_densities_ODF
    Density of an orientation distribution function (a sum of Peak profiles)
discretize_odf
    Split a phase into phases oriented along an ODF about a direction
L_eff
    Effective stiffness (``sim.L_eff``), orientation as a Rotation, phases as objects

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

import copy
import json
import warnings
from dataclasses import dataclass, field, fields
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Union

import numpy as np

from simcoon import _core

from simcoon.rotation import EULER_SEQ, Orientation, Rotation, as_rotation, euler_angles


# =============================================================================
# Data Classes
# =============================================================================

# The orientation coercion lives in simcoon.rotation, next to the Rotation class, because
# it serves every API that takes one (phases here, fibre directions in modular.py). The names
# this module has always published are re-exported above.


def _dataclass_eq(self, other):
    """Field-wise equality: the generated one compares Rotations by identity and raises
    on arrays, so a saved-then-loaded phase could never be compared to its source."""
    if other.__class__ is not self.__class__:
        return NotImplemented
    for f in fields(self):
        a, b = getattr(self, f.name), getattr(other, f.name)
        if isinstance(a, Rotation) and isinstance(b, Rotation):
            same = a.equals(b)
        elif isinstance(a, np.ndarray) or isinstance(b, np.ndarray):
            same = np.array_equal(a, b)
        else:
            same = a == b
        if not same:
            return False
    return True


def _coerce_fields(obj):
    """props given as a list, orientations in any accepted form, nested phases given as
    dicts (the JSON form)."""
    if not isinstance(obj.props, np.ndarray):
        obj.props = np.asarray([] if obj.props is None else obj.props, dtype=float).ravel()
    for name in ('material_orientation', 'geometry_orientation'):
        if hasattr(obj, name):
            setattr(obj, name, as_rotation(getattr(obj, name)))
    nested = getattr(obj, 'phases', None)
    if nested and any(isinstance(p, dict) for p in nested):
        cls, layout = _JSON_LAYOUTS[_json_kind_of(obj.umat_name)]
        obj.phases = [p if not isinstance(p, dict)
                      else _from_json_entry(p, cls, layout, None, f"phases of {obj.umat_name}")
                      for p in nested]


@dataclass(eq=False)
class Phase:
    """
    Generic phase for micromechanics homogenization.

    Mirrors the C++ phase_characteristics class.

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
    material_orientation : Rotation
        Orientation of the material frame (any form ``as_rotation`` accepts)
    nstatev : int
        Number of state variables
    props : np.ndarray
        Material properties array
    phases : list
        Sub-phases of a phase that is itself a mean-field model (MIHEN, MIMTN, MISCN:
        ellipsoids; MIPLN: layers).
        Empty for a homogeneous phase.
    """
    number: int = 0
    umat_name: str = "ELISO"
    save: int = 1
    concentration: float = 1.0
    material_orientation: Rotation = field(default_factory=Rotation.identity)
    nstatev: int = 1
    props: np.ndarray = field(default_factory=lambda: np.array([]))
    phases: List["Phase"] = field(default_factory=list)

    def __post_init__(self):
        _coerce_fields(self)

    __eq__ = _dataclass_eq


@dataclass(eq=False)
class Layer(Phase):
    """
    Layer phase for laminate homogenization.

    Mirrors the C++ layer class.
    Layers are oriented using geometry orientation angles.

    Additional Attributes
    ---------------------
    geometry_orientation : Rotation
        Orientation of the geometry (any form ``as_rotation`` accepts)
    layerup : int
        Index of layer above (0 by default, as the C++ layer)
    layerdown : int
        Index of layer below (0 by default, as the C++ layer)
    """
    geometry_orientation: Rotation = field(default_factory=Rotation.identity)
    layerup: int = 0
    layerdown: int = 0


@dataclass(eq=False)
class Ellipsoid(Phase):
    """
    Ellipsoidal inclusion for Eshelby-based homogenization.

    Mirrors the C++ ellipsoid class.

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
    geometry_orientation : Rotation
        Orientation of the geometry (any form ``as_rotation`` accepts)
    """
    coatingof: int = 0
    a1: float = 1.0
    a2: float = 1.0
    a3: float = 1.0
    geometry_orientation: Rotation = field(default_factory=Rotation.identity)

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


@dataclass(eq=False)
class Cylinder(Phase):
    """
    Cylindrical inclusion for micromechanics.

    Mirrors the C++ cylinder class.

    Additional Attributes
    ---------------------
    coatingof : int
        Index of phase this cylinder coats (0 if none)
    L : float
        Length parameter
    R : float
        Radius parameter
    geometry_orientation : Rotation
        Orientation of the geometry (any form ``as_rotation`` accepts)
    """
    coatingof: int = 0
    L: float = 1.0
    R: float = 1.0
    geometry_orientation: Rotation = field(default_factory=Rotation.identity)

    @property
    def aspect_ratio(self) -> float:
        """Length to radius ratio."""
        return self.L / self.R if self.R > 0 else float('inf')


@dataclass(eq=False)
class Section:
    """
    Section/yarn for textile composite homogenization.


    Attributes
    ----------
    number : int
        Section identification number
    name : str
        Section name
    umat_name : str
        Constitutive model name
    material_orientation : Rotation
        Orientation of the material frame (any form ``as_rotation`` accepts)
    nstatev : int
        Number of state variables
    props : np.ndarray
        Material properties array
    """
    number: int = 0
    name: str = "Section"
    umat_name: str = "ELISO"
    material_orientation: Rotation = field(default_factory=Rotation.identity)
    nstatev: int = 1
    props: np.ndarray = field(default_factory=lambda: np.array([]))

    def __post_init__(self):
        _coerce_fields(self)

    __eq__ = _dataclass_eq


# =============================================================================
# JSON I/O
# =============================================================================

def _props_to_dict(props: np.ndarray, prop_names: Optional[List[str]] = None) -> Dict[str, float]:
    """Convert props array to dict with named keys."""
    if prop_names and len(prop_names) == len(props):
        return {name: float(val) for name, val in zip(prop_names, props)}
    else:
        return {f'prop_{i}': float(val) for i, val in enumerate(props)}


def _props_from_json(props, prop_names: Optional[List[str]], context: str) -> np.ndarray:
    """Positional props out of a JSON entry.

    ``save_*_json`` writes ``props`` as a ``name -> value`` mapping when it is given
    ``prop_names``, but the C++ side reads properties **positionally** and the file
    records that order nowhere but in its own key order. Read back blindly, a file whose
    keys were reordered by hand or by a formatter silently describes another material:
    ``{"nu": 0.3, "E": 70000}`` is read as E = 0.3, nu = 70000, and ``L_iso`` then
    returns nonsense without an error. Pass ``prop_names`` to order explicitly.
    """
    if not isinstance(props, dict):
        return np.array(props, dtype=float)
    if prop_names:
        #either the caller's explicit order, or the one save_*_json recorded next to
        #the mapping it wrote — a file written by simcoon round-trips without warning
        missing = [n for n in prop_names if n not in props]
        if missing:
            raise ValueError(f"{context}: 'props' has no entry for {missing}")
        return np.array([float(props[n]) for n in prop_names])
    warnings.warn(
        f"{context}: 'props' is a name -> value mapping and no `prop_names` was given, "
        "so the values are read in the order the file lists them. The C++ side reads "
        "them positionally: reordering the keys changes the material silently.",
        UserWarning, stacklevel=4)
    return np.array(list(props.values()), dtype=float)


# The file layout of each kind: the dataclass and its entry keys, in file order. A plain
# name is a field written as is (an orientation as its {psi, theta, phi} dict); ('key',
# fields) is a nested object holding those fields; 'props' is the name -> value mapping
# followed by the 'prop_names' that record its order; 'phases' the sub-phases of a phase
# that is itself a mean-field model, written only when there are some. Defaults on
# reading are the dataclass defaults, and a nested field is also accepted flat at the top
# level.
_JSON_LAYOUTS = {
    'phases': (Phase, ('number', 'umat_name', 'save', 'concentration',
                       'material_orientation', 'nstatev', 'props', 'phases')),
    'layers': (Layer, ('number', 'umat_name', 'save', 'concentration',
                       'material_orientation', 'geometry_orientation', 'nstatev', 'props',
                       'layerup', 'layerdown', 'phases')),
    'ellipsoids': (Ellipsoid, ('number', 'coatingof', 'umat_name', 'save', 'concentration',
                               'material_orientation', ('semi_axes', ('a1', 'a2', 'a3')),
                               'geometry_orientation', 'nstatev', 'props', 'phases')),
    'cylinders': (Cylinder, ('number', 'coatingof', 'umat_name', 'save', 'concentration',
                             'material_orientation', ('geometry', ('L', 'R')),
                             'geometry_orientation', 'nstatev', 'props', 'phases')),
    'sections': (Section, ('number', 'name', 'umat_name', 'material_orientation',
                           'nstatev', 'props')),
}

def _kind_of(phase) -> str:
    """What to_phase_dict labels a phase with; the binding checks it against the geometry
    the model builds."""
    for cls, kind in ((Ellipsoid, 'ellipsoid'), (Cylinder, 'cylinder'), (Layer, 'layer')):
        if isinstance(phase, cls):
            return kind
    return 'phase'


def _json_kind_of(umat_name: str) -> str:
    """The layout of the sub-phases a mean-field model builds (sub_phase_shape in C++)."""
    if umat_name in ('MIHEN', 'MIMTN', 'MISCN'):
        return 'ellipsoids'
    if umat_name == 'MIPLN':
        return 'layers'
    raise ValueError(f"{umat_name} is not a mean-field model: it has no sub-phases")


def _to_json_entry(obj, layout, prop_names: Optional[List[str]]) -> Dict:
    entry = {}
    for col in layout:
        if col == 'props':
            props_data = _props_to_dict(obj.props, prop_names)
            entry['props'] = props_data
            entry['prop_names'] = list(props_data)
        elif col == 'phases':
            if obj.phases:
                _, sub_layout = _JSON_LAYOUTS[_json_kind_of(obj.umat_name)]
                entry['phases'] = [_to_json_entry(p, sub_layout, None) for p in obj.phases]
        elif isinstance(col, tuple):
            key, fields = col
            entry[key] = {f: getattr(obj, f) for f in fields}
        elif col.endswith('_orientation'):
            entry[col] = euler_angles(getattr(obj, col))
        else:
            entry[col] = getattr(obj, col)
    return entry


def _from_json_entry(entry: Dict, cls, layout, prop_names: Optional[List[str]], context: str):
    kwargs = {}
    for col in layout:
        if col == 'props':
            kwargs['props'] = _props_from_json(entry.get('props', []),
                                               prop_names or entry.get('prop_names'), context)
        elif col == 'phases':
            if entry.get('phases'):
                sub_cls, sub_layout = _JSON_LAYOUTS[_json_kind_of(entry.get('umat_name', 'ELISO'))]
                kwargs['phases'] = [_from_json_entry(p, sub_cls, sub_layout, None, context)
                                    for p in entry['phases']]
        elif isinstance(col, tuple):
            key, fields = col
            nested = entry.get(key, {})
            for f in fields:
                if f in nested:
                    kwargs[f] = nested[f]
                elif f in entry:
                    kwargs[f] = entry[f]
        elif col in entry:
            kwargs[col] = entry[col]   # orientation dicts are coerced by __post_init__
    return cls(**kwargs)


def _load_json(filepath: Union[str, Path], kind: str, prop_names: Optional[List[str]]) -> List:
    cls, layout = _JSON_LAYOUTS[kind]
    with open(filepath, 'r') as f:
        data = json.load(f)
    if kind not in data:
        raise ValueError(f"{filepath}: no '{kind}' entry (top-level keys: {sorted(data)})")
    return [_from_json_entry(entry, cls, layout, prop_names, str(filepath))
            for entry in data[kind]]


def _save_json(filepath: Union[str, Path], kind: str, items: List,
               prop_names: Optional[List[str]]):
    _, layout = _JSON_LAYOUTS[kind]
    payload = {kind: [_to_json_entry(item, layout, prop_names) for item in items]}
    with open(filepath, 'w') as f:   # opened once the payload exists: a failed save keeps the file
        json.dump(payload, f, indent=2)


def load_phases_json(filepath: Union[str, Path],
                     prop_names: Optional[List[str]] = None) -> List[Phase]:
    """Load phases from a JSON file: ``{"phases": [{"number": 0, "umat_name": "ELISO",
    "save": 1, "concentration": 0.8, "material_orientation": {"psi": 0, "theta": 0,
    "phi": 0}, "nstatev": 1, "props": {"E": 70000, "nu": 0.3, "alpha": 1e-5}}]}``.

    ``prop_names`` gives the order of a ``props`` mapping explicitly (see
    ``_props_from_json``); a file written by ``save_phases_json`` records it itself.
    """
    return _load_json(filepath, 'phases', prop_names)


def save_phases_json(filepath: Union[str, Path], phases: List[Phase],
                     prop_names: Optional[List[str]] = None):
    """Save phases to a JSON file (the layout ``load_phases_json`` reads), ``prop_names``
    naming the properties."""
    _save_json(filepath, 'phases', phases, prop_names)


def load_layers_json(filepath: Union[str, Path],
                     prop_names: Optional[List[str]] = None) -> List[Layer]:
    """Load layers from a JSON file for laminate homogenization: the phase entry plus
    ``"geometry_orientation": {"psi", "theta", "phi"}``, ``"layerup"`` and
    ``"layerdown"``, under the top-level key ``"layers"``."""
    return _load_json(filepath, 'layers', prop_names)


def save_layers_json(filepath: Union[str, Path], layers: List[Layer],
                     prop_names: Optional[List[str]] = None):
    """Save layers to a JSON file (the layout ``load_layers_json`` reads)."""
    _save_json(filepath, 'layers', layers, prop_names)


def load_ellipsoids_json(filepath: Union[str, Path],
                         prop_names: Optional[List[str]] = None) -> List[Ellipsoid]:
    """Load ellipsoids from a JSON file: the phase entry plus ``"coatingof"``,
    ``"semi_axes": {"a1", "a2", "a3"}`` and ``"geometry_orientation"``, under the
    top-level key ``"ellipsoids"``."""
    return _load_json(filepath, 'ellipsoids', prop_names)


def save_ellipsoids_json(filepath: Union[str, Path], ellipsoids: List[Ellipsoid],
                         prop_names: Optional[List[str]] = None):
    """Save ellipsoids to a JSON file (the layout ``load_ellipsoids_json`` reads)."""
    _save_json(filepath, 'ellipsoids', ellipsoids, prop_names)


def load_cylinders_json(filepath: Union[str, Path],
                        prop_names: Optional[List[str]] = None) -> List[Cylinder]:
    """Load cylinders from a JSON file: the phase entry plus ``"coatingof"``,
    ``"geometry": {"L", "R"}`` and ``"geometry_orientation"``, under the top-level key
    ``"cylinders"``."""
    return _load_json(filepath, 'cylinders', prop_names)


def save_cylinders_json(filepath: Union[str, Path], cylinders: List[Cylinder],
                        prop_names: Optional[List[str]] = None):
    """Save cylinders to a JSON file (the layout ``load_cylinders_json`` reads)."""
    _save_json(filepath, 'cylinders', cylinders, prop_names)


def load_sections_json(filepath: Union[str, Path],
                       prop_names: Optional[List[str]] = None) -> List[Section]:
    """Load sections from a JSON file for textile composites: ``"number"``, ``"name"``,
    ``"umat_name"``, ``"material_orientation"``, ``"nstatev"``, ``"props"``, under the
    top-level key ``"sections"``."""
    return _load_json(filepath, 'sections', prop_names)


def save_sections_json(filepath: Union[str, Path], sections: List[Section],
                       prop_names: Optional[List[str]] = None):
    """Save sections to a JSON file (the layout ``load_sections_json`` reads)."""
    _save_json(filepath, 'sections', sections, prop_names)


# =============================================================================
# Orientation distribution functions (ODF): the peaks, in memory
# =============================================================================

#: parameters each profile reads, for validation: (needs s_dev, needs width, len(params))
_PROFILES = {1: (False, False, 4), 2: (True, False, 0), 3: (True, False, 0), 4: (False, True, 0),
             5: (True, True, 1), 6: (False, True, 2), 7: (False, False, 0)}


@dataclass(eq=False)
class Peak:
    r"""One peak of an orientation (ODF) or parameter (PDF) distribution.

    ``method`` selects the profile, with :math:`d = x - \text{mean}`:

    1. standard-deviation kernel,
       :math:`|(a_1 \cos^{2p_1} d + a_2 \cos^{2p_2} d\,\sin^{2p_2} d) \cos d|`,
       ``params = [a1, a2, p1, p2]``
    2. hard cut-off, :math:`A \exp(-\tfrac{1}{2}(|d|/s)^2)`
    3. Gaussian, :math:`\frac{A}{s\sqrt{2\pi}} \exp(-\tfrac{1}{2}(d/s)^2)`
    4. Lorentzian, :math:`\frac{A\,w}{2\pi\,(d^2 + (w/2)^2)}`
    5. pseudo-Voigt, :math:`\eta\,L + (1-\eta)\,G`, ``params = [eta]``
    6. Pearson VII, :math:`M\,(1 + (d/w)^2/m)^{-m}`, ``params = [M, m]`` (``M = 0`` reads 1)
    7. uniform, 1

    with :math:`A` = ``ampl``, :math:`s` = ``s_dev``, :math:`w` = ``width``. For an ODF the
    angles (``mean``, ``s_dev``, ``width``) are degrees; the profile is evaluated in
    radians, which sets the Gaussian and Lorentzian normalisations.
    """
    number: int = 0
    method: int = 3
    mean: float = 0.0
    s_dev: float = 1.0
    width: float = 1.0
    ampl: float = 1.0
    params: np.ndarray = field(default_factory=lambda: np.array([]))

    def __post_init__(self):
        if not isinstance(self.params, np.ndarray):
            self.params = np.asarray(self.params, dtype=float).ravel()
        if self.method not in _PROFILES:
            raise ValueError(f"Peak.method = {self.method!r}: 1 to 7 expected (see the class docstring)")
        needs_s_dev, needs_width, n_params = _PROFILES[self.method]
        if needs_s_dev and not self.s_dev > 0.0:
            raise ValueError(f"Peak (method {self.method}): s_dev = {self.s_dev} must be > 0")
        if needs_width and not self.width > 0.0:
            raise ValueError(f"Peak (method {self.method}): width = {self.width} must be > 0")
        if self.params.size < n_params:
            raise ValueError(f"Peak (method {self.method}): params needs {n_params} values, got {self.params.size}")
        if self.method == 6 and not self.params[1] > 0.0:
            raise ValueError("Peak (method 6): the Pearson VII shape params[1] must be > 0")

    __eq__ = _dataclass_eq

    def density(self, x, periodic: bool = False, scale: float = 1.0) -> np.ndarray:
        """The profile at ``x``. ``scale`` converts ``mean``, ``s_dev`` and ``width`` to the
        unit of ``x`` (pi/180 for degrees against radians). ``periodic`` adds the images at
        plus and minus pi: a director distribution, as an ODF is (profiles 2 to 6)."""
        x = np.asarray(x, dtype=float)
        mean, s_dev, width = self.mean * scale, self.s_dev * scale, self.width * scale

        def profile(d):
            if self.method == 1:
                a1, a2, p1, p2 = self.params[:4]
                c = np.cos(d)
                with np.errstate(invalid='ignore'):   # a negative cosine to a fractional power
                    y = np.abs((a1 * c ** (2 * p1) + a2 * c ** (2 * p2) * np.sin(d) ** (2 * p2)) * c)
                y = np.where(np.abs(d - 0.5 * np.pi) < 1e-6, 0.0, y)
                return np.where(np.abs(d) < 1e-6, a1, y)
            if self.method == 2:
                return self.ampl * np.exp(-0.5 * (np.abs(d) / s_dev) ** 2)
            gauss = lambda: self.ampl / (s_dev * np.sqrt(2 * np.pi)) * np.exp(-0.5 * (d / s_dev) ** 2)
            lorentz = lambda: self.ampl * width / (2 * np.pi * (d ** 2 + (width / 2) ** 2))
            if self.method == 3:
                return gauss()
            if self.method == 4:
                return lorentz()
            if self.method == 5:
                eta = self.params[0]
                return eta * lorentz() + (1.0 - eta) * gauss()
            if self.method == 6:
                peak_max, shape = self.params[:2]
                return (peak_max if abs(peak_max) >= 1e-9 else 1.0) * (1.0 + (d / width) ** 2 / shape) ** (-shape)
            return np.ones_like(d)

        d = x - mean
        if periodic and 2 <= self.method <= 6:
            return profile(d) + profile(d - np.pi) + profile(d + np.pi)
        return profile(d)


_PEAK_FIELDS = ('number', 'method', 'mean', 's_dev', 'width', 'ampl', 'params')


def _as_peaks(peaks) -> List[Peak]:
    """Peak objects out of Peak objects or of their dicts (the JSON entries)."""
    out = []
    for pk in peaks:
        if isinstance(pk, dict):
            unknown = set(pk) - set(_PEAK_FIELDS)
            if unknown:
                raise ValueError(f"peak: unknown entries {sorted(unknown)}; expected {_PEAK_FIELDS}")
            if pk.get('method') is None:
                raise ValueError("peak: no 'method' entry (1 to 7, see Peak)")
            pk = Peak(**pk)
        out.append(pk)
    return out


def get_densities_ODF(x, peaks, radian: bool = False) -> np.ndarray:
    """Density of an orientation distribution function at the angles ``x``: the sum of its
    peaks (:class:`Peak` objects or their dicts), each a director distribution of period
    180 degrees. ``x`` lies in [0, 180] degrees, or [0, pi] with ``radian``; the angles of
    the peaks follow the same unit."""
    x = np.asarray(x, dtype=float)
    scale = 1.0 if radian else np.pi / 180.0
    x_rad = x * scale
    if x.size and (x_rad.min() < 0.0 or x_rad.max() > np.pi * (1 + 1e-12)):
        raise ValueError(f"get_densities_ODF: x must lie in [0, {'pi' if radian else '180'}], "
                         f"got [{x.min():g}, {x.max():g}]")
    return sum((pk.density(x_rad, periodic=True, scale=scale) for pk in _as_peaks(peaks)),
               np.zeros_like(x_rad))


def load_peaks_json(filepath: Union[str, Path]) -> List[Peak]:
    """Load the peaks of a distribution: ``{"peaks": [{"number": 0, "method": 3,
    "mean": 90, "s_dev": 10, "width": 0, "ampl": 1, "params": []}, ...]}``."""
    with open(filepath, 'r') as f:
        data = json.load(f)
    return [Peak(**{k: v for k, v in entry.items() if k in _PEAK_FIELDS})
            for entry in data.get('peaks', [])]


def save_peaks_json(filepath: Union[str, Path], peaks: List[Peak]):
    """Save the peaks of a distribution (the layout ``load_peaks_json`` reads)."""
    entries = [{f: (np.asarray(getattr(p, f), dtype=float).tolist() if f == 'params' else getattr(p, f))
                for f in _PEAK_FIELDS} for p in peaks]
    with open(filepath, 'w') as f:
        json.dump({'peaks': entries}, f, indent=2)


def discretize_odf(phases: List, num_phase: int, peaks: List, nphases: int,
                   axis=(0.0, 0.0, 1.0), angle_range=(0.0, 180.0),
                   rotate_material: bool = True) -> List:
    """Split one phase into ``nphases`` phases whose orientations follow an ODF.

    Phase ``num_phase`` of ``phases`` is replaced by ``nphases`` copies of itself, the
    k-th rotated by the angle ``alpha_k`` about the direction ``axis`` (a unit vector of
    the global frame): its geometry orientation becomes
    ``Rotation.from_rotvec(alpha_k * axis) * geometry_orientation``, and so does its
    material orientation when ``rotate_material`` is true. The base orientations of the
    phase are kept, so a tilted or pre-rotated inclusion is swept from where it stands.

    ``alpha_k = angle_min + k * d``, ``d = (angle_max - angle_min) / nphases``, degrees.
    Each copy takes the fraction of the parent's concentration given by the integral
    of the ODF density (``peaks``, see :class:`Peak`) over ``[alpha_k - d/2, alpha_k
    + d/2]`` by Simpson's rule, the fractions being normalised over the sweep. The
    density is that of a director distribution: periodic over 180 degrees, evaluated
    modulo 180. The default range, one half turn, is the right one when a half turn
    about ``axis`` brings the inclusion back onto itself, i.e. when ``axis`` is along
    or perpendicular to a principal axis of the inclusion (a fibre swept about a
    transverse direction, a disc about its normal). About any other direction the
    physically distinct orientations span a full turn: give ``angle_range=(0, 360)``,
    knowing that the density then repeats over the two half turns.

    Euler-angle sweeps map onto this: a sweep of ``psi`` or
    ``phi`` from an unrotated phase is ``axis=(0, 0, 1)``, a sweep of ``theta`` is
    ``axis=(1, 0, 0)``.

    Returns a new list of new objects, numbered by position; the caller's phases are
    left untouched. ``coatingof`` indices pointing past the swept phase are shifted; a
    coating OF the swept phase is refused (which of its copies would it coat?).
    """
    phases = phases_from_dicts([p for p in phases]) if any(isinstance(p, dict) for p in phases) else list(phases)
    if not 0 <= num_phase < len(phases):
        raise ValueError(f"num_phase = {num_phase} is outside the {len(phases)} phases given")
    if nphases < 1:
        raise ValueError("nphases must be >= 1")
    n = np.asarray(axis, dtype=float).ravel()
    if n.shape != (3,) or np.linalg.norm(n) < 1e-12:
        raise ValueError("axis must be a non-zero direction vector of 3 components")
    n /= np.linalg.norm(n)
    a_min, a_max = (float(a) for a in angle_range)
    if a_max <= a_min:
        raise ValueError("angle_range must be (min, max) with max > min, in degrees")

    parent = phases[num_phase]
    d = (a_max - a_min) / nphases
    alphas = a_min + d * np.arange(nphases)
    # Simpson over each bin, the density being periodic over 180 deg (director distribution)
    x = np.concatenate([alphas - d / 2, alphas, alphas + d / 2]) % 180.0
    rho = get_densities_ODF(x, peaks)
    weights = d / 6.0 * (rho[:nphases] + 4.0 * rho[nphases:2 * nphases] + rho[2 * nphases:])
    if weights.sum() <= 0.0:
        raise ValueError("the ODF density is zero over the whole angle_range")
    weights *= parent.concentration / weights.sum()

    swept = []
    for alpha, w in zip(alphas, weights):
        rot = Rotation.from_rotvec(np.deg2rad(alpha) * n)
        copy_k = copy.deepcopy(parent)
        copy_k.concentration = float(w)
        if hasattr(copy_k, 'geometry_orientation'):
            copy_k.geometry_orientation = rot * copy_k.geometry_orientation
        if rotate_material:
            copy_k.material_orientation = rot * copy_k.material_orientation
        swept.append(copy_k)
    others = [copy.copy(ph) for ph in phases]
    for ph in others:
        coated = getattr(ph, 'coatingof', 0)
        if coated == num_phase and num_phase != 0 and ph is not others[num_phase]:
            raise ValueError(f"phase {ph.number} coats phase {num_phase}, which is being discretised")
        if coated > num_phase:
            ph.coatingof = coated + nphases - 1
    out = others[:num_phase] + swept + others[num_phase + 1:]
    for i, ph in enumerate(out):
        ph.number = i
    return out


# =============================================================================
# Handing phases to the C++ side, in memory
# =============================================================================

def to_phase_dict(phase: Union[Phase, Layer, Ellipsoid, Cylinder],
                  number: Optional[int] = None) -> Dict:
    """The in-memory form of one phase, as the C++ side expects it (see ``sim.L_eff``).

    The keys are those of the JSON files, with one deliberate difference: ``props`` stays
    a plain sequence of numbers instead of the ``name -> value`` mapping ``save_*_json``
    writes, because the C++ side reads properties positionally. ``kind`` names the
    geometry, which the binding checks against the one the model builds, and ``phases``
    carries the sub-phases of a phase that is itself a mean-field model.

    Angles are left in degrees; the binding converts them to radians.
    """
    out = {
        'kind': _kind_of(phase),
        'number': phase.number if number is None else number,
        'umat_name': phase.umat_name,
        'save': phase.save,
        'concentration': phase.concentration,
        'material_orientation': euler_angles(phase.material_orientation),
        'nstatev': phase.nstatev,
        'props': np.asarray(phase.props, dtype=float),
    }
    geometry = getattr(phase, 'geometry_orientation', None)
    if geometry is not None:
        out['geometry_orientation'] = euler_angles(geometry)
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
    if phase.phases:
        out['phases'] = to_phase_dicts(phase.phases)
    return out


_CLASS_OF_KIND = {'phase': Phase, 'layer': Layer, 'ellipsoid': Ellipsoid, 'cylinder': Cylinder}


def phases_from_dicts(dicts: List[Dict]) -> List[Union[Phase, Layer, Ellipsoid, Cylinder]]:
    """The dataclasses back from the dicts ``to_phase_dicts`` made (or any dict in that
    form); ``kind`` picks the class, ``phases`` nests."""
    out = []
    for d in dicts:
        cls = _CLASS_OF_KIND[d.get('kind', 'phase')]
        _, layout = _JSON_LAYOUTS[d.get('kind', 'phase') + 's']
        out.append(_from_json_entry(d, cls, layout, None, "phases"))
    return out


def to_phase_dicts(phases: List[Union[Phase, Layer, Ellipsoid, Cylinder, Dict]]) -> List[Dict]:
    """The whole sub-phase list, ready to be passed as ``phases=`` to ``sim.L_eff`` or
    ``sim.solver.solve``; a dict in the list is passed through as is.

    The dataclasses are numbered by their **position** in the list, whatever each object
    carries. The mean-field schemes select the matrix with ``number == n_matrix`` and
    then index ``sub_phases[n_matrix]``, so the number has to be the position; the
    dataclass default is 0 for every phase, which would otherwise make every
    concentration tensor the identity and turn L_eff into a Voigt average.
    """
    return [p if isinstance(p, dict) else to_phase_dict(p, number=i)
            for i, p in enumerate(phases)]


# =============================================================================
# Effective stiffness
# =============================================================================

def L_eff(umat_name: str, props, nstatev: int, orientation: Orientation = None,
          phases: Optional[List] = None) -> np.ndarray:
    """Elastic stiffness tensor (6x6, Voigt) of a material in the global frame.

    ``orientation`` is the material frame, in any form :func:`as_rotation` accepts
    (a :class:`simcoon.Rotation`, ``(psi, theta, phi)`` in degrees, the JSON dict;
    ``None`` for the identity): the stiffness is ``R.apply_stiffness(L_local)``. For a
    mean-field model (MIHEN, MIMTN, MISCN, MIPLN) ``phases`` gives its sub-phases, as
    the dataclasses of this module or the dicts :func:`to_phase_dicts` makes of them.
    This is ``sim.L_eff``; ``sim._core.L_eff`` is the extension entry it wraps.
    """
    angles = None if orientation is None else euler_angles(orientation)
    dicts = None if phases is None else to_phase_dicts(phases)
    return np.asarray(_core.L_eff(str(umat_name), np.asarray(props, dtype=float).ravel(),
                                  int(nstatev), angles, dicts))


# =============================================================================
# Exports
# =============================================================================

__all__ = [
    # Orientations
    'EULER_SEQ',
    'as_rotation',
    'euler_angles',
    # Data classes
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
    # Handing phases to the C++ side, in memory, and back
    'to_phase_dict',
    'to_phase_dicts',
    'phases_from_dicts',
    # Orientation distribution functions
    'Peak',
    'get_densities_ODF',
    'load_peaks_json',
    'save_peaks_json',
    'discretize_odf',
    # Effective stiffness
    'L_eff',
]
