"""
Micromechanics data classes and JSON I/O for Simcoon.

Dataclasses and I/O functions describing the sub-phases of a mean-field model (phases,
layers, ellipsoids, ...), and their in-memory form for ``sim.L_eff`` and
``sim.solver.solve``. Nothing here calls the C++ extension, but the package that hosts
the module does import it: ``simcoon`` is not importable without ``simcoon._core``.

Orientations (material frame of a phase, geometry of an inclusion or a layer) are
``simcoon.Rotation`` objects; ``as_rotation`` also takes the Euler angles of the files.

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
    An orientation as a simcoon.Rotation, and back to the (psi, theta, phi) of the files
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

import json
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Union, Optional

import numpy as np

from simcoon import _core
from simcoon._core import get_densities_ODF
from simcoon.rotation import Rotation


# =============================================================================
# Data Classes
# =============================================================================

_ANGLES = ('psi', 'theta', 'phi')

#: The Euler convention of the C++ side. ``Rotation::from_euler(psi, theta, phi, "zxz")``
#: composes the three axis rotations as scipy's *extrinsic* ``'zxz'`` does, and the
#: solver applies it actively (material frame -> global frame): a phase at
#: ``(psi, theta, phi)`` responds with ``R.apply_stiffness(L_local)``,
#: ``R = Rotation.from_euler('zxz', [psi, theta, phi], degrees=True)``. Pinned by
#: test_micromechanics.py::TestOrientationConvention against the solver and L_eff.
EULER_SEQ = 'zxz'

Orientation = Union[Rotation, Dict[str, float], "Sequence[float]", None]


def as_rotation(value: Orientation) -> Rotation:
    """An orientation as a :class:`simcoon.Rotation`.

    ``value`` is a ``Rotation`` (returned as is), the Euler angles ``(psi, theta, phi)``
    in degrees as a 3-sequence or as the ``{"psi", "theta", "phi"}`` dict of the JSON
    files (missing angles are 0), or ``None`` for the identity. The angles are the
    ``'zxz'`` Euler angles the C++ side reads (see ``EULER_SEQ``).
    """
    if value is None:
        return Rotation.identity()
    if isinstance(value, Rotation):
        return value
    if isinstance(value, dict):
        unknown = set(value) - set(_ANGLES)
        if unknown:
            raise ValueError(f"orientation: unknown keys {sorted(unknown)}; expected {_ANGLES}")
        angles = [float(value.get(k, 0.0)) for k in _ANGLES]
    else:
        angles = np.asarray(value, dtype=float).ravel()
        if angles.size != 3:
            raise ValueError(f"orientation: 3 Euler angles (psi, theta, phi) in degrees "
                             f"expected, got {angles.size} values")
    return Rotation.from_euler(EULER_SEQ, angles, degrees=True)


def euler_angles(rotation: Orientation) -> Dict[str, float]:
    """The ``{"psi", "theta", "phi"}`` dict (degrees, ``EULER_SEQ``) of an orientation:
    the form of the JSON files and of the dicts the C++ binding reads.

    The decomposition is not unique when ``theta`` is 0 or 180 degrees (gimbal lock):
    scipy then puts the whole z rotation in ``psi`` and sets ``phi`` to 0, which is the
    same rotation as the angles that were given, written differently.
    """
    rot = as_rotation(rotation)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')   # scipy's "Gimbal lock detected" — see above
        psi, theta, phi = rot.as_euler(EULER_SEQ, degrees=True)
    return {'psi': float(psi) + 0.0, 'theta': float(theta) + 0.0, 'phi': float(phi) + 0.0}


def _coerce_fields(obj):
    """props given as a list, orientations in any accepted form, nested phases given as
    dicts (the JSON form)."""
    if isinstance(obj.props, list):
        obj.props = np.array(obj.props, dtype=float)
    for name in ('material_orientation', 'geometry_orientation'):
        if hasattr(obj, name):
            setattr(obj, name, as_rotation(getattr(obj, name)))
    nested = getattr(obj, 'phases', None)
    if nested and any(isinstance(p, dict) for p in nested):
        cls, layout = _JSON_LAYOUTS[_json_kind_of(obj.umat_name)]
        obj.phases = [p if not isinstance(p, dict)
                      else _from_json_entry(p, cls, layout, None, f"phases of {obj.umat_name}")
                      for p in nested]


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
    material_orientation : Rotation
        Orientation of the material frame (any form ``as_rotation`` accepts)
    nstatev : int
        Number of state variables
    props : np.ndarray
        Material properties array
    phases : list
        Sub-phases of a phase that is itself a mean-field model (MIHEN, MIMTN, MISCN:
        ellipsoids; MIPLN: layers), the way Nellipsoids<N>.dat chained through props[1].
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


@dataclass
class Layer(Phase):
    """
    Layer phase for laminate homogenization.

    Corresponds to Nlayers.dat format and C++ layer class.
    Layers are oriented using geometry orientation angles.

    Additional Attributes
    ---------------------
    geometry_orientation : Rotation
        Orientation of the geometry (any form ``as_rotation`` accepts)
    layerup : int
        Index of layer above (-1 if none)
    layerdown : int
        Index of layer below (-1 if none)
    """
    geometry_orientation: Rotation = field(default_factory=Rotation.identity)
    layerup: int = 0
    layerdown: int = 0


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


# =============================================================================
# JSON I/O
# =============================================================================

def _props_to_dict(props: np.ndarray, prop_names: List[str] = None) -> Dict[str, float]:
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
    return [_from_json_entry(entry, cls, layout, prop_names, str(filepath))
            for entry in data.get(kind, [])]


def _save_json(filepath: Union[str, Path], kind: str, items: List,
               prop_names: Optional[List[str]]):
    _, layout = _JSON_LAYOUTS[kind]
    with open(filepath, 'w') as f:
        json.dump({kind: [_to_json_entry(item, layout, prop_names) for item in items]},
                  f, indent=2)


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
                     prop_names: List[str] = None):
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
                     prop_names: List[str] = None):
    """Save layers to a JSON file (the layout ``load_layers_json`` reads)."""
    _save_json(filepath, 'layers', layers, prop_names)


def load_ellipsoids_json(filepath: Union[str, Path],
                         prop_names: Optional[List[str]] = None) -> List[Ellipsoid]:
    """Load ellipsoids from a JSON file: the phase entry plus ``"coatingof"``,
    ``"semi_axes": {"a1", "a2", "a3"}`` and ``"geometry_orientation"``, under the
    top-level key ``"ellipsoids"``."""
    return _load_json(filepath, 'ellipsoids', prop_names)


def save_ellipsoids_json(filepath: Union[str, Path], ellipsoids: List[Ellipsoid],
                         prop_names: List[str] = None):
    """Save ellipsoids to a JSON file (the layout ``load_ellipsoids_json`` reads)."""
    _save_json(filepath, 'ellipsoids', ellipsoids, prop_names)


def load_cylinders_json(filepath: Union[str, Path],
                        prop_names: Optional[List[str]] = None) -> List[Cylinder]:
    """Load cylinders from a JSON file: the phase entry plus ``"coatingof"``,
    ``"geometry": {"L", "R"}`` and ``"geometry_orientation"``, under the top-level key
    ``"cylinders"``."""
    return _load_json(filepath, 'cylinders', prop_names)


def save_cylinders_json(filepath: Union[str, Path], cylinders: List[Cylinder],
                        prop_names: List[str] = None):
    """Save cylinders to a JSON file (the layout ``load_cylinders_json`` reads)."""
    _save_json(filepath, 'cylinders', cylinders, prop_names)


def load_sections_json(filepath: Union[str, Path],
                       prop_names: Optional[List[str]] = None) -> List[Section]:
    """Load sections from a JSON file for textile composites: ``"number"``, ``"name"``,
    ``"umat_name"``, ``"material_orientation"``, ``"nstatev"``, ``"props"``, under the
    top-level key ``"sections"``."""
    return _load_json(filepath, 'sections', prop_names)


def save_sections_json(filepath: Union[str, Path], sections: List[Section],
                       prop_names: List[str] = None):
    """Save sections to a JSON file (the layout ``load_sections_json`` reads)."""
    _save_json(filepath, 'sections', sections, prop_names)


# =============================================================================
# Orientation distribution functions (ODF): the peaks, in memory
# =============================================================================

@dataclass
class Peak:
    """One peak of an orientation (ODF) or parameter (PDF) distribution.

    ``method`` selects the profile the C++ side evaluates: 1 standard deviation kernel
    (``params``), 2 hard cut-off, 3 Gaussian, 4 Lorentzian, 5 pseudo-Voigt (``params``),
    6 Pearson VII (``params``), 7 uniform. Angles (``mean``, ``s_dev``, ``width``) are
    degrees, as in the files.
    """
    number: int = 0
    method: int = 3
    mean: float = 0.0
    s_dev: float = 1.0
    width: float = 1.0
    ampl: float = 1.0
    params: np.ndarray = field(default_factory=lambda: np.array([]))

    def __post_init__(self):
        if isinstance(self.params, list):
            self.params = np.array(self.params, dtype=float)


_PEAK_FIELDS = ('number', 'method', 'mean', 's_dev', 'width', 'ampl', 'params')


def to_peak_dicts(peaks: List[Peak]) -> List[Dict]:
    """The peaks as ``sim.get_densities_ODF`` reads them."""
    return [{f: (np.asarray(getattr(p, f), dtype=float) if f == 'params' else getattr(p, f))
             for f in _PEAK_FIELDS} for p in peaks]


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

    The Euler-angle sweeps of the pre-2.0 files map onto this: a sweep of ``psi`` or
    ``phi`` from an unrotated phase is ``axis=(0, 0, 1)``, a sweep of ``theta`` is
    ``axis=(1, 0, 0)``.

    Returns a new list; the phases are renumbered by position, the others untouched.
    """
    import copy
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
    peak_dicts = [pk if isinstance(pk, dict) else to_peak_dicts([pk])[0] for pk in peaks]

    parent = phases[num_phase]
    d = (a_max - a_min) / nphases
    alphas = a_min + d * np.arange(nphases)
    # Simpson over each bin, the density being periodic over 180 deg (director distribution)
    x = np.concatenate([alphas - d / 2, alphas, alphas + d / 2]) % 180.0
    rho = np.asarray(get_densities_ODF(x, peak_dicts, False)).ravel()
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
    out = list(phases[:num_phase]) + swept + list(phases[num_phase + 1:])
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

    Angles are left in degrees, as in the files; the binding converts them to radians.
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
    'to_peak_dicts',
    'load_peaks_json',
    'save_peaks_json',
    'discretize_odf',
    # Effective stiffness
    'L_eff',
]
