"""
Orientation Distribution Function Module
========================================

Tools for working with Orientation Distribution Functions (ODF).

This module provides:
- ODF density evaluation
- ODF discretization for homogenization
- Peak file creation and management
- Visualization utilities

Examples
--------
>>> import numpy as np
>>> from simcoon.odf import get_densities, discretize, create_peaks_file, ODFPeak
>>>
>>> # Create ODF with two fiber families at +/-45 degrees
>>> peaks = [ODFPeak(angle=45, intensity=0.5, width=10),
...          ODFPeak(angle=-45, intensity=0.5, width=10)]
>>> create_peaks_file(peaks, "peaks.json")
>>>
>>> # Evaluate ODF density
>>> angles = np.linspace(-90, 90, 181)
>>> densities = get_densities(angles, "peaks.json", path=".")
"""

from dataclasses import dataclass, field
from typing import List, Optional, Dict, Any
import numpy as np
import json

import simcoon._core as _sim


# =============================================================================
# Data Classes
# =============================================================================

@dataclass
class ODFPeak:
    """
    Single ODF peak definition.

    Parameters
    ----------
    angle : float
        Peak position in degrees
    intensity : float
        Peak intensity/weight (should sum to 1 for normalized ODF)
    width : float
        Peak width (FWHM in degrees)
    distribution : str
        Distribution type: 'gaussian', 'lorentzian', 'pseudo-voigt'
    """
    angle: float
    intensity: float
    width: float
    distribution: str = 'gaussian'


@dataclass
class ODF:
    """
    Orientation Distribution Function.

    Parameters
    ----------
    peaks : list of ODFPeak
        List of peaks defining the ODF
    angle_min : float
        Minimum angle (degrees)
    angle_max : float
        Maximum angle (degrees)
    symmetric : bool
        If True, ODF is symmetric about 0
    """
    peaks: List[ODFPeak] = field(default_factory=list)
    angle_min: float = -90.0
    angle_max: float = 90.0
    symmetric: bool = False

    def add_peak(self, angle: float, intensity: float, width: float,
                 distribution: str = 'gaussian'):
        """Add a peak to the ODF."""
        self.peaks.append(ODFPeak(angle, intensity, width, distribution))

    def normalize(self):
        """Normalize peak intensities to sum to 1."""
        total = sum(p.intensity for p in self.peaks)
        if total > 0:
            for p in self.peaks:
                p.intensity /= total

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        return {
            'angle_min': self.angle_min,
            'angle_max': self.angle_max,
            'symmetric': self.symmetric,
            'peaks': [
                {
                    'angle': p.angle,
                    'intensity': p.intensity,
                    'width': p.width,
                    'distribution': p.distribution
                }
                for p in self.peaks
            ]
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'ODF':
        """Create ODF from dictionary."""
        odf = cls(
            angle_min=data.get('angle_min', -90.0),
            angle_max=data.get('angle_max', 90.0),
            symmetric=data.get('symmetric', False)
        )
        for p in data.get('peaks', []):
            odf.add_peak(
                p['angle'], p['intensity'], p['width'],
                p.get('distribution', 'gaussian')
            )
        return odf


# =============================================================================
# Core Functions (wrapping C++ bindings)
# =============================================================================

def get_densities(angles: np.ndarray, peaks_file: str,
                  path: str = ".", radian: bool = False) -> np.ndarray:
    """
    Compute ODF density values at given angles.

    Uses the C++ get_densities_ODF function.

    Parameters
    ----------
    angles : ndarray
        Angles at which to evaluate the ODF
    peaks_file : str
        Name of the peaks definition file (JSON format)
    path : str
        Directory containing the peaks file
    radian : bool
        If True, angles are in radians; otherwise degrees

    Returns
    -------
    ndarray
        ODF density values at each angle

    Examples
    --------
    >>> angles = np.linspace(-90, 90, 181)
    >>> densities = get_densities(angles, "peaks.json", path="data")
    """
    return _sim.get_densities_ODF(np.asarray(angles), path, peaks_file, radian)


def discretize(nphases: int, num_phase: int,
               angle_min: float, angle_max: float,
               umat_name: str, props: np.ndarray,
               peaks_file: str, rve_init_file: str, rve_output_file: str,
               path: str = ".", angle_mat: int = 0):
    """
    Discretize an ODF into distinct phases for homogenization.

    Uses the C++ ODF_discretization function.

    Parameters
    ----------
    nphases : int
        Number of phases in the discretized RVE
    num_phase : int
        Phase number to discretize (0-indexed)
    angle_min : float
        Minimum angle for discretization (degrees)
    angle_max : float
        Maximum angle for discretization (degrees)
    umat_name : str
        Material model name for the phases (e.g., 'MIHEN', 'MIMTN')
    props : ndarray
        Material properties array
    peaks_file : str
        Name of the ODF peaks file
    rve_init_file : str
        Input RVE JSON configuration file
    rve_output_file : str
        Output discretized RVE JSON file
    path : str
        Directory for input/output files
    angle_mat : int
        Material angle index to vary (0=psi, 1=theta, 2=phi)

    Examples
    --------
    >>> props = np.array([...])  # Material properties
    >>> discretize(
    ...     nphases=10, num_phase=0,
    ...     angle_min=-90, angle_max=90,
    ...     umat_name="MIHEN", props=props,
    ...     peaks_file="peaks.json",
    ...     rve_init_file="rve_init.json",
    ...     rve_output_file="rve_disc.json",
    ...     path="data"
    ... )
    """
    _sim.ODF_discretization(
        nphases, num_phase, angle_min, angle_max,
        umat_name, np.asarray(props), path, peaks_file,
        rve_init_file, rve_output_file, angle_mat
    )


# =============================================================================
# File I/O
# =============================================================================

def create_peaks_file(peaks: List[ODFPeak], filename: str,
                      angle_min: float = -90.0, angle_max: float = 90.0):
    """
    Create a peaks JSON file from ODFPeak objects.

    Parameters
    ----------
    peaks : list of ODFPeak
        List of peak definitions
    filename : str
        Output filename (JSON format)
    angle_min : float
        Minimum angle
    angle_max : float
        Maximum angle

    Examples
    --------
    >>> peaks = [
    ...     ODFPeak(angle=0, intensity=0.6, width=15),
    ...     ODFPeak(angle=90, intensity=0.4, width=20),
    ... ]
    >>> create_peaks_file(peaks, "my_odf.json")
    """
    data = {
        'angle_min': angle_min,
        'angle_max': angle_max,
        'peaks': [
            {
                'angle': p.angle,
                'intensity': p.intensity,
                'width': p.width,
                'distribution': p.distribution
            }
            for p in peaks
        ]
    }
    with open(filename, 'w') as f:
        json.dump(data, f, indent=2)


def load_peaks_file(filename: str) -> ODF:
    """
    Load an ODF from a peaks JSON file.

    Parameters
    ----------
    filename : str
        Path to the peaks file

    Returns
    -------
    ODF
        ODF object with loaded peaks
    """
    with open(filename, 'r') as f:
        data = json.load(f)
    return ODF.from_dict(data)


def save_odf(odf: ODF, filename: str):
    """
    Save an ODF to a JSON file.

    Parameters
    ----------
    odf : ODF
        ODF object to save
    filename : str
        Output filename
    """
    with open(filename, 'w') as f:
        json.dump(odf.to_dict(), f, indent=2)


# =============================================================================
# ODF Generation Utilities
# =============================================================================

def random_odf(n_peaks: int, angle_range: tuple = (-90, 90),
               width_range: tuple = (5, 20), seed: int = None) -> ODF:
    """
    Generate a random ODF with specified number of peaks.

    Parameters
    ----------
    n_peaks : int
        Number of peaks
    angle_range : tuple
        (min, max) angle range for peak positions
    width_range : tuple
        (min, max) range for peak widths
    seed : int, optional
        Random seed for reproducibility

    Returns
    -------
    ODF
        Random ODF object (normalized)
    """
    if seed is not None:
        np.random.seed(seed)

    odf = ODF(angle_min=angle_range[0], angle_max=angle_range[1])

    angles = np.random.uniform(angle_range[0], angle_range[1], n_peaks)
    widths = np.random.uniform(width_range[0], width_range[1], n_peaks)
    intensities = np.random.uniform(0.1, 1.0, n_peaks)

    for angle, width, intensity in zip(angles, widths, intensities):
        odf.add_peak(angle, intensity, width)

    odf.normalize()
    return odf


def uniform_odf(angle_min: float = -90, angle_max: float = 90) -> ODF:
    """
    Create a uniform (isotropic) ODF.

    Parameters
    ----------
    angle_min : float
        Minimum angle
    angle_max : float
        Maximum angle

    Returns
    -------
    ODF
        Uniform ODF (single very wide peak)
    """
    odf = ODF(angle_min=angle_min, angle_max=angle_max)
    center = (angle_min + angle_max) / 2
    width = (angle_max - angle_min) * 2  # Very wide to approximate uniform
    odf.add_peak(center, 1.0, width)
    return odf


def fiber_odf(fiber_angles: List[float], intensities: List[float] = None,
              width: float = 10.0) -> ODF:
    """
    Create an ODF for fiber-reinforced materials.

    Parameters
    ----------
    fiber_angles : list of float
        Fiber orientation angles (degrees)
    intensities : list of float, optional
        Volume fractions for each fiber family (default: equal)
    width : float
        Peak width for all fibers

    Returns
    -------
    ODF
        Fiber ODF

    Examples
    --------
    >>> # Cross-ply laminate (+/-45)
    >>> odf = fiber_odf([45, -45])
    >>>
    >>> # Quasi-isotropic (0/45/90/-45)
    >>> odf = fiber_odf([0, 45, 90, -45])
    """
    n = len(fiber_angles)
    if intensities is None:
        intensities = [1.0 / n] * n

    odf = ODF(angle_min=-90, angle_max=90)
    for angle, intensity in zip(fiber_angles, intensities):
        odf.add_peak(angle, intensity, width)

    return odf


# =============================================================================
# Visualization
# =============================================================================

def plot_odf(angles: np.ndarray, densities: np.ndarray,
             polar: bool = False, ax=None, **kwargs):
    """
    Plot an ODF.

    Parameters
    ----------
    angles : ndarray
        Angles (degrees)
    densities : ndarray
        ODF density values
    polar : bool
        If True, create polar plot
    ax : matplotlib.axes.Axes, optional
        Axes to plot on
    **kwargs
        Additional arguments passed to plot

    Returns
    -------
    ax : matplotlib.axes.Axes
        The axes object
    """
    import matplotlib.pyplot as plt

    if ax is None:
        if polar:
            fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
        else:
            fig, ax = plt.subplots()

    if polar:
        ax.plot(np.radians(angles), densities, **kwargs)
        ax.set_theta_zero_location('E')
        ax.set_theta_direction(1)
    else:
        ax.plot(angles, densities, **kwargs)
        ax.set_xlabel('Angle (degrees)')
        ax.set_ylabel('ODF Density')

    return ax


def plot_odf_polar_2d(angles: np.ndarray, densities: np.ndarray,
                      ax=None, fill: bool = True, **kwargs):
    """
    Create a 2D polar plot of ODF (rose diagram).

    Parameters
    ----------
    angles : ndarray
        Angles in degrees
    densities : ndarray
        ODF density values
    ax : matplotlib.axes.Axes, optional
        Polar axes to plot on
    fill : bool
        If True, fill the area under the curve
    **kwargs
        Additional arguments passed to plot/fill

    Returns
    -------
    ax : matplotlib.axes.Axes
    """
    import matplotlib.pyplot as plt

    if ax is None:
        fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})

    theta = np.radians(angles)

    # Close the curve
    theta = np.append(theta, theta[0])
    r = np.append(densities, densities[0])

    if fill:
        ax.fill(theta, r, alpha=0.3, **kwargs)
    ax.plot(theta, r, **kwargs)

    ax.set_theta_zero_location('E')
    ax.set_theta_direction(1)

    return ax


# =============================================================================
# Convenience Exports
# =============================================================================

__all__ = [
    # Data classes
    'ODFPeak',
    'ODF',
    # Core functions
    'get_densities',
    'discretize',
    # File I/O
    'create_peaks_file',
    'load_peaks_file',
    'save_odf',
    # Utilities
    'random_odf',
    'uniform_odf',
    'fiber_odf',
    # Visualization
    'plot_odf',
    'plot_odf_polar_2d',
]
