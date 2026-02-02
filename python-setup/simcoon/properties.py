"""
Elastic Properties Module
=========================

Tools for computing and analyzing elastic properties of materials.

This module wraps the C++ constitutive and recovery_props functions
to provide a convenient Python interface for:

- Building stiffness/compliance tensors for various symmetries
- Extracting engineering constants from tensors
- Computing effective properties of composites
- Analyzing directional elastic properties

Examples
--------
>>> import numpy as np
>>> from simcoon.properties import IsotropicMaterial, effective_properties
>>>
>>> # Create isotropic material and get properties
>>> mat = IsotropicMaterial(E=210000, nu=0.3)
>>> print(mat.L)  # Stiffness tensor
>>> print(mat.bulk_modulus, mat.shear_modulus)
>>>
>>> # Compute effective properties of a composite
>>> props = effective_properties("MIHEN", composite_props, nstatev=1)
>>> print(props)
"""

from dataclasses import dataclass
from typing import Optional, Tuple, Dict, Union
import numpy as np

import simcoon._core as _sim


# =============================================================================
# Material Classes
# =============================================================================

@dataclass
class IsotropicMaterial:
    """
    Isotropic elastic material.

    Parameters
    ----------
    E : float
        Young's modulus
    nu : float
        Poisson's ratio

    Attributes
    ----------
    L : ndarray
        6x6 stiffness tensor
    M : ndarray
        6x6 compliance tensor
    bulk_modulus : float
        Bulk modulus K
    shear_modulus : float
        Shear modulus G
    lame_lambda : float
        First Lame parameter
    """
    E: float
    nu: float

    def __post_init__(self):
        self._L = None
        self._M = None

    @property
    def L(self) -> np.ndarray:
        """Stiffness tensor (6x6)."""
        if self._L is None:
            self._L = _sim.L_iso(np.array([self.E, self.nu]), 'Enu')
        return self._L

    @property
    def M(self) -> np.ndarray:
        """Compliance tensor (6x6)."""
        if self._M is None:
            self._M = _sim.M_iso(np.array([self.E, self.nu]), 'Enu')
        return self._M

    @property
    def bulk_modulus(self) -> float:
        """Bulk modulus K."""
        return self.E / (3 * (1 - 2 * self.nu))

    @property
    def shear_modulus(self) -> float:
        """Shear modulus G."""
        return self.E / (2 * (1 + self.nu))

    @property
    def lame_lambda(self) -> float:
        """First Lame parameter."""
        return self.E * self.nu / ((1 + self.nu) * (1 - 2 * self.nu))

    @classmethod
    def from_stiffness(cls, L: np.ndarray) -> 'IsotropicMaterial':
        """Create from stiffness tensor using recovery_props."""
        props = _sim.L_iso_props(L)
        return cls(E=props[0], nu=props[1])

    @classmethod
    def from_compliance(cls, M: np.ndarray) -> 'IsotropicMaterial':
        """Create from compliance tensor using recovery_props."""
        props = _sim.M_iso_props(M)
        return cls(E=props[0], nu=props[1])

    @classmethod
    def from_bulk_shear(cls, K: float, G: float) -> 'IsotropicMaterial':
        """Create from bulk and shear moduli."""
        E = 9 * K * G / (3 * K + G)
        nu = (3 * K - 2 * G) / (2 * (3 * K + G))
        return cls(E=E, nu=nu)

    def to_dict(self) -> Dict[str, float]:
        """Return all properties as dictionary."""
        return {
            'E': self.E,
            'nu': self.nu,
            'K': self.bulk_modulus,
            'G': self.shear_modulus,
            'lambda': self.lame_lambda,
        }


@dataclass
class TransverselyIsotropicMaterial:
    """
    Transversely isotropic elastic material.

    Parameters
    ----------
    EL : float
        Longitudinal Young's modulus (along symmetry axis)
    ET : float
        Transverse Young's modulus
    nuTL : float
        Poisson's ratio for transverse strain under longitudinal stress
    nuTT : float
        Poisson's ratio in transverse plane
    GLT : float
        Longitudinal-transverse shear modulus
    axis : int
        Symmetry axis (1, 2, or 3)
    """
    EL: float
    ET: float
    nuTL: float
    nuTT: float
    GLT: float
    axis: int = 1

    def __post_init__(self):
        self._L = None
        self._M = None

    @property
    def L(self) -> np.ndarray:
        """Stiffness tensor (6x6)."""
        if self._L is None:
            props = np.array([self.EL, self.ET, self.nuTL, self.nuTT, self.GLT])
            self._L = _sim.L_isotrans(props, self.axis)
        return self._L

    @property
    def M(self) -> np.ndarray:
        """Compliance tensor (6x6)."""
        if self._M is None:
            props = np.array([self.EL, self.ET, self.nuTL, self.nuTT, self.GLT])
            self._M = _sim.M_isotrans(props, self.axis)
        return self._M

    @property
    def GTT(self) -> float:
        """Transverse shear modulus."""
        return self.ET / (2 * (1 + self.nuTT))

    @classmethod
    def from_stiffness(cls, L: np.ndarray, axis: int = 1) -> 'TransverselyIsotropicMaterial':
        """Create from stiffness tensor using recovery_props."""
        props = _sim.L_isotrans_props(L, axis)
        return cls(EL=props[0], ET=props[1], nuTL=props[2], nuTT=props[3], GLT=props[4], axis=axis)

    @classmethod
    def from_compliance(cls, M: np.ndarray, axis: int = 1) -> 'TransverselyIsotropicMaterial':
        """Create from compliance tensor using recovery_props."""
        props = _sim.M_isotrans_props(M, axis)
        return cls(EL=props[0], ET=props[1], nuTL=props[2], nuTT=props[3], GLT=props[4], axis=axis)

    def to_dict(self) -> Dict[str, float]:
        """Return all properties as dictionary."""
        return {
            'EL': self.EL,
            'ET': self.ET,
            'nuTL': self.nuTL,
            'nuTT': self.nuTT,
            'GLT': self.GLT,
            'GTT': self.GTT,
            'axis': self.axis,
        }


@dataclass
class OrthotropicMaterial:
    """
    Orthotropic elastic material.

    Parameters
    ----------
    E1, E2, E3 : float
        Young's moduli in principal directions
    nu12, nu13, nu23 : float
        Poisson's ratios
    G12, G13, G23 : float
        Shear moduli
    """
    E1: float
    E2: float
    E3: float
    nu12: float
    nu13: float
    nu23: float
    G12: float
    G13: float
    G23: float

    def __post_init__(self):
        self._L = None
        self._M = None

    @property
    def L(self) -> np.ndarray:
        """Stiffness tensor (6x6)."""
        if self._L is None:
            props = np.array([self.E1, self.E2, self.E3,
                              self.nu12, self.nu13, self.nu23,
                              self.G12, self.G13, self.G23])
            self._L = _sim.L_ortho(props, 'EnuG')
        return self._L

    @property
    def M(self) -> np.ndarray:
        """Compliance tensor (6x6)."""
        if self._M is None:
            props = np.array([self.E1, self.E2, self.E3,
                              self.nu12, self.nu13, self.nu23,
                              self.G12, self.G13, self.G23])
            self._M = _sim.M_ortho(props, 'EnuG')
        return self._M

    @classmethod
    def from_stiffness(cls, L: np.ndarray) -> 'OrthotropicMaterial':
        """Create from stiffness tensor using recovery_props."""
        props = _sim.L_ortho_props(L)
        return cls(
            E1=props[0], E2=props[1], E3=props[2],
            nu12=props[3], nu13=props[4], nu23=props[5],
            G12=props[6], G13=props[7], G23=props[8]
        )

    @classmethod
    def from_compliance(cls, M: np.ndarray) -> 'OrthotropicMaterial':
        """Create from compliance tensor using recovery_props."""
        props = _sim.M_ortho_props(M)
        return cls(
            E1=props[0], E2=props[1], E3=props[2],
            nu12=props[3], nu13=props[4], nu23=props[5],
            G12=props[6], G13=props[7], G23=props[8]
        )

    def to_dict(self) -> Dict[str, float]:
        """Return all properties as dictionary."""
        return {
            'E1': self.E1, 'E2': self.E2, 'E3': self.E3,
            'nu12': self.nu12, 'nu13': self.nu13, 'nu23': self.nu23,
            'G12': self.G12, 'G13': self.G13, 'G23': self.G23,
        }


@dataclass
class CubicMaterial:
    """
    Cubic elastic material.

    Parameters
    ----------
    E : float
        Young's modulus
    nu : float
        Poisson's ratio
    G : float
        Shear modulus (independent for cubic symmetry)
    """
    E: float
    nu: float
    G: float

    def __post_init__(self):
        self._L = None
        self._M = None

    @property
    def L(self) -> np.ndarray:
        """Stiffness tensor (6x6)."""
        if self._L is None:
            props = np.array([self.E, self.nu, self.G])
            self._L = _sim.L_cubic(props, 'EnuG')
        return self._L

    @property
    def M(self) -> np.ndarray:
        """Compliance tensor (6x6)."""
        if self._M is None:
            props = np.array([self.E, self.nu, self.G])
            self._M = _sim.M_cubic(props, 'EnuG')
        return self._M

    @property
    def zener_ratio(self) -> float:
        """Zener anisotropy ratio A = 2*G*(1+nu)/E."""
        return 2 * self.G * (1 + self.nu) / self.E

    @classmethod
    def from_stiffness(cls, L: np.ndarray) -> 'CubicMaterial':
        """Create from stiffness tensor using recovery_props."""
        props = _sim.L_cubic_props(L)
        return cls(E=props[0], nu=props[1], G=props[2])

    @classmethod
    def from_compliance(cls, M: np.ndarray) -> 'CubicMaterial':
        """Create from compliance tensor using recovery_props."""
        props = _sim.M_cubic_props(M)
        return cls(E=props[0], nu=props[1], G=props[2])

    def to_dict(self) -> Dict[str, float]:
        """Return all properties as dictionary."""
        return {
            'E': self.E,
            'nu': self.nu,
            'G': self.G,
            'A': self.zener_ratio,
        }


# =============================================================================
# Property Recovery Functions
# =============================================================================

def recover_isotropic(L: np.ndarray = None, M: np.ndarray = None) -> Dict[str, float]:
    """
    Recover isotropic elastic constants from stiffness or compliance tensor.

    Parameters
    ----------
    L : ndarray, optional
        6x6 stiffness tensor
    M : ndarray, optional
        6x6 compliance tensor

    Returns
    -------
    dict
        Dictionary with E (Young's modulus) and nu (Poisson's ratio)
    """
    if L is not None:
        props = _sim.L_iso_props(L)
    elif M is not None:
        props = _sim.M_iso_props(M)
    else:
        raise ValueError("Either L or M must be provided")

    return {'E': props[0], 'nu': props[1]}


def recover_transversely_isotropic(L: np.ndarray = None, M: np.ndarray = None,
                                    axis: int = 1) -> Dict[str, float]:
    """
    Recover transversely isotropic constants from stiffness or compliance tensor.

    Parameters
    ----------
    L : ndarray, optional
        6x6 stiffness tensor
    M : ndarray, optional
        6x6 compliance tensor
    axis : int
        Symmetry axis (1, 2, or 3)

    Returns
    -------
    dict
        Dictionary with EL, ET, nuTL, nuTT, GLT
    """
    if L is not None:
        props = _sim.L_isotrans_props(L, axis)
    elif M is not None:
        props = _sim.M_isotrans_props(M, axis)
    else:
        raise ValueError("Either L or M must be provided")

    return {
        'EL': props[0], 'ET': props[1],
        'nuTL': props[2], 'nuTT': props[3],
        'GLT': props[4]
    }


def recover_orthotropic(L: np.ndarray = None, M: np.ndarray = None) -> Dict[str, float]:
    """
    Recover orthotropic elastic constants from stiffness or compliance tensor.

    Parameters
    ----------
    L : ndarray, optional
        6x6 stiffness tensor
    M : ndarray, optional
        6x6 compliance tensor

    Returns
    -------
    dict
        Dictionary with E1, E2, E3, nu12, nu13, nu23, G12, G13, G23
    """
    if L is not None:
        props = _sim.L_ortho_props(L)
    elif M is not None:
        props = _sim.M_ortho_props(M)
    else:
        raise ValueError("Either L or M must be provided")

    return {
        'E1': props[0], 'E2': props[1], 'E3': props[2],
        'nu12': props[3], 'nu13': props[4], 'nu23': props[5],
        'G12': props[6], 'G13': props[7], 'G23': props[8]
    }


def recover_cubic(L: np.ndarray = None, M: np.ndarray = None) -> Dict[str, float]:
    """
    Recover cubic elastic constants from stiffness or compliance tensor.

    Parameters
    ----------
    L : ndarray, optional
        6x6 stiffness tensor
    M : ndarray, optional
        6x6 compliance tensor

    Returns
    -------
    dict
        Dictionary with E, nu, G
    """
    if L is not None:
        props = _sim.L_cubic_props(L)
    elif M is not None:
        props = _sim.M_cubic_props(M)
    else:
        raise ValueError("Either L or M must be provided")

    return {'E': props[0], 'nu': props[1], 'G': props[2]}


def recover_anisotropic(M: np.ndarray) -> Dict[str, float]:
    """
    Recover engineering constants from fully anisotropic compliance tensor.

    Parameters
    ----------
    M : ndarray
        6x6 compliance tensor

    Returns
    -------
    dict
        Dictionary with E1, E2, E3, nu12, nu13, nu23, G12, G13, G23,
        and coupling coefficients eta_ij
    """
    props = _sim.M_aniso_props(M)

    return {
        'E1': props[0], 'E2': props[1], 'E3': props[2],
        'nu12': props[3], 'nu13': props[4], 'nu23': props[5],
        'G12': props[6], 'G13': props[7], 'G23': props[8],
        'eta14': props[9], 'eta15': props[10], 'eta16': props[11],
        'eta24': props[12], 'eta25': props[13], 'eta26': props[14],
        'eta34': props[15], 'eta35': props[16], 'eta36': props[17],
        'eta45': props[18], 'eta46': props[19], 'eta56': props[20],
    }


# =============================================================================
# Effective Properties (Homogenization)
# =============================================================================

def effective_stiffness(umat_name: str, props: np.ndarray, nstatev: int = 1,
                        psi: float = 0., theta: float = 0., phi: float = 0.) -> np.ndarray:
    """
    Compute effective stiffness tensor for a composite material.

    Uses the C++ L_eff function for homogenization schemes like
    Mori-Tanaka, Self-Consistent, etc.

    Parameters
    ----------
    umat_name : str
        Homogenization scheme name (MIHEN, MIMTN, MISCN, etc.)
    props : ndarray
        Material properties array (phase properties, volume fractions, etc.)
    nstatev : int
        Number of state variables
    psi, theta, phi : float
        Euler angles for RVE orientation (degrees)

    Returns
    -------
    ndarray
        6x6 effective stiffness tensor
    """
    return _sim.L_eff(umat_name, props, nstatev, psi, theta, phi)


def effective_properties(umat_name: str, props: np.ndarray, nstatev: int = 1,
                         psi: float = 0., theta: float = 0., phi: float = 0.,
                         symmetry: str = 'auto') -> Dict[str, float]:
    """
    Compute effective elastic properties of a composite.

    Parameters
    ----------
    umat_name : str
        Homogenization scheme (MIHEN, MIMTN, MISCN, etc.)
    props : ndarray
        Material properties array
    nstatev : int
        Number of state variables
    psi, theta, phi : float
        Euler angles for RVE orientation
    symmetry : str
        Expected symmetry: 'iso', 'isotrans', 'ortho', 'cubic', or 'auto'
        If 'auto', attempts to detect symmetry.

    Returns
    -------
    dict
        Engineering constants appropriate for the symmetry
    """
    L_eff = effective_stiffness(umat_name, props, nstatev, psi, theta, phi)

    if symmetry == 'auto':
        # Try to detect symmetry by checking tensor structure
        # Simple heuristic based on diagonal elements
        diag = np.diag(L_eff)
        off_diag = L_eff[0, 1]

        # Check if isotropic (all diagonal elements equal)
        if np.allclose(diag[:3], diag[0], rtol=0.01) and np.allclose(diag[3:], diag[3], rtol=0.01):
            symmetry = 'iso'
        else:
            symmetry = 'ortho'

    if symmetry == 'iso':
        return recover_isotropic(L=L_eff)
    elif symmetry == 'isotrans':
        return recover_transversely_isotropic(L=L_eff, axis=1)
    elif symmetry == 'ortho':
        return recover_orthotropic(L=L_eff)
    elif symmetry == 'cubic':
        return recover_cubic(L=L_eff)
    else:
        # Return orthotropic as default
        return recover_orthotropic(L=L_eff)


# =============================================================================
# Directional Properties
# =============================================================================

def directional_modulus(L: np.ndarray, direction: np.ndarray) -> float:
    """
    Compute Young's modulus in a given direction.

    Parameters
    ----------
    L : ndarray
        6x6 stiffness tensor
    direction : ndarray
        3D direction vector (will be normalized)

    Returns
    -------
    float
        Young's modulus in the specified direction
    """
    d = np.asarray(direction, dtype=float)
    d = d / np.linalg.norm(d)

    M = np.linalg.inv(L)

    # Build the Voigt strain vector for uniaxial stress in direction d
    # For uniaxial stress sigma_d in direction d, strain is e_ij = M_ijkl * sigma * d_k * d_l
    # The compliance in direction d is: S_d = d_i d_j M_ijkl d_k d_l

    # Convert direction to Voigt form for double contraction
    d_voigt = np.array([
        d[0]**2, d[1]**2, d[2]**2,
        d[0]*d[1], d[0]*d[2], d[1]*d[2]
    ])

    # For proper Voigt notation with engineering shear
    d_voigt_full = np.array([
        d[0]**2, d[1]**2, d[2]**2,
        2*d[0]*d[1], 2*d[0]*d[2], 2*d[1]*d[2]
    ])

    S_d = d_voigt @ M @ d_voigt_full

    return 1.0 / S_d


def directional_modulus_surface(L: np.ndarray, n_theta: int = 36, n_phi: int = 18) -> Dict:
    """
    Compute Young's modulus for all directions (for 3D visualization).

    Parameters
    ----------
    L : ndarray
        6x6 stiffness tensor
    n_theta : int
        Number of azimuthal angles
    n_phi : int
        Number of polar angles

    Returns
    -------
    dict
        'theta': azimuthal angles
        'phi': polar angles
        'E': Young's modulus values (n_theta x n_phi array)
        'x', 'y', 'z': Cartesian coordinates for plotting
    """
    theta = np.linspace(0, 2*np.pi, n_theta)
    phi = np.linspace(0, np.pi, n_phi)

    E = np.zeros((n_theta, n_phi))

    for i, t in enumerate(theta):
        for j, p in enumerate(phi):
            direction = np.array([
                np.sin(p) * np.cos(t),
                np.sin(p) * np.sin(t),
                np.cos(p)
            ])
            E[i, j] = directional_modulus(L, direction)

    # Create surface coordinates (modulus as radius)
    THETA, PHI = np.meshgrid(theta, phi, indexing='ij')
    X = E * np.sin(PHI) * np.cos(THETA)
    Y = E * np.sin(PHI) * np.sin(THETA)
    Z = E * np.cos(PHI)

    return {
        'theta': theta,
        'phi': phi,
        'E': E,
        'x': X, 'y': Y, 'z': Z
    }


# =============================================================================
# Bulk Properties
# =============================================================================

def bulk_modulus(L: np.ndarray) -> float:
    """
    Compute bulk modulus from stiffness tensor (Voigt average).

    K = (L11 + L22 + L33 + 2*(L12 + L13 + L23)) / 9
    """
    return (L[0,0] + L[1,1] + L[2,2] + 2*(L[0,1] + L[0,2] + L[1,2])) / 9.0


def shear_modulus_voigt(L: np.ndarray) -> float:
    """
    Compute Voigt average shear modulus from stiffness tensor.

    G_V = (L11 + L22 + L33 - L12 - L13 - L23 + 3*(L44 + L55 + L66)) / 15
    """
    return (L[0,0] + L[1,1] + L[2,2] - L[0,1] - L[0,2] - L[1,2] +
            3*(L[3,3] + L[4,4] + L[5,5])) / 15.0


def shear_modulus_reuss(M: np.ndarray) -> float:
    """
    Compute Reuss average shear modulus from compliance tensor.

    1/G_R = (4*(M11 + M22 + M33) - 4*(M12 + M13 + M23) + 3*(M44 + M55 + M66)) / 15
    """
    inv_G = (4*(M[0,0] + M[1,1] + M[2,2]) - 4*(M[0,1] + M[0,2] + M[1,2]) +
             3*(M[3,3] + M[4,4] + M[5,5])) / 15.0
    return 1.0 / inv_G


def shear_modulus_hill(L: np.ndarray) -> float:
    """
    Compute Hill average shear modulus (Voigt-Reuss-Hill average).
    """
    M = np.linalg.inv(L)
    G_V = shear_modulus_voigt(L)
    G_R = shear_modulus_reuss(M)
    return (G_V + G_R) / 2.0


def universal_anisotropy_index(L: np.ndarray) -> float:
    """
    Compute universal elastic anisotropy index A^U.

    A^U = 5 * G_V/G_R + K_V/K_R - 6

    A^U = 0 for isotropic materials.

    Reference: Ranganathan & Ostoja-Starzewski (2008)
    """
    M = np.linalg.inv(L)

    K_V = bulk_modulus(L)
    G_V = shear_modulus_voigt(L)

    # Reuss averages
    K_R = 1.0 / (M[0,0] + M[1,1] + M[2,2] + 2*(M[0,1] + M[0,2] + M[1,2]))
    G_R = shear_modulus_reuss(M)

    return 5 * G_V / G_R + K_V / K_R - 6


# =============================================================================
# Convenience Exports
# =============================================================================

__all__ = [
    # Material classes
    'IsotropicMaterial',
    'TransverselyIsotropicMaterial',
    'OrthotropicMaterial',
    'CubicMaterial',
    # Recovery functions
    'recover_isotropic',
    'recover_transversely_isotropic',
    'recover_orthotropic',
    'recover_cubic',
    'recover_anisotropic',
    # Effective properties
    'effective_stiffness',
    'effective_properties',
    # Directional properties
    'directional_modulus',
    'directional_modulus_surface',
    # Bulk properties
    'bulk_modulus',
    'shear_modulus_voigt',
    'shear_modulus_reuss',
    'shear_modulus_hill',
    'universal_anisotropy_index',
]
