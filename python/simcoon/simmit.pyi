"""Type stubs for simmit C++ extension module

This file provides type hints for IDE autocomplete and type checking.
The actual implementation is in C++ (compiled from python/bindings/).

The simmit module provides high-performance computational algorithms for:
- Continuum mechanics
- Constitutive models
- Material behavior
- Stress/strain analysis
"""

from typing import Literal
import numpy as np
from numpy.typing import NDArray

# Type aliases
FloatArray = NDArray[np.float64]
Convention = Literal["Enu", "KG", "EnuG"]

# Continuum Mechanics - Constitutive models
def L_iso(props: FloatArray, convention: str) -> FloatArray:
    """Compute isotropic stiffness tensor.

    Args:
        props: Material properties array
               - If convention="Enu": [E, nu] (Young's modulus, Poisson's ratio)
               - If convention="KG": [K, G] (Bulk modulus, Shear modulus)
        convention: Material property convention ("Enu", "KG", etc.)

    Returns:
        6x6 stiffness matrix in Voigt notation
    """
    ...

def M_iso(props: FloatArray, convention: str) -> FloatArray:
    """Compute isotropic compliance tensor.

    Args:
        props: Material properties array (same format as L_iso)
        convention: Material property convention

    Returns:
        6x6 compliance matrix in Voigt notation
    """
    ...

# Add more function stubs as needed based on the actual API
# This is a starter template - expand with all exported functions

def __version__() -> str:
    """Return the version of the core extension module."""
    ...
