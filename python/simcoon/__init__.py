"""Simcoon: Simulation in Mechanics and Materials

A C++ library with Python bindings for multiphysics simulations,
focusing on constitutive models for heterogeneous materials.
"""

# C++ extension module (compiled from bindings/)
# Pure Python utility modules
from simcoon import ashby, clean_data, constant, data, parameter, simmit
from simcoon.__version__ import __version__

__all__ = [
    "__version__",
    "simmit",
    "ashby",
    "constant",
    "clean_data",
    "data",
    "parameter",
]
