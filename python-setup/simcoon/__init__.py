from simcoon._core import *
from simcoon.__version__ import __version__
from simcoon.rotation import Rotation  # override _CppRotation from star-import

# Backward compatibility alias - simmit was the legacy module name
from simcoon import _core as simmit


def identification(*args, **kwargs):
    """Removed in simcoon 2.0. Use scipy.optimize.differential_evolution instead."""
    raise NotImplementedError(
        "simcoon.identification() was removed in v2.0. "
        "Use scipy.optimize.differential_evolution with the simcoon "
        "Parameter key system instead. "
        "See docs/simulation/identification.rst and "
        "examples/heterogeneous/composite_parameter_identification.py"
    )


def calc_cost(*args, **kwargs):
    """Removed in simcoon 2.0. Use numpy or sklearn.metrics instead."""
    raise NotImplementedError(
        "simcoon.calc_cost() was removed in v2.0. "
        "Use np.mean((y_pred - y_exp)**2) or "
        "sklearn.metrics.mean_squared_error instead."
    )
