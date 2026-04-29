from simcoon._core import *
from simcoon.__version__ import __version__
from simcoon.rotation import Rotation  # override _CppRotation from star-import
from simcoon.identify import identification, calc_cost
from simcoon.parameter import (
    Parameter,
    read_parameters,
    copy_parameters,
    apply_parameters,
)
from simcoon.constant import (
    Constant,
    read_constants,
    copy_constants,
    apply_constants,
)

# Backward compatibility alias - simmit was the legacy module name
from simcoon import _core as simmit
