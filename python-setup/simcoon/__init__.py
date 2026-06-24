import os
import sys

if sys.platform == "win32":
    _pkg_dir = os.path.dirname(os.path.abspath(__file__))
    if hasattr(os, "add_dll_directory"):
        os.add_dll_directory(_pkg_dir)
    os.environ["PATH"] = _pkg_dir + os.pathsep + os.environ.get("PATH", "")

from simcoon._core import *
from simcoon.__version__ import __version__
from simcoon.rotation import Rotation  # override _CppRotation from star-import
from simcoon.tensor import Tensor2, Tensor4, dyadic, auto_dyadic, sym_dyadic, auto_sym_dyadic, double_contract  # unified tensor classes
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
