import os
import sys

if sys.platform == "win32":
    _pkg_dir = os.path.dirname(os.path.abspath(__file__))
    if hasattr(os, "add_dll_directory"):
        os.add_dll_directory(_pkg_dir)
    os.environ["PATH"] = _pkg_dir + os.pathsep + os.environ.get("PATH", "")

from simcoon._core import *
from simcoon import modular
# In-memory solver API (2.0): the subpackage replaces the pre-2.0 file-driven
# `solver` function from the star-import above (still reachable as
# simcoon._core.solver; parse legacy files with simcoon.solver.from_file).
# `import` (not `from ... import`) is required: the latter would return the
# existing function attribute without importing the subpackage.
import simcoon.solver
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


# Deprecated free projector functions (removed in simcoon 2.0). These shadow the same
# names brought in by the `from simcoon._core import *` above; see simcoon/_deprecated.py.
from simcoon._deprecated import Ireal, Ireal2, Ivol, Idev, Idev2
