from simcoon._core import *
from simcoon.__version__ import __version__
from simcoon.rotation import Rotation  # override _CppRotation from star-import
from simcoon.tensor import Tensor2, Tensor4, dyadic, auto_dyadic, double_contract  # unified tensor classes

# Backward compatibility alias - simmit was the legacy module name
from simcoon import _core as simmit
