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

# Backward compatibility alias - simmit was the legacy module name
from simcoon import _core as simmit
