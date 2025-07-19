from simcoon.__version__ import __version__

# Try to import simmit module - it may not be available during build
try:
    from simcoon.simmit import *
except ImportError:
    # simmit module not built yet or not available
    pass