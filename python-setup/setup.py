from setuptools import setup
import sys

# Determine the appropriate file extension based on the OS
if sys.platform == "win32":
    extension_suffix = ".pyd"
else:
    extension_suffix = ".so"

from simcoon.__version__ import __version__

setup(
    name="simcoon",
    version=__version__,
    description="Simulation in Mechanics and Materials: Interactive Tools",
    author="Yves Chemisky",
    author_email="yves.chemisky@gmail.com",
    packages=[
        "simcoon",
    ],
    package_data={
        "simcoon": [
            f"simmit{extension_suffix}"
        ],  # Include .pyd or .so depending on the platform
    },
    include_package_data=True,
    license="GPL",
)
