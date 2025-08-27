import sys

from setuptools import setup

# Determine the appropriate file extension based on the OS
if sys.platform == "win32":
    extension_suffix = ".pyd"
else:
    extension_suffix = ".so"

# Read version directly from file to avoid importing the package during build
with open("simcoon/__version__.py", "r") as f:
    exec(f.read())
# Now __version__ is available

setup(
    name="simcoon",
    version=__version__,
    description="Simulation in Mechanics and Materials: Interactive Tools",
    author="Yves Chemisky",
    author_email="yves.chemisky@gmail.com",
    packages=[
        "simcoon",
    ],
    package_data={"simcoon": ["*.so", "*.dylib", "*.dll", "*.pyd"]},
    include_package_data=True,
    license="GPL",
)
