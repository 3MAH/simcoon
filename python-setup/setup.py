from setuptools import setup

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
    package_data={"simcoon": ["simmit.*", "*.dll"]},
    include_package_data=True,
    license="GPL",
)
