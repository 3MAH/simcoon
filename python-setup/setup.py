from setuptools import setup, find_packages

setup(
    name="simcoon",
    version="1.9.3",
    description="Simulation in Mechanics and Materials: Interactive Tools",
    author="Yves Chemisky",
    author_email="yves.chemisky@gmail.com",
    # url=
    packages=[
        "simcoon",
    ],
    package_data={"simcoon": ["simmit.so"]},
    include_package_data=True,
    license="GPL",
)
