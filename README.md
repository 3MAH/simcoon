Simcoon
=========


[![Simcoon Logo](https://github.com/simcoon/simcoon/blob/master/docs/_static/simcoon_logo_small.png?raw=true)](https://github.com/simcoon/simcoon)

About
-----

Simcoon is a free, open-source library for the simulation of multiphysics systems. Its primarily objective was the simulation of heterogeneous materials, but now goes beyhond with the simulation of other physics, e.g. biological systems. It is developed with the aim to be a high-quality scientific library to facilitate the analysis of the complex, non-linear behavior of systems.
    It integrates tools to simulate the response of composite material response and thus integrates several algorithms for the analysis of heterogeneous materials.

Simcoon is a C++ library with emphasis on speed and ease-of-use. Its principle focus is to provide tools to facilitate the implementation of up-to-date constitutive model for materials in Finite Element Analysis Packages. This is done by providing a C++ API to generate user material subroutine based on a library of functions. Also, SMART+ provides tools to analyse the behavior of material, considering loading at the material point level. Such tools include a thermomechanical solver, a software to predict effective properties of composites, and a built-in identification software (using a combined genetic-gradient based algorithm)

Simcoon is mainly developed by faculty and researchers from University of Bordeaux and the I2M Laboratory (Institut de d'Ingénierie et de Mécanique). Fruitful contribution came from the LEM3 laboratory in Metz, France, TU Bergakademie Freiberg in Germany and the TIMC-IMAG laboratory in Grenoble, France. It is released under the GNU General Public License: GPL, version 3.

[![GitHub license](https://img.shields.io/badge/licence-GPL%203-blue.svg)](https://github.com/chemiskyy/simcoon/blob/master/LICENSE.txt)

Simcoon make use and therefore include the FTensor library (http://www.wlandry.net/Projects/FTensor) for convenience. FTensor is a library that handle complex tensor computations. FTensor is released under the GNU General Public License: GPL, version 2. You can get it there (but is is already included in simcoon): (https://bitbucket.org/wlandry/ftensor)

Documentation
--------------

Provider      | Status
--------      | ------
Read the Docs | [![Documentation Status](https://readthedocs.org/projects/simcoon/badge/?version=latest&style=plastic)](http://simcoon.readthedocs.io/en/latest)


Installation
------------

How to install Simcoon (In most UNIX-based systems):

1 - Make sure you have Boost (1.66 at least) installed and Armadillo (latest version is best) installed to use smartplus. Make sure that you have Boost.Python compiled against Python3 (Since Python 2 will be outdated). Using Homebrew and Pyenv for MacOS is a good idea to correctly handle the Python version to use!

2 - Unzip the file in a source location.

3 - Go to such folder "simcoon (i.e. using cd)"

4 - Execute the installation bash file : 

```bash
sh Install.sh
```

5 - Enjoy

Alternative
--------------------

You can also Cmake to compile the simcoon library, and then navigate to "simcoon-python-builder" to build the simmit library, which is the Python binding of Simcoon. You can run the tests using 'make test'.

If your OS is Ubuntu 18.04, you will soon be able to use a script "install-simcoon.sh" that contains all the commands to install the simcoon dependancies and simcoon. Stay tuned! You will then just have to navigate to the simcoon folder downloaded and run:

```bash
sudo sh install-simcoon.sh
```
on a terminal, with the sudo admin privileges

How to use simcoon
--------------------

Several possibilities 

1 - Use the python wrapper for simcoon

By doing, so, you can utilize simcoon in Ipython (jupyter) notebooks. Several examples are provided

2 - Use the executables provided by simcoon

For instance, a solver, an identification software, etc..

3 - Link simcoon with FEA Packages

An example is provided on how to use simcoon constitutive models as Umat libraries:

You can directly copy-paste "umat_singleM.cpp" (mechanical) or "umat_singleT.cpp" (thermomechanical) from 'pathtothefile'/simcoon/software to your Abaqus work directory and use it like a classical Umat.

Example : 
```bash
abaqus make library=umat_singleM.cpp
```

4- Build your own projects using the Simcoon library

Link with -lsimcoon
Have fun :)

Authors
-------
* [Yves Chemisky](https://github.com/chemiskyy)