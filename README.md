Simcoon
=========


[![Simcoon Logo](https://github.com/3MAH/3mah.github.io/blob/master/assets/images/logo_simcoon/simcoon_logo_text.png?raw=true)](https://github.com/3MAH/simcoon)

/assets/images/logo_simcoon/

About
-----

Simcoon is a free, open-source library for the simulation of multiphysics systems. Its primarily objective was the developement of constitutive models for the simulation of heterogeneous materials, but now goes beyond with tools to facilitate their full-field simulation. Together with [microgen](https://github.com/3MAH/microgen) for the CAD and meshing of heterogeneous materials and [fedoo](https://github.com/3MAH/fedoo) our Finite Element solver, we offer a comprehensive simulation set for the in-depth analysis of heterogeneous materials.

Simcoon is developed with the aim to be a high-quality scientific library to facilitate the analysis of the complex, non-linear behavior of systems. It integrates tools to simulate the response of composite material response and thus integrates several algorithms for the analysis of heterogeneous materials.

Simcoon integrates 
- a easy way to handle geometrical non-linearities : Use of Lagrangian measures, Eulerian measures and cumulative strains considering several spins : Jaumann, Green-Naghdi, Xi-Meyers-Bruhns logarithmic. With this last measure, cumulative strain correspond to a logarithmic strain measure and is the standard measure utilized for our constitutive laws. 

Simcoon is a C++ library with emphasis on speed and ease-of-use, that offers a python interface to facilitate its use. Its principle focus is to provide tools to facilitate the implementation of up-to-date constitutive model for materials in Finite Element Analysis Packages. This is done by providing a C++ API to generate user material subroutine based on a library of functions. Also, Simconnn provides tools to analyse the behavior of material, considering loading at the material point level. Such tools include a thermomechanical solver, a software to predict effective properties of composites, and a built-in identification software (using a combined genetic-gradient based algorithm)

Simcoon is mainly developed by faculty and researchers from University of Bordeaux and the I2M Laboratory (Institut de d'Ingénierie et de Mécanique). Fruitful contribution came from the LEM3 laboratory in Metz, France, TU Bergakademie Freiberg in Germany and the TIMC-IMAG laboratory in Grenoble, France. It is released under the GNU General Public License: GPL, version 3.

[![GitHub license](https://img.shields.io/badge/licence-GPL%203-blue.svg)](https://github.com/chemiskyy/simcoon/blob/master/LICENSE.txt)

Simcoon make use and therefore include the FTensor library (http://www.wlandry.net/Projects/FTensor) for convenience. FTensor is a library that handle complex tensor computations. FTensor is released under the GNU General Public License: GPL, version 2. You can get it there (but is is already included in simcoon): (https://bitbucket.org/wlandry/ftensor)

Documentation
--------------

Provider      | Status
--------      | ------
Read the Docs | [![Documentation Status](https://readthedocs.org/projects/simcoon/badge/?version=latest&style=plastic)](http://simcoon.readthedocs.io/en/latest)


Building doc :
requires doxygen, sphinx, breathe

```bash
cd doxdocs && doxygen && make html
```

open _build/index.html


Installation
------------

It is now possible to install simcoon directly with conda :
```bash
conda install -c conda-forge -c set3mah simcoon
```
In case there are any conflicts, it is preferable to do it in a new conda environment :
```bash
conda create --name scientific
```

OR

The easiest way to install simcoon is to create a *conda* environnement: You can utilize the Anaconda GUI or type:
(for the installation of an environment called "scientific")

```bash
conda create --name scientific
```

To activate the environment: 

```bash
conda activate scientific
```

The next step is to install the required packages:

```bash
conda install -c conda-forge armadillo boost cgal numpy
```

Next, after downloading the simcoon sources in the github repository of [Simcoon](https://github.com/3MAH/simcoon). Unzip the content in a folder.

The last step is to run the installation script:

```bash
sh Install.sh
```



Authors
-------
* [Yves Chemisky](https://github.com/chemiskyy)
