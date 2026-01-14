====================================================
Simcoon — Constitutive Modeling and Material Library
====================================================
.. image:: _static/simcoon_overview.pdf
   :alt: Simcoon overview
   :class: align-center
   :width: 700px

📌 Overview
===========

Simcoon is a free, open-source C++ library for simulating multiphysics systems. Its primary objective was the development of constitutive models for heterogeneous materials, and it now also provides tools for full-field simulation. Together with
[microgen](https://github.com/3MAH/microgen) for CAD and meshing of heterogeneous materials and
[fedoo](https://github.com/3MAH/fedoo), our finite-element solver, the project provides a comprehensive toolset for analyzing heterogeneous materials.

Simcoon aims to be a high-quality scientific library for analysing complex, nonlinear system behaviour. It emphasizes performance and ease of use and exposes a Python interface to simplify workflows. The principal focus is to provide a C++ API to generate user-material subroutines for finite-element packages. Simcoon also includes tools to analyse material-point behaviour under loading, such as a thermomechanical solver, a tool to predict effective properties of composites, and a built-in identification tool that combines genetic and gradient-based algorithms.

Simcoon supports geometric nonlinearity using Lagrangian and Eulerian measures, and cumulative strains with several objective rates (Jaumann, Green--Naghdi, and Xi--Meyers--Brühns logarithmic). The logarithmic cumulative-strain measure is the default used by the library's constitutive laws.

Simcoon makes use of the FTensor library (http://www.wlandry.net/Projects/FTensor) for tensor computations; FTensor is distributed under the GNU General Public License (GPL v2) and is included in this repository.

Simcoon is mainly developed by faculty and researchers from Université Grenoble Alpes (TIMC laboratory) and the University of Bordeaux (I2M Laboratory). Contributions also came from the LEM3 laboratory (Metz, France) and TU Bergakademie Freiberg (Germany). The project is released under the GNU General Public License (GPL v3).

[![GitHub license](https://img.shields.io/badge/license-GPL%203-blue.svg)](https://github.com/chemiskyy/simcoon/blob/master/LICENSE)

⚙️ Installation
===============

Install via conda:

.. code-block:: bash

   conda install -c conda-forge simcoon

Complete installation instructions are available on the
:doc:`installation` page.

The guide covers:
- Conda-based installation
- Building from source
- Optional dependencies

Key features
------------

- 📐 Finite strain mechanics
- 🧱 Modular constitutive models
- 🐍 Python bindings
- ⚡ High-performance C++ core

🚀 Getting started
------------------

New users should begin with:

- :doc:`simulation/solver`
- :doc:`examples/index`

📚 Documentation
----------------
.. toctree::
   :maxdepth: 2
   :caption: Contents:

   about
   installation
   external
   simulation/index
   examples/index
   continuum_mechanics/index
   homogenization/index   

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
