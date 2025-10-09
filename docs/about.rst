About
=====

Simcoon is a free, open-source library for the simulation of constitutive response of materials.    
Its primarily objective was the developement of constitutive models for the simulation of heterogeneous materials, but now goes beyond with tools
to facilitate their full-field simulation. Together with :ref:`microgen <https://microgen.readthedocs.io>` for the CAD and meshing of heterogeneous materials and :ref:`fedoo <https://3mah.github.io/fedoo-docs>` our Finite Element solver,
we offer a comprehensive simulation set for the in-depth analysis of heterogeneous materials.

Simcoon is developed with the aim to be a high-quality scientific library to facilitate the analysis of the complex, non-linear behavior of systems.
It integrates tools to simulate the response of composite material response and thus integrates several algorithms for the analysis of heterogeneous materials.
Simcoon integrates an easy way to handle geometrical non-linearities : Use of Lagrangian measures, Eulerian measures and cumulative strains considering several spins : Jaumann, Green-Naghdi, Xi-Meyers-Bruhns logarithmic. With this last measure, cumulative strain correspond to a logarithmic strain measure and is the standard measure utilized for our constitutive laws.

Some of its main features are:

- Built on top of **Armadillo**, a high-quality C++ linear algebra library.
- Uses **FTensor** for complex tensor operations.
- Integrates several algorithms for heterogeneous materials analysis.
- Provides both C++ and Python API to generate user material subroutines for Finite Element Analysis packages.

simcoon is a C++ library with emphasis on speed and ease-of-use. Its principle focus is to provide tools to facilitate
the implementation of up-to-date constitutive model for materials in Finite Element Analysis Packages. This is done by
providing a C++ and a Python API to generate user material subroutine based on a library of functions. Also, simcoon
provides tools to analyse the behavior of material, considering loading at the material point level.

Simcoon is mainly developed by contributors from the staff and students of TIMC laboratory (Université Grenoble Alpes, France),
the I2M laboratory (Université de Bordeaux, France). It is released under the GNU General Public License: GPL, version 3.
Several institutions have contributed to the development of simcoon:

* **Université de Grenoble Alpes**: University of Grenoble Alpes, France.
* **Université de Bordeaux**: University of Bordeaux, France.
* **TIMC Laboratory**: Recherche Translationnelle et Innovation en Médecine et Complexité, Grenoble, France.
* **I2M Laboratory**: Institut de Mécanique et d'Ingénierie, Bordeaux, France.
* **LEM3 Laboratory**: Laboratoire d’Étude des Microstructures et de Mécaniquedes Matériaux, Metz, Germany.
* **TU Bergakademie Freiberg**: School of Mines and Technical University in Freiberg, Germany.
* **CNRS**: National French Center for scientific research.

.. raw:: html

    <div style="text-align: center;">
        <img src="_static/UB.png" style="height:80px; margin-right:20px;">
        <img src="_static/UGA.png" style="height:80px; margin-right:20px;">
        <img src="_static/CNRS_logo_med.png" style="height:80px;">
    </div>
