The Ashby Module
================

The ``simcoon.ashby`` subpackage provides tools for working with material
property databases and creating Ashby-style material selection charts.  It
bridges the gap between material data (from curated datasets or online
databases like Materials Project) and simcoon's constitutive modelling
functions.

Highlights
----------

- **Built-in dataset** of ~60 common engineering materials (metals, ceramics,
  polymers, composites, foams) with elastic, thermal, and strength properties.
- **Ashby diagram plotting** with automatic grouping by material family,
  envelope drawing (ellipse, convex hull), and performance-index guide lines.
- **Materials Project integration** — fetch elastic properties for thousands
  of crystalline compounds directly from the Materials Project API.
- **Simcoon bridge** — convert any material to a 6x6 stiffness or compliance
  tensor ready for use with ``sim.L_iso()``, ``sim.solver()``, etc.

Installation
------------

The core data model (``Material``, ``load_builtin()``) requires only
**numpy**, which is already a simcoon dependency.

For plotting and database features, install the optional ``ashby`` extras:

.. code-block:: bash

   pip install 'simcoon[ashby]'

This pulls in ``matplotlib``, ``scipy``, and ``mp-api``.

Quick start
-----------

.. code-block:: python

   from simcoon.ashby import load_builtin, ashby_plot, to_stiffness

   # 1. Load the curated material dataset
   mats = load_builtin()                          # ~60 materials, no optional deps

   # 2. Create an Ashby diagram
   ax = ashby_plot(mats, "density", "E")           # E vs density, log-log

   # 3. Convert a material to a simcoon stiffness tensor
   steel = mats.filter(category="Metal")[0]
   L = to_stiffness(steel)                         # 6x6 ndarray in MPa

Package structure
-----------------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Module
     - Description
   * - :doc:`material`
     - ``Material`` dataclass, ``SymmetryType`` enum, ``MaterialCollection``
   * - :doc:`data`
     - Built-in dataset loader (``load_builtin``) and CSV importer (``load_csv``)
   * - :doc:`plotting`
     - ``ashby_plot()`` and legacy visualization helpers
   * - :doc:`bridge`
     - Convert materials to simcoon stiffness/compliance tensors
   * - :doc:`database`
     - Fetch materials from the Materials Project API

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   material
   data
   plotting
   bridge
   database
