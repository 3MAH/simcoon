Dataset Loaders
===============

Functions for loading material datasets, either from the bundled JSON file
or from user-supplied CSV files.

.. autofunction:: simcoon.ashby.data.load_builtin

.. autofunction:: simcoon.ashby.data.load_csv

Built-in dataset
----------------

The bundled ``materials.json`` contains ~60 curated engineering materials
drawn from standard references.  The dataset covers five material families:

.. list-table::
   :header-rows: 1
   :widths: 25 10 65

   * - Category
     - Count
     - Examples
   * - Metal
     - 24
     - Mild steel, AISI 304/316L, Al 6061-T6/7075-T6, Ti-6Al-4V, Inconel 718, Cu, Brass, W, …
   * - Ceramic
     - 10
     - Alumina, SiC, Si₃N₄, ZrO₂, B₄C, WC, soda-lime glass, fused silica, …
   * - Polymer
     - 15
     - HDPE, LDPE, PP, PVC, PMMA, Nylon 6.6, PET, PTFE, PC, ABS, epoxy, …
   * - Composite
     - 5
     - CFRP (UD), CFRP (woven), GFRP, Aramid/epoxy, plywood
   * - Foam
     - 3
     - PU foam, aluminium foam, cork

Every entry includes at least ``E``, ``nu``, and ``density``.  Many metals
also provide ``yield_strength``, ``tensile_strength``, and hardening
parameters.

Loading a CSV file
------------------

Any CSV file with a header row can be imported.  Columns whose names match
``Material`` fields are used automatically; others are ignored.

.. code-block:: python

   from simcoon.ashby import load_csv

   # Direct match — CSV headers = Material field names
   mats = load_csv("my_materials.csv")

   # With column renaming
   mats = load_csv("lab_data.csv", column_mapping={
       "Material": "name",
       "Youngs_GPa": "E",
       "Poisson": "nu",
       "Rho_kg_m3": "density",
   })
