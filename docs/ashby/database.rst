Database Integration
====================

The ``database`` module provides functions to fetch material data from
online databases.  Currently the `Materials Project <https://materialsproject.org>`_
is supported via the ``mp-api`` package.

.. autofunction:: simcoon.ashby.database.fetch_materials

Requirements
------------

Database queries require the ``mp-api`` package:

.. code-block:: bash

   pip install 'simcoon[ashby]'

You also need a Materials Project API key.  Sign up at
https://materialsproject.org and set the key either as an environment
variable:

.. code-block:: bash

   export MP_API_KEY="your_key_here"

or pass it directly:

.. code-block:: python

   from simcoon.ashby import fetch_materials
   mats = fetch_materials(elements=["Al"], api_key="your_key_here")

Example
-------

.. code-block:: python

   from simcoon.ashby import fetch_materials, ashby_plot

   # Fetch aluminium-containing compounds with elastic data
   al_mats = fetch_materials(elements=["Al", "O"], limit=50)
   print(f"Fetched {len(al_mats)} materials")

   # Inspect the first entry
   m = al_mats[0]
   print(m.name, m.source_id, m.E, m.nu)

   # Plot an Ashby diagram of the fetched data
   ashby_plot(al_mats, "density", "E", group_by="category")

Combining with the built-in dataset
------------------------------------

.. code-block:: python

   from simcoon.ashby import load_builtin, fetch_materials, MaterialCollection

   builtin = load_builtin()
   online  = fetch_materials(elements=["Ti"], limit=30)

   combined = MaterialCollection(list(builtin) + list(online))

Symmetry detection
------------------

When the Materials Project returns a full elastic tensor, the module
automatically calls ``simcoon.check_symetries()`` to determine the material
symmetry (isotropic, cubic, orthotropic, etc.).  This information is stored
in the ``symmetry`` field and used by the :doc:`bridge` to select the
correct simcoon tensor constructor.
