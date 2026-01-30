Heterogeneous materials simulation
-----------------------------------------------

Below are examples illustrating Simcoon's capabilities for simulating heterogeneous materials.

This gallery contains examples demonstrating:

- **Eshelby tensors** - Computing Eshelby and Hill interaction tensors for various inclusion shapes
- **Micromechanics** - Effective properties using Mori-Tanaka and self-consistent schemes, effect of volume fraction and aspect ratio

Phase Configuration
^^^^^^^^^^^^^^^^^^^

Phase configurations (ellipsoids, layers, cylinders) can be defined using the
Python ``simcoon.solver.micromechanics`` module with JSON format:

.. code-block:: python

   from simcoon.solver.micromechanics import (
       Ellipsoid, load_ellipsoids_json, save_ellipsoids_json,
   )

   # Load phases from JSON
   phases = load_ellipsoids_json('data/ellipsoids.json')

   # Modify programmatically
   phases[1].a1 = 20.0  # Change aspect ratio
   phases[1].concentration = 0.35

   # Save back
   save_ellipsoids_json('data/ellipsoids.json', phases)

See ``data/ellipsoids.json`` for the JSON format example.

.. note::
   The legacy ``.dat`` file formats (``Nellipsoids0.dat``, etc.) are deprecated.
   Use ``convert_legacy_ellipsoids()`` to migrate existing files to JSON.