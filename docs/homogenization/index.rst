Homogenization
==================================

This section covers homogenization methods for computing effective properties
of heterogeneous materials.

Phase Configuration
-------------------

Define inclusion phases (ellipsoids, layers, cylinders) using the Python
``simcoon.solver.micromechanics`` module:

.. code-block:: python

   from simcoon.solver.micromechanics import (
       Ellipsoid, load_ellipsoids_json, save_ellipsoids_json,
   )

   # Create fiber-reinforced composite
   matrix = Ellipsoid(number=0, concentration=0.7, a1=1, a2=1, a3=1,
                      props=[3500, 0.35])
   fiber = Ellipsoid(number=1, concentration=0.3, a1=50, a2=1, a3=1,
                     props=[72000, 0.22])

   save_ellipsoids_json('composite.json', [matrix, fiber])

See :doc:`../simulation/micromechanics` for complete documentation.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   doc_eshelby.rst