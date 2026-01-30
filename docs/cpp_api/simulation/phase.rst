=====
Phase
=====

This module provides classes and functions for managing material phases
and their properties.

Python Interface
----------------

For convenient phase configuration, use the Python ``simcoon.solver.micromechanics`` module:

- JSON-based configuration (self-documenting, easy to modify)
- Dataclasses for phases, layers, ellipsoids, cylinders, sections

See :doc:`../../simulation/micromechanics` for the Python API documentation.

.. code-block:: python

   from simcoon.solver.micromechanics import (
       Ellipsoid, load_ellipsoids_json, save_ellipsoids_json,
   )

Phase Characteristics
---------------------

.. doxygenfile:: simcoon/Simulation/Phase/phase_characteristics.hpp
   :project: simcoon
   :sections: briefdescription detaileddescription innerclass func typedef enum
