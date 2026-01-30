=====
Phase
=====

This module provides classes and functions for managing material phases
and their properties.

.. note::
   **Python Interface Recommended**

   For new projects, use the Python ``simcoon.solver.micromechanics`` module instead
   of the C++ file I/O functions. The Python interface provides:

   - JSON-based configuration (self-documenting, easy to modify)
   - Dataclasses for phases, layers, ellipsoids, cylinders, sections
   - Legacy ``.dat`` file conversion utilities

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
