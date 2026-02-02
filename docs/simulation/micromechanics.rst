====================================
Micromechanics I/O (Python Module)
====================================

Overview
--------

The ``simcoon.solver.micromechanics`` module provides:

- **Data classes** for phases, layers, ellipsoids, cylinders, and sections
- **JSON I/O functions** for loading and saving configurations
- **Standalone operation** without requiring ``simcoon._core`` to be built

Quick Start
-----------

.. code-block:: python

   from simcoon.solver.micromechanics import (
       Ellipsoid, MaterialOrientation,
       load_ellipsoids_json, save_ellipsoids_json,
   )
   import numpy as np

   # Create phases programmatically
   matrix = Ellipsoid(
       number=0,
       concentration=0.7,
       umat_name="ELISO",
       props=np.array([3500, 0.35, 60e-6]),  # E, nu, alpha
       a1=1, a2=1, a3=1,  # Spherical (matrix)
   )

   fibers = Ellipsoid(
       number=1,
       concentration=0.3,
       umat_name="ELISO",
       props=np.array([72000, 0.22, 5e-6]),
       a1=20, a2=1, a3=1,  # Prolate spheroid (aspect ratio 20)
   )

   # Save to JSON
   save_ellipsoids_json('composite.json', [matrix, fibers])

   # Load from JSON
   phases = load_ellipsoids_json('composite.json')
   for p in phases:
       print(f"Phase {p.number}: {p.shape_type}, c={p.concentration}")

Data Classes
------------

MaterialOrientation
^^^^^^^^^^^^^^^^^^^

Material orientation via Euler angles (z-x-z convention, degrees).

.. code-block:: python

   MaterialOrientation(psi=0.0, theta=0.0, phi=0.0)

GeometryOrientation
^^^^^^^^^^^^^^^^^^^

Geometry/phase orientation via Euler angles (degrees).

.. code-block:: python

   GeometryOrientation(psi=0.0, theta=0.0, phi=0.0)

Phase
^^^^^

Base class for all phase types.

.. code-block:: python

   Phase(
       number=0,                    # Phase ID
       umat_name="ELISO",           # Constitutive model
       save=1,                      # Save flag (1=save)
       concentration=1.0,           # Volume fraction
       material_orientation=MaterialOrientation(),
       nstatev=1,                   # Number of state variables
       props=np.array([]),          # Material properties
   )

Layer
^^^^^

Layer phase for laminate homogenization (extends Phase).

.. code-block:: python

   Layer(
       # ... Phase attributes ...
       geometry_orientation=GeometryOrientation(),
       layerup=-1,     # Index of layer above (-1 if none)
       layerdown=-1,   # Index of layer below (-1 if none)
   )

Ellipsoid
^^^^^^^^^

Ellipsoidal inclusion for Eshelby-based homogenization (extends Phase).

.. code-block:: python

   Ellipsoid(
       # ... Phase attributes ...
       coatingof=0,    # Index of phase this coats (0 if none)
       a1=1.0,         # First semi-axis
       a2=1.0,         # Second semi-axis
       a3=1.0,         # Third semi-axis
       geometry_orientation=GeometryOrientation(),
   )

   # Shape type property (computed from semi-axes):
   ell.shape_type  # "sphere", "prolate_spheroid", "oblate_spheroid", "general_ellipsoid"

Cylinder
^^^^^^^^

Cylindrical inclusion for micromechanics (extends Phase).

.. code-block:: python

   Cylinder(
       # ... Phase attributes ...
       coatingof=0,    # Index of phase this coats
       L=1.0,          # Length
       R=1.0,          # Radius
       geometry_orientation=GeometryOrientation(),
   )

   # Aspect ratio property:
   cyl.aspect_ratio  # L / R

Section
^^^^^^^

Section/yarn for textile composite homogenization.

.. code-block:: python

   Section(
       number=0,
       name="Section",
       umat_name="ELISO",
       material_orientation=MaterialOrientation(),
       nstatev=1,
       props=np.array([]),
   )

JSON I/O Functions
------------------

Loading Functions
^^^^^^^^^^^^^^^^^

.. code-block:: python

   # Load phases from JSON
   phases = load_phases_json('phases.json')
   layers = load_layers_json('layers.json')
   ellipsoids = load_ellipsoids_json('ellipsoids.json')
   cylinders = load_cylinders_json('cylinders.json')
   sections = load_sections_json('sections.json')

Saving Functions
^^^^^^^^^^^^^^^^

.. code-block:: python

   # Save with optional property names
   save_phases_json('phases.json', phases, prop_names=['E', 'nu', 'alpha'])
   save_layers_json('layers.json', layers)
   save_ellipsoids_json('ellipsoids.json', ellipsoids)
   save_cylinders_json('cylinders.json', cylinders)
   save_sections_json('sections.json', sections)

JSON Format Examples
--------------------

Ellipsoids JSON
^^^^^^^^^^^^^^^

.. code-block:: json

   {
     "ellipsoids": [
       {
         "number": 0,
         "coatingof": 0,
         "umat_name": "ELISO",
         "save": 1,
         "concentration": 0.7,
         "material_orientation": {"psi": 0, "theta": 0, "phi": 0},
         "semi_axes": {"a1": 1, "a2": 1, "a3": 1},
         "geometry_orientation": {"psi": 0, "theta": 0, "phi": 0},
         "nstatev": 1,
         "props": {"E": 3500, "nu": 0.35, "alpha": 6e-5}
       },
       {
         "number": 1,
         "coatingof": 0,
         "umat_name": "ELISO",
         "save": 1,
         "concentration": 0.3,
         "material_orientation": {"psi": 0, "theta": 0, "phi": 0},
         "semi_axes": {"a1": 20, "a2": 1, "a3": 1},
         "geometry_orientation": {"psi": 0, "theta": 0, "phi": 0},
         "nstatev": 1,
         "props": {"E": 72000, "nu": 0.22, "alpha": 5e-6}
       }
     ]
   }

Layers JSON
^^^^^^^^^^^

.. code-block:: json

   {
     "layers": [
       {
         "number": 0,
         "umat_name": "ELISO",
         "save": 1,
         "concentration": 0.5,
         "material_orientation": {"psi": 0, "theta": 0, "phi": 0},
         "geometry_orientation": {"psi": 0, "theta": 90, "phi": -90},
         "nstatev": 1,
         "props": {"E": 70000, "nu": 0.3, "alpha": 1e-5},
         "layerup": -1,
         "layerdown": 1
       }
     ]
   }

Cylinders JSON
^^^^^^^^^^^^^^

.. code-block:: json

   {
     "cylinders": [
       {
         "number": 0,
         "coatingof": 0,
         "umat_name": "ELISO",
         "save": 1,
         "concentration": 0.3,
         "material_orientation": {"psi": 0, "theta": 0, "phi": 0},
         "L": 50.0,
         "R": 1.0,
         "geometry_orientation": {"psi": 0, "theta": 0, "phi": 0},
         "nstatev": 1,
         "props": {"E": 72000, "nu": 0.22, "alpha": 5e-6}
       }
     ]
   }

Standalone Usage
----------------

The micromechanics module can be imported without building ``simcoon._core``:

.. code-block:: python

   # This works without the C++ extension module
   from simcoon.solver.micromechanics import (
       Ellipsoid, Layer,
       load_ellipsoids_json, save_ellipsoids_json,
   )

This is useful for:

- Pre-processing phase configurations before running simulations
- Post-processing results
- CI/CD pipelines that don't need the full simcoon build
- Teaching and documentation

See Also
--------

- :doc:`solver` - Full solver documentation
- :doc:`../homogenization/index` - Homogenization theory and examples
- :doc:`../examples/index` - Complete examples
