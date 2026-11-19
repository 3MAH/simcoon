Material Data Model
===================

The data model lives in ``simcoon.ashby.material`` and is imported directly
from the top-level ``simcoon.ashby`` namespace.

SymmetryType
------------

.. autoclass:: simcoon.ashby.material.SymmetryType
   :members:
   :undoc-members:

The enum values map directly to simcoon's constitutive-model names:

.. list-table::
   :header-rows: 1
   :widths: 35 20 45

   * - Member
     - Value
     - Description
   * - ``ISOTROPIC``
     - ``"ELISO"``
     - Isotropic elastic (2 constants: E, nu)
   * - ``CUBIC``
     - ``"ELCUB"``
     - Cubic elastic (3 constants: E, nu, G)
   * - ``TRANSVERSE_ISOTROPIC``
     - ``"ELIST"``
     - Transversely isotropic (5 constants)
   * - ``ORTHOTROPIC``
     - ``"ELORT"``
     - Orthotropic (9 constants)
   * - ``ANISOTROPIC``
     - ``"ELANI"``
     - Fully anisotropic (21 constants)
   * - ``UNKNOWN``
     - ``"UNKNOWN"``
     - Symmetry not determined

Material
--------

.. autoclass:: simcoon.ashby.material.Material
   :members:
   :undoc-members:

**Units convention**

.. list-table::
   :header-rows: 1
   :widths: 30 25 45

   * - Property
     - Unit
     - Notes
   * - ``E``, ``G``, ``K``
     - GPa
     - Converted to MPa (×1000) when passed to simcoon
   * - ``density``
     - kg/m³
     -
   * - ``CTE``
     - 1/K
     -
   * - ``yield_strength``, ``tensile_strength``
     - MPa
     -

**Hardening models**

The ``hardening_type`` field selects a hardening law and ``hardening_params``
provides the model-specific parameters:

- ``"linear"`` — :math:`\sigma = \sigma_y + H\,\varepsilon_p` with
  ``{"H": float}``
- ``"power_law"`` (Hollomon) — :math:`\sigma = K\,\varepsilon^n` with
  ``{"K": float, "n": float}``
- ``"johnson_cook"`` —
  :math:`\sigma = (A + B\,\varepsilon^n)(1 + C\,\ln\dot\varepsilon^*)(1 - T^{*m})`
  with ``{"A": float, "B": float, "n": float, "C": float, "m": float}``

**Creating materials**

.. code-block:: python

   from simcoon.ashby import Material, SymmetryType

   # From keyword arguments
   steel = Material(
       name="Mild steel",
       E=205, nu=0.29, G=79.5,
       density=7870,
       symmetry=SymmetryType.ISOTROPIC,
       category="Metal",
   )

   # From a dictionary (e.g. loaded from JSON)
   mat = Material.from_dict({
       "name": "Al 6061-T6",
       "E": 68.9, "nu": 0.33,
       "density": 2700,
       "symmetry": "ISOTROPIC",
   })

MaterialCollection
------------------

.. autoclass:: simcoon.ashby.material.MaterialCollection
   :members:
   :undoc-members:

**Filtering and grouping**

.. code-block:: python

   from simcoon.ashby import load_builtin

   mats = load_builtin()

   # Filter by any field
   metals = mats.filter(category="Metal")
   stiff  = mats.filter(symmetry=SymmetryType.ISOTROPIC)

   # Group into sub-collections
   by_cat = mats.group_by("category")
   # {'Metal': MaterialCollection(24), 'Ceramic': MaterialCollection(10), ...}

   # Extract arrays for plotting
   x, y, names = metals.get_property_arrays("density", "E")
