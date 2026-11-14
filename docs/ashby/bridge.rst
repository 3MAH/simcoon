Simcoon Bridge
==============

The bridge module converts ``Material`` objects into simcoon-compatible
stiffness and compliance tensors.  This lets you go directly from a material
database query to a simcoon simulation.

.. autofunction:: simcoon.ashby.bridge.to_stiffness
.. autofunction:: simcoon.ashby.bridge.to_compliance
.. autofunction:: simcoon.ashby.bridge.to_solver_props

Unit conversion
---------------

Material properties in the Ashby module are stored in **GPa** (elastic moduli)
and **MPa** (strengths), following common engineering-data conventions.
Simcoon's internal convention is **MPa** for all stresses and moduli.

The ``unit_factor`` parameter (default ``1e3``) handles the GPa → MPa
conversion automatically.  If your material data is already in MPa, pass
``unit_factor=1``.

Example
-------

.. code-block:: python

   from simcoon.ashby import load_builtin, to_stiffness, to_compliance, to_solver_props

   mats = load_builtin()
   al = [m for m in mats if "6061" in m.name][0]

   # 6x6 stiffness tensor (MPa)
   L = to_stiffness(al)

   # 6x6 compliance tensor (1/MPa)
   M = to_compliance(al)

   # Ready-made dict for sim.solver()
   sp = to_solver_props(al)
   # {'umat_name': 'ELISO', 'props': array([68900., 0.33]), 'nstatev': 1}

Supported symmetries
--------------------

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - SymmetryType
     - simcoon function
     - Required Material fields
   * - ``ISOTROPIC``
     - ``sim.L_iso`` / ``sim.M_iso``
     - ``E``, ``nu``
   * - ``CUBIC``
     - ``sim.L_cubic`` / ``sim.M_cubic``
     - ``E``, ``nu``, ``G``
   * - ``TRANSVERSE_ISOTROPIC``
     - (full tensor pass-through)
     - ``elastic_tensor`` (6x6)
   * - ``ORTHOTROPIC``
     - (full tensor pass-through)
     - ``elastic_tensor`` (6x6)

For transversely isotropic and orthotropic materials the bridge requires the
full ``elastic_tensor`` field to be populated (e.g. from Materials Project).
