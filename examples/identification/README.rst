Parameter identification examples
-----------------------------------------------

Below are examples illustrating Simcoon's parameter-identification workflow:
calibrating constitutive-model parameters from experimental data with the
identification API (``sim.identification`` and ``sim.calc_cost``) on top of the
material-point solver and homogenization tools.

This gallery contains examples demonstrating:

- **Chaboche cyclic plasticity** - Identifying 7 elasto-plastic parameters (``EPCHA`` UMAT) from cyclic uniaxial tests, with an NMSE-per-response cost

.. note::

   Two related hyperelastic identification examples live in other galleries:
   ``analysis/hyperelastic_parameter_identification.py`` (hand-rolled analytical
   stress + manual ``scipy``/``sklearn`` cost) and
   ``hyperelasticity/hyperelastic_umat_identification.py`` (the deployable
   ``MOORI`` UMAT driven through ``sim.solver`` with the identification API).
