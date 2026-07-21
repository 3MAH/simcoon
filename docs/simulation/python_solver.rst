In-memory Python solver
================================

The ``simcoon.solver`` package drives the C++ material-point solver directly
from Python: the loading path is defined with :class:`~simcoon.solver.Block`
and :class:`~simcoon.solver.StepMeca` / :class:`~simcoon.solver.StepThermomeca`
objects, and the results come back as numpy arrays — no ``path.txt``,
``output.dat`` or result files involved. Since simcoon 2.0 this package *is*
``sim.solver``; legacy loading files are parsed into loading objects with
:func:`~simcoon.solver.from_file` (the raw pre-2.0 file-to-file runner remains
available as the low-level binding ``simcoon._core.solver``).

Quick start
-----------

A uniaxial tension test on an elastic isotropic material:

.. code-block:: python

    import numpy as np
    from simcoon import solver

    step = solver.StepMeca(
        control=['strain'] + ['stress'] * 5,   # E11 driven, lateral stress-free
        value=[0.01, 0, 0, 0, 0, 0],           # targets, Voigt order [11,22,33,12,13,23]
        time=1.0, ninc=100,
    )
    res = solver.solve(step, "ELISO", [70000., 0.3, 1.E-5], nstatev=1)

    stress = res["Stress"]     # Cauchy stress history, shape (6, N)
    strain = res["Strain"]     # Green-Lagrange strain history, shape (6, N)

Results follow the fedoo ``DataSet`` conventions — components first, one
column per increment — so they interoperate directly with fedoo utilities
(e.g. ``fedoo.util.voigt_tensors.StressTensorList(res["Stress"])``).
Available fields include ``Stress`` (Cauchy), ``Kirchhoff``, ``PKII``,
``Strain``, ``LogStrain``, ``F``, ``R``, ``DR`` (``(3, 3, N)``),
``TangentMatrix`` (``(6, 6, N)``), ``Statev``, ``Wm``, ``Time``, ``Temp`` and,
for thermomechanical runs, ``Q``, ``r``, ``Wt`` and the coupled tangents
``dSdE``, ``dSdT``, ``drdE``, ``drdT``.

Loading control
---------------

* ``control`` sets each component to ``'strain'`` (kinematically driven) or
  ``'stress'`` (statically driven); mixed control is solved by Newton-Raphson.
* ``Block(control_type=...)`` selects the strain/stress measures:
  ``'small_strain'`` (default), ``'green_lagrange'`` (PKII control),
  ``'logarithmic'`` (Kirchhoff control), ``'biot'``, or the fully kinematic
  ``'F'`` / ``'gradU'`` (9 components of the deformation gradient).
* ``solve(corate=...)`` selects the objective rate for the finite-strain
  control types: ``'jaumann'``, ``'green_naghdi'``, ``'logarithmic'`` (default),
  ``'logarithmic_R'``, ``'truesdell'``, ``'logarithmic_F'``.
* ``mode='sinusoidal'`` interpolates the step sinusoidally instead of
  linearly; ``mode='tabular'`` follows a user table passed in memory:

.. code-block:: python

    t = np.linspace(0.01, 1.0, 100)
    e11 = 0.015 * np.sin(np.pi * t)
    step = solver.StepMeca(control=['strain'] + ['zero'] * 5,
                           mode='tabular',
                           tabular=np.column_stack([t, e11]))

* Cyclic loading repeats the steps of a block: ``Block(steps=[...], ncycle=10)``.
  Tabular steps cannot be cycled (their time column is absolute); unroll the
  cycles into explicit steps instead.
* A rotation rate can be superimposed on the mixed finite-strain control
  types through ``StepMeca(BC_w=...)`` (3x3 spin matrix).

Thermomechanical loading
------------------------

:class:`~simcoon.solver.StepThermomeca` activates the coupled heat equation
(block type 2, small strain), with ``thermal_control`` set to
``'temperature'`` (ramp to ``T_final``), ``'heat_flux'`` (prescribed ``Q``) or
``'convection'`` (0D convection with coefficient ``q_conv``):

.. code-block:: python

    step = solver.StepThermomeca(control=['stress'] * 6, value=[0.] * 6,
                                 T_final=340., ninc=50)
    # ELISO thermomechanical props: rho, c_p, E, nu, alpha
    res = solver.solve(step, "ELISO", [1.E-9, 1., 70000., 0.3, 1.E-5], nstatev=1)

Mechanical steps (:class:`~simcoon.solver.StepMeca`) also accept ``T_final``:
the temperature then ramps as an imposed condition of the mechanical problem
(thermal expansion without the heat equation).

Tabular thermomechanical steps support the three thermal controls: with
``thermal_control='temperature'`` the table carries a T column when
``tabular_T=True`` (constant temperature otherwise); with ``'heat_flux'`` the
thermal column is the prescribed flux Q; with ``'convection'`` there is no
thermal column and ``q_conv`` applies.

For finite-element couplers, the point-wise thermomechanical UMAT batch entry
``sim.umat_T(...)`` complements ``sim.umat(...)``; it returns
``(sigma, statev, Wm, Wt, r, dSdE, dSdT, drdE, drdT)``.

Legacy file formats
-------------------

Existing ``path.txt`` / ``material.dat`` inputs are parsed into loading
objects — same solve, one entry point:

.. code-block:: python

    blocks, T_init = solver.from_file("data", "path.txt")
    material = solver.material_from_file("data", "material.dat")
    res = solver.solve(blocks, T_init=T_init, **material)

JSON configuration
------------------

Materials and loading paths round-trip through JSON
(:func:`~simcoon.solver.save_material_json`,
:func:`~simcoon.solver.save_path_json`,
:func:`~simcoon.solver.load_simulation_json`):

.. code-block:: python

    solver.save_material_json("material.json", "ELISO", [70000., 0.3, 1.E-5], 1)
    solver.save_path_json("path.json", [solver.Block(steps=[step])], T_init=293.15)
    res = solver.solve(**solver.load_simulation_json("material.json", "path.json"))

Results can be persisted with ``res.save("run.npz")`` /
``SolverResults.load("run.npz")``, or flattened with ``res.to_dataframe()``.

Solver parameters
-----------------

``solve()`` exposes the numeric controls of the adaptive Newton loop as
keyword arguments: ``precision`` (default 1e-6), ``maxiter``/``miniter``,
``div_tnew_dt``/``mul_tnew_dt`` (time-step cut/growth factors), ``inforce``,
``lambda_solver`` (penalty stiffness of strain-driven components), plus
``tangent_mode`` (``'none'``, ``'continuum'``, ``'algorithmic'`` — default)
and ``solver_type``.

API reference
-------------

.. autoclass:: simcoon.solver.StepMeca
   :members:

.. autoclass:: simcoon.solver.StepThermomeca
   :members:

.. autoclass:: simcoon.solver.Block
   :members:

.. autofunction:: simcoon.solver.solve

.. autofunction:: simcoon.solver.from_file

.. autofunction:: simcoon.solver.material_from_file

.. autoclass:: simcoon.solver.SolverResults
   :members:
