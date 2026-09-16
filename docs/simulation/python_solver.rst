In-memory Python solver
================================

The ``simcoon.solver`` package drives the C++ material-point solver directly
from Python: the loading path is defined with :class:`~simcoon.solver.Block`
and :class:`~simcoon.solver.StepMeca` / :class:`~simcoon.solver.StepThermomeca`
objects, and the results come back as numpy arrays — no ``path.txt``,
``output.dat`` or result files involved. Since simcoon 2.0 this package *is*
``sim.solver``, JSON is the only file format it reads, and the pre-2.0 text
inputs are converted once with ``scripts/legacy_to_json.py``.

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
    strain = res["Strain"]     # logarithmic strain history, shape (6, N)

Results follow the fedoo ``DataSet`` conventions — components first, one
column per increment — so they interoperate directly with fedoo utilities
(e.g. ``fedoo.util.voigt_tensors.StressTensorList(res["Stress"])``).
Available fields include ``Stress`` (Cauchy), ``Kirchhoff``, ``PKII``,
``Strain`` (logarithmic, alias ``LogStrain``), ``GreenLagrange``, ``F``, ``R``,
``DR`` (``(3, 3, N)``),
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
  control types: ``'jaumann'``, ``'green_naghdi'``, ``'logarithmic'``,
  ``'logarithmic_R'`` (default), ``'truesdell'``, ``'logarithmic_F'``.
  ``'logarithmic_R'`` transports by the exact polar rotation increment
  :math:`\Delta\mathbf{R} = \mathbf{R}_1\mathbf{R}_0^T`, for which the
  tangent transport is exact — including with rotated internal-variable
  history (plasticity at finite rotation).
* ``mode='sinusoidal'`` interpolates the step sinusoidally instead of
  linearly; ``mode='tabular'`` follows a user table passed in memory:

.. code-block:: python

    t = np.linspace(0.01, 1.0, 100)
    e11 = 0.015 * np.sin(np.pi * t)
    step = solver.StepMeca(control=['strain'] + ['zero'] * 5,
                           mode='tabular',
                           tabular=np.column_stack([t, e11]))

  In a saved path the table is the one input that is not JSON: ``save_path_json``
  writes it as ``<stem>_tab<k>.csv`` next to the JSON (``#`` header naming the
  columns, one row per increment) and the step's ``"tabular"`` entry holds that
  filename; ``load_path_json`` reads it back, comma- or whitespace-separated.
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

Mean-field composites
---------------------

The mean-field models (MIHEN, MIMTN, MISCN, MIPLN) take their sub-phases in memory,
as :class:`~simcoon.solver.micromechanics.Ellipsoid` or
:class:`~simcoon.solver.micromechanics.Layer` objects passed to ``solve(phases=...)``
or ``sim.L_eff(..., phases=...)`` (see :doc:`solver`, ``phases``). Every orientation
in those objects is a :class:`simcoon.Rotation`: the material frame of a phase
(``material_orientation``) and the geometry of an inclusion or a layer
(``geometry_orientation``). Any of these forms is accepted and coerced:

.. code-block:: python

    from simcoon.solver.micromechanics import Ellipsoid, as_rotation, euler_angles

    fibre = Ellipsoid(umat_name="ELISO", concentration=0.2, nstatev=1,
                      props=[50000., 0.3, 0.], a1=50.,
                      geometry_orientation=sim.Rotation.from_rotvec([0, 0, np.pi / 4]))
    fibre.geometry_orientation = (45., 0., 0.)                       # Euler angles, degrees
    fibre.geometry_orientation = {"psi": 45., "theta": 0., "phi": 0.}  # the JSON form
    euler_angles(fibre.geometry_orientation)   # {'psi': 45.0, 'theta': 0.0, 'phi': 0.0}

The Euler angles are the ``'zxz'`` sequence the C++ side reads, in degrees:
``as_rotation((psi, theta, phi))`` is
``Rotation.from_euler('zxz', [psi, theta, phi], degrees=True)`` (scipy's extrinsic
``zxz``), applied actively. A phase at that orientation responds with
``R.apply_stiffness(L_local)``, which the test suite pins against the solver and
``L_eff``. The JSON files and the dicts handed to the extension keep the angles;
``euler_angles`` writes them back, with the usual caveat that the decomposition is not
unique at ``theta = 0`` (the z angles merge into ``psi``). ``solve(orientation=...)``
and ``sim.L_eff(umat_name, props, nstatev, orientation=..., phases=...)`` take the
same forms for the frame of the material or of the whole RVE; ``L_eff`` also takes the
phase objects directly.

An orientation distribution splits one phase into phases rotated about a direction,
with concentrations following the ODF:

.. code-block:: python

    from simcoon.solver.micromechanics import Peak, discretize_odf

    peaks = [Peak(method=3, mean=90., s_dev=10.)]     # Gaussian, degrees
    phases = discretize_odf([matrix, fibre], num_phase=1, peaks=peaks, nphases=18,
                            axis=(0., 0., 1.), angle_range=(0., 180.))
    L = sim.L_eff("MIMTN", props, nstatev, phases=phases)

The k-th copy is rotated by ``alpha_k`` about ``axis`` on top of the orientations the
phase already has (``Rotation.from_rotvec(alpha_k * axis) * orientation``), the
material frame following unless ``rotate_material=False``; ``alpha_k`` runs from
``angle_min`` in ``nphases`` equal steps and each copy takes the Simpson integral of
the density over its step, normalised to the parent's concentration. The density is a
director distribution, periodic over 180 degrees: the default half turn is right when
a half turn about ``axis`` maps the inclusion onto itself (``axis`` along or normal to
a principal axis of the inclusion); about any other direction give
``angle_range=(0., 360.)``. The pre-2.0 Euler sweeps are the cases ``axis=(0, 0, 1)``
(psi or phi) and ``axis=(1, 0, 0)`` (theta) from an unrotated phase. Peak profiles
(``method``: 1 standard-deviation kernel, 2 hard cut-off, 3 Gaussian, 4 Lorentzian,
5 pseudo-Voigt, 6 Pearson VII, 7 uniform) are evaluated by ``sim.get_densities_ODF``.

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

Legacy text inputs
------------------

JSON is the only format simcoon reads or writes since 2.0: nothing in the package
parses ``path.txt``, ``material.dat``, ``tab_file_<n>.txt`` or ``N<kind><n>.dat``.
A pre-2.0 ``data`` directory is converted once with the migration script shipped in
the repository (not in the package), the tables its mode-3 steps referenced being
rewritten as ``<stem>_tab<k>.csv`` files next to the path JSON:

.. code-block:: bash

    python scripts/legacy_to_json.py data   # writes path.json, material.json, ellipsoids<N>.json, ...

.. code-block:: python

    res = solver.solve(**solver.load_simulation_json("data/material.json", "data/path.json"))

Solver parameters
-----------------

``solve()`` exposes the numeric controls of the adaptive Newton loop as
keyword arguments: ``precision`` (default 1e-6), ``maxiter``/``miniter``,
``div_tnew_dt``/``mul_tnew_dt`` (time-step cut/growth factors), ``inforce``,
``lambda_solver`` (penalty stiffness of strain-driven components), plus
``tangent_mode`` (``'none'``, ``'continuum'``, ``'algorithmic'`` — default)
and ``solver_type``.

Constitutive laws written in Python
-----------------------------------

``umat_name`` also accepts a law object implementing :class:`simcoon.PythonUMAT`
(numpy, PyTorch, ...): it is registered under the ``PYEXT`` name for the duration
of the call and integrated by the C++ solver like a built-in kernel; ``props`` and
``nstatev`` default to the object's attributes. See :doc:`python_umat`.

API reference
-------------

.. autoclass:: simcoon.solver.StepMeca
   :members:

.. autoclass:: simcoon.solver.StepThermomeca
   :members:

.. autoclass:: simcoon.solver.Block
   :members:

.. autofunction:: simcoon.solver.solve

.. autoclass:: simcoon.solver.SolverResults
   :members:
