Read the results
================================

``sim.solver.solve`` writes no file. It returns a :class:`~simcoon.solver.SolverResults`
holding the whole history in memory, one column per recorded increment, and
Python picks the measures it wants from it. There is nothing to configure
beforehand: every measure below is recorded at every increment, the tangent
being the only optional one (``solve(record_tangent=False)`` skips it).

.. code-block:: python

    res = sim.solver.solve(blocks, umat_name, props, nstatev, T_init=T_init)

    e11, e22, e33, e12, e13, e23 = res["Strain"]   # strain integrated with the objective rate, (6, N)
    s11, s22, s33, s12, s13, s23 = res["Stress"]   # Cauchy stress, (6, N)
    time, T = res["Time"], res["Temp"]             # (N,)
    Wm, Wm_r, Wm_ir, Wm_d = res["Wm"]              # mechanical energies, (4, N)

    print(res)          # SolverResults(mechanical, 100 increments, status=0, fields=[...])
    res.keys()          # the names available in this run

Layout
------

Arrays follow the fedoo ``DataSet`` conventions: components first, one column
per increment. A 6-component history is a ``(6, N)`` array in Voigt order
[11, 22, 33, 12, 13, 23], a tensor history is ``(3, 3, N)``, a scalar history
``(N,)``. ``fedoo.util.voigt_tensors.StressTensorList(res["Stress"])`` works
directly on the returned arrays, and ``res.get_data(key)`` is the fedoo-style
alias of ``res[key]``.

Increment bookkeeping
---------------------

.. list-table::
   :header-rows: 1
   :widths: 15 15 70

   * - Key
     - Shape
     - Content
   * - ``Block``
     - (N,)
     - Block number of the increment (from 0)
   * - ``Cycle``
     - (N,)
     - Cycle of the block (from 0)
   * - ``Step``
     - (N,)
     - Step within the block (from 0)
   * - ``Inc``
     - (N,)
     - Increment within the step (from 0)
   * - ``Time``
     - (N,)
     - Simulation time at the end of the increment
   * - ``Temp``
     - (N,)
     - Temperature :math:`T` at the end of the increment

Strain and stress measures
--------------------------

Every measure is recorded, whatever the control type of the block, so a single
run gives every conjugate pair:

.. list-table::
   :header-rows: 1
   :widths: 18 12 70

   * - Key
     - Shape
     - Content
   * - ``Strain``
     - (6, N)
     - Eulerian strain integrated along the path with the objective rate of the run (``corate``). With the logarithmic rates (``"logarithmic"``, ``"logarithmic_R"``, the default) it is the logarithmic strain :math:`\boldsymbol{\varepsilon} = \ln \mathbf{V}`; with ``"jaumann"`` or ``"green_naghdi"`` it departs from it under large rotations (simple shear). ``LogStrain`` is the same array
   * - ``GreenLagrange``
     - (6, N)
     - Green-Lagrange strain :math:`\mathbf{E} = \frac{1}{2}(\mathbf{F}^T\mathbf{F} - \mathbf{I})`
   * - ``F``
     - (3, 3, N)
     - Deformation gradient :math:`\mathbf{F}`
   * - ``Stress``
     - (6, N)
     - Cauchy stress :math:`\boldsymbol{\sigma}`
   * - ``Kirchhoff``
     - (6, N)
     - Kirchhoff stress :math:`\boldsymbol{\tau} = J\boldsymbol{\sigma}`
   * - ``PKII``
     - (6, N)
     - 2nd Piola-Kirchhoff stress :math:`\mathbf{S}`
   * - ``R``
     - (3, 3, N)
     - Rotation accumulated by the objective rate of the run; it is the :math:`\mathbf{R}` of the polar decomposition :math:`\mathbf{F} = \mathbf{R}\mathbf{U}` for ``"green_naghdi"`` and ``"logarithmic_R"``, and differs from it for the other rates
   * - ``DR``
     - (3, 3, N)
     - Rotation increment :math:`\Delta\mathbf{R}` of the objective rate over the increment

.. note::

   In small deformations (``control_type="small_strain"``) all strain measures
   reduce to the infinitesimal strain and all stress measures to the Cauchy
   stress. Shear strain
   components are engineering shears :math:`\gamma_{ij} = 2\varepsilon_{ij}`.

Energies, state variables, tangent
----------------------------------

.. list-table::
   :header-rows: 1
   :widths: 18 12 70

   * - Key
     - Shape
     - Content
   * - ``Wm``
     - (4, N)
     - Mechanical energies :math:`[W_m, W_m^r, W_m^{ir}, W_m^d]`: total, stored (recoverable), irrecoverable stored, dissipated
   * - ``Statev``
     - (nstatev, N)
     - The internal state variables of the constitutive model, in the order the model defines them (see :doc:`umat_catalog`)
   * - ``TangentMatrix``
     - (6, 6, N)
     - Tangent operator :math:`\mathbf{L}_t` of the mechanical problem (``record_tangent=True``, the default)

Thermomechanical runs
---------------------

A block of :class:`~simcoon.solver.StepThermomeca` steps adds the thermal
histories and, with ``record_tangent``, the coupled tangents in place of
``TangentMatrix``:

.. list-table::
   :header-rows: 1
   :widths: 18 12 70

   * - Key
     - Shape
     - Content
   * - ``Q``
     - (N,)
     - Heat flux prescribed on the RVE
   * - ``r``
     - (N,)
     - Heat source
   * - ``Wt``
     - (3, N)
     - Thermal energies :math:`[W_t, W_t^r, W_t^{ir}]`
   * - ``dSdE``
     - (6, 6, N)
     - :math:`\partial\boldsymbol{\sigma}/\partial\boldsymbol{\varepsilon}`
   * - ``dSdT``
     - (6, N)
     - :math:`\partial\boldsymbol{\sigma}/\partial T`
   * - ``drdE``
     - (6, N)
     - :math:`\partial r/\partial\boldsymbol{\varepsilon}`
   * - ``drdT``
     - (N,)
     - :math:`\partial r/\partial T`

``Q`` and ``r`` are of opposite signs: the solver enforces the 0D heat balance
of the RVE, :math:`Q = -r`, at every increment. ``r`` is the heat rate the
constitutive model produces (thermomechanical couplings and dissipation, net of
the heat stored as :math:`\rho c_p \dot T`) and ``Q`` the flux exchanged with
the surroundings that balances it. An adiabatic step (``thermal_control="heat_flux"``,
``Q = 0``) therefore reads ``r = 0``, the temperature rise absorbing the sources.

Completion status
-----------------

``res.status`` is 0 when the whole path ran, 1 when the solver aborted early
(Newton loop not converging at the minimal increment). By default ``solve``
raises a ``RuntimeError`` on abort; ``solve(raise_on_abort=False)`` returns the
partial history instead, so it can be inspected up to the failing increment.

Saving and tabulating
---------------------

.. code-block:: python

    res.save("run.npz")                              # compressed numpy archive
    res = sim.solver.SolverResults.load("run.npz")   # same object back

    df = res.to_dataframe()                          # pandas DataFrame
    df[["Time", "Strain_11", "Stress_11"]].to_csv("run.csv", index=False)

``to_dataframe`` flattens the scalar histories and every 2-D history (the
``LogStrain`` alias excepted): the 6-component ones become ``<key>_11`` ...
``<key>_23`` columns, ``Wm`` and
``Statev`` become ``Wm_0`` ... and ``Statev_0`` ...; the ``(3, 3, N)`` and
``(6, 6, N)`` histories are left out of the table and read from ``res`` directly.

.. note::

   Pre-2.0, an ``output.dat`` file in ``data/`` selected the measures written
   to a result text file. Nothing reads it any more: everything it could select
   is in ``SolverResults``, and the C++ solver writes no file.
