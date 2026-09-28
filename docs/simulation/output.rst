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
    lnV = res["LogStrain"]                          # ln V from F, whatever the rate
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
     - Eulerian strain integrated along the path with the objective rate of the run (``corate``): the strain the constitutive law sees. With the logarithmic rates (``"logarithmic"``, ``"logarithmic_R"``, the default, ``"logarithmic_F"``) it is the logarithmic strain :math:`\ln \mathbf{V}` up to the integration error; with ``"jaumann"`` or ``"green_naghdi"`` it is the integral of :math:`\mathbf{D}` along the spin of the rate: path dependent, not a function of :math:`\mathbf{F}`, equal to :math:`\ln\mathbf{V}` only on paths whose principal axes do not rotate (in simple shear it departs at third order in :math:`\gamma`); with ``"truesdell"`` it is the **Almansi strain** :math:`\mathbf{e}_A = \frac{1}{2}(\mathbf{I} - \mathbf{b}^{-1})`, exactly, not a logarithmic strain
   * - ``LogStrain``
     - (6, N)
     - Logarithmic strain :math:`\ln \mathbf{V} = \frac{1}{2}\ln(\mathbf{F}\mathbf{F}^T)`, computed from :math:`\mathbf{F}`: exact whatever the rate. Use it, not ``Strain``, when the logarithmic strain is wanted under ``"truesdell"``, ``"jaumann"`` or ``"green_naghdi"``
   * - ``GreenLagrange``
     - (6, N)
     - Green-Lagrange strain :math:`\mathbf{E} = \frac{1}{2}(\mathbf{F}^T\mathbf{F} - \mathbf{I})`, from :math:`\mathbf{F}`: exact whatever the rate
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
     - Rotation accumulated by the objective rate of the run; it is the :math:`\mathbf{R}` of the polar decomposition :math:`\mathbf{F} = \mathbf{R}\mathbf{U}` for ``"green_naghdi"`` and ``"logarithmic_R"``, differs from it for ``"jaumann"`` and ``"logarithmic"``, and for the convected rates ``"truesdell"`` and ``"logarithmic_F"`` it is not a rotation: it accumulates the frame increments :math:`\Delta\mathbf{F}`, i.e. it is :math:`\mathbf{F}` itself (to round-off: runs start at :math:`\mathbf{F} = \mathbf{I}`)
   * - ``DR``
     - (3, 3, N)
     - Frame increment of the objective rate over the increment: the rotation :math:`\Delta\mathbf{R}`, or :math:`\Delta\mathbf{F} = \mathbf{F}_1\mathbf{F}_0^{-1}` for ``"truesdell"`` and ``"logarithmic_F"``

.. note::

   Internally the finite-strain solver carries the Kirchhoff stress
   :math:`\boldsymbol{\tau}`: the constitutive kernels and the energies
   :math:`W_m` work on it, and the prescribed PKII and Biot stresses are
   derived from it. ``Kirchhoff`` is that stress as integrated; ``Stress`` is
   the Cauchy stress :math:`\boldsymbol{\sigma} = \boldsymbol{\tau}/J` formed
   at output and stays the default measure. See
   :ref:`stress-measure-tangent-rate` for the kernel contract and the
   ``sim.umat`` boundary used by fedoo.

.. note::

   In small deformations (``control_type="small_strain"``) the rate plays no
   role: all strain measures reduce to the infinitesimal strain and all stress
   measures to the Cauchy stress. Shear strain
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
     - Mechanical energies :math:`[W_m, W_m^r, W_m^{ir}, W_m^d]`: total, stored (recoverable), irrecoverable stored, dissipated, per reference volume. :math:`W_m` is the work :math:`\int \mathbf{P} : \mathrm{d}\mathbf{F}` whatever the rate (see :ref:`stress-measure-tangent-rate`)
   * - ``Statev``
     - (nstatev, N)
     - The internal state variables of the constitutive model, in the order the model defines them (see :doc:`umat_catalog`). Under finite strain, the tensorial ones of the kernels fed the logarithmic strain are components in the frame that follows the material (the material axes rotated with the body), not in the lab frame
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

``to_dataframe`` flattens the scalar histories and every 2-D history: the 6-component ones become ``<key>_11`` ...
``<key>_23`` columns, ``Wm`` and
``Statev`` become ``Wm_0`` ... and ``Statev_0`` ...; the ``(3, 3, N)`` and
``(6, 6, N)`` histories are left out of the table and read from ``res`` directly.

.. note::

   Pre-2.0, an ``output.dat`` file in ``data/`` selected the measures written
   to a result text file. Nothing reads it any more: everything it could select
   is in ``SolverResults``, and the C++ solver writes no file.
