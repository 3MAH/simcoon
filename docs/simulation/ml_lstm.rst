======================================================
Recurrent neural network constitutive models
======================================================

:mod:`simcoon.ml` implements recurrent constitutive cells: networks that map a
strain history to the stress, whose internal state plays the role of the internal
variables. Two cells share one interface, one UMAT wrapper and one fedoo law:

* :class:`simcoon.ml.StressLSTM`, the stress LSTM of Danoun, Prulière and
  Chemisky [Danoun2022]_, [Danoun2024]_ (reference PyTorch implementation:
  *StressLSTM*): a gated network, the most expressive of the two;
* :class:`simcoon.ml.LMSC`, the linearized minimal state cell of Bonatti and
  Mohr [BonattiMohr2022]_: stationary, self-consistent and rate-independent by
  construction, with a closed-form algorithmic tangent and a state two orders of
  magnitude smaller.

A cell is trained in PyTorch and then used as a constitutive law by the simcoon
solver, by ``sim.umat`` and by finite-element codes.

.. code-block:: bash

   conda install -c conda-forge pytorch scikit-learn   # conda environment (one OpenMP runtime)
   pip install simcoon[ml]                              # all-PyPI setup

(both optional at inference; never ``pip install torch`` inside a conda environment, see
:doc:`../installation`).

Architecture
============

StressLSTM
----------

A stacked ``torch.nn.LSTM`` (2 layers, 64 units by default) followed by a linear
head; inputs and outputs are standardised per component
(:math:`\hat{x}_c = (x_c - \mu_c)/s_c`). At each time step

.. math::

   \mathrm{LSTM}(\boldsymbol{\varepsilon}_t, \mathbf{h}_{t-1}, \mathbf{c}_{t-1})
   \mapsto (\boldsymbol{\sigma}_t, \mathbf{h}_t, \mathbf{c}_t).

The strain components fed to the network (``components``) are any subset of the
Voigt vector ``xx yy zz xy xz yz`` (engineering shear), so a 2D model has 3 inputs
and a 3D one 6. Optional input features: ``"dstrain"`` (strain increment),
``"dtime"``, ``"temperature"`` — recommended as soon as the time step may vary
(solver step cuts, adaptive FE increments), since the pure StressLSTM model
implicitly assumes the sampling of its training sequences.

A gated transition function is written
:math:`\boldsymbol{\chi}' = \mathbf{a}(\boldsymbol{\chi}, \Delta\boldsymbol{\varepsilon})
\odot \boldsymbol{\chi} + \mathbf{b}(\boldsymbol{\chi}, \Delta\boldsymbol{\varepsilon})`.
Leaving the state unchanged at a zero increment requires
:math:`\boldsymbol{\chi} \odot (\mathbf{1} - \mathbf{a}(\boldsymbol{\chi}, \mathbf{0}))
= \mathbf{b}(\boldsymbol{\chi}, \mathbf{0})`, and nothing in an LSTM or a GRU forces
that identity [BonattiMohr2022]_. Such a cell therefore loses information on zero-norm
increments and its response depends on how finely a path is discretised. The committed-state
rule of the law (``commit_tol`` below) *bounds* that dependence; the LMSC removes it.

LMSC (linearized minimal state cell)
------------------------------------

For a state :math:`\boldsymbol{\chi}` and a strain increment
:math:`\Delta\boldsymbol{\varepsilon}` of norm :math:`\nu` and direction
:math:`\mathbf{n}` [BonattiMohr2022]_, Eqs. (20)-(26):

.. math::

   \mathbf{l}_0 &= [\boldsymbol{\chi} ; \mathbf{n}], \qquad
   \mathbf{l}_i = \tanh(\mathbf{W}^\alpha_i \mathbf{l}_{i-1} + \mathbf{b}^\alpha_i)
                 \odot \tanh(\mathbf{W}^\beta_i \mathbf{l}_{i-1} + \mathbf{b}^\beta_i) \\
   \boldsymbol{\alpha} &= \exp(\mathbf{W}_A \mathbf{l}_d + \mathbf{b}_A), \qquad
   \boldsymbol{\beta} = \tanh(\mathbf{W}_B \mathbf{l}_d + \mathbf{b}_B) \\
   \boldsymbol{\chi}' &= \exp(-\nu\boldsymbol{\alpha}) \odot
                          (\boldsymbol{\chi} - \boldsymbol{\beta}) + \boldsymbol{\beta},
   \qquad \boldsymbol{\sigma} = \mathbf{W}_\sigma \boldsymbol{\chi}'

Each state component relaxes exponentially towards its target
:math:`\beta_i` at the rate :math:`\alpha_i`, along the arc length of the strain
path. What that buys, by construction and not by training:

* **stationarity is exact**: :math:`\nu = 0` gives :math:`\boldsymbol{\chi}' =
  \boldsymbol{\chi}` bit for bit, so a held strain gives a held stress;
* **self-consistency is exact at frozen coefficients**, since
  :math:`\exp(-\nu_1\boldsymbol{\alpha})\exp(-\nu_2\boldsymbol{\alpha}) =
  \exp(-(\nu_1+\nu_2)\boldsymbol{\alpha})`: splitting an increment changes nothing
  as long as :math:`\boldsymbol{\alpha}` and :math:`\boldsymbol{\beta}` do not vary
  along it, and the residual error decreases under refinement;
* **rate independence is structural**: no time increment enters the equations;
* **the state is minimal**: ``n_state = 6`` is the theoretical minimum for 3D von Mises
  plasticity and the size the authors deploy (``depth=3``, ``width=25``, 3.5 k parameters,
  against 175 k for a 64x2 LSTM).

The total strain is *not* an input, only the state and the increment direction, and the
output map is linear and unbiased so that :math:`\boldsymbol{\chi} = \mathbf{0}` is
exactly the stress-free state. The authors report that a minimal state is only reachable
when training on long sequences of small increments; short sequences of large increments
need excess state variables.

Workflow
========

1. Data
-------

Random non-proportional strain paths (four linear segments of 25 increments
towards uniformly drawn targets, as in the reference datasets) integrated with
any simcoon model, or the StressLSTM CSV format:

.. code-block:: python

   from simcoon import ml

   targets, ninc = ml.random_strain_paths(1000, n_segments=4, n_sub=25, seed=0)
   ds = ml.generate_dataset("EPICP", [70000., 0.3, 1e-5, 300., 1000., 0.3], 8,
                            targets=targets, ninc=ninc, mode="3D")
   train_ds, test_ds = ml.split_dataset(ds, test_size=0.3, seed=0)

   ds2d = ml.load_csv("dataset/train_dataset.csv")          # StressLSTM format

``mode`` follows the ``ndi`` convention: ``"3D"`` and ``"plane_strain"`` are
strain-driven (``ndi = 3``), ``"plane_stress"`` drives the out-of-plane stress to
zero (``ndi = 2``) and ``"uniaxial"`` all stresses but the axial one (``ndi = 1``),
so the data satisfy the constraint the model will be used with.

2. Training (PyTorch or scikit-learn API)
-----------------------------------------

The loss is :func:`simcoon.ml.torch_cost`, the differentiable twin of the
identification cost :func:`simcoon.identify.calc_cost` — same three-level
weights (``w_test``, ``w_response``, ``w_point``), same metric names (``mse``,
``nmse``, ``nmse_per_response``, ``rmse``, ``mae``, ``mape``, ``wmape``); the
padding mask is a null point weight.

.. code-block:: python

   model = ml.StressLSTM(hidden_size=64, num_layers=2)          # 3D, features=("strain",)
   train_losses, val_losses = ml.train(model, train_ds, test_ds, epochs=2000,
                                       batch_size=64, lr=1e-3, loss="mse")
   model.save("lstm_epicp.pt")

or, through scikit-learn (``GridSearchCV``, ``cross_val_score``, pipelines):

.. code-block:: python

   reg = ml.StressLSTMRegressor(hidden_size=64, epochs=2000, loss="nmse", random_state=0)
   reg.fit(train_ds.x.numpy(), train_ds.y.numpy())
   reg.score(test_ds.x.numpy(), test_ds.y.numpy())              # R²

Nothing in a loss on the stress constrains the *derivative* the solver consumes: measured
on a J2 + Voce surrogate, a gated model within 17 % of the reference stress is off by 70 %
on the tangent, and does not recover the elastic operator in the elastic regime. Adding a
tangent term to the loss was tried and abandoned — it reaches its target (median tangent
error 15 % to 4.5 %) but costs stress accuracy at a fixed budget, runs six times slower per
epoch, and needs a second derivative of the recurrent backward that the Apple MPS backend
does not provide. The LMSC answers the same need by construction, its tangent being a
closed form of the update rather than a fitted by-product.

:func:`simcoon.ml.sequence_tangent` remains available as a **diagnostic**: it returns the
block-diagonal ``d sigma_t / d eps_t`` of a cell at frozen state for a whole batch of
sequences, to be compared with the operator the solver recorded
(``generate_dataset(..., record=("tangent",))``), which is the analytic total derivative
:math:`\mathbf{D}^\varepsilon = \mathbf{L} - \sum_j \boldsymbol{\kappa}^j \mathbf{P}^j_\varepsilon`
of the theory manual — the internal variables are eliminated through the consistency
condition, not held fixed.

3. Evaluation
-------------

:func:`simcoon.ml.evaluate` de-standardises the predictions and scores them
with :func:`simcoon.identify.calc_cost`, globally, per component and per
sequence:

.. code-block:: python

   report = ml.evaluate(model, test_ds, metrics=("mse", "nmse", "r2", "mape", "wmape"),
                        per_component=True, per_sequence=True)

``wmape`` (:math:`\sum|e| / \sum|y|`) is the percentage error to prefer for
stress histories that cross zero; ``mape`` follows the scikit-learn definition.

The tangent of the LMSC, in closed form
---------------------------------------

Differentiating the update above at frozen incoming state, with
:math:`\mathbf{e} = \exp(-\nu\boldsymbol{\alpha})`,
:math:`\mathbf{u} = \boldsymbol{\chi}-\boldsymbol{\beta}` and
:math:`\varphi(t) = (1-e^{-t})/t`:

.. math::

   \frac{\partial\boldsymbol{\chi}'}{\partial\boldsymbol{\varepsilon}}
   = -\mathrm{diag}(\mathbf{e}\odot\mathbf{u})
     \Big[\boldsymbol{\alpha}\,\mathbf{n}^T
     + \frac{\partial\boldsymbol{\alpha}}{\partial\mathbf{n}}
       (\mathbf{I}-\mathbf{n}\mathbf{n}^T)\Big]
   + \mathrm{diag}\big(\boldsymbol{\alpha}\odot\varphi(\nu\boldsymbol{\alpha})\big)
     \frac{\partial\boldsymbol{\beta}}{\partial\mathbf{n}}
     (\mathbf{I}-\mathbf{n}\mathbf{n}^T)

and the operator the Newton loop consumes is
:math:`\mathbf{W}_\sigma\,\partial\boldsymbol{\chi}'/\partial\boldsymbol{\varepsilon}`.
This is the **total** derivative, a partial derivative plus a state-sensitivity term, the
same structure as the analytical tangents of the classical kernels
(:math:`\mathbf{D}^\varepsilon = \mathbf{L} - \sum_j \kappa^j \mathbf{P}^j_\varepsilon`).
The three Jacobians involved are those of a feed-forward network with respect to six
inputs, so ``torch.func.jacrev`` gives them exactly, on CPU as on MPS, second derivatives
included — nothing is differentiated through the recurrence.

The :math:`\varphi` form keeps the expression regular as :math:`\nu \to 0`, and the
limit stays direction-dependent, which is the correct behaviour of a rate-independent law
(elastic stiffness on unloading, plastic tangent on loading). At a strictly zero increment
there is no direction and the law returns the elastic operator, the way a classical UMAT
answers a zero trial increment.

:func:`simcoon.ml.sequence_tangent` computes the same block-diagonal operator by autograd
over a whole batch of sequences, which is how the closed form is checked and how a gated
cell's tangent is measured.

4. A cell as a constitutive law
-------------------------------

:class:`simcoon.ml.RecurrentLaw` (exported as ``LSTMLaw`` too) wraps any trained cell as a
:class:`simcoon.PythonUMAT` (see :doc:`python_umat`):

.. code-block:: python

   law = ml.LSTMLaw(model)                     # or ml.LSTMLaw.load("lstm_epicp.pt")
   res = solve(StepMeca(control=["strain"] + ["stress"] * 5,
                        value=[0.02, 0, 0, 0, 0, 0], ninc=100), law)

* The cell's flat state lives in ``statev`` (``nstatev = state_size + n_components``,
  the last block being the committed strain): the solver's rollback rewinds the network
  correctly on Newton retrials and step cuts. ``state_size`` is
  ``2 * num_layers * hidden_size`` for the LSTM (256 by default) and ``n_state`` for the
  LMSC (6 to 20).
* The tangent is the cell's closed form when it has one (the LMSC), and the autograd
  Jacobian of the step at frozen state otherwise — in both cases the algorithmic tangent
  of the model, which drives the Newton loop under mixed stress/strain control.
  ``tangent_mode = 0`` returns the elastic operator identified at the origin.
* ``ndi``: a model trained on the reduced components (``mode="plane_stress"``,
  ``"uniaxial"``) is used as is. A 3D model used with ``ndi = 2`` or ``1`` is
  **condensed**: a local Newton iteration finds the strain of the stress-free
  directions (:math:`\sigma_S = 0`), the tangent
  returned is :math:`\mathbf{L}_{FF} - \mathbf{L}_{FS}\mathbf{L}_{SS}^{-1}\mathbf{L}_{SF}`
  (the algebra of ``el_pred(ndi=2)`` and of fedoo's ``get_H_plane_stress``) and
  the state is advanced once, with the converged strain.
* Energies: ``Wm`` is accumulated by the trapezoidal rule; ``Wm_r``, ``Wm_ir``,
  ``Wm_d`` are zero (no free energy in the base model).

5. Finite-element coupling (fedoo)
----------------------------------

:meth:`simcoon.ml.LSTMLaw.step_batch` evaluates all Gauss points at once
(``strain (6, N)``, ``h``/``c (state_size, N)`` → ``stress (6, N)``,
``Lt (6, 6, N)``, new state) on CPU or GPU. fedoo ships the corresponding law,
``fedoo.constitutivelaw.LSTMLaw(law)`` (a ``Mechanical3D``): it keeps ``h``/``c`` in
``assembly.sv`` (``'LSTM_h'``, ``'LSTM_c'``, read from ``sv_start`` and written as new
arrays, so the assembly rollback rewinds the network), calls ``step_batch`` in
``update``, passes ``ndi = 2`` for ``2Dstress`` and declares the corotational box
tangent like the ``Simcoon`` law:

.. code-block:: python

   import fedoo as fd
   material = fd.constitutivelaw.LSTMLaw(ml.LSTMLaw.load("lstm_epicp.pt"))
   wf = fd.weakform.StressEquilibrium(material)

The serial ``sim.umat("PYEXT", ...)`` route also works but is slower at mesh scale.

Limits and options
==================

* **Time discretisation**: a recurrent surrogate of a rate-independent law must be
  *self-consistent* (the same strain path described with finer steps gives the same
  response) and *stationary* (a held strain leaves stress and state unchanged)
  [BonattiMohr2022]_. The **LMSC satisfies both by construction** and needs nothing else
  (``commit_tol`` defaults to 0 for it). A **gated cell satisfies neither**, and two
  mechanisms bound the departure: **committed-state inference** in the law — the state
  advances only when the strain moved by more than ``commit_tol`` (default 0.1 x the
  median training increment, recorded by ``fit_scalers``) since the last committed strain,
  smaller moves getting a trial stress from the committed state, exactly the trial-state
  contract of an implicit UMAT — and, at training time, resampling the paths at several
  increment sizes (:func:`simcoon.ml.data.resample_time`, or simply regenerating the same
  targets at several ``ninc``), as in the LSTM-GNN paper [GuevaraGarban2026]_. Without
  them a gated model trained on a fixed discretisation drifts as soon as the solver cuts
  or refines its increments (measured on J2 + Voce cycles: a refinement by 16 turned a
  45 % strain error into 55 %; the committed-state rule alone brings the refined case to
  5 %). Below the trained step range the accuracy still degrades: subdivide such inputs
  rather than extrapolate.
  Measured on J2 + Voce at an identical budget, running each model against *itself* on the
  same path described with 50 to 800 increments (spread of the answer, in per cent):

  .. list-table::
     :header-rows: 1
     :widths: 34 22 22

     * - model
       - fixed strain path
       - stress-driven cycle
     * - ``MODUL`` J2 + Voce (reference law)
       - 0.00
       - 0.00
     * - LMSC, 20 state variables
       - 0.01
       - 0.06
     * - StressLSTM (committed-state rule on)
       - 0.09
       - 65.4

  The committed-state rule makes a gated cell self-consistent under *strain* control; under
  *stress* control, where the solver chooses the strain increment, it still spreads by more
  than half the amplitude. Below the trained increment range every cell degrades, the LMSC
  included: subdivide such inputs rather than extrapolate. Note also that scoring a model
  against the reference law under refinement mixes fit quality with discretisation
  dependence and is *not* a self-consistency test.

* **Amplification of the solver residual** (gated cells). The recurrent state is
  *expanding* along a loading path, unlike the internal variables of a classical model,
  which contract towards the yield surface. Feeding a strain-driven run's own stress
  history back as a stress-driven one (inverting the model's own map, which must return
  the same strain) measures it: the reference kernel closes the loop to ``1e-12``, the
  LSTM to ``2e-4`` at the default Newton tolerance, because each increment multiplies the
  residual seeded by the solver by about 1.3, i.e. ``1e9`` over 100 increments. The error
  is strictly proportional to that tolerance (``precision=1e-12`` closes the loop to
  ``6e-12``), so it is not a modelling error: **drive a gated recurrent law with a tighter
  Newton tolerance than a classical one** (``solve(..., precision=1e-9)``, and the
  equivalent in the FE solver). Committed-state inference cuts the amplification by a
  further factor of about 700. What this buys is reproducibility — the same path,
  replayed, restarted or reached after a step cut, gives the same answer; it does not
  improve the agreement with a reference law, which is set by the training error alone.
* **The elastic-plastic knee is smoothed** (LMSC). The linearization behind Eq. (25)
  neglects the second derivatives of :math:`\boldsymbol{\alpha}` and
  :math:`\boldsymbol{\beta}`, and no increment is small enough for that at the onset of
  plastic flow, where the rates of the physical state variables are discontinuous
  ([BonattiMohr2022]_, Section 7.2). Expect a rounded yield knee, and better accuracy on
  smooth paths than on random walks.
* **Objectivity**: small-strain models. Under NLGEOM they receive the corotational
  logarithmic strain (log_R route) like any small-strain kernel, but the recurrent
  state cannot be rotated by ``DR``: use them for moderate rotations.
* **Thermodynamic consistency** ([Danoun2022]_): ``StressLSTM(psi_head=True)`` adds
  a free-energy head and ``train(..., dissipation_weight=w)`` penalises a negative
  discrete dissipation :math:`\boldsymbol{\sigma}_t\!:\!\Delta\boldsymbol{\varepsilon}_t - \Delta\psi_t < 0`.
  This option is provided as is and not yet validated.

.. [Danoun2022] A. Danoun, E. Prulière, Y. Chemisky, *Thermodynamically consistent
   Recurrent Neural Networks to predict non linear behaviors of dissipative materials
   subjected to non-proportional loading paths*, Mechanics of Materials 173 (2022) 104436.
.. [Danoun2024] A. Danoun, E. Prulière, Y. Chemisky, *FE-LSTM: A hybrid approach to
   accelerate multiscale simulations of architectured materials using Recurrent Neural
   Networks and Finite Element Analysis*, CMAME 429 (2024) 117192.
.. [BonattiMohr2022] C. Bonatti, D. Mohr, *On the importance of self-consistency in
   recurrent neural network models representing elasto-plastic solids*, J. Mech. Phys.
   Solids 158 (2022) 104697.
.. [GuevaraGarban2026] M. R. Guevara Garban, E. Prulière, Y. Chemisky, *Non-linear
   mechanical field reconstruction coupling recurrent neural networks with
   physics-informed graph neural networks* (LSTM-GNN, revised manuscript).

Example
=======

* :ref:`sphx_glr_examples_ml_plot_lstm_epicp.py` — train an LSTM on EPICP paths
  and run it in the solver under strain and stress control.

API reference
=============

.. autoclass:: simcoon.ml.cells.StateModel
   :members:

.. autoclass:: simcoon.ml.StressLSTM
   :members:

.. autoclass:: simcoon.ml.LMSC
   :members:

.. autoclass:: simcoon.ml.RecurrentLaw
   :members:

.. autoclass:: simcoon.ml.StressLSTMRegressor
   :members:

.. autofunction:: simcoon.ml.random_strain_paths

.. autofunction:: simcoon.ml.generate_dataset

.. autofunction:: simcoon.ml.load_csv

.. autofunction:: simcoon.ml.split_dataset

.. autoclass:: simcoon.ml.SequenceDataset
   :members:

.. autofunction:: simcoon.ml.torch_cost

.. autofunction:: simcoon.ml.sequence_tangent

.. autofunction:: simcoon.ml.train

.. autofunction:: simcoon.ml.evaluate
