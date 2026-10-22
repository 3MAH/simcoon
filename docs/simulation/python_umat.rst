=====================================
Constitutive laws in Python (PYEXT)
=====================================

A constitutive law can be written in Python and driven by the C++ material-point
solver exactly like a built-in kernel. The C++ side keeps the Newton loop, the
step cutting, the start-of-increment rollback and the finite-strain kinematics;
the Python side only integrates **one increment at one material point**. The law
is served under the UMAT name ``PYEXT`` and reaches the C++ code through a
process-wide callback registered by the bindings.

Typical uses: a research model prototyped in numpy, a machine-learning model
(PyTorch, JAX), or any law for which a C++ port is not worth
the effort yet.

Quick start
===========

.. code-block:: python

   import numpy as np
   import simcoon as sim
   from simcoon.solver import StepMeca, solve

   class LinearElastic(sim.PythonUMAT):
       nstatev = 0                      # no internal variable

       def __init__(self, E, nu):
           self.L = sim.L_iso([E, nu], "Enu")

       def integrate(self, *, Etot, DEtot, sigma, statev, Wm, **kw):
           stress = self.L @ (Etot + DEtot)
           Wm = Wm.copy()
           Wm[0] += 0.5 * (sigma + stress) @ DEtot     # cumulative work
           return stress, self.L, statev, Wm            # sigma, Lt, statev, Wm

   step = StepMeca(control=["strain"] + ["stress"] * 5, value=[0.01, 0, 0, 0, 0, 0], ninc=50)
   res = solve(step, LinearElastic(70000., 0.3))       # the law object replaces the name
   res["Stress"][0, -1]                                 # 700.0

``solve`` registers the object under ``PYEXT`` for the duration of the call;
``props`` and ``nstatev`` are taken from the object (``props`` may stay empty).

The contract
============

:class:`simcoon.PythonUMAT` mirrors the C++ small-strain UMAT convention. The
``integrate`` method receives keyword arguments and returns a tuple:

.. code-block:: python

   def integrate(self, *, Etot, DEtot, sigma, DR, props, statev, T, DT, Time, DTime,
                 Wm, ndi, nshr, start, tangent_mode, **_):
       ...
       return sigma, Lt, statev, Wm            # optionally a 5th value: L (elastic operator)

* **Voigt convention**: order ``11 22 33 12 13 23``, **engineering** shear strains,
  quantities in the material (local) frame. ``Etot`` is the strain at the beginning
  of the increment, ``DEtot`` the increment, ``sigma`` the stress at the beginning of
  the increment; ``DR`` is the rotation increment (identity in small strain).
* **State**: the law must be a pure function of its arguments except through
  ``statev`` and ``Wm``. The solver calls the same increment several times (Newton
  iterations, step cuts) with ``statev``, ``sigma`` and ``Wm`` **reset to their
  start-of-increment values** — nothing may be cached in the Python object between
  calls. A recurrent model keeps its hidden state inside ``statev``.
* **start** is ``True`` on the first call of a block (``Time == 0``); ``statev``
  arrives zero-filled — initialise it there if needed. The solver primes the
  tangent at the start of **every** block with a zero-increment call
  (``DTime == 0``, ``DEtot == 0``), and finite-element couplers do the same at
  initialisation; :class:`simcoon.PythonUMAT` treats that probe as a pure tangent
  query: stress, ``statev`` and ``Wm`` are returned exactly as received and only
  the tangent comes from the evaluation, so a history-dependent law does not count
  it as a loading step (a bare callable registered instead of a ``PythonUMAT`` must
  do it itself).
* **ndi** follows the classical convention: 3 = 3D (plane strain is 3D with a null
  out-of-plane strain), 2 = plane stress, 1 = uniaxial. ``nshr`` is the number of
  shear components kept, in the order 12, 13, 23 (Abaqus). With ``ndi < 3`` every
  component outside the ``(ndi, nshr)`` state is stress-free: ``(1, 0)`` truss,
  ``(1, 2)`` beam fibre (11, 12, 13), ``(2, 1)`` plane stress.
  ``sim.umat(..., ndi=1, nshr=2)`` selects a state other than the default of ``ndi``.
  The law itself is always called 3D (``ndi = 3``): ``sim.umat`` solves the
  stress-free strain components around it and condenses its tangent, so a Python
  law never condenses anything. The solved strain comes back in ``Detot`` (in
  place; the caller carries it with ``etot += Detot``). Under finite strain,
  ``condense='spatial_face'`` (default) frees the stress the law returns and
  ``condense='material_face'`` the second Piola-Kirchhoff stress, i.e. the traction on
  the free material faces (see :func:`simcoon.umat`).
  The two conditions differ by the stress measure and the frame they are written
  in. ``'spatial_face'`` cancels the law's own stress on the free components:
  :math:`\hat{\boldsymbol{\tau}}_C` in the rotated frame under the rate form,
  :math:`\mathbf{T}_C` under the total one. The face it frees has the normal
  :math:`\mathbf{R}\mathbf{e}_3`, whereas the material face has the current
  normal :math:`\mathbf{F}^{-T}\mathbf{e}_3`. The two coincide when the free
  directions are principal for :math:`\mathbf{U}` (truss, plane stress without
  transverse shear, isotropic material). Otherwise ``'spatial_face'`` leaves a
  traction on the material face: :math:`S_{33} \simeq -\mu\gamma_t^2` for a shell
  under a transverse shear :math:`\gamma_t`, i.e. of second order in the shears
  that couple kept and free directions. It is of first order for anisotropy axes
  tilted out of the free planes. ``'material_face'`` cancels
  :math:`\mathbf{S}\mathbf{e}_3`, the pull-back of
  :math:`\boldsymbol{\sigma}\mathbf{n} = \mathbf{0}`, i.e. a traction-free
  material face.
  Which face is free depends on the element:

  - **plate or shell in plane stress** (``ndi = 2``): the top and bottom faces
    are material surfaces, and ``'material_face'`` is the traction-free
    condition;
  - **beam fibre** (``ndi = 1``, ``nshr = 1`` or ``2``):
    :math:`\sigma_{22} = \sigma_{33} = \sigma_{23} = 0` is a hypothesis on the
    section, not a free face, since the kept shears 12 and 13 load the lateral
    faces. The condition consistent with the beam kinematics frees the stress
    work-conjugate to the section deformation the element allows: the first
    Piola-Kirchhoff stress in the section frame :math:`\boldsymbol{\Lambda}`,
    :math:`\operatorname{sym}(\boldsymbol{\Lambda}^T\mathbf{P})_{\alpha\beta} = 0`.
    It is neither option. ``'spatial_face'`` poses the condition in the frame
    :math:`\mathbf{R}` of :math:`\mathbf{F}`, which is the section frame only
    for a shear-free rod in extension and bending. Transverse shear
    (Timoshenko, Simo-Reissner) and torsion at a fibre off the axis tilt the
    section away from :math:`\mathbf{R}` by a first-order angle. All three
    conditions agree to first order in shear and torsion and differ at second
    order, the order the rigid-section hypothesis already neglects (Poynting
    and Wagner effects). ``'spatial_face'``, the cheaper option, is the
    default for beams;
  - **truss** (``ndi = 1``, ``nshr = 0``): the free directions are principal
    for an isotropic material, and the two options coincide.

  ``'spatial_face'`` is the default and the cheaper option. For a shell under
  transverse shear (``ndi = 2``, ``nshr = 3``), ``'material_face'`` frees the
  material face exactly where ``'spatial_face'`` leaves
  :math:`S_{33} \simeq -\mu\gamma_t^2`. It costs 1.3 to 2.9 times more per
  point, depending on the law and the form, needs ``F0`` and ``F1``, and returns
  the material tangent :math:`\partial\mathbf{S}/\partial\mathbf{E}`
  (``tangent_output='material'``). For plane stress without transverse shear
  on an isotropic material, the two give the same point.

  ``sim.condensation_serves(umat_name, ndi, nshr=None, kinematic_form=None)``
  tells a caller (a finite-element code) whether that condensed point is
  validated for a law, under either ``condense`` option: the condition holds at
  convergence and the tangent is the Schur complement of the law's own, exact
  whenever the 3D tangent is (``tangent_mode=3`` for ``EPHIL``, ``EPKCP``,
  ``EPCHA``). It is ``True`` for ``ndi = 3`` and ``False`` for a law not
  validated in that form; a caller selects its path on it rather than on the
  version number.
* **tangent_mode**: 0 = return the elastic operator as ``Lt``, 1 = continuum tangent,
  2 = algorithmic (consistent) tangent (default of the solver).
* **Wm** is the accumulated ``(Wm, Wm_r, Wm_ir, Wm_d)`` at the beginning of the
  increment and must be returned accumulated.
* **Return**: ``sigma (6,)``, ``Lt (6,6)``, ``statev (nstatev,)``, ``Wm (4,)`` and
  optionally ``L (6,6)`` (defaults to ``Lt``). Lists, ``float32`` arrays and CPU
  torch tensors are converted; shapes are checked and a non-finite value is an
  error.
* **Finite strain**: under NLGEOM the caller feeds the corotational logarithmic
  strain and expects the Kirchhoff stress — the same convention as ``ELISO`` /
  ``EPICP`` (``PYEXT`` belongs to the Kirchhoff-box set). The law never sees
  the Cauchy stress: the solver and ``sim.umat`` form
  :math:`\boldsymbol{\sigma} = \boldsymbol{\tau}/J` at their own boundaries
  (see :ref:`stress-measure-tangent-rate`). On the finite solver route the
  law runs in the frame that follows the material: it receives ``DR`` equal to
  the identity and its tensorial history stays in that frame, with nothing to
  rotate. Through ``sim.umat`` it receives the caller's ``DR`` and must rotate
  its history with it. ``PYEXT`` is refused under corates 4 and 5, which need
  the history's variance.

Step cuts and errors
--------------------

Raise :class:`simcoon.StepCut` to ask the solver for a smaller increment (the
trial is discarded and the increment retried; the effective factor is the
solver's ``div_tnew_dt``). Any other exception aborts the solve and is re-raised
**unchanged** in Python (type and traceback preserved), including
``KeyboardInterrupt``. With ``inforce=0`` the solver aborts (``status = 1``, a
``RuntimeError`` unless ``raise_on_abort=False``) when the increment falls below
``Dn_mini``; with the default ``inforce=1`` it forces the minimal increment.

The batch entry points ``sim.umat`` and ``sim.umat_T`` cannot subdivide an
increment. When any kernel asks for a smaller one -- a Python law raising
:class:`simcoon.StepCut`, or a built-in kernel (the modular engine on a
non-finite or runaway return mapping, the SMR* SMA and LLDM0 damage laws on a
failed local iteration) -- they raise :class:`simcoon.StepCut` carrying the
smallest requested ``ratio``. The input arrays are left untouched: a
finite-element coupler discards the call and retries with a smaller increment.

Batch entry point and explicit registration
===========================================

``sim.umat("PYEXT", ...)`` (the per-Gauss-point batch call used by finite-element
couplers) works with a registered law:

.. code-block:: python

   with sim.registered(law):
       stress, statev, Wm, Lt = sim.umat("PYEXT", etot, Detot, F0, F1, sigma, DR,
                                         props, statev, time, dtime, Wm, ndi=3)

The points are integrated **serially** on the calling thread (the callback
re-enters the interpreter, so it never runs inside the parallel region;
``n_threads`` is ignored). :func:`simcoon.registered` restores the previously
registered law on exit; :func:`simcoon.pyumat.register` /
:func:`simcoon.pyumat.unregister` are the explicit forms.

Performance
===========

Per call, the bridge acquires the GIL (released by the solver around the C++
loop), builds small numpy copies of the inputs and copies the outputs back —
about 5–10 µs, negligible against any non-trivial law. A numpy J2 law costs a few
tens of µs per call, a small LSTM step with its autograd tangent about a
millisecond. For finite-element scale, batch the evaluation on the Python side
(a batched step of its own) rather than calling ``sim.umat``
point by point.

Examples
========

* :ref:`sphx_glr_examples_ml_pyumat_numpy_j2.py` — J2 plasticity in numpy
  compared with ``EPICP``.

API reference
=============

.. autoclass:: simcoon.PythonUMAT
   :members:

.. autoclass:: simcoon.StepCut

.. autofunction:: simcoon.registered

.. autofunction:: simcoon.pyumat.register

.. autofunction:: simcoon.pyumat.unregister

.. autofunction:: simcoon.pyumat.current
