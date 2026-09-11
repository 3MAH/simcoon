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
  out-of-plane strain), 2 = plane stress (the law must condense, cf. ``el_pred``),
  1 = uniaxial. ``nshr`` is the number of shear components.
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
  ``EPICP`` (``PYEXT`` belongs to the Kirchhoff-box set). Internal tensorial history
  is not rotated by the solver; rotate it with ``DR`` in the law if needed.

Step cuts and errors
--------------------

Raise :class:`simcoon.StepCut` to ask the solver for a smaller increment (the
trial is discarded and the increment retried; the effective factor is the
solver's ``div_tnew_dt``). Any other exception aborts the solve and is re-raised
**unchanged** in Python (type and traceback preserved), including
``KeyboardInterrupt``. With ``inforce=0`` the solver aborts (``status = 1``, a
``RuntimeError`` unless ``raise_on_abort=False``) when the increment falls below
``Dn_mini``; with the default ``inforce=1`` it forces the minimal increment.

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
