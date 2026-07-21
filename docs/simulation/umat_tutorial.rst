=====================================
Hands-on: writing a UMAT with simcoon
=====================================

This tutorial builds a complete constitutive subroutine (UMAT) from scratch:
**J2 plasticity with linear isotropic hardening**, written with the typed
:doc:`Tensor2/Tensor4 API <../continuum_mechanics/functions/doc_tensor>` so
the code reads like the equations. The model is chosen because *everything*
has a closed form — the return mapping needs no local iteration and the
continuum tangent can be derived by hand, so every line of the implementation
can be checked analytically.

The complete companion code lives in ``examples/umat_tutorial/``
(``umat_tutorial_J2.hpp/.cpp``) and is **continuously validated** by the
``Ttutorial_umat`` test against the in-tree EPICP reference, the shared
tangent assembly, and a finite-difference Jacobian.

The model
=========

Additive small-strain decomposition, linear isotropic elasticity, von Mises
yield with a linear hardening law:

.. math::

   \boldsymbol{\sigma} = \mathbf{L} : (\boldsymbol{\varepsilon}
       - \boldsymbol{\varepsilon}^{th} - \boldsymbol{\varepsilon}^p),
   \qquad
   f = \sigma_{eq} - (\sigma_Y + H\,p) \le 0,
   \qquad
   \sigma_{eq} = \sqrt{\tfrac{3}{2}\,\mathbf{s}:\mathbf{s}},

with associated flow
:math:`\dot{\boldsymbol{\varepsilon}}^p = \dot{p}\,\boldsymbol{\Lambda}`,
:math:`\boldsymbol{\Lambda} = \tfrac{3}{2}\,\mathbf{s}/\sigma_{eq}`.

Props (5): ``E, nu, alpha, sigma_Y, H``. State (8): ``T_init, p, EP(6)`` —
deliberately the EPICP layout, so the two are directly comparable (EPICP with
``m = 1`` and ``k = H`` is the same model).

Step 1 — anatomy of a UMAT
==========================

Every simcoon UMAT has the same signature (see any kernel or the tutorial
header): strain state in (``Etot``, ``DEtot``, Voigt order
``[11, 22, 33, 12, 13, 23]``, MPa), stress and tangent out (``sigma``,
``Lt``), material properties (``props``), persistent state (``statev``),
the rotation increment ``DR`` for objectivity, cumulative work outputs
(``Wm, Wm_r, Wm_ir, Wm_d``) and the ``tangent_mode`` selector.

Three contract points that every UMAT must honor:

1. **start**: on the first increment, initialize the state and RESET the
   cumulative work accumulators.
2. **objectivity**: stored tensorial state co-rotates with ``DR``. With the
   typed API this is one line — ``EP.rotate(...)`` picks the strain rotation
   kernel (factor-2 shear convention) from the tensor's own type.
3. **cumulative work**: ``Wm += increment`` — the solver reports path
   integrals, never overwrite them.

Step 2 — trial state, typed
===========================

.. literalinclude:: ../../examples/umat_tutorial/umat_tutorial_J2.cpp
   :language: cpp
   :start-at: // ------------------------------------------------------------------ 4.
   :end-before: // ------------------------------------------------------------------ 5.

Because ``L_el`` is a *stiffness-typed* ``tensor4`` and the operand is a
*strain-typed* ``tensor2``, the contraction produces a *stress-typed* result
— the shear-convention bookkeeping (engineering factor 2) is carried by the
types, not by the programmer.

Step 3 — radial return, closed form
===================================

The trial direction is preserved by the return (the correction is radial in
deviatoric space), and with linear hardening the consistency condition is
linear in :math:`\Delta p`:

.. math::

   f(\Delta p) = f^{trial} - (3\mu + H)\,\Delta p = 0
   \quad\Longrightarrow\quad
   \Delta p = \frac{f^{trial}}{3\mu + H}.

.. literalinclude:: ../../examples/umat_tutorial/umat_tutorial_J2.cpp
   :language: cpp
   :start-at: if (f_trial > 0.) {
   :end-at: sigma = sig.to_arma_voigt();

Step 4 — the continuum tangent, derived and built
=================================================

Differentiating the converged stress w.r.t. the total strain with the flow
direction frozen gives the classical rank-one update

.. math::

   \mathbf{L}_t = \mathbf{L}
     - \frac{(\mathbf{L}:\boldsymbol{\Lambda}) \otimes (\mathbf{L}:\boldsymbol{\Lambda})}
            {\boldsymbol{\Lambda}:\mathbf{L}:\boldsymbol{\Lambda} + H}.

For isotropic :math:`\mathbf{L}` and the deviatoric von Mises direction the
pieces are analytic — :math:`\mathbf{L}:\boldsymbol{\Lambda} =
2\mu\boldsymbol{\Lambda}` and :math:`\boldsymbol{\Lambda}:\mathbf{L}:
\boldsymbol{\Lambda} = 3\mu` — hence the closed form

.. math::

   \mathbf{L}_t = \mathbf{L} - \frac{4\mu^2}{3\mu + H}\,
       \boldsymbol{\Lambda}\otimes\boldsymbol{\Lambda}.

.. literalinclude:: ../../examples/umat_tutorial/umat_tutorial_J2.cpp
   :language: cpp
   :start-at: if (Dp <= simcoon::iota || tangent_mode == tangent_none) {
   :end-before: // ------------------------------------------------------------------ 8.

.. warning::
   **The classic Voigt trap.** The closed form above holds for the
   *tensorial* :math:`\boldsymbol{\Lambda}`. In Voigt components the outer
   product must be taken with the **stress-Voigt** normal (shear *not*
   doubled): :math:`\boldsymbol{\kappa} = \mathbf{L}:\boldsymbol{\Lambda}`
   lands in stress Voigt automatically, but ``eta_stress()`` returns the
   *engineering* (doubled-shear) normal — using it directly in
   :math:`\boldsymbol{\Lambda}\boldsymbol{\Lambda}^T` overweights the shear
   block of the correction by up to a factor 4. The typed ``auto_dyadic``
   route sidesteps the trap (the type carries the convention); note also
   that ``sym_dyadic`` is the :math:`(ik)(jl)`-symmetrized product — a
   *different* operator from the plain dyadic needed here.

For the **algorithmic** (Simo–Hughes) operator the tutorial does what the
production kernels do — hand the physical inputs to the shared dispatch
``compute_tangent_operator`` with a lazy flow-Hessian provider; for J2 the
resulting operator is the exact Jacobian of the discrete radial-return map.

Step 5 — validation
===================

``test/Libraries/Umat/Ttutorial_umat.cpp`` asserts, on every test run:

1. **Physics**: the tutorial matches the EPICP reference to the CCP
   tolerance (< 1e-8 relative) over a cyclic tension/shear path.
2. **Tangent**: the hand-built continuum operator equals both the
   independent closed form (with the stress-Voigt normal!) and
   ``assemble_continuum_tangent`` to machine precision.
3. **Consistency**: the algorithmic operator equals the central
   finite-difference Jacobian of the discrete stress update (< 1e-5).

Where to go next
================

- Swap the yield criterion: replace ``flow_normal`` by a Hill or DFA
  gradient (``criteria.hpp``) — the *generic* :math:`\kappa`-form of the
  tangent in the code keeps working; only the closed-form shortcut is
  J2-specific.
- Add kinematic hardening: shift the criterion to
  :math:`\boldsymbol{\sigma} - \mathbf{X}` and consult the kept EPCHA kernel
  for the classic Armstrong–Frederick treatment.
- Or skip hand-written kernels entirely: the same model is one line of
  :class:`simcoon.modular.ModularMaterial` — see :doc:`umat_catalog`.
