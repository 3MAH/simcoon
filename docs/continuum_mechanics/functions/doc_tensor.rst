Typed Tensors (Tensor2 / Tensor4)
=================================

``simcoon.Tensor2`` and ``simcoon.Tensor4`` are unified typed wrappers for 2nd-
and 4th-order tensors. A single object transparently represents either one tensor
or a batch (scipy ``Rotation`` style): shape ``(6,)`` / ``(6,6)`` is a single
tensor, ``(N,6)`` / ``(N,6,6)`` a batch. ``Tensor4`` stores its data internally
in the **Kelvin-Mandel** convention, so the double contraction, inverse and
composition of 4th-order tensors become ordinary 6x6 linear algebra; the
engineering Voigt form is recovered on demand via ``.mat`` / ``.voigt``.

Objects are built through typed factories, e.g. ``Tensor2.stress(v)``,
``Tensor2.strain(v)``, ``Tensor4.stiffness(m)``, ``Tensor4.compliance(m)``,
``Tensor4.from_voigt(v, type_str)``.

Batch operations
----------------

Batch ``contract``, ``rotate``, ``push_forward``, ``pull_back`` and ``inverse``
dispatch to vectorised kernels over the ``(N, ...)`` data.

.. tip:: **Shared-tangent contraction is ~20x faster — pass a single Tensor4.**

   When the same 4th-order tensor is contracted against many 2nd-order tensors
   (a shared elastic tangent applied at every integration point), pass a
   **single** ``Tensor4`` and contract it with a **batch** ``Tensor2``::

       L   = Tensor4.stiffness(L6x6)     # single (6,6)
       eps = Tensor2.strain(eps_Nx6)     # batch  (N,6)
       sig = L.contract(eps)             # one BLAS GEMM  ->  ~20x

   This collapses to a single ``L @ X`` matrix-matrix product (BLAS ``dgemm``):
   measured ~3.7 ns/point vs ~74 ns/point for the per-point path (N = 1e5).
   The C++ ``batch_contract`` implements the same fast path when the tensor4
   cube has a single slice (``N4 == 1``).

   **Anti-pattern:** materialising the shared tangent as a *tiled* ``(N,6,6)``
   batch of identical slices — that forces ``N`` independent 6x6 matrix-vector
   products and the GEMM speed-up is lost. Tile the tangent **only** when it
   genuinely differs per point (a distinct consistent tangent, e.g. plasticity);
   there a GEMM collapse is not possible and per-point evaluation is correct.

.. note:: A dedicated symmetric ``symtensor2`` (Kelvin-Mandel 6-vector storage)
   was benchmarked as a way to speed the general contraction path and
   **rejected**: it yields only ~1.1x on a single contraction and ~1.07x on
   realistic distinct-tangent batch work, and Armadillo fixed-size containers
   negate the 6-vs-9-double storage saving. The real ~20x lever is the
   shared-tangent GEMM path above, which needs no new type; the C++ ``tensor2``
   stays the general 3x3 type (the only one able to hold a non-symmetric tensor
   such as ``F``, ``L``, ``R``). The Python ``Tensor2`` stores the 6 Voigt
   components for the symmetric types (``"stress"``, ``"strain"``,
   ``"symmetric"``) and the 9 components of any 3x3 for the type ``"none"`` (no
   Voigt convention); a non-symmetric matrix given to a
   symmetric type is refused rather than silently symmetrised.

Basis: the reference system of the components
---------------------------------------------

Objectivity implies that the physical laws are independent of the choice of reference frame. 
The components of a tensor only mean something together with the basis they are
written in. By default that basis is the fixed orthonormal laboratory basis
:math:`\mathbf{e}_i`: ``t.basis`` is ``None`` and nothing is stored. A
``simcoon.Basis`` makes another choice explicit. It holds the basis vectors
:math:`\mathbf{g}_i` through the matrix :math:`\mathbf{A}` whose columns are
their laboratory frame components, and comes in two kinds:

* **orthonormal**, built from a rotation, :math:`\mathbf{g}_i = \mathbf{R}\,\mathbf{e}_i`
  (a material frame, a corotational frame). Its metric is the identity, so every
  formula is the one of the lab.
* **natural**, any three independent vectors, typically the convected basis
  :math:`\mathbf{g}_i = \mathbf{F}\,\mathbf{G}_i` (``Basis.from_F``). Its metric
  :math:`g_{ij} = \mathbf{g}_i \cdot \mathbf{g}_j`, i.e.
  :math:`\mathbf{g} = \mathbf{A}^T \mathbf{A}`, enters the invariants.
  'Natural' definition come from the convected basis of a material point.

The type tag gives the variance of the components: a stress is contravariant,
:math:`\boldsymbol{\sigma} = \sigma^{ij}\,\mathbf{g}_i \otimes \mathbf{g}_j`, a
strain covariant, :math:`\boldsymbol{\varepsilon} = \varepsilon_{ij}\,\mathbf{g}^i
\otimes \mathbf{g}^j`; a stiffness is contravariant on its four indices, a
compliance covariant, the concentration tensors mixed. The lab components follow:

.. math::

   \boldsymbol{\sigma}_{lab} = \mathbf{A}\,\hat{\boldsymbol{\sigma}}\,\mathbf{A}^T,
   \qquad
   \boldsymbol{\varepsilon}_{lab} = \mathbf{A}^{-T}\,\hat{\boldsymbol{\varepsilon}}\,\mathbf{A}^{-1}.

.. list-table:: Default variance of each type (``t.variance``)
   :header-rows: 1
   :widths: 30 22 48

   * - Type
     - Components
     - Transport by :math:`\mathbf{F}` (``push_forward``)
   * - ``Tensor2`` ``"stress"``
     - :math:`\sigma^{ij}` (contravariant)
     - :math:`\mathbf{F}\,\boldsymbol{\sigma}\,\mathbf{F}^T`, weight :math:`1/J` with ``metric=True``
   * - ``Tensor2`` ``"strain"``
     - :math:`\varepsilon_{ij}` (covariant)
     - :math:`\mathbf{F}^{-T}\,\boldsymbol{\varepsilon}\,\mathbf{F}^{-1}`
   * - ``Tensor2`` ``"symmetric"``, ``"none"`` (any 3x3, no Voigt convention)
     - lab or orthonormal-frame components only
     - none: no variance, hence no natural basis and no transport
   * - ``Tensor4`` ``"stiffness"``, ``"generic"``
     - :math:`C^{ijkl}`
     - :math:`\mathbf{F}` on the four indices, weight :math:`1/J`
   * - ``Tensor4`` ``"compliance"``
     - :math:`M_{ijkl}`
     - :math:`\mathbf{F}^{-T}` on the four indices, weight :math:`J`
   * - ``Tensor4`` concentrations
     - mixed, :math:`A_{ij}{}^{kl}` / :math:`B^{ij}{}_{kl}`
     - pair-wise: :math:`\mathbf{F}^{-T}` on the covariant pair, :math:`\mathbf{F}` on the contravariant one

The variance is a tag of the components, ``t.variance``, carried next to the type
and the basis: ``"contravariant"`` or ``"covariant"`` for a ``Tensor2``, a pair
``(output, input)`` for a ``Tensor4``. The table gives its default; it is read,
not the type, by the transports, by ``to_basis`` and by the metric invariants.
Two things the default is **not**:

* It is not the only possible representation. Variance is a choice made for the
  components, not a property of the tensor: with the metric any tensor can be
  given contravariant or covariant components (``det`` and ``eigvals`` use the
  mixed ones :math:`T^i{}_j = T^{ik} g_{kj}` internally). The default is the
  one under which the conjugate pairs :math:`\mathbf{S} : \mathbf{E}` and
  :math:`\boldsymbol{\tau} : \mathbf{e}` contract without a metric and the
  convected components stay constant under transport. ``t.to_variance(v)``
  gives the same tensor with components of the other variance: a retag in the
  lab or an orthonormal basis (the numbers do not change), a contraction with
  the metric in a natural basis, :math:`\sigma_{ij} = g_{ik}\,\sigma^{kl}\,g_{lj}`
  (or its inverse), pair by pair for a ``Tensor4``. The Voigt shear factors stay
  those of the type. Mixed components of a single pair are not representable in
  6 components (they are not symmetric), hence a tag per pair. The variance is
  part of what a transport means: a covariantly tagged stiffness pushes forward
  with :math:`\mathbf{F}^{-T}` on its four indices, which is another tensor than
  the push-forward of the contravariant one. A ``"symmetric"`` ``Tensor2`` has no
  default variance: ``to_variance`` declares it before any transport or natural
  basis; type ``"none"`` has none.
* It is not a statement about every strain measure. The covariant transport is
  that of the Green-Lagrange / Almansi pair. A logarithmic strain is a tensor
  function of a stretch, not the transport of anything: ``push_forward`` on a
  strain-typed :math:`\ln \mathbf{V}` is accepted by the type tag and has no
  meaning.

Transport in index and Mandel notation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**The index products.** The transports of the table act with one kernel
:math:`\mathbf{K}` on every index: :math:`\mathbf{K} = \mathbf{F}` on a
contravariant index, :math:`\mathbf{K} = \mathbf{F}^{-T}` on a covariant one. For a
symmetric second-order tensor and a fourth-order tensor with both minor
symmetries,

.. math::

   t'_{ij} = K_{ia} K_{jb}\, t_{ab},
   \qquad
   C'_{ijkl} = K_{ia} K_{jb} K_{kc} K_{ld}\, C_{abcd}.

:math:`\mathbf{F}` is any deformation gradient: neither symmetric nor orthogonal.
The transported tensors keep their symmetries (:math:`t'` is symmetric because
:math:`\mathbf{K}\,\mathbf{t}\,\mathbf{K}^T` is, whatever :math:`\mathbf{K}`; the minor
symmetries of :math:`C'` follow from those of :math:`C` in the same way).

**Why a 6 x 6 matrix does the job.** The first transport is the linear map
:math:`\mathbf{X} \mapsto \mathbf{K}\mathbf{X}\mathbf{K}^T` restricted to the symmetric
second-order tensors, a six-dimensional space. In its Kelvin-Mandel orthonormal
basis (:math:`e_{(ii)} = \mathbf{e}_i \otimes \mathbf{e}_i` and
:math:`e_{(ij)} = (\mathbf{e}_i \otimes \mathbf{e}_j + \mathbf{e}_j \otimes \mathbf{e}_i)/\sqrt2`
for :math:`i<j`) that map has the matrix

.. math::

   (\mathbb{Q}_{\mathbf{K}})_{IJ} = \frac{c_I}{c_J}\left( K_{ik} K_{jl} + (1-\delta_{kl})\, K_{il} K_{jk} \right),
   \qquad I = (ij),\ J = (kl),\ c = (1,1,1,\sqrt2,\sqrt2,\sqrt2).

A fourth-order tensor with both minor symmetries is a bilinear form on that same
space, and transporting both of its arguments by the same map is a congruence.
The two index products above therefore read, exactly,

.. math::

   \{t'\} = \mathbb{Q}_{\mathbf{K}}\, \{t\},
   \qquad
   \mathbb{C}' = \mathbb{Q}_{\mathbf{K}}\, \mathbb{C}\, \mathbb{Q}_{\mathbf{K}}^T .

No approximation enters: the Mandel components are the tensor components in an
orthonormal basis, so the matrix product is the index product. simcoon never
forms the :math:`3^4` components (``mandel_congruence`` in ``tensor.cpp``; checked
against the explicit four-index contraction to :math:`5 \times 10^{-16}` on random
non-symmetric :math:`\mathbf{F}`).

**What** :math:`\mathbb{Q}_{\mathbf{K}}` **is not.** It is not orthogonal for a general
:math:`\mathbf{F}`: :math:`\mathbb{Q}_{\mathbf{K}}^{-1} = \mathbb{Q}_{\mathbf{K}^{-1}} \neq \mathbb{Q}_{\mathbf{K}}^T`.
That is why ``pull_back`` applies the congruence with the inverse kernel
(:math:`\mathbf{F}^{-1}` for a stiffness, :math:`\mathbf{F}^{T}` for a compliance) and never
the transpose of :math:`\mathbb{Q}_{\mathbf{F}}`. :math:`\mathbb{Q}_{\mathbf{K}}` is orthogonal
exactly when :math:`\mathbf{K}` is a rotation, in which case the congruence is a
rotation of components: ``rotate`` and ``push_forward`` share the formula, not the
property.

**Where the congruence stops.** It needs one kernel on every index and both minor
symmetries. A tensor with mixed variance (the concentration tensors
:math:`A_{ij}{}^{kl}`, :math:`B^{ij}{}_{kl}`) transports as
:math:`\mathbb{Q}_{\mathbf{F}^{-T}}\, \mathbb{A}\, \mathbb{Q}_{\mathbf{F}}^T`, two different
kernels: not a congruence, and ``push_forward`` refuses those types (as before). A
non-symmetric second-order tensor has no Mandel vector; ``Tensor2`` transports it
as a 3 x 3 product directly.

**The dyadic products of** ``contimech`` follow the same rule:
:math:`C_{ijkl} = a_{ij} b_{kl}` is the outer product of the two Voigt vectors, and
:math:`\partial \mathbf{S}^{-1} / \partial \mathbf{S} =
-\tfrac12 (S^{-1}_{ik} S^{-1}_{jl} + S^{-1}_{il} S^{-1}_{jk})` is the symmetric dyadic
operator of :math:`\mathbf{S}^{-1}`, negated.

Two-point tensors have no basis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The deformation gradient :math:`\mathbf{F} = F^i{}_J\,\mathbf{g}_i \otimes
\mathbf{G}^J` has one index in each configuration, so there is no single basis
to attach it to; in the convected pair :math:`(\mathbf{G}_i, \mathbf{g}_i =
\mathbf{F}\,\mathbf{G}_i)` its mixed components are :math:`\delta^i_J`: it *is*
the relation between the two bases. This is why ``F`` is always a plain array of
lab components in the API, consumed by ``Basis.from_F`` and by the transports,
never wrapped. The same holds for the rotation :math:`\mathbf{R}` of the polar
decomposition, the rotation increment :math:`\Delta\mathbf{R}` and the first
Piola-Kirchhoff stress. They can be held in a ``Tensor2`` of type ``"none"`` (no
Voigt convention: 9 components, lab-lab). Such a tensor has no variance, hence no
transport and no natural basis (re-expressing its components would need the
variance of each leg); it does accept an orthonormal ``Basis``, both legs being
read in that frame, since :math:`\mathbf{Q}^T \mathbf{X}\,\mathbf{Q}` is the same
for every variance. It offers the trace, determinant, norm and eigenvalues, but
no Voigt vector.

Three operations, three meanings:

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - Call
     - Meaning
   * - ``t.to_basis(b)``
     - The same tensor, its components re-expressed in ``b`` (``None`` = lab).
   * - ``t.with_basis(b)``
     - The same components in another basis: a transport, hence another tensor.
       No arithmetic, no copy.
   * - ``t.rotate(R, active=False)``
     - The same tensor in the frame turned by ``R``; the result carries that basis.

Tensors written in different bases cannot be combined: ``+``, ``-``, ``@``,
``%`` and the dyadic products raise ``ValueError("Mixed basis")``. A ``name``
given to a basis is a label for ``repr`` and error messages; two bases are the
same when their vectors are.

.. code-block:: python

   import numpy as np
   import simcoon as sim

   material = sim.Basis(rotation=sim.Rotation.from_euler("z", 30, degrees=True),
                        name="material")
   sigma = sim.Tensor2.stress(np.array([100., 0., 0., 0., 0., 0.]))
   sigma_m = sigma.to_basis(material)     # material-frame components
   sigma + sigma_m                        # ValueError: Mixed basis
   sigma_m.to_basis(None)                 # lab components again

One basis can serve a whole batch (it is stored once and shared by reference by
every tensor derived from it), or a batch of N bases can give one per tensor.

A basis also knows its reciprocal vectors :math:`\mathbf{g}^i` (``reciprocal``,
the columns of :math:`\mathbf{A}^{-T}`, on which covariant components live), its
right stretch :math:`\mathbf{U} = \sqrt{\mathbf{g}}` (``stretch``) and the
orthonormal basis closest to it, the rotation of the polar decomposition
:math:`\mathbf{A} = \mathbf{R}\,\mathbf{U}` (``polar``). For a convected basis
these are the right stretch and the rotation of :math:`\mathbf{F}`; ``polar`` is
the rotated orthonormal frame a convected tensor is often read in.

Push-forward and pull-back carry the basis along
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A push-forward is a transport, not a change of basis: the result is another
tensor, living in the current configuration. But in the convected basis its
components are those the original tensor had in the reference basis:

.. math::

   \tau^{ij} = S^{IJ}, \qquad e_{ij} = E_{IJ}, \qquad
   (\mathcal{L}_v \boldsymbol{\tau})^{ij} = \mathbb{C}^{IJKL}\,\dot{E}_{KL},

with :math:`\boldsymbol{\tau}` the Kirchhoff stress, :math:`\mathbf{S}` the
second Piola-Kirchhoff stress, :math:`\mathbf{e}` the Almansi and
:math:`\mathbf{E}` the Green-Lagrange strain. Pushing forward carries the basis
:math:`\mathbf{G}_i` to :math:`\mathbf{g}_i = \mathbf{F}\,\mathbf{G}_i` and keeps
the components; the familiar :math:`\mathbf{F}\,\mathbf{S}\,\mathbf{F}^T` is that
statement followed by "write the result in the lab". Accordingly:

* on a lab tensor, ``push_forward(F)`` returns lab components, as before;
* on a tensor that has its own basis, it keeps the components and returns the
  convected basis :math:`\mathbf{F}\,\mathbf{A}` (``pull_back``:
  :math:`\mathbf{F}^{-1}\mathbf{A}`). Only the weight :math:`1/J` of
  ``metric=True`` (Kirchhoff to Cauchy, a density) touches the numbers.

.. code-block:: python

   F = np.array([[1., 0.5, 0.], [0., 1., 0.], [0., 0., 1.]])    # simple shear
   S = sim.Tensor2.stress(np.array([0., 80., 0., 0., 0., 0.]))  # PK2 components
   tau = S.with_basis(sim.Basis.from_F(F))    # Kirchhoff stress, same numbers
   tau.to_basis(None)                         # == S.push_forward(F, metric=False)
   tau.trace()                                # 100.0 = g_22 * 80, not 80

The same holds for an active rotation: a lab tensor is rotated in the lab, a
tensor with its own basis keeps its components while its basis vectors turn
(:math:`\mathbf{A} \to \mathbf{Q}\,\mathbf{A}`), which leaves the metric and
every invariant unchanged.

Invariants with the metric
~~~~~~~~~~~~~~~~~~~~~~~~~~

Contractions between dual variances need no metric and act on the components in
any basis: :math:`\boldsymbol{\sigma} : \boldsymbol{\varepsilon} = \sigma^{ij}
\varepsilon_{ij}`, :math:`\mathbb{L} : \boldsymbol{\varepsilon}`, the inverse of
a stiffness. Contracting two indices of the same variance does. For
contravariant components (covariant ones exchange :math:`\mathbf{g}` and
:math:`\mathbf{g}^{-1}`):

.. math::

   \mathrm{tr}\,\mathbf{T} = g_{ij}\,T^{ij}, \qquad
   \mathrm{dev}\,\hat{\mathbf{T}} = \hat{\mathbf{T}} - \tfrac13\,(\mathrm{tr}\,\mathbf{T})\,\mathbf{g}^{-1}, \qquad
   \mathbf{a} : \mathbf{b} = \mathrm{tr}(\hat{\mathbf{a}}\,\mathbf{g}\,\hat{\mathbf{b}}\,\mathbf{g}),

.. math::

   \sigma_{eq} = \sqrt{\tfrac32\,\mathbf{s} : \mathbf{s}}, \qquad
   \det\mathbf{T} = \det\hat{\mathbf{T}}\,\det\mathbf{g}, \qquad
   \det(\hat{\mathbf{T}} - \lambda\,\mathbf{g}^{-1}) = 0 .

``trace``, ``dev``, ``mises``, ``norm``, ``det``, ``eigvals`` and ``%`` use these
forms in a natural basis and the usual ones otherwise. The identity tensor is the
metric itself: ``Tensor2.identity("stress", basis=b)`` has components
:math:`g^{ij}`, ``Tensor2.identity("strain", basis=b)`` has :math:`g_{ij}`. For
fourth-order tensors, ``Tensor4.identity``, ``volumetric`` and ``deviatoric``
take the same ``basis`` argument; the stiffness-type identity becomes
:math:`\tfrac12 (g^{ik} g^{jl} + g^{il} g^{jk})` and the volumetric projector
:math:`\tfrac13\,g^{ij} g^{kl}`.

.. note::

   The Kelvin-Mandel form is an isometry only in an orthonormal basis. In a
   natural basis ``.mandel`` and ``.voigt`` still return the stored components,
   but ``eye(6)`` is no longer the stiffness-type identity and the plain
   Frobenius norm of the 6x6 is not the norm of the tensor.

Cost: with ``basis=None`` nothing changes. A basis costs one reference per
tensor object, not per point; ``with_basis`` and the transports of a tensor with
its own basis are free; ``to_basis`` costs one transport (the same kernel as
``push_forward``); the metric of a natural basis is computed on first use and
cached, and its invariants cost about one 3x3 product per point more than in the
lab.

API reference
-------------

.. autoclass:: simcoon.Basis
   :members:

.. autoclass:: simcoon.Tensor2
   :members:
   :inherited-members:

.. autoclass:: simcoon.Tensor4
   :members:
   :inherited-members:

.. autofunction:: simcoon.dyadic
.. autofunction:: simcoon.auto_dyadic
.. autofunction:: simcoon.sym_dyadic
.. autofunction:: simcoon.auto_sym_dyadic
.. autofunction:: simcoon.double_contract
