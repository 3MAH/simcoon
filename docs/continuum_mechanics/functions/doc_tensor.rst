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
   shared-tangent GEMM path above, which needs no new type; ``Tensor2`` stays
   the general 3x3 type (it remains the only one able to hold a non-symmetric
   tensor such as ``F``, ``L``, ``R``).

Basis: the reference system of the components
---------------------------------------------

The components of a tensor only mean something together with the basis they are
written in. By default that basis is the fixed orthonormal lab basis
:math:`\mathbf{e}_i`: ``t.basis`` is ``None`` and nothing is stored. A
``simcoon.Basis`` makes another choice explicit. It holds the basis vectors
:math:`\mathbf{g}_i` through the matrix :math:`\mathbf{A}` whose columns are
their lab components, and comes in two kinds:

* **orthonormal**, built from a rotation, :math:`\mathbf{g}_i = \mathbf{R}\,\mathbf{e}_i`
  (a material frame, a corotational frame). Its metric is the identity, so every
  formula is the one of the lab.
* **natural**, any three independent vectors, typically the convected basis
  :math:`\mathbf{g}_i = \mathbf{F}\,\mathbf{G}_i` (``Basis.from_F``). Its metric
  :math:`g_{ij} = \mathbf{g}_i \cdot \mathbf{g}_j`, i.e.
  :math:`\mathbf{g} = \mathbf{A}^T \mathbf{A}`, enters the invariants.

The type tag gives the variance of the components: a stress is contravariant,
:math:`\boldsymbol{\sigma} = \sigma^{ij}\,\mathbf{g}_i \otimes \mathbf{g}_j`, a
strain covariant, :math:`\boldsymbol{\varepsilon} = \varepsilon_{ij}\,\mathbf{g}^i
\otimes \mathbf{g}^j`; a stiffness is contravariant on its four indices, a
compliance covariant, the concentration tensors mixed. The lab components follow:

.. math::

   \boldsymbol{\sigma}_{lab} = \mathbf{A}\,\hat{\boldsymbol{\sigma}}\,\mathbf{A}^T,
   \qquad
   \boldsymbol{\varepsilon}_{lab} = \mathbf{A}^{-T}\,\hat{\boldsymbol{\varepsilon}}\,\mathbf{A}^{-1}.

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
