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

API reference
-------------

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
