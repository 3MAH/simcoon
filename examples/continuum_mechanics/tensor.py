"""
Tensor2 and Tensor4 Examples
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``Tensor2`` and ``Tensor4`` classes wrap simcoon's C++ tensor objects.
They carry a **type tag** (stress/strain for ``Tensor2``, stiffness/compliance
for ``Tensor4``), which automatically selects the correct Voigt convention,
rotation rule, and push-forward/pull-back formula.

A single class handles **both** individual tensors and **batches** of N tensors
(like ``scipy.spatial.transform.Rotation``).  Shape determines which mode:

* ``(6,)`` or ``(3,3)`` → single ``Tensor2``
* ``(N,6)`` or ``(N,3,3)`` → batch of N ``Tensor2``
* ``(6,6)`` → single ``Tensor4``
* ``(N,6,6)`` → batch of N ``Tensor4``

Batch operations loop in C++ over the same functions used by single tensors,
so results are always bit-identical.
"""

import numpy as np
import simcoon as sim

# %%
# 1. Building typed tensors
# -------------------------
# The factory methods ``stress()``, ``strain()``, ``stiffness()``, etc.
# embed the Voigt convention into the object.

sigma = sim.Tensor2.stress(np.array([100, 200, 300, 30, 20, 40], dtype=float))
print("sigma.vtype:", sigma.vtype)
print("sigma.voigt:", sigma.voigt)
print("sigma.mat:\n", sigma.mat)

eps = sim.Tensor2.strain(np.array([0.01, 0.02, 0.03, 0.01, 0.006, 0.008]))
print("\neps.vtype:", eps.vtype)
print("eps.voigt:", eps.voigt)

# %%
# 2. Stiffness and compliance
# ----------------------------
# ``Tensor4`` wraps a 6x6 Voigt matrix with its type.  Inversion automatically
# flips the type (stiffness <-> compliance).

L = sim.Tensor4.stiffness(sim.L_iso([70000, 0.3], "Enu"))
print("L.type:", L.type)

M = L.inverse()
print("M.type:", M.type)  # compliance

L_back = M.inverse()
print("Roundtrip match:", np.allclose(L.mat, L_back.mat))

# %%
# 3. Contraction: sigma = L : epsilon
# ------------------------------------
# The ``@`` operator contracts a ``Tensor4`` with a ``Tensor2``.
# The output type is inferred: stiffness @ strain -> stress.

sigma = L @ eps
print("sigma.vtype:", sigma.vtype)  # stress
print("sigma.voigt:", np.array_str(sigma.voigt, precision=4))

# %%
# 4. Scalar invariants
# ---------------------
# ``trace()`` and ``mises()`` are type-aware.

print("Trace:", sigma.trace())
print("Von Mises:", sigma.mises())

# %%
# 5. Rotation via simcoon.Rotation
# ---------------------------------
# Rotation dispatches on the type tag: stress uses QS, strain uses QE,
# stiffness uses QS on both sides, compliance uses QE, etc.
# The rotation roundtrip is exact.

rot = sim.Rotation.from_euler("ZXZ", [0.3, 0.5, 0.7])

sigma_rot = sigma.rotate(rot)
sigma_back = sigma_rot.rotate(rot.inv())
print("Rotation roundtrip error:", np.max(np.abs(sigma.mat - sigma_back.mat)))

# Isotropic stiffness is invariant under any rotation
L_rot = L.rotate(rot)
print("Isotropic L unchanged:", np.allclose(L.mat, L_rot.mat, atol=1e-8))

# %%
# 6. Push-forward and pull-back
# -------------------------------
# Stress (contravariant): ``F * sigma * F^T``.
# Strain (covariant): ``F^{-T} * eps * F^{-1}``.
# The type tag selects the correct formula automatically.

F = np.array([[1.1, 0.1, 0.05],
              [0.02, 0.95, 0.03],
              [0.01, 0.04, 1.05]])

eps_push = eps.push_forward(F)
eps_back = eps_push.pull_back(F)
print("Push/pull roundtrip error:", np.max(np.abs(eps.mat - eps_back.mat)))

# %%
# 7. Material-frame equivalence
# --------------------------------
# ``L : eps`` in the material frame equals ``(rotate L) : (rotate eps)``
# in the global frame.  This is the fundamental rotation identity.

L_cubic = sim.Tensor4.stiffness(sim.L_cubic([200000, 0.3, 80000], "EnuG"))
eps_test = sim.Tensor2.strain(np.array([0.01, -0.003, -0.003, 0.005, 0.002, 0.001]))

sigma_mat = L_cubic @ eps_test
sigma_global = L_cubic.rotate(rot) @ eps_test.rotate(rot)
sigma_check = sigma_mat.rotate(rot)

print("Frame equivalence error:", np.max(np.abs(sigma_global.mat - sigma_check.mat)))

# %%
# 8. Arithmetic and dyadic product
# ----------------------------------
# Tensors support ``+``, ``-``, scalar ``*``, and the dyadic product.

t1 = sim.Tensor2.stress(np.array([100, 0, 0, 0, 0, 0], dtype=float))
t2 = sim.Tensor2.stress(np.array([0, 100, 0, 0, 0, 0], dtype=float))

t3 = t1 + t2
print("t1 + t2:", t3.voigt)
print("2 * t1:", (2 * t1).voigt)

C = sim.dyadic(t1, t2)
print("dyadic type:", C.type)
print("C[0,1] =", C.mat[0, 1])

# %%
# 9. Building isotropic stiffness from projectors
# --------------------------------------------------
# Using the volumetric and deviatoric identity tensors.

E, nu = 70000.0, 0.3
K = E / (3.0 * (1.0 - 2.0 * nu))
mu = E / (2.0 * (1.0 + nu))

L_proj = 3.0 * K * sim.Tensor4.volumetric() + 2.0 * mu * sim.Tensor4.deviatoric()
print("L from projectors matches L_iso:",
      np.allclose(L_proj.mat, sim.L_iso([E, nu], "Enu"), atol=1e-8))

# %%
# =====================================================================
# Batch operations — same class, automatic from shape
# =====================================================================

# %%
# 10. Creating a batch from an (N, 6) array
# -------------------------------------------
# Passing a 2D array to ``stress()`` creates a batch.  ``.single`` tells you
# which mode you are in.

N = 100
rng = np.random.default_rng(42)

eps_batch = sim.Tensor2.strain(rng.standard_normal((N, 6)) * 0.01)
print("eps_batch.single:", eps_batch.single)
print("len(eps_batch):", len(eps_batch))
print("eps_batch.voigt.shape:", eps_batch.voigt.shape)  # (100, 6)

# %%
# 11. Indexing and iteration
# ---------------------------
# Batch tensors support ``len``, ``[]``, ``for t in batch``, ``reversed``.

first = eps_batch[0]
print("first.single:", first.single)
print("first.voigt:", np.array_str(first.voigt, precision=6))

sub = eps_batch[10:15]
print("sub-batch length:", len(sub))

for i, t in enumerate(eps_batch[:3]):
    print(f"  [{i}] mises = {t.mises():.6f}")

# %%
# 12. Batch contraction: sigma = L : eps for 100 Gauss points
# --------------------------------------------------------------
# A single ``Tensor4`` is broadcast across the batch of ``Tensor2``.

L = sim.Tensor4.stiffness(sim.L_iso([70000, 0.3], "Enu"))
sigma_batch = L @ eps_batch
print("sigma_batch.vtype:", sigma_batch.vtype)  # stress
print("sigma_batch shape:", sigma_batch.voigt.shape)  # (100, 6)

# Or with a batch of stiffness tensors
L_batch = sim.Tensor4.from_tensor(L, N)
sigma_batch2 = L_batch @ eps_batch
print("Batch L @ batch eps match:", np.allclose(sigma_batch.voigt, sigma_batch2.voigt))

# %%
# 13. Batch von Mises and trace
# --------------------------------
# Scalar invariants return ``(N,)`` arrays for batches.

vm = sigma_batch.mises()
tr = sigma_batch.trace()
print("Von Mises: min={:.2f}, max={:.2f}".format(vm.min(), vm.max()))
print("Trace: min={:.2f}, max={:.2f}".format(tr.min(), tr.max()))

# %%
# 14. Batch rotation with scipy
# --------------------------------
# Rotate N tensors by N different rotations (one per Gauss point).
# The C++ loop calls the exact same rotation code as single tensors.

from scipy.spatial.transform import Rotation as ScipyRotation

rotations = ScipyRotation.random(N, random_state=rng)

sigma_rotated = sigma_batch.rotate(rotations)
sigma_restored = sigma_rotated.rotate(rotations.inv())
print("Rotation roundtrip error:", np.max(np.abs(sigma_batch.voigt - sigma_restored.voigt)))

# Single rotation broadcast to all tensors
single_rot = ScipyRotation.from_euler("z", 45, degrees=True)
sigma_rot_bc = sigma_batch.rotate(single_rot)
print("Broadcast rotation shape:", sigma_rot_bc.voigt.shape)

# %%
# 15. Batch push-forward / pull-back
# -------------------------------------
# Works with a single F (broadcast) or (N, 3, 3) per-point.

F_single = np.eye(3) + 0.05 * rng.standard_normal((3, 3))
pushed = eps_batch.push_forward(F_single)
pulled = pushed.pull_back(F_single)
print("Push/pull roundtrip error:", np.max(np.abs(eps_batch.voigt - pulled.voigt)))

# Per-point deformation gradients
F_batch = np.eye(3) + 0.05 * rng.standard_normal((N, 3, 3))
pushed_pp = eps_batch.push_forward(F_batch)
pulled_pp = pushed_pp.pull_back(F_batch)
print("Per-point push/pull error:", np.max(np.abs(eps_batch.voigt - pulled_pp.voigt)))

# %%
# 16. Batch Tensor4 inverse and rotation
# -----------------------------------------

L_batch = sim.Tensor4.from_tensor(
    sim.Tensor4.stiffness(sim.L_iso([70000, 0.3], "Enu")), N
)
M_batch = L_batch.inverse()
print("M_batch.type:", M_batch.type)  # compliance
print("L -> M -> L roundtrip:", np.allclose(L_batch.voigt, M_batch.inverse().voigt, atol=1e-8))

L_rotated = L_batch.rotate(rotations)
L_restored = L_rotated.rotate(rotations.inv())
print("L rotation roundtrip error:", np.max(np.abs(L_batch.voigt - L_restored.voigt)))

# %%
# 17. Concatenation and construction from lists
# ------------------------------------------------

part1 = sim.Tensor2.stress(rng.standard_normal((50, 6)))
part2 = sim.Tensor2.stress(rng.standard_normal((50, 6)))
combined = sim.Tensor2.concatenate([part1, part2])
print("Combined length:", len(combined))

singles = [sim.Tensor2.stress(np.array([float(i)] * 6)) for i in range(5)]
batch_from_list = sim.Tensor2(singles)
print("From list:", len(batch_from_list), "tensors")

# %%
# 18. numpy interop
# -------------------
# ``np.asarray(batch)`` returns the Voigt data.

arr = np.asarray(sigma_batch)
print("np.asarray shape:", arr.shape)  # (100, 6)
print("dtype:", arr.dtype)

# %%
# 19. Arithmetic on batches
# ----------------------------
# Element-wise ``+``, ``-``, scalar ``*``, per-point ``*``.

double = sigma_batch * 2.0
diff = double - sigma_batch
print("2*sigma - sigma == sigma:", np.allclose(diff.voigt, sigma_batch.voigt))

factors = rng.uniform(0.5, 1.5, N)
scaled = sigma_batch * factors
print("Per-point scaling shape:", scaled.voigt.shape)

# %%
# 20. Material-frame equivalence — batch version
# --------------------------------------------------
# Same identity as the single case, but across 100 orientations at once.

L_mat = sim.Tensor4.stiffness(sim.L_cubic([200000, 0.3, 80000], "EnuG"))
L_batch = sim.Tensor4.from_tensor(L_mat, N)
eps_batch = sim.Tensor2.strain(rng.standard_normal((N, 6)) * 0.01)

sigma_mat = L_batch @ eps_batch                          # material frame
L_global = L_batch.rotate(rotations)                     # rotate L
eps_global = eps_batch.rotate(rotations)                  # rotate eps
sigma_global = L_global @ eps_global                      # global frame
sigma_check = sigma_mat.rotate(rotations)                 # rotate result

print("Batch frame equivalence error:",
      np.max(np.abs(sigma_global.voigt - sigma_check.voigt)))
