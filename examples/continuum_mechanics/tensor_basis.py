"""
The basis of a tensor: frames, convected bases and the metric
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The components of a ``Tensor2`` / ``Tensor4`` only mean something together with
the basis they are written in. By default this is the fixed orthonormal lab
basis: ``t.basis`` is ``None`` and nothing is stored. A ``simcoon.Basis`` makes
another choice explicit:

* an **orthonormal** basis, built from a rotation (a material frame, a
  corotational frame). Its metric is the identity;
* a **natural** basis, any three independent vectors, typically the convected
  basis ``g_i = F G_i`` of a deformation gradient. Its metric ``g = A^T A``
  enters the invariants.

Three operations, three meanings:

* ``t.to_basis(b)``: the same tensor, components re-expressed in ``b``;
* ``t.with_basis(b)``: the same components in another basis, a transport;
* ``t.rotate(R, active=False)``: the same tensor in the frame turned by ``R``.

The type tag gives the variance of the components: stress contravariant, strain
covariant, stiffness and compliance likewise on four indices, concentration
tensors mixed.
"""

import copy
import pickle

import numpy as np
import simcoon as sim

np.set_printoptions(precision=4, suppress=True)

# %%
# 1. The lab basis is the default
# ---------------------------------
# Nothing changes for a tensor that never leaves the lab: ``basis`` is ``None``.

sigma = sim.Tensor2.stress(np.array([100.0, 50.0, 0.0, 30.0, 0.0, 0.0]))
eps = sim.Tensor2.strain(np.array([0.010, -0.003, -0.003, 0.005, 0.002, 0.001]))
L_iso = sim.Tensor4.stiffness(sim.L_iso([70000.0, 0.3], "Enu"))

print("sigma.basis:", sigma.basis)
print("(L_iso @ eps).basis:", (L_iso @ eps).basis)

# %%
# 2. An orthonormal basis and its description
# ----------------------------------------------
# A ``Basis`` holds the basis vectors as the columns of ``matrix`` (their lab
# components). The ``name`` is a label for printing; two bases are the same
# when their vectors are.

rot = sim.Rotation.from_euler("zxz", [30.0, 20.0, 10.0], degrees=True)
material = sim.Basis(rotation=rot, name="material")

print(material)
print("orthonormal:", material.orthonormal, "| single:", material.single,
      "| name:", material.name)
print("rotation is kept:", material.rotation.equals(rot))
print("matrix (columns = basis vectors):\n", material.matrix)
print("metric is the identity:", np.allclose(material.metric, np.eye(3)))
print("same vectors, other name -> equal:",
      material.equals(sim.Basis(rotation=rot, name="other")))

# %%
# 3. ``to_basis``: the same tensor in another basis
# ----------------------------------------------------
# The components change, the tensor does not: its invariants are those of the
# lab, and ``to_basis(None)`` gives the lab components back.

sigma_m = sigma.to_basis(material)
print(sigma_m)
print("material components:", sigma_m.voigt)
print("same Mises:", np.isclose(sigma_m.mises(), sigma.mises()),
      "| same trace:", np.isclose(sigma_m.trace(), sigma.trace()))
print("back to the lab:", np.allclose(sigma_m.to_basis(None).voigt, sigma.voigt))

# Between two frames directly
other = sim.Basis(rotation=sim.Rotation.from_euler("z", 45.0, degrees=True))
print("material -> other == lab -> other:",
      np.allclose(sigma_m.to_basis(other).voigt, sigma.to_basis(other).voigt))

# %%
# 4. Tensors in different bases cannot be mixed
# ------------------------------------------------
# Adding material-frame components to lab components is meaningless; it is
# refused instead of returning numbers.

for label, operation in [("sigma + sigma_m", lambda: sigma + sigma_m),
                         ("L_iso @ eps_m", lambda: L_iso @ eps.to_basis(material)),
                         ("sigma % sigma_m", lambda: sigma % sigma_m)]:
    try:
        operation()
    except ValueError as err:
        print(f"{label}: {str(err).split(' (')[0]}")

print("comparison across bases is False:", sigma == sigma_m)
print("in a common basis it works:", (sigma_m + sigma.to_basis(material)).basis)

# %%
# 5. A passive rotation ships its frame
# ----------------------------------------
# ``rotate(R, active=False)`` returns the same tensor seen from the frame
# turned by ``R``. The result now remembers that frame; a second passive
# rotation is read in the current frame, so the frames compose.

sigma_p = sigma.rotate(rot, active=False)
print("basis after a passive rotation:", sigma_p.basis)
print("same numbers as to_basis:", np.allclose(sigma_p.voigt, sigma_m.voigt))

rot2 = sim.Rotation.from_euler("x", 35.0, degrees=True)
sigma_pp = sigma_p.rotate(rot2, active=False)
print("composed frame = rot * rot2:", sigma_pp.basis.rotation.equals(rot * rot2))
print("still the same tensor:", np.allclose(sigma_pp.to_basis(None).voigt, sigma.voigt))

# %%
# 6. An active rotation is a transport
# ---------------------------------------
# ``rotate(R)`` gives another tensor, ``Q X Q^T``. A lab tensor stays written
# in the lab. A tensor with its own basis keeps its components and has its
# basis vectors turned: the same result, without touching the numbers.

sigma_rot = sigma.rotate(rot2)
print("lab tensor, active rotation -> basis:", sigma_rot.basis)

turned = sigma_m.rotate(rot2)
print("framed tensor: components kept:", np.array_equal(turned.voigt, sigma_m.voigt))
print("its basis is rot2 * rot:", turned.basis.rotation.equals(rot2 * rot))
print("same tensor as the lab rotation:",
      np.allclose(turned.to_basis(None).voigt, sigma_rot.voigt))

# %%
# 7. A constitutive law written in its material frame
# -------------------------------------------------------
# A cubic stiffness is known in the crystal frame. Tag it with that basis,
# bring the strain to the same basis, contract, and take the result back.

L_crystal = sim.Tensor4.stiffness(
    sim.L_cubic([200000.0, 0.3, 80000.0], "EnuG")).with_basis(material)
sigma_crystal = L_crystal @ eps.to_basis(material)
print("stress is in the crystal frame:", sigma_crystal.basis)

L_lab = L_crystal.to_basis(None)                      # the same tensor, lab components
print("same stress as the lab computation:",
      np.allclose(sigma_crystal.to_basis(None).voigt, (L_lab @ eps).voigt))
print("compliance keeps the basis:", L_crystal.inverse().basis is material)

# %%
# 8. Batches: one basis for all, or one per tensor
# ---------------------------------------------------
# A single basis serves a whole batch: it is stored once and shared by
# reference by every tensor derived from it. A batch of N rotations gives one
# basis per tensor (one orientation per grain or per Gauss point).

N = 1000
rng = np.random.default_rng(42)
sigma_batch = sim.Tensor2.stress(rng.standard_normal((N, 6)) * 100.0)

shared = sigma_batch.to_basis(material)
print("one object for", N, "tensors:", shared.basis is material,
      "| also after indexing:", shared[17].basis is material,
      "| and after arithmetic:", (2.0 * shared - shared).basis is material)

grains = sim.Basis(rotation=sim.Rotation.random(N, random_state=1), name="grains")
per_grain = sigma_batch.to_basis(grains)
print(per_grain)
print("basis of tensor 17:", per_grain[17].basis)
print("round trip error:", np.max(np.abs(per_grain.to_basis(None).voigt - sigma_batch.voigt)))

# Stacking: equal bases stay shared, different ones become a per-tensor basis
a = sigma.to_basis(material)
b = sigma.to_basis(other)
print("from_list, same basis:", sim.Tensor2.from_list([a, a]).basis)
print("from_list, two bases:", sim.Tensor2.from_list([a, b]).basis)
print("concatenate:", sim.Tensor2.concatenate([shared[:3], shared[10:12]]).basis)

# %%
# 9. A natural basis and its metric
# ------------------------------------
# Simple shear ``F = I + gamma e1 (x) e2``. The convected basis ``g_i = F e_i``
# is not orthonormal: ``g_2 = (gamma, 1, 0)``. Its metric is the right
# Cauchy-Green tensor ``C = F^T F``. ``from_F`` keeps a reference to ``F``.

gamma = 0.5
F = np.array([[1.0, gamma, 0.0],
              [0.0, 1.0, 0.0],
              [0.0, 0.0, 1.0]])
convected = sim.Basis.from_F(F, name="convected")

print(convected, "| orthonormal:", convected.orthonormal)
print("metric g_ij = C:\n", convected.metric)
print("inverse metric g^ij:\n", convected.inverse_metric)
print("volume det(A) = J:", convected.det)
print("reciprocal vectors g^i (columns of A^-T):\n", convected.reciprocal)
print("right stretch U = sqrt(g):\n", convected.stretch)
print("closest orthonormal basis, R = A U^-1:", convected.polar)
print(convected.polar.matrix)
print("any three vectors work too:", sim.Basis(vectors=np.diag([2.0, 1.0, 0.5])).metric.diagonal())

# %%
# 10. ``with_basis``: transport by keeping the components
# -----------------------------------------------------------
# The Kirchhoff stress has, in the convected basis, the components of the
# second Piola-Kirchhoff stress: ``tau^ij = S^IJ``. Likewise the Almansi strain
# has the components of the Green-Lagrange strain: ``e_ij = E_IJ``.

S = sim.Tensor2.stress(np.array([0.0, 80.0, 0.0, 0.0, 0.0, 0.0]))      # PK2 = 80 e2 (x) e2
tau = S.with_basis(convected)
print("same numbers:", np.array_equal(tau.voigt, S.voigt), "| basis:", tau.basis)
print("lab components of tau:\n", tau.to_basis(None).mat)
print("== F S F^T:", np.allclose(tau.to_basis(None).mat, F @ S.mat @ F.T))

E = sim.Tensor2.strain(0.5 * (F.T @ F - np.eye(3)))                    # Green-Lagrange
almansi = E.with_basis(convected).to_basis(None)
print("Almansi from the convected components:\n", almansi.mat)
print("== 1/2 (I - b^-1):",
      np.allclose(almansi.mat, 0.5 * (np.eye(3) - np.linalg.inv(F @ F.T))))

# %%
# 11. Invariants with the metric
# ---------------------------------
# In a natural basis the trace is ``g_ij tau^ij``, not the sum of the diagonal
# components. Deviator, norm, von Mises, determinant and eigenvalues follow.

tau_lab = tau.to_basis(None)
print("trace:      ", tau.trace(), "| naive diagonal sum:", tau.voigt[:3].sum(),
      "| lab:", tau_lab.trace())
print("von Mises:  ", tau.mises(), "| lab:", tau_lab.mises())
print("norm:       ", tau.norm(), "| lab:", tau_lab.norm())

# Determinant and principal values, on a full stress in a basis that also changes volume:
# det T = det(T^ij) det(g), and the eigenvalues solve det(T^ij - lambda g^ij) = 0.
F_gen = np.array([[1.20, 0.15, -0.05],
                  [0.10, 0.90, 0.20],
                  [0.00, -0.10, 1.10]])
general = sim.Basis.from_F(F_gen)
full = sim.Tensor2.stress(np.array([100.0, 50.0, 20.0, 30.0, -10.0, 5.0]))
full_g = full.to_basis(general)
print("det:        ", round(full_g.det(), 3), "| lab:", round(full.det(), 3),
      "| naive det of the components:", round(np.linalg.det(full_g.mat), 3))
print("eigenvalues:", full_g.eigvals(), "| lab:", np.linalg.eigvalsh(full.mat))

s = tau.dev()                                         # tau - tr/3 g^-1, same basis
print("deviator is traceless:", np.isclose(s.trace(), 0.0), "| basis kept:", s.basis is convected)
print("deviator, lab components match:",
      np.allclose(s.to_basis(None).voigt, tau_lab.dev().voigt))

# %%
# 12. Contractions: dual variances need no metric
# --------------------------------------------------
# ``sigma : eps = sigma^ij eps_ij`` acts on the components in any basis. Two
# tensors of the same variance need the metric on both indices.

e = eps.to_basis(convected)                           # covariant components
t = sigma.to_basis(convected)                         # contravariant components
print(f"sigma : eps, lab:         {sigma % eps:.6f}")
print(f"sigma : eps, convected:   {t % e:.6f}  (plain sum of products)")
print(f"sigma : sigma, lab:       {sigma % sigma:.3f}")
print(f"sigma : sigma, convected: {t % t:.3f}  (with the metric)")
print(f"naive sum of squares:     {np.sum(t.mat * t.mat):.3f}  (not an invariant)")
print("double_contract:", sim.double_contract(t, e))

# %%
# 13. The identity tensor is the metric
# ----------------------------------------
# Written with contravariant components the identity is ``g^ij``; with
# covariant components it is ``g_ij``.

print("identity, stress-typed  == g^ij:",
      np.allclose(sim.Tensor2.identity("stress", basis=convected).mat, convected.inverse_metric))
print("identity, strain-typed  == g_ij:",
      np.allclose(sim.Tensor2.identity("strain", basis=convected).mat, convected.metric))
print("in an orthonormal basis:", sim.Tensor2.identity("stress", basis=material).voigt)

# %%
# 14. Fourth-order tensors in a natural basis
# ----------------------------------------------
# The convected components of the spatial stiffness are the lab components of
# its pull-back. Contraction and inverse act on the components; the identity
# and the projectors are built from the metric.

L_c = L_iso.to_basis(convected)
print("L in the convected basis == L.pull_back(F):",
      np.allclose(L_c.mat, L_iso.pull_back(F, metric=False).mat))
print("L : eps there, back in the lab == lab result:",
      np.allclose((L_c @ e).to_basis(None).voigt, (L_iso @ eps).voigt))
print("inverse is the compliance in that basis:",
      np.allclose(L_c.inverse().to_basis(None).mat, L_iso.inverse().mat))

I4 = sim.Tensor4.identity("stiffness", basis=convected)
ginv = convected.inverse_metric
print("identity 1/2 (g^ik g^jl + g^il g^jk) raises both indices:",
      np.allclose((I4 @ e).mat, ginv @ e.mat @ ginv))
print("it is not eye(6) any more:", not np.allclose(I4.mandel, np.eye(6)))
print("the mixed (concentration) identity still is:",
      np.allclose(sim.Tensor4.identity("strain_concentration", basis=convected).mat, np.eye(6)))

K, mu = 70000.0 / (3 * (1 - 2 * 0.3)), 70000.0 / (2 * (1 + 0.3))
P_vol = sim.Tensor4.volumetric("stiffness", basis=convected)
P_dev = sim.Tensor4.deviatoric("stiffness", basis=convected)
print("L = 3K P_vol + 2 mu P_dev with the metric projectors:",
      np.allclose(L_c.mat, (3 * K * P_vol + 2 * mu * P_dev).mat))

# A concentration tensor (mixed variance) changes basis too
A = sim.Tensor4.strain_concentration(sim.A_R(F_gen))
print("concentration tensor, same result through the convected basis:",
      np.allclose((A.to_basis(convected) @ e).to_basis(None).voigt, (A @ eps).voigt))

try:
    L_c @ t                                           # stiffness with a stress: not dual
except ValueError as err:
    print("refused:", err)

# %%
# 15. Push-forward and pull-back carry the basis along
# -------------------------------------------------------
# A push-forward is a transport: another tensor, in the current configuration.
# In the convected basis it has the components the original had in the
# reference basis. On a lab tensor ``push_forward`` returns lab components, as
# always. On a tensor that has its own basis it keeps the components and
# convects the basis: nothing is computed until lab components are asked for.

reference = sim.Basis(rotation=sim.Rotation.identity(), name="reference")
S_ref = S.with_basis(reference)

tau_lazy = S_ref.push_forward(F, metric=False)
print("components untouched:", np.array_equal(tau_lazy.voigt, S.voigt))
print("basis is now F:", np.allclose(tau_lazy.basis.matrix, F))
print("lab components == eager push-forward:",
      np.allclose(tau_lazy.to_basis(None).voigt, S.push_forward(F, metric=False).voigt))

# metric=True adds the Piola weight 1/J (Kirchhoff -> Cauchy); use a volume change to see it
F_vol = F @ np.diag([1.1, 1.0, 1.0])
cauchy_lazy = S_ref.push_forward(F_vol)               # metric=True is the default
print("only the 1/J weight touches the numbers:",
      np.allclose(cauchy_lazy.voigt, S.voigt / np.linalg.det(F_vol)))
print("Cauchy, lab == eager:",
      np.allclose(cauchy_lazy.to_basis(None).voigt, S.push_forward(F_vol).voigt))

back = cauchy_lazy.pull_back(F_vol)
print("pull_back returns to the reference basis and components:",
      np.allclose(back.basis.matrix, np.eye(3)), np.allclose(back.voigt, S.voigt))

# The same holds for a stiffness: the lazy transport costs nothing
C_ref = L_iso.with_basis(reference)
c_lazy = C_ref.push_forward(F, metric=False)
print("Tensor4, lab == eager:",
      np.allclose(c_lazy.to_basis(None).mat, L_iso.push_forward(F, metric=False).mat))

# %%
# 16. Objectivity
# ------------------
# A superposed rigid rotation turns the convected basis and leaves the
# components, the metric and every invariant unchanged.

Q = sim.Rotation.from_euler("zxz", [70.0, 25.0, -40.0], degrees=True)
tau_Q = tau.rotate(Q)
print("components unchanged:", np.array_equal(tau_Q.voigt, tau.voigt))
print("basis is Q F:", np.allclose(tau_Q.basis.matrix, Q.as_matrix() @ F))
print("metric unchanged:", np.allclose(tau_Q.basis.metric, tau.basis.metric))
print("Mises unchanged:", np.isclose(tau_Q.mises(), tau.mises()))
print("lab components are those of Q tau Q^T:",
      np.allclose(tau_Q.to_basis(None).voigt, tau_lab.rotate(Q).voigt))

# %%
# 17. A convected basis per Gauss point
# ----------------------------------------
# With one deformation gradient per point the basis is a batch; it refers to
# the caller's ``F`` array without copying it.

F_batch = np.eye(3) + 0.15 * rng.standard_normal((N, 3, 3))
conv_batch = sim.Basis.from_F(F_batch)
S_batch = sim.Tensor2.stress(rng.standard_normal((N, 6)) * 50.0)
tau_batch = S_batch.with_basis(conv_batch)

print(tau_batch)
print("no arithmetic: components are those of S:", np.array_equal(tau_batch.voigt, S_batch.voigt))
print("Mises per point matches the lab:",
      np.allclose(tau_batch.mises(), tau_batch.to_basis(None).mises()))
print("lab components == batch push-forward:",
      np.allclose(tau_batch.to_basis(None).voigt,
                  S_batch.push_forward(F_batch, metric=False).voigt))
print("metric of point 3:\n", tau_batch[3].basis.metric)

# %%
# 18. Housekeeping
# -------------------
# ``np.asarray`` returns the components in the tensor's own basis. A basis
# survives copying and pickling. Tensor types without a variance cannot be
# written in a natural basis, and a name alone is not a basis.

print("np.asarray(tau):", np.asarray(tau))
clone = pickle.loads(pickle.dumps(sigma_m))
print("pickled:", clone, "| equal:", clone == sigma_m, "| deepcopy equal:",
      copy.deepcopy(sigma_m) == sigma_m)

symmetric = sim.Tensor2.from_voigt(sigma.voigt, "symmetric")
print("a symmetric tensor in an orthonormal basis is fine:", symmetric.to_basis(material).basis)
for label, operation in [("symmetric in a natural basis", lambda: symmetric.with_basis(convected)),
                         ("a name instead of a basis", lambda: sigma.with_basis("material"))]:
    try:
        operation()
    except (ValueError, TypeError) as err:
        print(f"{label}: {type(err).__name__}")

# A two-point or non-symmetric tensor (F, R, PK1) takes the type "none" -- no
# Voigt convention, 9 components: lab-lab components, rotation of both legs,
# invariants, but no Voigt vector, no variance and no transport.
F_t = sim.Tensor2.from_mat(F_gen, "none")
print(F_t, "| stored exactly:", np.array_equal(F_t.mat, F_gen), "| det J =", round(F_t.det(), 6))
print("rotated with both legs:", np.allclose(F_t.rotate(Q).mat, Q.as_matrix() @ F_gen @ Q.as_matrix().T))
for label, operation in [("a symmetric type refuses F", lambda: sim.Tensor2.stress(F_gen)),
                         ("F has no convected basis", lambda: F_t.with_basis(convected)),
                         ("F is not transported", lambda: F_t.push_forward(F_gen))]:
    try:
        operation()
    except ValueError as err:
        print(f"{label}: ValueError")

# %%
# 19. Bases that come from elsewhere: ply frames, shells, successive transports
# ------------------------------------------------------------------------------
# A ``Basis`` is just the matrix of its vectors, so every case goes through the
# same three operations. A composite ply knows its stress in the fibre frame;
# pushing it forward convects that frame (the first vector becomes the stretched
# fibre ``F a0``). ``F`` and the rotations of ``rotate`` are always given in lab
# components; a deformation gradient known in ply axes is brought to the lab
# first, ``F = R F_hat R^T``.

ply = sim.Basis(rotation=sim.Rotation.from_euler("z", 30.0, degrees=True), name="ply")
S_ply = sim.Tensor2.stress(np.array([500.0, 20.0, 0.0, 10.0, 0.0, 0.0])).with_basis(ply)  # PK2, ply axes
F_ply = np.array([[1.1, 0.4, 0.0], [0.0, 0.95, 0.0], [0.0, 0.0, 1.0]])

tau_ply = S_ply.push_forward(F_ply, metric=False)
print("convected ply basis = F R:", np.allclose(tau_ply.basis.matrix, F_ply @ ply.matrix))
fibre = tau_ply.basis.matrix[:, 0]
print("convected fibre F a0:", fibre, "| stretch:", round(np.linalg.norm(fibre), 4))
print("fibre component kept: tau^11 = S^11 =", tau_ply.voigt[0])
print("lab components:", tau_ply.to_basis(None).voigt)

# The orthonormal frame closest to the convected one: its polar rotation (the ply
# frame turned by the rotation of F)
print("in the polar (rotated orthonormal) ply frame:", tau_ply.to_basis(tau_ply.basis.polar).voigt)

# Successive transports compose, and the total pull-back restores the ply frame
F1 = np.array([[1.05, 0.1, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
F2 = np.array([[1.0, 0.0, 0.2], [0.0, 1.1, 0.0], [0.0, 0.0, 0.9]])
two_steps = S_ply.push_forward(F1, metric=False).push_forward(F2, metric=False)
print("F2 after F1 == F2 F1:", np.allclose(two_steps.basis.matrix, F2 @ F1 @ ply.matrix),
      "| components kept:", np.array_equal(two_steps.voigt, S_ply.voigt))
print("pull-back by F2 F1 restores the ply frame:",
      np.allclose(two_steps.pull_back(F2 @ F1, metric=False).basis.matrix, ply.matrix))

# A shell: the covariant surface basis (a_1, a_2, n) is a natural basis; a
# membrane stress has contravariant components N^{ab} in it.
a1 = np.array([1.0, 0.2, 0.0])
a2 = np.array([-0.1, 0.8, 0.3])
normal = np.cross(a1, a2)
normal /= np.linalg.norm(normal)
shell = sim.Basis(vectors=np.column_stack([a1, a2, normal]), name="shell")
N_shell = sim.Tensor2.stress(np.array([120.0, 40.0, 0.0, 15.0, 0.0, 0.0])).with_basis(shell)
print("surface metric a_ab:\n", shell.metric)
print("membrane stress: trace", round(N_shell.trace(), 4), "| Mises", round(N_shell.mises(), 4),
      "| lab Mises", round(N_shell.to_basis(None).mises(), 4))
lamina = sim.Basis(rotation=sim.Rotation.from_matrix(np.linalg.qr(shell.matrix)[0]), name="lamina")
print("the same stress in the orthonormal lamina frame:", N_shell.to_basis(lamina).voigt)

# %%
# 20. The variance tag: raising and lowering indices
# -----------------------------------------------------
# The variance of the components is a tag defaulted from the type (stress
# contravariant, strain covariant, stiffness (contra, contra), compliance (co, co),
# concentration tensors mixed). ``to_variance`` gives the same tensor with the
# other variance: nothing changes in the lab, the metric acts in a natural basis.

print("defaults:", sigma.variance, "|", eps.variance, "|", L_iso.variance, "|",
      sim.Tensor4.strain_concentration(np.eye(6)).variance)

tau_low = tau.to_variance("covariant")                  # tau_ij = g_ik tau^kl g_lj
print(tau_low)
print("lowered with the metric:", np.allclose(tau_low.mat, convected.metric @ tau.mat @ convected.metric))
print("same tensor: same Mises", np.isclose(tau_low.mises(), tau.mises()),
      "| same lab components", np.allclose(tau_low.to_basis(None).voigt, tau.to_basis(None).voigt))
print("raised back:", np.allclose(tau_low.to_variance("contravariant").mat, tau.mat))

# Contractions read the tags: dual variances contract as they are, equal ones
# through the metric; the value is the lab one whatever the tags
s_conv, e_conv = sigma.to_basis(convected), eps.to_basis(convected)
print("sigma : eps, lab:", round(sigma % eps, 6),
      "| sigma^sharp : eps_flat:", round(s_conv % e_conv, 6),
      "| sigma_flat : eps^sharp:", round(s_conv.to_variance("covariant") % e_conv.to_variance("contravariant"), 6),
      "| sigma_flat : eps_flat (metric):", round(s_conv.to_variance("covariant") % e_conv, 6))

# A stiffness with lowered indices contracts a contravariant strain
L_low = L_iso.to_basis(convected).to_variance(("covariant", "covariant"))
print(L_low)
print("L_low : e^sharp back in the lab == L : eps:",
      np.allclose((L_low @ e_conv.to_variance("contravariant")).to_basis(None).voigt, (L_iso @ eps).voigt))
try:
    L_low @ e_conv
except ValueError as err:
    print("not dual:", str(err).split(" (")[0])

# The variance is part of what a transport means: lowering the indices of a
# lab tensor and pushing it forward is not the push-forward of the original.
pushed = sigma.push_forward(F_gen, metric=False)
pushed_low = sigma.to_variance("covariant").push_forward(F_gen, metric=False)
print("F sigma F^T vs F^-T sigma F^-1 differ:", not np.allclose(pushed.mat, pushed_low.mat))

# A "symmetric" tensor has no default variance: declare it first
sym = sim.Tensor2.from_voigt(sigma.voigt, "symmetric")
print("symmetric, undeclared:", sym.variance, "| declared:", sym.to_variance("contravariant").variance)
