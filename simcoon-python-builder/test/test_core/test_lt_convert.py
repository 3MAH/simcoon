"""Value contract of ``Lt_convert``.

The only other test touching this function checks that it does not leak memory
(test_lt_convert_memory.py); nothing pinned what it *returns*. These tests are the
regression gate for changes inside the converter and its per-point loop.

No golden values are stored. Everything here is a property that must hold whatever
the implementation does internally:

* the batched (6, 6, N) form must equal the single-point (6, 6) form, point by point;
* the map's inverse pairs must round-trip;
* an unknown converter key must raise rather than silently apply a different one.
"""

import numpy as np
import pytest

import simcoon as sim

# Every key of the converter map (objective_rates.cpp). Split by what they consume so
# the tests can drive each with a physically sensible input.
BOX_TO_MATERIAL = [          # consume the box tangent d(tau_hat)/d(De) -> material dS/dE
    "Dsigma_LieDD_2_DSDE",
    "DsigmaDe_2_DSDE",
    "DsigmaDe_JaumannDD_2_DSDE",
    "DsigmaDe_GreenNaghdiDD_2_DSDE",
]
MATERIAL_TO_SPATIAL = [      # consume dS/dE -> a spatial tangent
    "DSDE_2_Dsigma_GreenNaghdiDD",
    "DSDE_2_Dsigma_JaumannDD",
    "DSDE_2_Dsigma_LieDD",
    "DSDE_2_Dsigma_logarithmicDD",
]
SPATIAL_TO_SPATIAL = [       # cross maps between spatial rates
    "Dsigma_LieDD_Dsigma_JaumannDD",
    "Dsigma_LieDD_Dsigma_GreenNaghdiDD",
    "Dsigma_LieDD_Dsigma_logarithmicDD",
]
ALL_KEYS = BOX_TO_MATERIAL + MATERIAL_TO_SPATIAL + SPATIAL_TO_SPATIAL


def _batch(n, seed=12345):
    """A deliberately awkward batch: anisotropic Lt, compressive J, non-coaxial stress."""
    rng = np.random.default_rng(seed)
    # symmetric (major-symmetric) 6x6 tangents, well scaled
    A = rng.normal(size=(n, 6, 6))
    Lt = np.empty((6, 6, n), order="F")
    for k in range(n):
        s = 0.5 * (A[k] + A[k].T) + 6.0 * np.eye(6)
        Lt[:, :, k] = s
    # deformation gradients with det spread over [0.2, 3]. det F must stay POSITIVE:
    # a random perturbation of the identity can flip orientation, and an F with
    # det < 0 is not a deformation at all.
    F = np.empty((3, 3, n), order="F")
    k = 0
    while k < n:
        M = np.eye(3) + 0.35 * rng.normal(size=(3, 3))
        J = np.linalg.det(M)
        if J <= 0.05:
            continue                                # orientation flipped or near-singular
        if J < 0.2 or J > 3.0:                      # rescale into range, keep it non-symmetric
            M = M * (1.0 / J) ** (1.0 / 3.0)
        F[:, :, k] = M
        k += 1
    # stress in Voigt, with shear -> not coaxial with anything
    stress = np.asfortranarray(rng.normal(scale=40.0, size=(6, n)))
    return Lt, F, stress


@pytest.mark.parametrize("key", ALL_KEYS)
def test_batched_equals_pointwise(key):
    """The batch loop must do exactly what the single-point call does.

    This is the invariant that any change to the per-point loop (allocation,
    ordering, parallelisation) has to preserve, and it needs no reference data.
    """
    n = 24
    Lt, F, stress = _batch(n)

    batched = np.asarray(sim.Lt_convert(Lt, F, stress, key))
    assert batched.shape == (6, 6, n)

    for k in range(n):
        single = np.asarray(sim.Lt_convert(
            np.asfortranarray(Lt[:, :, k]),
            np.asfortranarray(F[:, :, k]),
            np.asfortranarray(stress[:, k]),
            key))
        np.testing.assert_allclose(
            batched[:, :, k], single, rtol=1e-13, atol=1e-11,
            err_msg=f"{key}: batch slice {k} differs from the single-point call")


def test_lie_round_trip_carries_exactly_one_factor_of_J():
    """`DSDE_2_Dsigma_LieDD` then `Dsigma_LieDD_2_DSDE` returns Lt / J, not Lt.

    That is deliberate and is the convention split documented at the top of
    Lt_convert: the INVERSE keys consume the no-J Kirchhoff box, while the FORWARD
    keys produce the Cauchy (1/J) spatial tangent. Composing one of each therefore
    leaves exactly one 1/J behind. Pinning it here means a future change that
    "tidies up" the J placement on either side cannot pass unnoticed -- which
    matters because fedoo chains exactly these two keys in its UL path.
    """
    Lt, F, stress = _batch(16, seed=7)
    there = sim.Lt_convert(Lt, F, stress, "DSDE_2_Dsigma_LieDD")
    back = np.asarray(sim.Lt_convert(np.asfortranarray(there), F, stress,
                                     "Dsigma_LieDD_2_DSDE"))
    for k in range(Lt.shape[2]):
        J = np.linalg.det(F[:, :, k])
        np.testing.assert_allclose(
            back[:, :, k], Lt[:, :, k] / J, rtol=1e-9,
            atol=1e-9 * np.abs(Lt[:, :, k]).max(),
            err_msg=f"round trip is not exactly Lt/J at slice {k} (det F = {J:.4f})")


@pytest.mark.parametrize("key", ALL_KEYS)
def test_undeformed_unstressed_is_a_no_op(key):
    """At F = I with zero stress every conversion degenerates to the identity.

    J = 1, the strain-concentration operators reduce to the identity and every
    stress-dependent correction vanishes, so the tangent must come back untouched.
    """
    n = 4
    rng = np.random.default_rng(3)
    Lt = np.empty((6, 6, n), order="F")
    for k in range(n):
        s = rng.normal(size=(6, 6))
        Lt[:, :, k] = 0.5 * (s + s.T) + 6.0 * np.eye(6)
    F = np.asfortranarray(np.tile(np.eye(3)[:, :, None], (1, 1, n)))
    stress = np.zeros((6, n), order="F")

    out = np.asarray(sim.Lt_convert(Lt, F, stress, key))
    np.testing.assert_allclose(out, Lt, rtol=1e-12, atol=1e-12,
                               err_msg=f"{key} is not the identity at F = I, sigma = 0")


def test_unknown_key_raises():
    """An unrecognised converter key must fail loudly.

    It used to index a non-const std::map with operator[], which default-inserts 0 --
    so a typo silently performed the *first* conversion in the map and returned a
    plausible-looking but wrong tangent. fedoo passes this key from a table, so a
    rename would have surfaced as poor Newton convergence, not as an error.
    """
    Lt, F, stress = _batch(4)
    with pytest.raises(Exception):
        sim.Lt_convert(Lt, F, stress, "NotAConverterKey")


def test_unknown_key_is_not_silently_the_first_entry():
    """Sharper form of the above: whatever happens, it must not equal key 0's answer."""
    Lt, F, stress = _batch(4)
    first = np.asarray(sim.Lt_convert(Lt, F, stress, ALL_KEYS[0]))
    try:
        bogus = np.asarray(sim.Lt_convert(Lt, F, stress, "NotAConverterKey"))
    except Exception:
        return                                    # raised: the correct behaviour
    assert not np.allclose(bogus, first), (
        "an unknown key silently fell through to the first converter in the map")
