"""Batch Lt_convert runs its points in parallel: it must match the point-by-point conversion,
for every converter, on distinct per-point data (a wrong slice or an index race shows)."""

import numpy as np
import pytest

import simcoon as sim

_KEYS = ["Dsigma_LieDD_2_DSDE", "DsigmaDe_2_DSDE", "DsigmaDe_JaumannDD_2_DSDE",
         "Dsigma_LieDD_Dsigma_JaumannDD", "Dsigma_LieDD_Dsigma_GreenNaghdiDD",
         "Dsigma_LieDD_Dsigma_logarithmicDD", "DsigmaDe_GreenNaghdiDD_2_DSDE",
         "DSDE_2_Dsigma_GreenNaghdiDD", "DSDE_2_Dsigma_JaumannDD", "DSDE_2_Dsigma_LieDD",
         "DSDE_2_Dsigma_logarithmicDD"]


def _batch(n, seed=3):
    rng = np.random.default_rng(seed)
    L = np.asarray(sim.L_iso([70000., 0.3], "Enu"))
    Lt = np.asfortranarray(L[:, :, None] * (1. + 0.1 * rng.random(n)))
    F = np.asfortranarray(np.eye(3)[:, :, None] + 0.05 * rng.standard_normal((3, 3, n)))
    stress = np.asfortranarray(100. * rng.standard_normal((6, n)))
    return Lt, F, stress


@pytest.mark.parametrize("key", _KEYS)
def test_batch_matches_single_points(key):
    n = 300                                   # past the parallel cutoff
    Lt, F, stress = _batch(n)
    batch = np.asarray(sim.Lt_convert(Lt, F, stress, key))
    for i in (0, 1, 149, n - 1):
        single = np.asarray(sim.Lt_convert(np.ascontiguousarray(Lt[:, :, i]),
                                           np.ascontiguousarray(F[:, :, i]),
                                           np.ascontiguousarray(stress[:, i]), key))
        np.testing.assert_allclose(batch[:, :, i], single, rtol=1e-12, atol=1e-9)


def test_batch_rejects_mismatched_point_counts():
    Lt, F, stress = _batch(10)
    with pytest.raises(ValueError, match="one entry per point"):
        sim.Lt_convert(Lt, np.asfortranarray(F[:, :, :9]), stress, "DsigmaDe_2_DSDE")
