import pytest
import numpy as np
import numpy.typing as npt
from simcoon import _core as sim


@pytest.fixture(scope="session")
def L_gt() -> npt.NDArray[np.float64]:
    L_gt = np.array(
        [
            [1.34615385, 0.57692308, 0.57692308, 0.0, 0.0, 0.0],
            [0.57692308, 1.34615385, 0.57692308, 0.0, 0.0, 0.0],
            [0.57692308, 0.57692308, 1.34615385, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.38461538, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.38461538, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.38461538],
        ]
    )
    return L_gt


def test_L_iso(L_gt):
    """Test L_iso with E=1., nu=0.3."""
    L = sim.L_iso([1.0, 0.3], "Enu")
    assert np.allclose(L, L_gt)


def test_isochoric_invariants():

    V = np.random.rand(3, 3) + np.eye(3, 3)
    V = 0.5 * (V + V.transpose())
    b = V * V
    lambda_bar = sim.isochoric_invariants(b)
    b_vec = sim.t2v_strain(b)
    lambda_bar_vec = sim.isochoric_invariants(b_vec)
    assert np.allclose(lambda_bar, lambda_bar_vec)


# ---------------------------------------------------------------------------
# Shared fixture: small deformation gradient pair for objective_rate tests
# ---------------------------------------------------------------------------
@pytest.fixture(scope="module")
def deformation_pair():
    """Return (F0, F1, DTime) for a small deformation increment."""
    F0 = np.eye(3)
    F1 = np.eye(3)
    F1[0, 0] = 1.01
    F1[0, 1] = 0.005
    F1[1, 0] = -0.005
    F1[1, 1] = 1.02
    F1[2, 2] = 0.98
    DTime = 0.001
    return F0, F1, DTime


ALL_RATE_NAMES = ["jaumann", "green_naghdi", "logarithmic", "logarithmic_R", "truesdell", "logarithmic_F"]


def test_objective_rate_all_rates_D_consistent(deformation_pair):
    """All 6 objective rates should produce the same D tensor."""
    F0, F1, DTime = deformation_pair
    D_ref = None
    for name in ALL_RATE_NAMES:
        result = sim.objective_rate(name, F0, F1, DTime)
        D = result[0]
        assert D.shape == (3, 3), f"{name}: D shape {D.shape}"
        # D should be symmetric
        assert np.allclose(D, D.T, atol=1e-9), f"{name}: D not symmetric"
        if D_ref is None:
            D_ref = D
        else:
            assert np.allclose(D, D_ref, atol=1e-9), f"{name}: D differs from jaumann"


def test_objective_rate_truesdell(deformation_pair):
    """Test truesdell rate: shapes and return_de."""
    F0, F1, DTime = deformation_pair
    D, DR, Omega = sim.objective_rate("truesdell", F0, F1, DTime)
    assert D.shape == (3, 3)
    assert DR.shape == (3, 3)
    assert Omega.shape == (3, 3)
    # D symmetric
    assert np.allclose(D, D.T, atol=1e-9)

    # return_de mode
    de, D2, DR2, Omega2 = sim.objective_rate("truesdell", F0, F1, DTime, return_de=True)
    assert de.shape[0] == 6
    assert np.allclose(D, D2, atol=1e-12)


def test_objective_rate_log_F(deformation_pair):
    """Test logarithmic_F rate and its alias log_F."""
    F0, F1, DTime = deformation_pair
    D1, DR1, Omega1 = sim.objective_rate("logarithmic_F", F0, F1, DTime)
    D2, DR2, Omega2 = sim.objective_rate("log_F", F0, F1, DTime)

    assert D1.shape == (3, 3)
    assert np.allclose(D1, D1.T, atol=1e-9)
    # alias should give identical results
    assert np.allclose(D1, D2, atol=1e-12)
    assert np.allclose(DR1, DR2, atol=1e-12)
    assert np.allclose(Omega1, Omega2, atol=1e-12)

    # return_de mode
    de, D3, DR3, Omega3 = sim.objective_rate("log_F", F0, F1, DTime, return_de=True)
    assert de.shape[0] == 6
    assert np.allclose(D1, D3, atol=1e-12)


def test_objective_rate_batch(deformation_pair):
    """Test 3D batch mode (F0=2D, F1=3D) for truesdell and log_F."""
    F0, F1, DTime = deformation_pair
    n_points = 4
    # Stack n_points copies of F1 into a 3×3×n cube
    F1_batch = np.stack([F1] * n_points, axis=-1)
    assert F1_batch.shape == (3, 3, n_points)

    for name in ["truesdell", "logarithmic_F"]:
        D_batch, DR_batch, Omega_batch = sim.objective_rate(name, F0, F1_batch, DTime)
        assert D_batch.shape == (3, 3, n_points), f"{name}: D_batch shape {D_batch.shape}"
        assert DR_batch.shape == (3, 3, n_points)
        assert Omega_batch.shape == (3, 3, n_points)

        # Single-point reference
        D_single, DR_single, Omega_single = sim.objective_rate(name, F0, F1, DTime)
        for i in range(n_points):
            assert np.allclose(D_batch[:, :, i], D_single, atol=1e-12), f"{name}: batch D[{i}] mismatch"


def test_objective_rate_return_de_small_increment():
    """For very small deformation, de ≈ DTime * t2v_strain(D) for all rates."""
    F0 = np.eye(3)
    eps = 1e-6
    F1 = np.eye(3) + eps * np.array([
        [1.0, 0.1, 0.0],
        [0.1, -0.5, 0.0],
        [0.0, 0.0, -0.5],
    ])
    DTime = 1e-3

    for name in ALL_RATE_NAMES:
        result_de = sim.objective_rate(name, F0, F1, DTime, return_de=True)
        de = np.asarray(result_de[0]).ravel()
        D = result_de[1]
        assert de.shape == (6,), f"{name}: de shape {de.shape}"
        # For small increments, de ≈ DTime * t2v_strain(D)
        de_approx = np.asarray(DTime * sim.t2v_strain(D)).ravel()
        assert np.allclose(de, de_approx, rtol=1e-3), (
            f"{name}: de not close to DTime*t2v_strain(D), "
            f"max diff = {np.max(np.abs(de - de_approx))}"
        )
