"""
Tests for UMAT material models using the Python solver.
"""

import pytest
import numpy as np
from pathlib import Path

from simcoon.solver import (
    StateVariablesM,
    StepMeca,
    Block,
    Solver,
)
from simcoon.solver.micromechanics import (
    Layer,
    Ellipsoid,
    MaterialOrientation,
    GeometryOrientation,
    save_layers_json,
    save_ellipsoids_json,
)


# =============================================================================
# Test Data Directory
# =============================================================================

TEST_DATA_DIR = Path(__file__).parent / "data"


# =============================================================================
# Helper Functions
# =============================================================================

def _has_simcoon_core() -> bool:
    """Check if simcoon._core is available."""
    try:
        from simcoon import _core as scc
        _ = scc.L_iso([1.0, 0.3], "Enu")
        return True
    except (ImportError, AttributeError, Exception):
        return False


def run_uniaxial_test(umat_name: str, props: np.ndarray, nstatev: int,
                      strain_max: float = 0.02, n_increments: int = 100,
                      control_type: str = 'small_strain') -> list:
    """
    Run a uniaxial tension test with loading and unloading.

    Parameters
    ----------
    umat_name : str
        Name of the UMAT (e.g., 'ELISO', 'ELIST')
    props : np.ndarray
        Material properties array
    nstatev : int
        Number of state variables
    strain_max : float
        Maximum strain in direction 1
    n_increments : int
        Number of increments per step
    control_type : str
        Control type ('small_strain', 'logarithmic', etc.)

    Returns
    -------
    list
        History of state variables
    """
    # Loading step: strain to max
    step_load = StepMeca(
        DEtot_end=np.array([strain_max, 0, 0, 0, 0, 0]),
        Dsigma_end=np.zeros(6),
        control=['strain', 'stress', 'stress', 'stress', 'stress', 'stress'],
        Dn_init=1,
        Dn_inc=n_increments,
        time=1.0
    )

    # Unloading step: strain back to zero
    step_unload = StepMeca(
        DEtot_end=np.array([-strain_max, 0, 0, 0, 0, 0]),
        Dsigma_end=np.zeros(6),
        control=['strain', 'stress', 'stress', 'stress', 'stress', 'stress'],
        Dn_init=1,
        Dn_inc=n_increments,
        time=1.0
    )

    block = Block(
        steps=[step_load, step_unload],
        umat_name=umat_name,
        props=props,
        nstatev=nstatev,
        control_type=control_type
    )

    solver = Solver(blocks=[block], max_iter=20, tol=1e-9)
    return solver.solve()


# =============================================================================
# ELISO Tests - Isotropic Elasticity
# =============================================================================

@pytest.mark.skipif(not _has_simcoon_core(), reason="simcoon._core not available")
class TestELISO:
    """Tests for ELISO (isotropic elastic) material model."""

    def test_uniaxial_tension(self):
        """Test uniaxial tension with ELISO."""
        E = 70000.0  # Young's modulus (MPa)
        nu = 0.3     # Poisson's ratio
        props = np.array([E, nu])

        history = run_uniaxial_test('ELISO', props, nstatev=1, strain_max=0.01)

        # At max strain (middle of history)
        mid_idx = len(history) // 2
        mid = history[mid_idx]

        # Check axial stress: sigma = E * epsilon for uniaxial
        expected_sigma = E * 0.01
        assert np.isclose(mid.sigma[0], expected_sigma, rtol=0.01), \
            f"Expected sigma_11={expected_sigma}, got {mid.sigma[0]}"

        # Check lateral strain: eps_lat = -nu * eps_axial
        expected_lat = -nu * 0.01
        assert np.isclose(mid.Etot[1], expected_lat, rtol=0.01), \
            f"Expected eps_22={expected_lat}, got {mid.Etot[1]}"

        # After unloading, stress should be near zero
        final = history[-1]
        assert np.allclose(final.sigma, 0, atol=1.0), \
            f"Expected zero stress after unload, got {final.sigma}"

    def test_pure_shear(self):
        """Test pure shear with ELISO."""
        E = 70000.0
        nu = 0.3
        G = E / (2 * (1 + nu))  # Shear modulus
        props = np.array([E, nu])

        # Shear strain
        gamma = 0.01
        step = StepMeca(
            DEtot_end=np.array([0, 0, 0, gamma, 0, 0]),
            Dsigma_end=np.zeros(6),
            control=['stress', 'stress', 'stress', 'strain', 'stress', 'stress'],
            Dn_init=1,
            Dn_inc=10
        )

        block = Block(
            steps=[step],
            umat_name='ELISO',
            props=props,
            nstatev=1
        )

        solver = Solver(blocks=[block])
        history = solver.solve()

        final = history[-1]
        # Shear stress = G * gamma (factor of 2 due to engineering shear strain)
        expected_tau = G * gamma
        assert np.isclose(final.sigma[3], expected_tau, rtol=0.01), \
            f"Expected tau_12={expected_tau}, got {final.sigma[3]}"

    def test_hydrostatic_compression(self):
        """Test hydrostatic compression with ELISO."""
        E = 70000.0
        nu = 0.3
        K = E / (3 * (1 - 2 * nu))  # Bulk modulus
        props = np.array([E, nu])

        # Apply equal strain in all directions
        eps = -0.001
        step = StepMeca(
            DEtot_end=np.array([eps, eps, eps, 0, 0, 0]),
            control=['strain'] * 6,
            Dn_init=1,
            Dn_inc=10
        )

        block = Block(
            steps=[step],
            umat_name='ELISO',
            props=props,
            nstatev=1
        )

        solver = Solver(blocks=[block])
        history = solver.solve()

        final = history[-1]
        # Hydrostatic pressure = K * volumetric_strain
        vol_strain = 3 * eps
        expected_p = K * vol_strain
        actual_p = np.mean(final.sigma[:3])
        assert np.isclose(actual_p, expected_p, rtol=0.02), \
            f"Expected pressure={expected_p}, got {actual_p}"


# =============================================================================
# ELIST Tests - Transversely Isotropic Elasticity
# =============================================================================

@pytest.mark.skipif(not _has_simcoon_core(), reason="simcoon._core not available")
class TestELIST:
    """Tests for ELIST (transversely isotropic elastic) material model."""

    def test_longitudinal_tension(self):
        """Test tension along fiber direction."""
        # Transversely isotropic properties
        # E_L, E_T, nu_TL, nu_TT, G_LT
        E_L = 150000.0   # Longitudinal modulus
        E_T = 10000.0    # Transverse modulus
        nu_TL = 0.3      # Poisson's ratio (transverse-longitudinal)
        nu_TT = 0.4      # Poisson's ratio (transverse-transverse)
        G_LT = 5000.0    # Shear modulus
        props = np.array([E_L, E_T, nu_TL, nu_TT, G_LT])

        # Tension along direction 1 (fiber direction)
        strain = 0.01
        step = StepMeca(
            DEtot_end=np.array([strain, 0, 0, 0, 0, 0]),
            Dsigma_end=np.zeros(6),
            control=['strain', 'stress', 'stress', 'stress', 'stress', 'stress'],
            Dn_init=1,
            Dn_inc=10
        )

        block = Block(
            steps=[step],
            umat_name='ELIST',
            props=props,
            nstatev=1
        )

        solver = Solver(blocks=[block])
        history = solver.solve()

        final = history[-1]
        # Stress should be approximately E_L * strain
        expected_sigma = E_L * strain
        assert np.isclose(final.sigma[0], expected_sigma, rtol=0.05), \
            f"Expected sigma_11={expected_sigma}, got {final.sigma[0]}"

    def test_transverse_tension(self):
        """Test tension perpendicular to fiber direction."""
        E_L = 150000.0
        E_T = 10000.0
        nu_TL = 0.3
        nu_TT = 0.4
        G_LT = 5000.0
        props = np.array([E_L, E_T, nu_TL, nu_TT, G_LT])

        # Tension along direction 2 (transverse)
        strain = 0.01
        step = StepMeca(
            DEtot_end=np.array([0, strain, 0, 0, 0, 0]),
            Dsigma_end=np.zeros(6),
            control=['stress', 'strain', 'stress', 'stress', 'stress', 'stress'],
            Dn_init=1,
            Dn_inc=10
        )

        block = Block(
            steps=[step],
            umat_name='ELIST',
            props=props,
            nstatev=1
        )

        solver = Solver(blocks=[block])
        history = solver.solve()

        final = history[-1]
        # Transverse stress response
        # The stress will be lower than E_T * strain due to coupling
        assert final.sigma[1] > 0, "Transverse stress should be positive"
        assert final.sigma[1] < E_T * strain * 1.5, "Stress seems too high"


# =============================================================================
# ELORT Tests - Orthotropic Elasticity
# =============================================================================

@pytest.mark.skipif(not _has_simcoon_core(), reason="simcoon._core not available")
class TestELORT:
    """Tests for ELORT (orthotropic elastic) material model."""

    def test_uniaxial_1(self):
        """Test uniaxial tension along direction 1."""
        # Orthotropic properties: E1, E2, E3, nu12, nu13, nu23, G12, G13, G23
        E1, E2, E3 = 150000.0, 10000.0, 10000.0
        nu12, nu13, nu23 = 0.3, 0.3, 0.4
        G12, G13, G23 = 5000.0, 5000.0, 3500.0
        props = np.array([E1, E2, E3, nu12, nu13, nu23, G12, G13, G23])

        strain = 0.01
        step = StepMeca(
            DEtot_end=np.array([strain, 0, 0, 0, 0, 0]),
            Dsigma_end=np.zeros(6),
            control=['strain', 'stress', 'stress', 'stress', 'stress', 'stress'],
            Dn_init=1,
            Dn_inc=10
        )

        block = Block(
            steps=[step],
            umat_name='ELORT',
            props=props,
            nstatev=1
        )

        solver = Solver(blocks=[block])
        history = solver.solve()

        final = history[-1]
        # Check approximate stress
        assert np.isclose(final.sigma[0], E1 * strain, rtol=0.1), \
            f"Sigma_11 should be approximately E1 * strain"


# =============================================================================
# EPICP Tests - Isotropic Elastoplasticity (Chaboche)
# =============================================================================

@pytest.mark.skipif(not _has_simcoon_core(), reason="simcoon._core not available")
class TestEPICP:
    """Tests for EPICP (isotropic plasticity with Chaboche hardening)."""

    def test_yield_detection(self):
        """Test that material yields at expected stress level."""
        # EPICP props: E, nu, alpha, sigmaY, k, m, kx[1-3], Dx[1-3]
        E = 70000.0
        nu = 0.3
        alpha = 1e-5  # thermal expansion (not used here)
        sigma_y = 300.0  # yield stress
        k = 0.0  # isotropic hardening
        m = 1.0  # hardening exponent
        kx1, kx2, kx3 = 0.0, 0.0, 0.0  # kinematic hardening
        Dx1, Dx2, Dx3 = 0.0, 0.0, 0.0  # dynamic recovery

        props = np.array([E, nu, alpha, sigma_y, k, m, kx1, Dx1, kx2, Dx2, kx3, Dx3])

        # Apply strain beyond yield
        strain = 0.02  # Well beyond yield strain (~sigma_y/E = 0.0043)
        step = StepMeca(
            DEtot_end=np.array([strain, 0, 0, 0, 0, 0]),
            Dsigma_end=np.zeros(6),
            control=['strain', 'stress', 'stress', 'stress', 'stress', 'stress'],
            Dn_init=1,
            Dn_inc=100
        )

        block = Block(
            steps=[step],
            umat_name='EPICP',
            props=props,
            nstatev=14  # EPICP uses 14 state variables
        )

        solver = Solver(blocks=[block], max_iter=20, tol=1e-8)
        history = solver.solve()

        final = history[-1]
        # Stress should be near yield stress (no hardening)
        assert np.isclose(final.sigma[0], sigma_y, rtol=0.05), \
            f"Expected stress near {sigma_y}, got {final.sigma[0]}"

    def test_plastic_loading_unloading(self):
        """Test plastic loading followed by elastic unloading."""
        E = 70000.0
        nu = 0.3
        sigma_y = 300.0

        props = np.array([E, nu, 0.0, sigma_y, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

        history = run_uniaxial_test('EPICP', props, nstatev=14,
                                    strain_max=0.02, n_increments=100)

        # After unloading, should have residual plastic strain
        final = history[-1]

        # Stress should be near zero
        assert np.isclose(final.sigma[0], 0, atol=10), \
            f"Expected near-zero stress after unload, got {final.sigma[0]}"

        # Should have residual strain (plastic deformation)
        elastic_strain = sigma_y / E
        expected_plastic = 0.02 - elastic_strain
        # Residual strain = applied strain - 2*elastic recovery
        residual = final.Etot[0]
        assert residual > 0, f"Should have positive residual strain, got {residual}"


# =============================================================================
# MIPLN Tests - Periodic Layers Homogenization
# =============================================================================

@pytest.mark.skipif(not _has_simcoon_core(), reason="simcoon._core not available")
class TestMIPLN:
    """Tests for MIPLN (periodic layers homogenization)."""

    def test_laminate_effective_stiffness(self, tmp_path, monkeypatch):
        """Test effective stiffness of 50/50 laminate using L_eff."""
        from simcoon.properties import effective_stiffness

        # Create 50/50 laminate with two isotropic materials
        layers = [
            Layer(
                number=0,
                umat_name='ELISO',
                concentration=0.5,
                props=np.array([70000, 0.3, 0]),  # E=70000, nu=0.3
                material_orientation=MaterialOrientation(0, 0, 0),
                geometry_orientation=GeometryOrientation(0, 90, -90)
            ),
            Layer(
                number=1,
                umat_name='ELISO',
                concentration=0.5,
                props=np.array([3000, 0.4, 0]),  # E=3000, nu=0.4
                material_orientation=MaterialOrientation(0, 0, 0),
                geometry_orientation=GeometryOrientation(0, 90, -90)
            ),
        ]

        # Save to JSON (L_eff reads from current working directory)
        save_layers_json(str(tmp_path / "layers0.json"), layers)

        # Change to temp directory so L_eff can find the JSON files
        monkeypatch.chdir(tmp_path)

        # Compute effective stiffness
        # MIPLN expects: nphases, method (0=Voigt, 1=Reuss, 2=Periodic)
        nphases = 2
        method = 2  # Periodic homogenization
        props = np.array([nphases, method])

        L_eff = effective_stiffness('MIPLN', props, nstatev=0)

        # Check that L_eff is a valid 6x6 stiffness matrix
        assert L_eff.shape == (6, 6), f"Expected 6x6 matrix, got {L_eff.shape}"

        # Stiffness should be positive definite (all eigenvalues > 0)
        eigenvalues = np.linalg.eigvals(L_eff)
        assert np.all(eigenvalues > 0), "Stiffness matrix should be positive definite"

        # Effective modulus should be between the two constituents
        # For laminate in fiber direction (direction 3 for layers)
        E_soft = 3000
        E_stiff = 70000
        # L_33 gives effective modulus in layer stacking direction
        assert L_eff[2, 2] > E_soft, "Effective stiffness should exceed soft phase"
        assert L_eff[2, 2] < E_stiff, "Effective stiffness should be less than stiff phase"


# =============================================================================
# MIMTN Tests - Mori-Tanaka Homogenization
# =============================================================================

@pytest.mark.skipif(not _has_simcoon_core(), reason="simcoon._core not available")
class TestMIMTN:
    """Tests for MIMTN (Mori-Tanaka homogenization)."""

    def test_spherical_inclusions_effective_stiffness(self, tmp_path, monkeypatch):
        """Test Mori-Tanaka with spherical inclusions."""
        from simcoon.properties import effective_stiffness

        # Create composite with spherical inclusions (30% volume fraction)
        ellipsoids = [
            Ellipsoid(
                number=0,
                coatingof=0,
                umat_name='ELISO',
                concentration=0.7,
                props=np.array([3000, 0.4, 0]),  # Matrix: E=3000, nu=0.4
                material_orientation=MaterialOrientation(0, 0, 0),
                geometry_orientation=GeometryOrientation(0, 0, 0),
                a1=1, a2=1, a3=1  # Sphere
            ),
            Ellipsoid(
                number=1,
                coatingof=0,
                umat_name='ELISO',
                concentration=0.3,
                props=np.array([70000, 0.3, 0]),  # Inclusions: E=70000, nu=0.3
                material_orientation=MaterialOrientation(0, 0, 0),
                geometry_orientation=GeometryOrientation(0, 0, 0),
                a1=1, a2=1, a3=1  # Sphere
            ),
        ]

        # Save to JSON (L_eff reads from current working directory)
        save_ellipsoids_json(str(tmp_path / "ellipsoids0.json"), ellipsoids)

        # Change to temp directory so L_eff can find the JSON files
        monkeypatch.chdir(tmp_path)

        # Compute effective stiffness
        # MIMTN expects: nphases
        nphases = 2
        props = np.array([nphases])

        L_eff = effective_stiffness('MIMTN', props, nstatev=0)

        # Check that L_eff is a valid 6x6 stiffness matrix
        assert L_eff.shape == (6, 6), f"Expected 6x6 matrix, got {L_eff.shape}"

        # For spherical inclusions, result should be isotropic
        # Check that diagonal terms are approximately equal
        diag = np.diag(L_eff)
        assert np.isclose(diag[0], diag[1], rtol=0.01), "Should be isotropic (L11=L22)"
        assert np.isclose(diag[0], diag[2], rtol=0.01), "Should be isotropic (L11=L33)"

        # Effective modulus should be between matrix and inclusion
        E_matrix = 3000
        E_inclusion = 70000
        # Extract effective E from L_eff (approximate for isotropic)
        # For isotropic: L11 = E(1-nu)/((1+nu)(1-2nu))
        assert L_eff[0, 0] > E_matrix, "Effective stiffness should exceed matrix"
        assert L_eff[0, 0] < E_inclusion, "Effective stiffness should be less than inclusion"


# =============================================================================
# Comparison Tests (against analytical solutions)
# =============================================================================

@pytest.mark.skipif(not _has_simcoon_core(), reason="simcoon._core not available")
class TestAnalyticalComparison:
    """Tests comparing results against analytical solutions."""

    def test_eliso_linear_response(self):
        """Compare ELISO results against analytical linear elastic solution."""
        E = 70000.0
        nu = 0.3
        props = np.array([E, nu])

        # Run test
        history = run_uniaxial_test('ELISO', props, nstatev=1,
                                    strain_max=0.02, n_increments=100)

        # Extract strain-stress curve
        strains = np.array([h.Etot[0] for h in history])
        stresses = np.array([h.sigma[0] for h in history])

        # Analytical solution for linear elastic uniaxial
        expected_stresses = E * strains

        # Check within tolerance
        max_error = np.max(np.abs(stresses - expected_stresses))
        assert max_error < 10.0, f"Max error {max_error} exceeds tolerance"


# =============================================================================
# Run Tests
# =============================================================================

if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])
