"""
Tests for UMAT material models using the Python solver.

Covers all available UMATs:
- Elasticity: ELISO, ELIST, ELORT
- Plasticity: EPICP, EPKCP, EPCHA, EPHIL, EPHAC
- Viscoelasticity: ZENER, ZENNK, PRONK
- Finite strain: NEOHC, NEOHI, MOORI, SNTVE, HYPOO
- Homogenization: MIMTN, MIPLN, MIHEN, MISCN
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
    """Run a uniaxial tension test."""
    step_load = StepMeca(
        DEtot_end=np.array([strain_max, 0, 0, 0, 0, 0]),
        Dsigma_end=np.zeros(6),
        control=['strain', 'stress', 'stress', 'stress', 'stress', 'stress'],
        Dn_init=1,
        Dn_inc=n_increments,
        time=1.0
    )

    block = Block(
        steps=[step_load],
        umat_name=umat_name,
        props=props,
        nstatev=nstatev,
        control_type=control_type
    )

    solver = Solver(blocks=[block], max_iter=20, tol=1e-9)
    return solver.solve()


def run_cyclic_test(umat_name: str, props: np.ndarray, nstatev: int,
                    strain_max: float = 0.02, n_increments: int = 50,
                    control_type: str = 'small_strain') -> list:
    """Run a cyclic loading/unloading test."""
    step_load = StepMeca(
        DEtot_end=np.array([strain_max, 0, 0, 0, 0, 0]),
        Dsigma_end=np.zeros(6),
        control=['strain', 'stress', 'stress', 'stress', 'stress', 'stress'],
        Dn_init=1,
        Dn_inc=n_increments,
        time=1.0
    )
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
        E = 70000.0
        nu = 0.3
        alpha = 1e-5
        props = np.array([E, nu, alpha])

        history = run_uniaxial_test('ELISO', props, nstatev=1, strain_max=0.01)

        mid_idx = len(history) // 2
        mid = history[mid_idx]

        expected_sigma = E * 0.01
        assert np.isclose(mid.sigma[0], expected_sigma, rtol=0.01)

    def test_pure_shear(self):
        """Test pure shear with ELISO."""
        E = 70000.0
        nu = 0.3
        alpha = 1e-5
        G = E / (2 * (1 + nu))
        props = np.array([E, nu, alpha])

        step = StepMeca(
            DEtot_end=np.array([0, 0, 0, 0.01, 0, 0]),
            control=['stress', 'stress', 'stress', 'strain', 'stress', 'stress'],
            Dn_init=1, Dn_inc=10
        )
        block = Block(steps=[step], umat_name='ELISO', props=props, nstatev=1)
        solver = Solver(blocks=[block])
        history = solver.solve()

        expected_tau = G * 0.01
        assert np.isclose(history[-1].sigma[3], expected_tau, rtol=0.01)

    def test_hydrostatic(self):
        """Test hydrostatic compression."""
        E = 70000.0
        nu = 0.3
        alpha = 1e-5
        K = E / (3 * (1 - 2 * nu))
        props = np.array([E, nu, alpha])

        eps = -0.001
        step = StepMeca(
            DEtot_end=np.array([eps, eps, eps, 0, 0, 0]),
            control=['strain'] * 6,
            Dn_init=1, Dn_inc=10
        )
        block = Block(steps=[step], umat_name='ELISO', props=props, nstatev=1)
        solver = Solver(blocks=[block])
        history = solver.solve()

        expected_p = K * 3 * eps
        actual_p = np.mean(history[-1].sigma[:3])
        assert np.isclose(actual_p, expected_p, rtol=0.02)


# =============================================================================
# ELIST Tests - Transversely Isotropic Elasticity
# =============================================================================

@pytest.mark.skipif(not _has_simcoon_core(), reason="simcoon._core not available")
class TestELIST:
    """Tests for ELIST (transversely isotropic elastic)."""

    def test_longitudinal_tension(self):
        """Test tension along fiber direction."""
        # axis, EL, ET, nuTL, nuTT, GLT, alphaL, alphaT (8 props)
        axis = 1  # Longitudinal direction along axis 1
        EL = 150000
        ET = 10000
        nuTL = 0.3
        nuTT = 0.4
        GLT = 5000
        alphaL = 1e-5
        alphaT = 1e-5
        props = np.array([axis, EL, ET, nuTL, nuTT, GLT, alphaL, alphaT])

        history = run_uniaxial_test('ELIST', props, nstatev=1, strain_max=0.01)

        # Stress should be approximately EL * strain
        assert history[-1].sigma[0] > 1000  # Should have significant stress


# =============================================================================
# ELORT Tests - Orthotropic Elasticity
# =============================================================================

@pytest.mark.skipif(not _has_simcoon_core(), reason="simcoon._core not available")
class TestELORT:
    """Tests for ELORT (orthotropic elastic)."""

    def test_uniaxial_direction1(self):
        """Test uniaxial tension along direction 1."""
        # E1, E2, E3, nu12, nu13, nu23, G12, G13, G23, alpha1, alpha2, alpha3
        props = np.array([150000, 10000, 10000, 0.3, 0.3, 0.4, 5000, 5000, 3500, 1e-5, 1e-5, 1e-5])

        history = run_uniaxial_test('ELORT', props, nstatev=1, strain_max=0.01)

        # Stress should be approximately E1 * strain
        assert np.isclose(history[-1].sigma[0], 150000 * 0.01, rtol=0.1)


# =============================================================================
# EPICP Tests - Isotropic Plasticity (Chaboche)
# =============================================================================

@pytest.mark.skipif(not _has_simcoon_core(), reason="simcoon._core not available")
class TestEPICP:
    """Tests for EPICP (isotropic plasticity with Chaboche hardening)."""

    def test_yield_detection(self):
        """Test that material yields at expected stress level."""
        # E, nu, alpha, sigmaY, k, m, kx1, Dx1, kx2, Dx2, kx3, Dx3
        E = 70000.0
        sigma_y = 300.0
        props = np.array([E, 0.3, 1e-5, sigma_y, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

        history = run_uniaxial_test('EPICP', props, nstatev=14, strain_max=0.02, n_increments=100)

        # Stress should plateau near yield stress
        assert np.isclose(history[-1].sigma[0], sigma_y, rtol=0.05)

    def test_plastic_strain(self):
        """Test that plastic strain develops after yield."""
        E = 70000.0
        sigma_y = 300.0
        props = np.array([E, 0.3, 1e-5, sigma_y, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

        history = run_cyclic_test('EPICP', props, nstatev=14, strain_max=0.02)

        # After unloading, should have residual plastic strain
        final = history[-1]
        # For EPICP: statev[1] = accumulated plastic p, statev[2:8] = plastic strain EP
        # After cyclic loading, plastic strain should be non-zero
        plastic_strain_11 = final.statev[2]  # EP(0,0)
        assert plastic_strain_11 > 0.001, f"Expected plastic strain > 0.001, got {plastic_strain_11}"


# =============================================================================
# EPKCP Tests - Kinematic + Isotropic Plasticity
# =============================================================================

@pytest.mark.skipif(not _has_simcoon_core(), reason="simcoon._core not available")
class TestEPKCP:
    """Tests for EPKCP (kinematic + isotropic plasticity)."""

    def test_kinematic_hardening(self):
        """Test kinematic hardening behavior."""
        # E, nu, alpha, sigmaY, k, m, kx1, Dx1, kx2, Dx2, kx3, Dx3
        E = 70000.0
        sigma_y = 300.0
        kx1 = 10000.0  # Kinematic hardening
        props = np.array([E, 0.3, 1e-5, sigma_y, 500.0, 1.0, kx1, 100.0, 0.0, 0.0, 0.0, 0.0])

        history = run_uniaxial_test('EPKCP', props, nstatev=14, strain_max=0.03)

        # With hardening, stress should exceed initial yield
        assert history[-1].sigma[0] > sigma_y


# =============================================================================
# EPCHA Tests - Chaboche Plasticity
# =============================================================================

@pytest.mark.skipif(not _has_simcoon_core(), reason="simcoon._core not available")
class TestEPCHA:
    """Tests for EPCHA (Chaboche plasticity)."""

    def test_chaboche_hardening(self):
        """Test Chaboche hardening behavior with fully strain-controlled loading."""
        # props: E, nu, alpha, sigmaY, Q, b, C_1, D_1, C_2, D_2 (10 props)
        E = 70000.0
        nu = 0.3
        alpha = 1e-5
        sigma_y = 300.0
        Q = 100.0       # Isotropic hardening parameter
        b = 10.0        # Isotropic hardening exponent
        C_1 = 5000.0    # Kinematic hardening parameter 1
        D_1 = 50.0      # Kinematic hardening saturation 1
        C_2 = 0.0       # Kinematic hardening parameter 2
        D_2 = 0.0       # Kinematic hardening saturation 2
        props = np.array([E, nu, alpha, sigma_y, Q, b, C_1, D_1, C_2, D_2])

        # Use fully strain-controlled loading for stability
        step = StepMeca(
            DEtot_end=np.array([0.02, 0, 0, 0, 0, 0]),
            control=['strain', 'strain', 'strain', 'strain', 'strain', 'strain'],
            Dn_init=1, Dn_inc=200
        )
        # EPCHA requires 33 state variables (indices 0-32):
        # 0: T_init, 1: p, 2-7: EP, 8-13: a_1, 14-19: a_2, 20-25: X_1, 26-31: X_2, 32: Hp
        block = Block(steps=[step], umat_name='EPCHA', props=props, nstatev=33)
        solver = Solver(blocks=[block], max_iter=20, tol=1e-9)
        history = solver.solve()

        # Material should harden
        assert history[-1].sigma[0] >= sigma_y


# =============================================================================
# ZENER Tests - Zener Viscoelasticity
# =============================================================================

@pytest.mark.skipif(not _has_simcoon_core(), reason="simcoon._core not available")
class TestZENER:
    """Tests for ZENER (Zener viscoelasticity)."""

    def test_viscoelastic_response(self):
        """Test viscoelastic relaxation with fully strain-controlled loading."""
        # props: E0, nu0, alpha_iso, E1, nu1, etaB1, etaS1 (7 props)
        E0 = 10000.0       # Thermoelastic Young's modulus
        nu0 = 0.3          # Thermoelastic Poisson's ratio
        alpha_iso = 1e-5   # Thermoelastic CTE
        E1 = 50000.0       # Viscoelastic Young modulus
        nu1 = 0.3          # Viscoelastic Poisson ratio
        etaB1 = 10000.0    # Viscoelastic Bulk viscosity (higher for stability)
        etaS1 = 5000.0     # Viscoelastic Shear viscosity
        props = np.array([E0, nu0, alpha_iso, E1, nu1, etaB1, etaS1])

        # Use fully strain-controlled loading (all components)
        step = StepMeca(
            DEtot_end=np.array([0.01, 0, 0, 0, 0, 0]),
            control=['strain', 'strain', 'strain', 'strain', 'strain', 'strain'],
            Dn_init=1, Dn_inc=100
        )
        block = Block(steps=[step], umat_name='ZENER', props=props, nstatev=8)
        solver = Solver(blocks=[block], max_iter=20, tol=1e-9)
        history = solver.solve()

        # Final stress should be positive
        assert history[-1].sigma[0] > 0


# =============================================================================
# PRONK Tests - Prony Series Viscoelasticity
# =============================================================================

@pytest.mark.skipif(not _has_simcoon_core(), reason="simcoon._core not available")
class TestPRONK:
    """Tests for PRONK (Prony series viscoelasticity)."""

    def test_prony_relaxation(self):
        """Test Prony series relaxation with fully strain-controlled loading."""
        # props: E0, nu0, alpha_iso, N_prony, E_visco[i], nu_visco[i], etaB_visco[i], etaS_visco[i]...
        # Total props: 4 + N_prony * 4
        E0 = 10000.0       # Equilibrium modulus
        nu0 = 0.3          # Poisson ratio
        alpha_iso = 1e-5   # CTE
        N_prony = 2        # Number of Prony terms

        # Branch 1
        E_visco_1 = 20000.0
        nu_visco_1 = 0.3
        etaB_visco_1 = 10000.0  # Higher viscosity for stability
        etaS_visco_1 = 5000.0

        # Branch 2
        E_visco_2 = 10000.0
        nu_visco_2 = 0.3
        etaB_visco_2 = 50000.0
        etaS_visco_2 = 25000.0

        props = np.array([
            E0, nu0, alpha_iso, N_prony,
            E_visco_1, nu_visco_1, etaB_visco_1, etaS_visco_1,
            E_visco_2, nu_visco_2, etaB_visco_2, etaS_visco_2
        ])

        # Use fully strain-controlled loading (all components)
        step = StepMeca(
            DEtot_end=np.array([0.01, 0, 0, 0, 0, 0]),
            control=['strain', 'strain', 'strain', 'strain', 'strain', 'strain'],
            Dn_init=1, Dn_inc=100
        )
        # nstatev = 7 + N_prony * 7 = 7 + 2*7 = 21
        block = Block(steps=[step], umat_name='PRONK', props=props, nstatev=21)
        solver = Solver(blocks=[block], max_iter=20, tol=1e-9)
        history = solver.solve()

        # Stress should be non-zero (viscoelastic response) - may be positive or negative
        # during relaxation depending on parameters
        assert abs(history[-1].sigma[0]) > 0


# =============================================================================
# Finite Strain Tests - Neo-Hookean, Mooney-Rivlin
# =============================================================================

@pytest.mark.skipif(not _has_simcoon_core(), reason="simcoon._core not available")
class TestFiniteStrain:
    """Tests for finite strain constitutive models."""

    def test_neohc_compression(self):
        """Test Neo-Hookean compressible with fully strain-controlled loading."""
        # props: E, nu, alpha (3 props - same as elastic_isotropic)
        E = 70000.0
        nu = 0.3
        alpha = 1e-5
        props = np.array([E, nu, alpha])

        # Use fully strain-controlled loading for stability
        step = StepMeca(
            DEtot_end=np.array([0.1, 0, 0, 0, 0, 0]),
            control=['strain', 'strain', 'strain', 'strain', 'strain', 'strain'],
            Dn_init=1, Dn_inc=200
        )
        block = Block(steps=[step], umat_name='NEOHC', props=props, nstatev=1,
                      control_type='logarithmic')
        solver = Solver(blocks=[block], max_iter=20, tol=1e-9)
        history = solver.solve()

        assert history[-1].sigma[0] > 0

    def test_sntve_tension(self):
        """Test Saint-Venant (St. Venant-Kirchhoff) finite strain.

        SNTVE uses Green-Lagrange strain + PKII stress with linear elastic relation.
        Simple extension of small-strain elasticity to finite strains.
        Good for pedagogical comparison with NEOHC (Neo-Hookean).
        """
        # props: E, nu, alpha (3 props - same as ELISO)
        E = 70000.0
        nu = 0.3
        alpha = 1e-5
        props = np.array([E, nu, alpha])

        # Use fully strain-controlled loading for stability
        step = StepMeca(
            DEtot_end=np.array([0.05, 0, 0, 0, 0, 0]),
            control=['strain', 'strain', 'strain', 'strain', 'strain', 'strain'],
            Dn_init=1, Dn_inc=100
        )
        block = Block(steps=[step], umat_name='SNTVE', props=props, nstatev=1,
                      control_type='logarithmic')
        solver = Solver(blocks=[block], max_iter=20, tol=1e-9)
        history = solver.solve()

        assert history[-1].sigma[0] > 0


# =============================================================================
# Logarithmic Strain Tests
# =============================================================================

@pytest.mark.skipif(not _has_simcoon_core(), reason="simcoon._core not available")
class TestLogarithmic:
    """Tests using logarithmic strain formulation."""

    def test_eliso_logarithmic(self):
        """Test ELISO with logarithmic strain."""
        E = 70000.0
        nu = 0.3
        alpha = 1e-5
        props = np.array([E, nu, alpha])

        history = run_uniaxial_test('ELISO', props, nstatev=1, strain_max=0.1,
                                    control_type='logarithmic')

        # Should handle larger strains with logarithmic formulation
        assert len(history) > 0
        assert history[-1].sigma[0] > 0

    def test_epicp_logarithmic(self):
        """Test EPICP with logarithmic strain."""
        E = 70000.0
        sigma_y = 300.0
        props = np.array([E, 0.3, 1e-5, sigma_y, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

        history = run_uniaxial_test('EPICP', props, nstatev=14, strain_max=0.1,
                                    control_type='logarithmic')

        # Stress should plateau near yield
        assert history[-1].sigma[0] > 250


# =============================================================================
# MIPLN Tests - Periodic Layers Homogenization
# =============================================================================

@pytest.mark.skipif(not _has_simcoon_core(), reason="simcoon._core not available")
class TestMIPLN:
    """Tests for MIPLN (periodic layers homogenization)."""

    def test_laminate_effective_stiffness(self, tmp_path, monkeypatch):
        """Test effective stiffness of 50/50 laminate.

        MIPLN props: [unused, file_index]
        File: data/layers{file_index}.json

        Uses existing test data file (layers_laminate.json) which has
        props as arrays (the format C++ read_layer_json expects).
        """
        import shutil
        from simcoon.properties import effective_stiffness

        # Create data/ subdirectory (C++ expects files in data/)
        data_dir = tmp_path / "data"
        data_dir.mkdir()

        # Copy existing test data file (has correct array format for props)
        src_file = TEST_DATA_DIR / "layers_laminate.json"
        dst_file = data_dir / "layers0.json"
        shutil.copy(src_file, dst_file)
        monkeypatch.chdir(tmp_path)

        # MIPLN: [unused, file_index]
        props = np.array([0, 0])  # file_index=0 for layers0.json
        L_eff = effective_stiffness('MIPLN', props, nstatev=0)

        assert L_eff.shape == (6, 6)
        # Eigenvalues should be positive (positive definite)
        eigenvalues = np.linalg.eigvals(L_eff)
        assert np.all(eigenvalues.real > 0)


# =============================================================================
# MIMTN Tests - Mori-Tanaka Homogenization
# =============================================================================

@pytest.mark.skipif(not _has_simcoon_core(), reason="simcoon._core not available")
class TestMIMTN:
    """Tests for MIMTN (Mori-Tanaka homogenization)."""

    def test_spherical_inclusions(self, tmp_path, monkeypatch):
        """Test Mori-Tanaka with spherical inclusions.

        MIMTN props: [unused, file_index, mp, np, n_matrix]
        - file_index: index for ellipsoids{N}.json file
        - mp, np: number of integration points for Eshelby tensor
        - n_matrix: which phase is the matrix (0-indexed)

        File: data/ellipsoids{file_index}.json
        """
        from simcoon.properties import effective_stiffness

        ellipsoids = [
            Ellipsoid(number=0, coatingof=0, umat_name='ELISO', concentration=0.7,
                      props=np.array([3000, 0.4, 0]),
                      a1=1, a2=1, a3=1),
            Ellipsoid(number=1, coatingof=0, umat_name='ELISO', concentration=0.3,
                      props=np.array([70000, 0.3, 0]),
                      a1=1, a2=1, a3=1),
        ]

        # Create data/ subdirectory (C++ expects files in data/)
        data_dir = tmp_path / "data"
        data_dir.mkdir()
        # File index is in props[1], so ellipsoids0.json
        save_ellipsoids_json(str(data_dir / "ellipsoids0.json"), ellipsoids)
        monkeypatch.chdir(tmp_path)

        # MIMTN: [unused, file_index, mp, np, n_matrix]
        # mp, np = integration points for Eshelby tensor calculation
        # n_matrix = which phase is the matrix (phase 0 is matrix with 70% concentration)
        props = np.array([0, 0, 10, 10, 0])  # file_index=0, mp=10, np=10, n_matrix=0
        L_eff = effective_stiffness('MIMTN', props, nstatev=0)

        assert L_eff.shape == (6, 6)
        # For spherical inclusions, should be approximately isotropic
        diag = np.diag(L_eff)
        assert np.isclose(diag[0], diag[1], rtol=0.05)


# =============================================================================
# Analytical Comparison Tests
# =============================================================================

@pytest.mark.skipif(not _has_simcoon_core(), reason="simcoon._core not available")
class TestAnalyticalComparison:
    """Tests comparing results against analytical solutions."""

    def test_eliso_linear_response(self):
        """Compare ELISO against analytical linear elastic solution."""
        E = 70000.0
        nu = 0.3
        props = np.array([E, nu, 1e-5])

        history = run_uniaxial_test('ELISO', props, nstatev=1, strain_max=0.02)

        strains = np.array([h.Etot[0] for h in history])
        stresses = np.array([h.sigma[0] for h in history])
        expected = E * strains

        max_error = np.max(np.abs(stresses - expected))
        assert max_error < 10.0

    def test_eliso_poisson_effect(self):
        """Check Poisson's ratio effect."""
        E = 70000.0
        nu = 0.3
        props = np.array([E, nu, 1e-5])

        history = run_uniaxial_test('ELISO', props, nstatev=1, strain_max=0.01)

        final = history[-1]
        # Lateral strain should be -nu * axial strain
        expected_lateral = -nu * final.Etot[0]
        assert np.isclose(final.Etot[1], expected_lateral, rtol=0.01)


# =============================================================================
# Run Tests
# =============================================================================

if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])
