"""
Tests for the simcoon.solver module.

Tests the Python 0D solver implementation with various material models
and loading conditions.
"""

import pytest
import numpy as np
import numpy.typing as npt

from simcoon.solver import (
    StateVariables, StateVariablesM, StateVariablesT,
    Step, StepMeca, StepThermomeca,
    Block, Solver,
    Lt_2_K, Lth_2_K,
    CONTROL_TYPES, CORATE_TYPES,
    HistoryPoint,
)


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


# =============================================================================
# State Variables Tests
# =============================================================================

class TestStateVariables:
    """Tests for StateVariables classes."""

    def test_state_variables_default_init(self):
        """Test default initialization of StateVariables."""
        sv = StateVariables()
        assert sv.Etot.shape == (6,)
        assert sv.DEtot.shape == (6,)
        assert sv.sigma.shape == (6,)
        assert sv.F0.shape == (3, 3)
        assert sv.F1.shape == (3, 3)
        assert np.allclose(sv.F0, np.eye(3))
        assert np.allclose(sv.F1, np.eye(3))
        assert sv.T == 293.15
        assert sv.DT == 0.0

    def test_state_variables_m_default_init(self):
        """Test default initialization of StateVariablesM."""
        sv = StateVariablesM(nstatev=5)
        assert sv.Wm.shape == (4,)
        assert sv.Lt.shape == (6, 6)
        assert sv.L.shape == (6, 6)
        assert sv.statev.shape == (5,)

    def test_state_variables_t_default_init(self):
        """Test default initialization of StateVariablesT."""
        sv = StateVariablesT(nstatev=3)
        assert sv.Wm.shape == (4,)
        assert sv.Wt.shape == (4,)
        assert sv.dSdE.shape == (6, 6)
        assert sv.dSdT.shape == (6, 1)
        assert sv.drdE.shape == (1, 6)
        assert sv.Q == 0.0
        assert sv.r == 0.0

    def test_state_variables_copy(self):
        """Test deep copy of StateVariables."""
        sv1 = StateVariablesM(nstatev=2)
        sv1.Etot[0] = 0.01
        sv1.sigma[0] = 100.0

        sv2 = sv1.copy()
        sv2.Etot[0] = 0.02
        sv2.sigma[0] = 200.0

        # Original should be unchanged
        assert sv1.Etot[0] == 0.01
        assert sv1.sigma[0] == 100.0
        # Copy should have new values
        assert sv2.Etot[0] == 0.02
        assert sv2.sigma[0] == 200.0

    def test_to_start_and_set_start(self):
        """Test to_start and set_start methods (in-place operations)."""
        sv = StateVariablesM(nstatev=2)
        sv.sigma[:] = [100.0, 50.0, 50.0, 0.0, 0.0, 0.0]
        sv.statev[:] = [0.001, 0.002]

        # Save to start (in-place copy)
        sv.to_start()
        assert np.allclose(sv.sigma_start, sv.sigma)
        assert np.allclose(sv.statev_start, sv.statev)

        # Modify current (in-place)
        sv.sigma[:] = [200.0, 100.0, 100.0, 0.0, 0.0, 0.0]
        sv.statev[:] = [0.003, 0.004]

        # Restore from start (in-place copy)
        sv.set_start(1)
        assert np.allclose(sv.sigma, np.array([100.0, 50.0, 50.0, 0.0, 0.0, 0.0]))
        assert np.allclose(sv.statev, np.array([0.001, 0.002]))

    def test_copy_to_in_place(self):
        """Test copy_to method for in-place copying between objects."""
        sv1 = StateVariablesM(nstatev=2)
        sv1.Etot[:] = [0.01, -0.003, -0.003, 0, 0, 0]
        sv1.sigma[:] = [100.0, 50.0, 50.0, 0.0, 0.0, 0.0]
        sv1.Lt[:] = np.eye(6) * 200000.0

        sv2 = StateVariablesM(nstatev=2)
        sv1.copy_to(sv2)

        # Values should be equal
        assert np.allclose(sv2.Etot, sv1.Etot)
        assert np.allclose(sv2.sigma, sv1.sigma)
        assert np.allclose(sv2.Lt, sv1.Lt)

        # But arrays should be different objects (not aliased)
        assert sv1.Etot is not sv2.Etot
        assert sv1.sigma is not sv2.sigma
        assert sv1.Lt is not sv2.Lt

        # Modifying sv1 should not affect sv2
        sv1.Etot[0] = 0.02
        assert sv2.Etot[0] == 0.01


# =============================================================================
# Step Tests
# =============================================================================

class TestStep:
    """Tests for Step classes."""

    def test_step_default(self):
        """Test default Step initialization."""
        step = Step()
        assert step.Dn_init == 1
        assert step.Dn_mini == 1
        assert step.Dn_inc == 100
        assert step.control == ['strain'] * 6

    def test_step_meca(self):
        """Test StepMeca initialization."""
        step = StepMeca(
            DEtot_end=np.array([0.01, -0.003, -0.003, 0, 0, 0]),
            control=['strain', 'stress', 'stress', 'stress', 'stress', 'stress'],
            Dn_init=10
        )
        assert step.DEtot_end[0] == 0.01
        assert step.control[0] == 'strain'
        assert step.control[1] == 'stress'
        assert step.Dn_init == 10

    def test_get_cBC_meca(self):
        """Test cBC_meca array generation."""
        step = StepMeca(
            control=['strain', 'stress', 'stress', 'stress', 'stress', 'stress']
        )
        cBC = step.get_cBC_meca()
        assert np.array_equal(cBC, np.array([0, 1, 1, 1, 1, 1]))

    def test_step_thermomeca(self):
        """Test StepThermomeca initialization."""
        step = StepThermomeca(
            DT_end=50.0,
            thermal_control='temperature'
        )
        assert step.DT_end == 50.0
        assert step.get_cBC_T() == 0  # temperature controlled

        step2 = StepThermomeca(
            thermal_control='heat_flux'
        )
        assert step2.get_cBC_T() == 1  # heat flux controlled


# =============================================================================
# Block Tests
# =============================================================================

class TestBlock:
    """Tests for Block class."""

    def test_block_default(self):
        """Test default Block initialization."""
        block = Block()
        assert block.steps == []
        assert block.nstatev == 0
        assert block.umat_name == "ELISO"
        assert block.control_type == 'small_strain'
        assert block.corate_type == 'jaumann'

    def test_block_with_steps(self):
        """Test Block with steps."""
        step1 = StepMeca(DEtot_end=np.array([0.01, 0, 0, 0, 0, 0]))
        step2 = StepMeca(DEtot_end=np.array([0.02, 0, 0, 0, 0, 0]))

        block = Block(
            steps=[step1],
            umat_name="ELISO",
            props=np.array([210000.0, 0.3]),
            nstatev=1
        )
        block.add_step(step2)

        assert len(block.steps) == 2
        assert block.props[0] == 210000.0

    def test_control_type_int(self):
        """Test control type integer conversion."""
        block = Block(control_type='small_strain')
        assert block.get_control_type_int() == 1

        block2 = Block(control_type='logarithmic')
        assert block2.get_control_type_int() == 3

    def test_corate_type_int(self):
        """Test corate type integer conversion."""
        block = Block(corate_type='jaumann')
        assert block.get_corate_type_int() == 0

        block2 = Block(corate_type='green_naghdi')
        assert block2.get_corate_type_int() == 1


# =============================================================================
# Jacobian Helper Tests
# =============================================================================

class TestJacobianHelpers:
    """Tests for Jacobian helper functions."""

    def test_Lt_2_K_full_stress_control(self):
        """Test Lt_2_K with full stress control."""
        Lt = np.eye(6) * 200000.0  # Simple diagonal tangent
        cBC_meca = np.ones(6, dtype=int)
        lambda_solver = 10000.0

        K = Lt_2_K(Lt, cBC_meca, lambda_solver)
        assert np.allclose(K, Lt)

    def test_Lt_2_K_full_strain_control(self):
        """Test Lt_2_K with full strain control."""
        Lt = np.eye(6) * 200000.0
        cBC_meca = np.zeros(6, dtype=int)
        lambda_solver = 10000.0

        K = Lt_2_K(Lt, cBC_meca, lambda_solver)
        expected = np.eye(6) * lambda_solver
        assert np.allclose(K, expected)

    def test_Lt_2_K_mixed_control(self):
        """Test Lt_2_K with mixed control (uniaxial)."""
        Lt = np.array([
            [1.3461, 0.5769, 0.5769, 0, 0, 0],
            [0.5769, 1.3461, 0.5769, 0, 0, 0],
            [0.5769, 0.5769, 1.3461, 0, 0, 0],
            [0, 0, 0, 0.3846, 0, 0],
            [0, 0, 0, 0, 0.3846, 0],
            [0, 0, 0, 0, 0, 0.3846],
        ]) * 1e5
        cBC_meca = np.array([0, 1, 1, 1, 1, 1])  # Strain on 11, stress on rest
        lambda_solver = 10000.0

        K = Lt_2_K(Lt, cBC_meca, lambda_solver)

        # First row should be lambda on diagonal
        assert K[0, 0] == lambda_solver
        assert K[0, 1] == 0.0
        # Other rows should be from Lt
        assert np.allclose(K[1, :], Lt[1, :])
        assert np.allclose(K[2, :], Lt[2, :])

    def test_Lth_2_K(self):
        """Test Lth_2_K for thermomechanical Jacobian."""
        dSdE = np.eye(6) * 200000.0
        dSdT = np.ones((6, 1)) * -100.0
        dQdE = np.ones((1, 6)) * 50.0
        dQdT = 1000.0
        cBC_meca = np.array([0, 1, 1, 1, 1, 1])
        cBC_T = 1  # Heat flux controlled
        lambda_solver = 10000.0

        K = Lth_2_K(dSdE, dSdT, dQdE, dQdT, cBC_meca, cBC_T, lambda_solver)

        assert K.shape == (7, 7)
        # Check first row (strain controlled)
        assert K[0, 0] == lambda_solver
        # Check thermal row
        assert np.allclose(K[6, 0:6], dQdE.flatten())
        assert K[6, 6] == dQdT


# =============================================================================
# Solver Tests
# =============================================================================

class TestSolver:
    """Tests for Solver class."""

    def test_solver_initialization(self):
        """Test Solver initialization."""
        solver = Solver(max_iter=20, tol=1e-10)
        assert solver.max_iter == 20
        assert solver.tol == 1e-10
        assert solver.blocks == []
        assert solver.history == []

    @pytest.mark.skipif(
        not _has_simcoon_core(),
        reason="simcoon._core not available"
    )
    def test_solver_eliso_strain_controlled(self):
        """Test solver with ELISO under pure strain control."""
        # Material properties for ELISO (E, nu, alpha)
        E = 210000.0
        nu = 0.3
        alpha = 0.0  # CTE - no thermal expansion
        props = np.array([E, nu, alpha])

        # Pure strain-controlled step (all components)
        strain_11 = 0.001
        step = StepMeca(
            DEtot_end=np.array([strain_11, 0, 0, 0, 0, 0]),
            control=['strain'] * 6,
            Dn_init=1
        )

        block = Block(
            steps=[step],
            umat_name="ELISO",
            props=props,
            nstatev=1
        )

        solver = Solver(blocks=[block])
        history = solver.solve()

        # Check final strain
        final = history[-1]
        assert np.isclose(final.Etot[0], strain_11, rtol=1e-6)

        # Check stress (isotropic elasticity)
        # For pure strain control, sigma = L * eps
        lam = E * nu / ((1 + nu) * (1 - 2 * nu))
        mu = E / (2 * (1 + nu))
        expected_sigma_11 = (lam + 2 * mu) * strain_11
        assert np.isclose(final.sigma[0], expected_sigma_11, rtol=1e-3)

    @pytest.mark.skipif(
        not _has_simcoon_core(),
        reason="simcoon._core not available"
    )
    def test_solver_eliso_uniaxial(self):
        """Test solver with ELISO under uniaxial tension."""
        # Material properties for ELISO (E, nu, alpha)
        E = 210000.0
        nu = 0.3
        alpha = 0.0  # CTE - no thermal expansion
        props = np.array([E, nu, alpha])

        # Uniaxial tension: strain on 11, stress-free on others
        strain_11 = 0.01
        step = StepMeca(
            DEtot_end=np.array([strain_11, 0, 0, 0, 0, 0]),
            Dsigma_end=np.zeros(6),
            control=['strain', 'stress', 'stress', 'stress', 'stress', 'stress'],
            Dn_init=10,
            Dn_inc=100
        )

        block = Block(
            steps=[step],
            umat_name="ELISO",
            props=props,
            nstatev=1
        )

        solver = Solver(blocks=[block], max_iter=10, tol=1e-9)
        history = solver.solve()

        # Check final state
        final = history[-1]

        # Axial strain should match target
        assert np.isclose(final.Etot[0], strain_11, rtol=1e-6)

        # Axial stress should be E * strain for uniaxial
        expected_sigma = E * strain_11
        assert np.isclose(final.sigma[0], expected_sigma, rtol=1e-2)

        # Lateral stresses should be near zero (uniaxial)
        assert np.allclose(final.sigma[1:], 0.0, atol=1.0)

        # Lateral strains should be -nu * axial
        expected_lateral = -nu * strain_11
        assert np.isclose(final.Etot[1], expected_lateral, rtol=1e-2)
        assert np.isclose(final.Etot[2], expected_lateral, rtol=1e-2)


# =============================================================================
# HistoryPoint Tests
# =============================================================================

class TestHistoryPoint:
    """Tests for HistoryPoint class."""

    def test_history_point_from_state(self):
        """Test creating HistoryPoint from StateVariablesM."""
        sv = StateVariablesM(nstatev=3)
        sv.Etot[:] = [0.01, -0.003, -0.003, 0, 0, 0]
        sv.sigma[:] = [2100.0, 0, 0, 0, 0, 0]
        sv.Wm[:] = [10.5, 5.0, 3.0, 2.5]
        sv.statev[:] = [0.001, 0.002, 0.003]
        sv.T = 350.0

        hp = HistoryPoint.from_state(sv)

        np.testing.assert_array_equal(hp.Etot, sv.Etot)
        np.testing.assert_array_equal(hp.sigma, sv.sigma)
        np.testing.assert_array_equal(hp.Wm, sv.Wm)
        np.testing.assert_array_equal(hp.statev, sv.statev)
        assert hp.T == 350.0

    def test_history_point_independence(self):
        """Test that HistoryPoint is independent of source state."""
        sv = StateVariablesM(nstatev=2)
        sv.Etot[0] = 0.01
        sv.sigma[0] = 100.0

        hp = HistoryPoint.from_state(sv)

        # Modify original
        sv.Etot[0] = 0.02
        sv.sigma[0] = 200.0

        # HistoryPoint should be unchanged
        assert hp.Etot[0] == 0.01
        assert hp.sigma[0] == 100.0


# =============================================================================
# Control Type Mapping Tests
# =============================================================================

class TestControlMappings:
    """Tests for control type and corate mappings."""

    def test_control_types_dict(self):
        """Test CONTROL_TYPES mapping."""
        assert CONTROL_TYPES['small_strain'] == 1
        assert CONTROL_TYPES['green_lagrange'] == 2
        assert CONTROL_TYPES['logarithmic'] == 3

    def test_corate_types_dict(self):
        """Test CORATE_TYPES mapping."""
        assert CORATE_TYPES['jaumann'] == 0
        assert CORATE_TYPES['green_naghdi'] == 1
        assert CORATE_TYPES['logarithmic'] == 2


# =============================================================================
# Advanced Solver Tests
# =============================================================================

class TestSolverAdvanced:
    """Advanced tests for Solver class."""

    @pytest.mark.skipif(
        not _has_simcoon_core(),
        reason="simcoon._core not available"
    )
    def test_solver_multiple_steps(self):
        """Test solver with multiple steps in a block."""
        E = 210000.0
        nu = 0.3
        props = np.array([E, nu, 0.0])

        # Loading step
        step1 = StepMeca(
            DEtot_end=np.array([0.005, 0, 0, 0, 0, 0]),
            control=['strain'] * 6,
            Dn_init=5
        )
        # Unloading step
        step2 = StepMeca(
            DEtot_end=np.array([-0.005, 0, 0, 0, 0, 0]),
            control=['strain'] * 6,
            Dn_init=5
        )

        block = Block(
            steps=[step1, step2],
            umat_name="ELISO",
            props=props,
            nstatev=1
        )

        solver = Solver(blocks=[block])
        history = solver.solve()

        # Should have history points from both steps (at least 2)
        assert len(history) >= 2

        # Final strain should be back to zero
        final = history[-1]
        assert np.isclose(final.Etot[0], 0.0, atol=1e-10)

    @pytest.mark.skipif(
        not _has_simcoon_core(),
        reason="simcoon._core not available"
    )
    def test_solver_cyclic_loading(self):
        """Test solver with cyclic loading (ncycle > 1)."""
        E = 210000.0
        nu = 0.3
        props = np.array([E, nu, 0.0])

        step1 = StepMeca(
            DEtot_end=np.array([0.002, 0, 0, 0, 0, 0]),
            control=['strain'] * 6,
            Dn_init=2
        )
        step2 = StepMeca(
            DEtot_end=np.array([-0.002, 0, 0, 0, 0, 0]),
            control=['strain'] * 6,
            Dn_init=2
        )

        block = Block(
            steps=[step1, step2],
            umat_name="ELISO",
            props=props,
            nstatev=1,
            ncycle=3  # 3 cycles
        )

        solver = Solver(blocks=[block])
        history = solver.solve()

        # 3 cycles * 2 steps * 2 increments each + 1 initial = 13 points
        assert len(history) == 13

    @pytest.mark.skipif(
        not _has_simcoon_core(),
        reason="simcoon._core not available"
    )
    def test_solver_pure_shear(self):
        """Test solver under pure shear loading."""
        E = 210000.0
        nu = 0.3
        props = np.array([E, nu, 0.0])

        shear_strain = 0.01
        step = StepMeca(
            DEtot_end=np.array([0, 0, 0, shear_strain, 0, 0]),
            control=['strain'] * 6,
            Dn_init=5
        )

        block = Block(
            steps=[step],
            umat_name="ELISO",
            props=props,
            nstatev=1
        )

        solver = Solver(blocks=[block])
        history = solver.solve()

        final = history[-1]

        # Shear modulus G = E / (2 * (1 + nu))
        G = E / (2 * (1 + nu))

        # Shear stress = 2 * G * shear_strain (factor 2 for engineering strain)
        expected_shear_stress = G * shear_strain
        assert np.isclose(final.sigma[3], expected_shear_stress, rtol=1e-3)

    @pytest.mark.skipif(
        not _has_simcoon_core(),
        reason="simcoon._core not available"
    )
    def test_solver_with_initial_state(self):
        """Test solver with provided initial state."""
        E = 210000.0
        nu = 0.3
        props = np.array([E, nu, 0.0])

        # Create initial state with pre-strain
        sv_init = StateVariablesM(nstatev=1)
        sv_init.Etot[0] = 0.001  # Pre-existing strain

        step = StepMeca(
            DEtot_end=np.array([0.001, 0, 0, 0, 0, 0]),
            control=['strain'] * 6,
            Dn_init=5
        )

        block = Block(
            steps=[step],
            umat_name="ELISO",
            props=props,
            nstatev=1
        )

        solver = Solver(blocks=[block])
        history = solver.solve(sv_init=sv_init)

        # Final strain should include pre-strain
        final = history[-1]
        assert np.isclose(final.Etot[0], 0.002, rtol=1e-6)

    @pytest.mark.skipif(
        not _has_simcoon_core(),
        reason="simcoon._core not available"
    )
    def test_solver_hydrostatic_compression(self):
        """Test solver under hydrostatic compression."""
        E = 210000.0
        nu = 0.3
        props = np.array([E, nu, 0.0])

        vol_strain = -0.003  # compression
        step = StepMeca(
            DEtot_end=np.array([vol_strain, vol_strain, vol_strain, 0, 0, 0]),
            control=['strain'] * 6,
            Dn_init=5
        )

        block = Block(
            steps=[step],
            umat_name="ELISO",
            props=props,
            nstatev=1
        )

        solver = Solver(blocks=[block])
        history = solver.solve()

        final = history[-1]

        # Bulk modulus K = E / (3 * (1 - 2*nu))
        K = E / (3 * (1 - 2*nu))

        # Pressure = K * volumetric_strain (for small strain)
        volumetric_strain = 3 * vol_strain
        expected_pressure = K * volumetric_strain

        # All normal stresses should be equal (hydrostatic)
        assert np.isclose(final.sigma[0], final.sigma[1], rtol=1e-6)
        assert np.isclose(final.sigma[1], final.sigma[2], rtol=1e-6)

        # Mean stress should match pressure
        mean_stress = (final.sigma[0] + final.sigma[1] + final.sigma[2]) / 3
        assert np.isclose(mean_stress, expected_pressure, rtol=1e-2)


# =============================================================================
# Solver Parameter Tests
# =============================================================================

class TestSolverParameters:
    """Tests for Solver parameter handling."""

    def test_solver_custom_tolerance(self):
        """Test solver with custom tolerance."""
        solver = Solver(tol=1e-12)
        assert solver.tol == 1e-12

    def test_solver_custom_max_iter(self):
        """Test solver with custom max iterations."""
        solver = Solver(max_iter=50)
        assert solver.max_iter == 50

    def test_solver_custom_lambda(self):
        """Test solver with custom lambda_solver."""
        solver = Solver(lambda_solver=50000.0)
        assert solver.lambda_solver == 50000.0


# =============================================================================
# Run Tests
# =============================================================================

if __name__ == "__main__":
    pytest.main([__file__, "-v"])
