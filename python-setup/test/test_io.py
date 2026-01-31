"""
Tests for the simcoon.solver.io module.

Tests JSON I/O for materials, paths, and simulation configurations.
"""

import pytest
import json
import tempfile
import os
import numpy as np
from pathlib import Path

from simcoon.solver import (
    Block, StepMeca, StepThermomeca,
)
from simcoon.solver.io import (
    load_material_json,
    save_material_json,
    load_path_json,
    save_path_json,
    load_simulation_json,
    # Micromechanics (re-exported)
    Phase,
    Section,
    load_phases_json,
    save_phases_json,
    load_sections_json,
    save_sections_json,
)


# =============================================================================
# Material JSON Tests
# =============================================================================

class TestMaterialJSON:
    """Tests for material JSON I/O."""

    def test_save_and_load_material_dict_props(self):
        """Test saving and loading material with dict properties."""
        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = os.path.join(tmpdir, 'material.json')

            # Save with dict properties
            save_material_json(
                filepath,
                name='ELISO',
                props={'E': 70000.0, 'nu': 0.3, 'alpha': 1e-5},
                nstatev=1,
                psi=0, theta=0, phi=0
            )

            # Load and verify
            mat = load_material_json(filepath)
            assert mat['name'] == 'ELISO'
            assert mat['nstatev'] == 1
            assert len(mat['props']) == 3
            assert mat['orientation']['psi'] == 0.0

    def test_save_and_load_material_array_props(self):
        """Test saving and loading material with array properties."""
        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = os.path.join(tmpdir, 'material.json')

            # Save with array properties
            props = np.array([210000.0, 0.3, 1e-5])
            save_material_json(
                filepath,
                name='ELISO',
                props=props,
                nstatev=1,
                prop_names=['E', 'nu', 'alpha']
            )

            # Load and verify
            mat = load_material_json(filepath)
            assert mat['name'] == 'ELISO'
            np.testing.assert_array_almost_equal(mat['props'], props)

    def test_save_material_with_orientation(self):
        """Test material with non-zero orientation."""
        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = os.path.join(tmpdir, 'material.json')

            save_material_json(
                filepath,
                name='ELORT',
                props=np.array([1000, 500, 0.3, 0.2, 200]),
                nstatev=1,
                psi=45.0, theta=30.0, phi=60.0
            )

            mat = load_material_json(filepath)
            assert mat['orientation']['psi'] == 45.0
            assert mat['orientation']['theta'] == 30.0
            assert mat['orientation']['phi'] == 60.0

    def test_load_material_defaults(self):
        """Test loading material with missing optional fields."""
        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = os.path.join(tmpdir, 'material.json')

            # Write minimal JSON
            with open(filepath, 'w') as f:
                json.dump({'props': [70000, 0.3]}, f)

            mat = load_material_json(filepath)
            assert mat['name'] == 'ELISO'  # default
            assert mat['nstatev'] == 1  # default
            assert mat['orientation']['psi'] == 0.0  # default


# =============================================================================
# Path JSON Tests
# =============================================================================

class TestPathJSON:
    """Tests for path JSON I/O."""

    def test_save_and_load_path_mechanical(self):
        """Test saving and loading mechanical path."""
        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = os.path.join(tmpdir, 'path.json')

            # Create blocks
            step = StepMeca(
                DEtot_end=np.array([0.01, 0, 0, 0, 0, 0]),
                Dsigma_end=np.zeros(6),
                control=['strain', 'stress', 'stress', 'stress', 'stress', 'stress'],
                Dn_init=10,
                time=1.0
            )
            block = Block(
                steps=[step],
                umat_name='ELISO',
                props=np.array([70000, 0.3, 1e-5]),
                nstatev=1,
                control_type='small_strain'
            )

            # Save
            save_path_json(filepath, [block], initial_temperature=300.0)

            # Load and verify
            path = load_path_json(filepath)
            assert path['initial_temperature'] == 300.0
            assert len(path['blocks']) == 1

            loaded_block = path['blocks'][0]
            assert loaded_block.umat_name == 'ELISO'
            assert len(loaded_block.steps) == 1

    def test_save_and_load_path_thermomechanical(self):
        """Test saving and loading thermomechanical path."""
        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = os.path.join(tmpdir, 'path.json')

            step = StepThermomeca(
                DEtot_end=np.array([0.005, 0, 0, 0, 0, 0]),
                control=['strain', 'stress', 'stress', 'stress', 'stress', 'stress'],
                DT_end=50.0,
                thermal_control='temperature',
                Dn_init=20
            )
            block = Block(
                steps=[step],
                umat_name='ELISO',
                umat_type='thermomechanical',
                props=np.array([70000, 0.3, 1e-5]),
                nstatev=1
            )

            save_path_json(filepath, [block])
            path = load_path_json(filepath)

            assert len(path['blocks']) == 1
            loaded_step = path['blocks'][0].steps[0]
            assert isinstance(loaded_step, StepThermomeca)
            assert loaded_step.DT_end == 50.0

    def test_load_path_multiple_blocks(self):
        """Test loading path with multiple blocks."""
        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = os.path.join(tmpdir, 'path.json')

            # Write JSON directly
            path_data = {
                'initial_temperature': 293.15,
                'blocks': [
                    {
                        'type': 'mechanical',
                        'umat_name': 'ELISO',
                        'props': [70000, 0.3, 0],
                        'nstatev': 1,
                        'control_type': 'small_strain',
                        'steps': [
                            {
                                'DEtot': [0.01, 0, 0, 0, 0, 0],
                                'control': ['strain'] * 6,
                                'Dn_init': 10
                            }
                        ]
                    },
                    {
                        'type': 'mechanical',
                        'umat_name': 'ELISO',
                        'props': [70000, 0.3, 0],
                        'nstatev': 1,
                        'control_type': 'small_strain',
                        'steps': [
                            {
                                'DEtot': [-0.01, 0, 0, 0, 0, 0],
                                'control': ['strain'] * 6,
                                'Dn_init': 10
                            }
                        ]
                    }
                ]
            }

            with open(filepath, 'w') as f:
                json.dump(path_data, f)

            path = load_path_json(filepath)
            assert len(path['blocks']) == 2

    def test_load_path_cyclic(self):
        """Test loading path with cyclic loading (ncycle > 1)."""
        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = os.path.join(tmpdir, 'path.json')

            path_data = {
                'blocks': [
                    {
                        'type': 'mechanical',
                        'umat_name': 'EPICP',
                        'props': [70000, 0.3, 0, 300, 1000, 0.5],
                        'nstatev': 8,
                        'ncycle': 5,
                        'steps': [
                            {'DEtot': [0.01, 0, 0, 0, 0, 0], 'control': ['strain'] * 6},
                            {'DEtot': [-0.01, 0, 0, 0, 0, 0], 'control': ['strain'] * 6},
                        ]
                    }
                ]
            }

            with open(filepath, 'w') as f:
                json.dump(path_data, f)

            path = load_path_json(filepath)
            assert path['blocks'][0].ncycle == 5
            assert len(path['blocks'][0].steps) == 2


# =============================================================================
# Simulation JSON Tests
# =============================================================================

class TestSimulationJSON:
    """Tests for combined simulation JSON I/O."""

    def test_load_simulation_json(self):
        """Test loading complete simulation configuration."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create material file
            mat_file = os.path.join(tmpdir, 'material.json')
            save_material_json(mat_file, 'ELISO', {'E': 70000, 'nu': 0.3, 'alpha': 0})

            # Create path file
            path_file = os.path.join(tmpdir, 'path.json')
            step = StepMeca(
                DEtot_end=np.array([0.01, 0, 0, 0, 0, 0]),
                control=['strain', 'stress', 'stress', 'stress', 'stress', 'stress']
            )
            block = Block(steps=[step], umat_name='ELISO', nstatev=1)
            save_path_json(path_file, [block])

            # Load combined
            sim = load_simulation_json(mat_file, path_file)

            assert 'material' in sim
            assert 'blocks' in sim
            assert sim['material']['name'] == 'ELISO'
            assert len(sim['blocks']) == 1


# =============================================================================
# Phase JSON Tests
# =============================================================================

class TestPhaseJSON:
    """Tests for Phase JSON I/O."""

    def test_phase_creation(self):
        """Test Phase dataclass creation."""
        phase = Phase(
            number=0,
            umat_name='ELISO',
            props=np.array([70000, 0.3]),
            concentration=0.6
        )
        assert phase.number == 0
        assert phase.concentration == 0.6

    def test_save_and_load_phases(self):
        """Test saving and loading phases."""
        phases = [
            Phase(number=0, umat_name='ELISO', props=np.array([70000, 0.3]),
                  concentration=0.6, nstatev=1),
            Phase(number=1, umat_name='ELISO', props=np.array([400000, 0.2]),
                  concentration=0.4, nstatev=1),
        ]

        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = os.path.join(tmpdir, 'phases.json')
            save_phases_json(filepath, phases)

            loaded = load_phases_json(filepath)
            assert len(loaded) == 2
            assert loaded[0].concentration == 0.6
            assert loaded[1].concentration == 0.4


# =============================================================================
# Section JSON Tests
# =============================================================================

class TestSectionJSON:
    """Tests for Section JSON I/O."""

    def test_section_creation(self):
        """Test Section dataclass creation."""
        section = Section(
            number=0,
            name='yarn_0',
            umat_name='ELISO',
            props=np.array([70000, 0.3]),
            nstatev=1
        )
        assert section.number == 0
        assert section.name == 'yarn_0'

    def test_save_and_load_sections(self):
        """Test saving and loading sections."""
        sections = [
            Section(number=0, name='yarn_0', umat_name='ELISO',
                    props=np.array([70000, 0.3]), nstatev=1),
            Section(number=1, name='yarn_1', umat_name='ELISO',
                    props=np.array([400000, 0.2]), nstatev=1),
        ]

        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = os.path.join(tmpdir, 'sections.json')
            save_sections_json(filepath, sections)

            loaded = load_sections_json(filepath)
            assert len(loaded) == 2
            assert loaded[0].name == 'yarn_0'
            assert loaded[1].name == 'yarn_1'


# =============================================================================
# Error Handling Tests
# =============================================================================

class TestIOErrors:
    """Tests for I/O error handling."""

    def test_load_nonexistent_file(self):
        """Test loading non-existent file raises error."""
        with pytest.raises(FileNotFoundError):
            load_material_json('/nonexistent/path/material.json')

    def test_load_invalid_json(self):
        """Test loading invalid JSON raises error."""
        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = os.path.join(tmpdir, 'invalid.json')
            with open(filepath, 'w') as f:
                f.write('not valid json {{{')

            with pytest.raises(json.JSONDecodeError):
                load_material_json(filepath)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
