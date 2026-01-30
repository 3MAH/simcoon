"""
Tests for the simcoon.solver.micromechanics module.

Tests the JSON-based I/O for phase configurations (ellipsoids, layers, cylinders).
"""

import pytest
import json
import tempfile
import os
import numpy as np
from pathlib import Path

from simcoon.solver.micromechanics import (
    MaterialOrientation,
    GeometryOrientation,
    Phase,
    Layer,
    Ellipsoid,
    Cylinder,
    load_layers_json,
    save_layers_json,
    load_ellipsoids_json,
    save_ellipsoids_json,
    load_cylinders_json,
    save_cylinders_json,
)


# =============================================================================
# MaterialOrientation Tests
# =============================================================================

class TestMaterialOrientation:
    """Tests for MaterialOrientation dataclass."""

    def test_default_values(self):
        """Test default orientation is identity (no rotation)."""
        orient = MaterialOrientation()
        assert orient.psi == 0.0
        assert orient.theta == 0.0
        assert orient.phi == 0.0

    def test_custom_values(self):
        """Test custom Euler angles."""
        orient = MaterialOrientation(psi=45.0, theta=30.0, phi=60.0)
        assert orient.psi == 45.0
        assert orient.theta == 30.0
        assert orient.phi == 60.0


# =============================================================================
# GeometryOrientation Tests
# =============================================================================

class TestGeometryOrientation:
    """Tests for GeometryOrientation dataclass."""

    def test_default_values(self):
        """Test default orientation."""
        orient = GeometryOrientation()
        assert orient.psi == 0.0
        assert orient.theta == 0.0
        assert orient.phi == 0.0

    def test_layer_orientation(self):
        """Test typical layer orientation (horizontal)."""
        orient = GeometryOrientation(psi=0, theta=90, phi=-90)
        assert orient.theta == 90.0
        assert orient.phi == -90.0


# =============================================================================
# Layer Tests
# =============================================================================

class TestLayer:
    """Tests for Layer dataclass."""

    def test_default_layer(self):
        """Test default layer initialization."""
        layer = Layer()
        assert layer.number == 0
        assert layer.umat_name == 'ELISO'
        assert layer.concentration == 1.0
        assert layer.save == 1

    def test_layer_with_props_array(self):
        """Test layer with properties as numpy array."""
        layer = Layer(
            number=0,
            umat_name='ELISO',
            concentration=0.5,
            props=np.array([70000, 0.3])
        )
        assert layer.concentration == 0.5


# =============================================================================
# Ellipsoid Tests
# =============================================================================

class TestEllipsoid:
    """Tests for Ellipsoid dataclass."""

    def test_default_ellipsoid(self):
        """Test default ellipsoid (sphere)."""
        ell = Ellipsoid()
        assert ell.a1 == 1.0
        assert ell.a2 == 1.0
        assert ell.a3 == 1.0

    def test_fiber_ellipsoid(self):
        """Test fiber-like ellipsoid (high aspect ratio)."""
        ell = Ellipsoid(
            number=1,
            concentration=0.3,
            a1=50, a2=1, a3=1,
            props=np.array([400000, 0.2])
        )
        assert ell.a1 / ell.a2 == 50.0

    def test_coated_ellipsoid(self):
        """Test coated ellipsoid (core-shell)."""
        core = Ellipsoid(number=0, coatingof=0)
        shell = Ellipsoid(number=1, coatingof=0)
        assert shell.coatingof == 0


# =============================================================================
# Cylinder Tests
# =============================================================================

class TestCylinder:
    """Tests for Cylinder dataclass."""

    def test_default_cylinder(self):
        """Test default cylinder."""
        cyl = Cylinder()
        assert cyl.L == 1.0
        assert cyl.R == 1.0

    def test_fiber_cylinder(self):
        """Test fiber-like cylinder."""
        cyl = Cylinder(
            number=0,
            concentration=0.3,
            L=100.0,
            R=1.0
        )
        assert cyl.L / cyl.R == 100.0


# =============================================================================
# JSON I/O Tests
# =============================================================================

class TestLayersJSON:
    """Tests for layers JSON I/O."""

    def test_save_and_load_layers(self):
        """Test saving and loading layers to/from JSON."""
        layers = [
            Layer(
                number=0,
                umat_name='ELISO',
                concentration=0.5,
                props=np.array([70000, 0.3]),
                material_orientation=MaterialOrientation(0, 0, 0),
                geometry_orientation=GeometryOrientation(0, 90, -90)
            ),
            Layer(
                number=1,
                umat_name='ELISO',
                concentration=0.5,
                props=np.array([150000, 0.25]),
                material_orientation=MaterialOrientation(0, 0, 0),
                geometry_orientation=GeometryOrientation(0, 90, -90)
            ),
        ]

        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = os.path.join(tmpdir, 'layers.json')
            save_layers_json(filepath, layers)

            assert os.path.exists(filepath)

            loaded = load_layers_json(filepath)
            assert len(loaded) == 2
            assert loaded[0].number == 0
            assert loaded[0].concentration == 0.5
            assert loaded[1].umat_name == 'ELISO'


class TestEllipsoidsJSON:
    """Tests for ellipsoids JSON I/O."""

    def test_save_and_load_ellipsoids(self):
        """Test saving and loading ellipsoids to/from JSON."""
        ellipsoids = [
            Ellipsoid(
                number=0,
                coatingof=0,
                umat_name='ELISO',
                concentration=0.7,
                props=np.array([70000, 0.3]),
                a1=1, a2=1, a3=1
            ),
            Ellipsoid(
                number=1,
                coatingof=0,
                umat_name='ELISO',
                concentration=0.3,
                props=np.array([400000, 0.2]),
                a1=10, a2=1, a3=1
            ),
        ]

        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = os.path.join(tmpdir, 'ellipsoids.json')
            save_ellipsoids_json(filepath, ellipsoids)

            loaded = load_ellipsoids_json(filepath)
            assert len(loaded) == 2
            assert loaded[0].a1 == 1.0
            assert loaded[1].a1 == 10.0


class TestCylindersJSON:
    """Tests for cylinders JSON I/O."""

    def test_save_and_load_cylinders(self):
        """Test saving and loading cylinders to/from JSON."""
        cylinders = [
            Cylinder(
                number=0,
                concentration=0.7,
                L=1.0,
                R=1.0,
                props=np.array([70000, 0.3])
            ),
            Cylinder(
                number=1,
                concentration=0.3,
                L=50.0,
                R=1.0,
                props=np.array([400000, 0.2])
            ),
        ]

        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = os.path.join(tmpdir, 'cylinders.json')
            save_cylinders_json(filepath, cylinders)

            loaded = load_cylinders_json(filepath)
            assert len(loaded) == 2
            assert loaded[0].L == 1.0
            assert loaded[1].L == 50.0


# =============================================================================
# Integration Tests
# =============================================================================

class TestMicromechanicsIntegration:
    """Integration tests for micromechanics module."""

    def test_composite_definition(self):
        """Test defining a complete composite microstructure."""
        matrix = Ellipsoid(
            number=0,
            umat_name='ELISO',
            concentration=0.6,
            props=np.array([3500, 0.35]),
            a1=1, a2=1, a3=1
        )

        fibers = Ellipsoid(
            number=1,
            umat_name='ELISO',
            concentration=0.4,
            props=np.array([230000, 0.2]),
            a1=100, a2=1, a3=1,
            geometry_orientation=GeometryOrientation(0, 0, 0)
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = os.path.join(tmpdir, 'composite.json')
            save_ellipsoids_json(filepath, [matrix, fibers])

            loaded = load_ellipsoids_json(filepath)
            assert sum(e.concentration for e in loaded) == pytest.approx(1.0)

    def test_laminate_definition(self):
        """Test defining a laminate structure."""
        layers = [
            Layer(number=0, concentration=0.25, geometry_orientation=GeometryOrientation(0, 90, -90)),
            Layer(number=1, concentration=0.25, geometry_orientation=GeometryOrientation(90, 90, -90)),
            Layer(number=2, concentration=0.25, geometry_orientation=GeometryOrientation(90, 90, -90)),
            Layer(number=3, concentration=0.25, geometry_orientation=GeometryOrientation(0, 90, -90)),
        ]

        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = os.path.join(tmpdir, 'laminate.json')
            save_layers_json(filepath, layers)

            loaded = load_layers_json(filepath)
            assert len(loaded) == 4
            assert sum(l.concentration for l in loaded) == pytest.approx(1.0)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
