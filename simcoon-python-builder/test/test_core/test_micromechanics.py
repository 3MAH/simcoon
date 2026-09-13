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
    Section,
    load_phases_json,
    save_phases_json,
    load_layers_json,
    save_layers_json,
    load_ellipsoids_json,
    save_ellipsoids_json,
    load_cylinders_json,
    save_cylinders_json,
    load_sections_json,
    save_sections_json,
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


# =============================================================================
# Phase Tests
# =============================================================================

class TestPhase:
    """Tests for Phase dataclass."""

    def test_default_phase(self):
        """Test default phase initialization."""
        phase = Phase()
        assert phase.number == 0
        assert phase.umat_name == 'ELISO'
        assert phase.concentration == 1.0
        assert phase.save == 1

    def test_phase_with_properties(self):
        """Test phase with material properties."""
        phase = Phase(
            number=1,
            umat_name='EPICP',
            props=np.array([70000, 0.3, 0, 300, 1000, 0.5]),
            concentration=0.4,
            nstatev=8
        )
        assert phase.number == 1
        assert phase.concentration == 0.4
        assert phase.nstatev == 8

    def test_phase_with_orientation(self):
        """Test phase with material orientation."""
        orient = MaterialOrientation(psi=45.0, theta=30.0, phi=0.0)
        phase = Phase(
            number=0,
            umat_name='ELORT',
            material_orientation=orient
        )
        assert phase.material_orientation.psi == 45.0


class TestPhasesJSON:
    """Tests for phases JSON I/O."""

    def test_save_and_load_phases(self):
        """Test saving and loading phases to/from JSON."""
        phases = [
            Phase(
                number=0,
                umat_name='ELISO',
                props=np.array([70000, 0.3]),
                concentration=0.6,
                nstatev=1
            ),
            Phase(
                number=1,
                umat_name='ELISO',
                props=np.array([400000, 0.2]),
                concentration=0.4,
                nstatev=1
            ),
        ]

        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = os.path.join(tmpdir, 'phases.json')
            save_phases_json(filepath, phases)

            assert os.path.exists(filepath)

            loaded = load_phases_json(filepath)
            assert len(loaded) == 2
            assert loaded[0].concentration == 0.6
            assert loaded[1].concentration == 0.4


# =============================================================================
# Section Tests
# =============================================================================

class TestSection:
    """Tests for Section dataclass."""

    def test_default_section(self):
        """Test default section initialization."""
        section = Section()
        assert section.number == 0
        assert section.name == 'Section'
        assert section.umat_name == 'ELISO'

    def test_section_with_properties(self):
        """Test section with material properties."""
        section = Section(
            number=0,
            name='yarn_weft',
            umat_name='ELISO',
            props=np.array([70000, 0.3]),
            nstatev=1
        )
        assert section.name == 'yarn_weft'
        assert len(section.props) == 2


class TestSectionsJSON:
    """Tests for sections JSON I/O."""

    def test_save_and_load_sections(self):
        """Test saving and loading sections to/from JSON."""
        sections = [
            Section(
                number=0,
                name='yarn_0',
                umat_name='ELISO',
                props=np.array([70000, 0.3]),
                nstatev=1
            ),
            Section(
                number=1,
                name='yarn_1',
                umat_name='ELISO',
                props=np.array([400000, 0.2]),
                nstatev=1
            ),
        ]

        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = os.path.join(tmpdir, 'sections.json')
            save_sections_json(filepath, sections)

            assert os.path.exists(filepath)

            loaded = load_sections_json(filepath)
            assert len(loaded) == 2
            assert loaded[0].name == 'yarn_0'
            assert loaded[1].name == 'yarn_1'


# =============================================================================
# Edge Cases and Error Handling
# =============================================================================

class TestMicromechanicsEdgeCases:
    """Tests for edge cases and error handling."""

    def test_empty_props_array(self):
        """Test geometry with empty properties array."""
        ell = Ellipsoid(number=0, props=np.array([]))
        assert len(ell.props) == 0

    def test_large_aspect_ratio(self):
        """Test ellipsoid with very large aspect ratio (fiber)."""
        ell = Ellipsoid(
            number=0,
            a1=1000.0,
            a2=1.0,
            a3=1.0
        )
        assert ell.a1 / ell.a2 == 1000.0

    def test_thin_layer(self):
        """Test very thin layer."""
        layer = Layer(
            number=0,
            concentration=0.001
        )
        assert layer.concentration == 0.001

    def test_json_file_not_found(self):
        """Test loading non-existent JSON file."""
        with pytest.raises(FileNotFoundError):
            load_ellipsoids_json('/nonexistent/path.json')


# The legacy .dat readers arrived with the JSON-only migration, when the parsing
# of those files moved out of the C++ src/Simulation/Phase/read.cpp.
from simcoon.solver.micromechanics import (
    convert_dat_to_json,
    kind_from_dat_name,
    load_cylinders_dat,
    load_ellipsoids_dat,
    load_layers_dat,
    load_phases_dat,
    load_sections_dat,
)

# Verbatim copies of the testBin fixtures, ragged tabs included.
PHASES_DAT = (
    "Number\tumat\tsave \tc\tpsi_mat\ttheta_mat \tphi_mat\tnprops\tnstatev\tprops\n"
    "0\tELISO\t1\t0.8\t0\t0\t\t0\t3\t1\t70000\t0.4\t0\n"
    "1\tELISO\t1\t0.2\t0\t0\t\t0\t3\t1\t3000\t0.4\t0\n"
    "\n\n\n"
)
LAYERS_DAT = (
    "Number\tumat\tsave\tc\tpsi_mat\ttheta_mat\tphi_mat\tpsi_geom\ttheta_geom\tphi_geom"
    "\tnprops\tnstatev\tprops\n"
    "0\tELISO\t1\t0.8\t0\t0\t0\t0\t90\t-90\t3\t1\t3000\t0.4\t0\n"
    "1\tELISO\t1\t0.2\t0\t0\t0\t0\t90\t-90\t3\t1\t70000\t0.3\t0\n"
)
ELLIPSOIDS_DAT = (
    "Number\tCoatingof\tumat\tsave\tc\tpsi_mat\ttheta_mat\tphi_mat\ta1\ta2\ta3"
    "\tpsi_geom\ttheta_geom\tphi_geom\tnprops\tnstatev\tprops\n"
    "0\t0\tELISO\t1\t0.8\t0\t0\t0\t1\t1\t1\t0\t0\t0\t3\t1\t3000\t0.4\t0\n"
    "1\t0\tELISO\t1\t0.2\t0\t0\t0\t50\t1\t1\t0\t0\t0\t3\t1\t70000\t0.3\t0\n"
)
CYLINDERS_DAT = (
    "Number\tCoatingof\tumat\tsave\tc\tpsi_mat\ttheta_mat\tphi_mat\tL\tR"
    "\tpsi_geom\ttheta_geom\tphi_geom\tnprops\tnstatev\tprops\n"
    "0\t0\tELISO\t1\t0.8\t0\t0\t0\t1\t1\t0\t0\t0\t3\t1\t3000\t0.4\t0\n"
    "1\t0\tELISO\t1\t0.2\t0\t0\t0\t50\t1\t0\t0\t0\t3\t1\t70000\t0.3\t0\n"
)
SECTIONS_DAT = (
    "Number\tSection_name\tumat\tpsi_mat\ttheta_mat \tphi_mat\tnprops\tnstatev\tprops\n"
    "0\tYarn0\t\tELISO\t0\t0\t\t0\t3\t1\t70000\t0.4\t0\n"
    "1\tYarn1\t\tELISO\t0\t0\t\t0\t3\t1\t3000\t0.4\t0\n"
)


def _write_dat(tmp_path, name, content):
    path = tmp_path / name
    path.write_text(content)
    return path


class TestLegacyDatReaders:
    """The historical .dat formats, parsed in Python instead of C++."""

    def test_phases(self, tmp_path):
        phases = load_phases_dat(_write_dat(tmp_path, "Nphases0.dat", PHASES_DAT))
        assert len(phases) == 2
        assert phases[0].umat_name == "ELISO"
        assert phases[0].save == 1
        assert phases[0].concentration == 0.8
        assert phases[0].nstatev == 1
        np.testing.assert_allclose(phases[0].props, [70000, 0.4, 0])
        np.testing.assert_allclose(phases[1].props, [3000, 0.4, 0])

    def test_layers_carry_geometry_orientation(self, tmp_path):
        layers = load_layers_dat(_write_dat(tmp_path, "Nlayers0.dat", LAYERS_DAT))
        assert len(layers) == 2
        assert layers[0].geometry_orientation.theta == 90
        assert layers[0].geometry_orientation.phi == -90
        assert layers[0].material_orientation.psi == 0

    def test_ellipsoids_semi_axes_and_shape(self, tmp_path):
        ells = load_ellipsoids_dat(_write_dat(tmp_path, "Nellipsoids0.dat", ELLIPSOIDS_DAT))
        assert [e.number for e in ells] == [0, 1]
        assert ells[0].shape_type == "sphere"
        assert (ells[1].a1, ells[1].a2, ells[1].a3) == (50, 1, 1)
        assert ells[1].shape_type == "prolate_spheroid"
        assert ells[1].coatingof == 0

    def test_cylinders_length_and_radius(self, tmp_path):
        cyls = load_cylinders_dat(_write_dat(tmp_path, "Ncylinders0.dat", CYLINDERS_DAT))
        assert (cyls[0].L, cyls[0].R) == (1, 1)
        assert cyls[1].aspect_ratio == 50

    def test_sections_keep_their_name(self, tmp_path):
        secs = load_sections_dat(_write_dat(tmp_path, "Nsections0.dat", SECTIONS_DAT))
        assert [s.name for s in secs] == ["Yarn0", "Yarn1"]
        assert secs[0].umat_name == "ELISO"
        np.testing.assert_allclose(secs[0].props, [70000, 0.4, 0])

    def test_trailing_blank_lines_are_ignored(self, tmp_path):
        # PHASES_DAT ends with three empty lines, as the shipped fixture does.
        assert len(load_phases_dat(_write_dat(tmp_path, "Nphases0.dat", PHASES_DAT))) == 2

    def test_row_contradicting_its_nprops_is_rejected(self, tmp_path):
        truncated = PHASES_DAT.replace("\t3\t1\t70000\t0.4\t0", "\t3\t1\t70000\t0.4")
        with pytest.raises(ValueError, match="announces"):
            load_phases_dat(_write_dat(tmp_path, "Nphases0.dat", truncated))

    def test_non_integer_nprops_is_rejected(self, tmp_path):
        broken = PHASES_DAT.replace("\t3\t1\t70000", "\tthree\t1\t70000")
        with pytest.raises(ValueError, match="nprops"):
            load_phases_dat(_write_dat(tmp_path, "Nphases0.dat", broken))

    def test_empty_file_is_rejected(self, tmp_path):
        with pytest.raises(ValueError, match="header"):
            load_phases_dat(_write_dat(tmp_path, "Nphases0.dat", ""))

    def test_missing_file(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            load_phases_dat(tmp_path / "absent.dat")

    def test_abaqus_deck_is_named_for_what_it_is(self, tmp_path):
        # testBin ships Nsections1.dat in this form; it is not tabular data.
        deck = "** ==================\n*Material, name=ELISO-0\n*Depvar\n     5\n"
        with pytest.raises(ValueError, match="Abaqus"):
            load_sections_dat(_write_dat(tmp_path, "Nsections1.dat", deck))

    def test_kind_inferred_from_name(self):
        assert kind_from_dat_name("Nellipsoids0.dat") == "ellipsoids"
        assert kind_from_dat_name("Nphases12.dat") == "phases"
        with pytest.raises(ValueError):
            kind_from_dat_name("whatever.dat")


class TestDatToJsonConversion:
    """The one-way door: a legacy file in, the JSON we now write out."""

    def test_ellipsoids_round_trip_through_json(self, tmp_path):
        dat = _write_dat(tmp_path, "Nellipsoids0.dat", ELLIPSOIDS_DAT)
        out = convert_dat_to_json(dat)
        assert out.name == "ellipsoids0.json"

        from_dat = load_ellipsoids_dat(dat)
        from_json = load_ellipsoids_json(out)
        assert len(from_json) == len(from_dat)
        for a, b in zip(from_dat, from_json):
            assert (a.number, a.umat_name, a.coatingof) == (b.number, b.umat_name, b.coatingof)
            assert (a.a1, a.a2, a.a3) == (b.a1, b.a2, b.a3)
            assert a.concentration == b.concentration
            assert a.geometry_orientation.theta == b.geometry_orientation.theta
            np.testing.assert_allclose(a.props, b.props)

    def test_layers_conversion_to_an_explicit_path(self, tmp_path):
        dat = _write_dat(tmp_path, "Nlayers0.dat", LAYERS_DAT)
        out = convert_dat_to_json(dat, tmp_path / "chosen.json")
        assert out.name == "chosen.json"
        layers = load_layers_json(out)
        assert [lay.geometry_orientation.phi for lay in layers] == [-90, -90]

    def test_unknown_kind_is_refused(self, tmp_path):
        dat = _write_dat(tmp_path, "mystery.dat", PHASES_DAT)
        with pytest.raises(ValueError):
            convert_dat_to_json(dat)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
