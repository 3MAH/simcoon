"""
Tests for the simcoon.solver.micromechanics module.

Tests the JSON-based I/O for phase configurations (ellipsoids, layers, cylinders).
"""

import pytest
import tempfile
import os
import numpy as np

import simcoon as sim
from simcoon.solver import StepMeca, solve
from simcoon.solver.micromechanics import (
    EULER_SEQ,
    Peak,
    as_rotation,
    discretize_odf,
    euler_angles,
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
    to_phase_dicts,
)


# =============================================================================
# Orientation tests: simcoon.Rotation in, Euler angles of the files out
# =============================================================================

class TestOrientations:
    """Orientations are Rotation objects; the files keep (psi, theta, phi) in degrees."""

    def test_default_is_identity(self):
        assert as_rotation(None).is_identity()
        assert euler_angles(None) == {"psi": 0.0, "theta": 0.0, "phi": 0.0}

    def test_accepted_forms_agree(self):
        by_seq = as_rotation((45.0, 30.0, 60.0))
        by_dict = as_rotation({"psi": 45.0, "theta": 30.0, "phi": 60.0})
        by_rot = as_rotation(sim.Rotation.from_euler(EULER_SEQ, [45, 30, 60], degrees=True))
        assert by_seq.equals(by_dict) and by_seq.equals(by_rot)
        assert as_rotation(by_rot) is by_rot

    def test_euler_round_trip(self):
        angles = euler_angles((45.0, 30.0, 60.0))
        assert angles == pytest.approx({"psi": 45.0, "theta": 30.0, "phi": 60.0})
        # gimbal lock (theta = 0): the same rotation, the z angles merged into psi
        merged = euler_angles((10.0, 0.0, 20.0))
        assert merged["theta"] == pytest.approx(0.0)
        assert as_rotation(merged).equals(as_rotation((10.0, 0.0, 20.0)))

    def test_bad_forms_are_rejected(self):
        with pytest.raises(ValueError):
            as_rotation((1.0, 2.0))
        with pytest.raises(ValueError):
            as_rotation({"alpha": 1.0})

    def test_dataclass_fields_are_rotations(self):
        ell = Ellipsoid(number=1, geometry_orientation=(0, 90, -90),
                        material_orientation={"psi": 45})
        assert isinstance(ell.geometry_orientation, sim.Rotation)
        assert isinstance(ell.material_orientation, sim.Rotation)
        d = to_phase_dicts([ell])[0]
        assert d["geometry_orientation"] == pytest.approx({"psi": 0.0, "theta": 90.0, "phi": -90.0})
        assert d["material_orientation"] == pytest.approx({"psi": 45.0, "theta": 0.0, "phi": 0.0})


class TestOrientationConvention:
    """Pins EULER_SEQ against the C++ side: the material route of the solver and the
    geometry route of the mean-field schemes."""

    ELIST = [3, 230000., 15000., 0.02, 0.4, 50000., 0., 0.]

    def test_solver_material_frame(self):
        psi, theta, phi = 30.0, 40.0, 50.0
        L = sim.L_isotrans(self.ELIST[1:6], 3)
        e = np.array([0.01, 0.004, -0.002, 0.003, 0.001, 0.002])
        for orientation in ((psi, theta, phi), as_rotation((psi, theta, phi))):
            res = solve(StepMeca(control="strain", value=e, ninc=1), "ELIST", self.ELIST, 1,
                        orientation=orientation)
            R = sim.Rotation.from_euler(EULER_SEQ, [psi, theta, phi], degrees=True)
            np.testing.assert_allclose(res["Stress"][:, -1], R.apply_stiffness(L) @ e,
                                       rtol=1e-10, atol=1e-8)

    def test_l_eff_geometry_frame(self):
        # rotating the fibre of an isotropic-matrix composite rotates its stiffness
        props = np.array([20.0, 20.0, 0.0])   # mp, np, matrix index
        def composite(rot):
            matrix = Ellipsoid(umat_name="ELISO", concentration=0.8, nstatev=1,
                               props=[5000.0, 0.3, 0.0])
            fibre = Ellipsoid(umat_name="ELISO", concentration=0.2, nstatev=1,
                              props=[50000.0, 0.3, 0.0], a1=50.0, geometry_orientation=rot)
            return to_phase_dicts([matrix, fibre])
        R = sim.Rotation.from_euler(EULER_SEQ, [30, 40, 50], degrees=True)
        L0 = np.asarray(sim.L_eff("MIMTN", props, 10000, phases=composite(None)))
        L_R = np.asarray(sim.L_eff("MIMTN", props, 10000, phases=composite(R)))
        np.testing.assert_allclose(L_R, R.apply_stiffness(L0), rtol=1e-9, atol=1e-6)
        # the RVE frame itself, given as a Rotation or as its Euler angles
        L_rve = sim.L_eff("MIMTN", props, 10000, orientation=R, phases=composite(None))
        np.testing.assert_allclose(L_rve, R.apply_stiffness(L0), rtol=1e-9, atol=1e-6)
        L_deg = sim.L_eff("MIMTN", props, 10000, orientation=(30, 40, 50), phases=composite(None))
        np.testing.assert_allclose(L_deg, L_rve, rtol=1e-12)


class TestDiscretizeODF:
    """A phase split along an ODF about a direction."""

    @staticmethod
    def composite(fibre_rot=None):
        matrix = Ellipsoid(umat_name="ELISO", concentration=0.7, nstatev=1,
                           props=[5000.0, 0.3, 0.0])
        fibre = Ellipsoid(umat_name="ELISO", concentration=0.3, nstatev=1,
                          props=[50000.0, 0.3, 0.0], a1=50.0, geometry_orientation=fibre_rot)
        return [matrix, fibre]

    def test_uniform_odf_shares_the_concentration_evenly(self):
        phases = discretize_odf(self.composite(), 1, [Peak(method=7)], 6)
        assert len(phases) == 7 and [p.number for p in phases] == list(range(7))
        assert phases[0].concentration == pytest.approx(0.7)
        np.testing.assert_allclose([p.concentration for p in phases[1:]], 0.3 / 6, rtol=1e-12)

    def test_orientations_sweep_about_the_axis_from_the_base(self):
        base = sim.Rotation.from_euler(EULER_SEQ, [10, 20, 30], degrees=True)
        n = np.array([0.0, 1.0, 0.0])
        phases = discretize_odf(self.composite(base), 1, [Peak(method=7)], 4, axis=n,
                                angle_range=(0.0, 180.0), rotate_material=False)
        for k, ph in enumerate(phases[1:]):
            expected = sim.Rotation.from_rotvec(np.deg2rad(45.0 * k) * n) * base
            assert ph.geometry_orientation.equals(expected, tol=1e-12)
            assert ph.material_orientation.is_identity()
        phases = discretize_odf(self.composite(base), 1, [Peak(method=7)], 4, axis=n)
        assert phases[3].material_orientation.equals(phases[3].geometry_orientation * base.inv(), tol=1e-12)

    def test_malformed_peaks_are_refused(self):
        x = np.linspace(0.0, 90.0, 5)
        for peaks, message in (([{"mean": 30.0}], "no 'method'"), ([{"method": 0}], "1 to 7"),
                               ([{"method": 3, "methd": 1}], "unknown entries"),
                               ([{"method": 3, "s_dev": 0.0}], "s_dev"),
                               ([{"method": 4, "width": 0.0}], "width"),
                               ([{"method": 5, "params": []}], "params needs 1"),
                               ([{"method": 6, "params": [1.0, 0.0]}], "shape")):
            with pytest.raises(ValueError, match=message):
                sim.get_densities_ODF(x, peaks)
        with pytest.raises(ValueError, match="must lie in"):
            sim.get_densities_ODF([-1.0, 10.0], [Peak()])

    @pytest.mark.parametrize("peak, expected", [
    ({'method': 1, 'mean': 40.0, 'params': [1.0, 0.5, 2.0, 1.0]},
     [0.3566651602447377, 0.9945303330076852, 0.1876587073661953, 0.005979608223535716, 0.3566651602447375]),
    ({'method': 2, 'mean': 170.0, 's_dev': 12.0, 'ampl': 2.0},
     [1.4132965557154336, 0.0009331060403447629, 4.467279520345654e-10, 0.15911901743645537, 1.4132965557154336]),
    ({'method': 3, 'mean': 5.0, 's_dev': 10.0, 'ampl': 1.5},
     [3.0257786004732883, 0.020489728789630575, 7.019228677806676e-16, 0.0005065783521005087, 3.0257786004732883]),
    ({'method': 4, 'mean': 90.0, 'width': 15.0, 'ampl': 0.7},
     [0.024790974022486378, 0.04109480521872039, 1.7080960442387314, 0.04109480521872037, 0.024790974022486378]),
    ({'method': 5, 'mean': 120.0, 's_dev': 8.0, 'width': 20.0, 'ampl': 1.2, 'params': [0.35]},
     [0.02683521260168469, 0.02012112892597051, 0.08369086662031242, 0.16247010119366517, 0.027312595059575475]),
    ({'method': 6, 'mean': 60.0, 'width': 25.0, 'params': [0.0, 1.7]},
     [0.09256584316280778, 0.5094880999984644, 0.3592136259508243, 0.05391559545149874, 0.09199049533876534]),
    ({'method': 7},
     [1.0, 1.0, 1.0, 1.0, 1.0]),
    ])
    def test_densities_match_the_cpp_implementation_they_replace(self, peak, expected):
        """Values of simcoon::get_densities_ODF (1.x to 2.0) at x = 0, 37, 90, 143, 180 deg."""
        x = np.array([0.0, 37.0, 90.0, 143.0, 180.0])
        np.testing.assert_allclose(sim.get_densities_ODF(x, [peak]), expected, rtol=1e-12, atol=1e-15)
        np.testing.assert_allclose(sim.get_densities_ODF(x, [Peak(**peak)]), expected, rtol=1e-12, atol=1e-15)
        rad = {k: (np.deg2rad(v) if k in ("mean", "s_dev", "width") else v) for k, v in peak.items()}
        np.testing.assert_allclose(sim.get_densities_ODF(np.deg2rad(x), [rad], radian=True), expected,
                                   rtol=1e-12, atol=1e-15)

    def test_density_is_periodic_over_a_half_turn(self):
        pk = Peak(method=3, mean=175.0, s_dev=10.0)
        np.testing.assert_allclose(sim.get_densities_ODF([0.0], [pk]), sim.get_densities_ODF([180.0], [pk]), rtol=1e-12)
        assert pk.density(np.deg2rad(5.0), periodic=True, scale=np.pi / 180) > 100 * pk.density(np.deg2rad(5.0), scale=np.pi / 180)

    def test_gaussian_peak_weights_follow_the_density(self):
        # a narrow Gaussian at 90 deg: the phases near 90 deg carry the mass, the sum is the parent
        peaks = [Peak(method=3, mean=90.0, s_dev=10.0, ampl=1.0)]
        phases = discretize_odf(self.composite(), 1, peaks, 18)
        w = np.array([p.concentration for p in phases[1:]])
        assert w.sum() == pytest.approx(0.3, rel=1e-12)
        assert w.argmax() == 9 and w[0] < 1e-3 * w[9]

    def test_uniform_sweep_about_z_is_transversely_isotropic(self):
        # a 4th-order tensor rotated about z carries harmonics up to 4 alpha: 8 equally
        # spaced angles over a half turn average them out exactly
        props = np.array([20.0, 20.0, 0.0])
        phases = discretize_odf(self.composite(), 1, [Peak(method=7)], 8)
        L = np.asarray(sim.L_eff("MIMTN", props, 10000, phases=to_phase_dicts(phases)))
        assert L[0, 0] == pytest.approx(L[1, 1], rel=1e-9)
        assert L[0, 2] == pytest.approx(L[1, 2], rel=1e-9)
        assert L[4, 4] == pytest.approx(L[5, 5], rel=1e-9)
        assert L[3, 3] == pytest.approx(0.5 * (L[0, 0] - L[0, 1]), rel=1e-9)
        assert abs(L[0, 3]) < 1e-9 * L[0, 0]


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
                material_orientation=(0, 0, 0),
                geometry_orientation=(0, 90, -90)
            ),
            Layer(
                number=1,
                umat_name='ELISO',
                concentration=0.5,
                props=np.array([150000, 0.25]),
                material_orientation=(0, 0, 0),
                geometry_orientation=(0, 90, -90)
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
            geometry_orientation=(0, 0, 0)
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = os.path.join(tmpdir, 'composite.json')
            save_ellipsoids_json(filepath, [matrix, fibers])

            loaded = load_ellipsoids_json(filepath)
            assert sum(e.concentration for e in loaded) == pytest.approx(1.0)

    def test_laminate_definition(self):
        """Test defining a laminate structure."""
        layers = [
            Layer(number=0, concentration=0.25, geometry_orientation=(0, 90, -90)),
            Layer(number=1, concentration=0.25, geometry_orientation=(90, 90, -90)),
            Layer(number=2, concentration=0.25, geometry_orientation=(90, 90, -90)),
            Layer(number=3, concentration=0.25, geometry_orientation=(0, 90, -90)),
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
        phase = Phase(
            number=0,
            umat_name='ELORT',
            material_orientation=(45.0, 30.0, 0.0)
        )
        assert euler_angles(phase.material_orientation)["psi"] == pytest.approx(45.0)


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


class TestNestedPhases:
    """A sub-phase that is itself a mean-field model carries its own sub-phases."""

    def _composite(self):
        inner = [Ellipsoid(number=0, concentration=0.8, props=[5000.0, 0.3, 0.0]),
                 Ellipsoid(number=1, concentration=0.2, a1=50.0, props=[50000.0, 0.3, 0.0],
                           geometry_orientation=(45.0, 0.0, 0.0))]
        outer = [Ellipsoid(number=0, umat_name="MIMTN", concentration=0.8, nstatev=1000,
                           props=[20.0, 20.0, 0.0], phases=inner),
                 Ellipsoid(number=1, concentration=0.2, a1=50.0, props=[50000.0, 0.3, 0.0])]
        return outer

    def test_json_round_trip(self, tmp_path):
        path = tmp_path / "nested.json"
        save_ellipsoids_json(path, self._composite())
        loaded = load_ellipsoids_json(path)
        assert [p.umat_name for p in loaded] == ["MIMTN", "ELISO"]
        inner = loaded[0].phases
        assert [type(p) for p in inner] == [Ellipsoid, Ellipsoid]
        assert inner[1].a1 == 50.0
        assert euler_angles(inner[1].geometry_orientation)["psi"] == pytest.approx(45.0)
        np.testing.assert_array_equal(inner[0].props, [5000.0, 0.3, 0.0])
        assert loaded[1].phases == []

    def test_to_phase_dict_carries_kind_and_nesting(self):
        dicts = to_phase_dicts(self._composite())
        assert dicts[0]["kind"] == "ellipsoid"
        assert [d["number"] for d in dicts[0]["phases"]] == [0, 1]
        assert "phases" not in dicts[1]
        assert to_phase_dicts([Layer()])[0]["kind"] == "layer"
        assert to_phase_dicts([Phase()])[0]["kind"] == "phase"

    def test_dict_form_is_coerced_back(self):
        outer = Ellipsoid(umat_name="MIMTN", props=[2, 1, 20, 20, 0],
                          phases=[{"umat_name": "ELISO", "concentration": 0.5, "props": [1, 0.3, 0]},
                                  {"umat_name": "ELISO", "concentration": 0.5, "props": [2, 0.3, 0]}])
        assert all(isinstance(p, Ellipsoid) for p in outer.phases)
        assert outer.phases[1].props[0] == 2
