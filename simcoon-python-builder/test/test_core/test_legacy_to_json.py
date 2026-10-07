"""The legacy text formats (path.txt, material.dat, N<kind>.dat, ...) are parsed by
scripts/legacy_to_json.py only, a migration tool kept outside the package: these tests
load the script and exercise its parsers and the one-way conversion to JSON."""

import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

import simcoon as sim
from simcoon import solver as slv
from simcoon.solver import StepMeca, StepThermomeca, solve
from simcoon.solver.micromechanics import (euler_angles, load_ellipsoids_json,
                                           load_layers_json)
from test_solver_run import (ELISO_PROPS, ELISO_T_PROPS, _FILE_ORDER,
                             _UNIAXIAL, assert_ran_and_responded)

_SCRIPT = Path(__file__).resolve().parents[3] / "scripts" / "legacy_to_json.py"
_spec = importlib.util.spec_from_file_location("legacy_to_json", _SCRIPT)
legacy = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(legacy)


def _meca_state_lines(flags, values):
    """Emit the #prescribed_mechanical_state lines in path-file order."""
    toks = []
    for k in _FILE_ORDER:
        toks.append(f"{flags[k]} {values[k]}")
    return "{}\n{} {}\n{} {} {}".format(*toks)


def _step_text(flags, values, time=1.0, ninc=100, mode=1, T=290.0, BC_w=None):
    txt = f"""#Mode
{mode}
#Dn_init 1.
#Dn_mini 1.
#Dn_inc {1.0/ninc}
#time
{time}
#prescribed_mechanical_state
{_meca_state_lines(flags, values)}
"""
    if BC_w is not None:
        BC_w = np.asarray(BC_w)
        txt += "#Rotation\n" + "\n".join(
            " ".join(str(BC_w[i, j]) for j in range(3)) for i in range(3)
        ) + "\n"
    txt += f"#prescribed_temperature_state\nT {T}\n"
    return txt


def _thermo_step_text(flags, values, thermal_token, time=1.0, ninc=50):
    return f"""#Mode
1
#Dn_init 1.
#Dn_mini 1.
#Dn_inc {1.0/ninc}
#time
{time}
#prescribed_mechanical_state
{_meca_state_lines(flags, values)}
#prescribed_thermal_state
{thermal_token}
"""


def _path_text(steps_text, control_type=1, loading_type=1, ncycle=1, T_init=290.0):
    return f"""#Initial_temperature
{T_init}
#Number_of_blocks
1

#Block
1
#Loading_type
{loading_type}
#Control_type(NLGEOM)
{control_type}
#Repeat
{ncycle}
#Steps
{len(steps_text)}

""" + "\n".join(steps_text)


def write_path_file(tmp_path, path_text, extra_files=None):
    """Drop a legacy path.txt (and its tabular files) in tmp_path/data."""
    data = tmp_path / "data"
    data.mkdir(exist_ok=True)
    (data / "path.txt").write_text(path_text)
    for name, content in (extra_files or {}).items():
        (data / name).write_text(content)
    return str(data)




def test_from_file_mechanical(tmp_path):
    # parse a legacy path.txt into Blocks and run the programme it describes
    flags = ["E"] + ["S"] * 5
    steps = [_step_text(flags, [0.02, 0, 0, 0, 0, 0], ninc=50),
             _step_text(flags, [0.0, 0, 0, 0, 0, 0], ninc=50)]
    data = write_path_file(tmp_path, _path_text(steps, ncycle=2))

    blocks, T_init = legacy.from_file(data, "path.txt")
    assert T_init == 290.0
    assert len(blocks) == 1 and blocks[0].ncycle == 2 and len(blocks[0].steps) == 2

    res = solve(blocks, "ELISO", ELISO_PROPS, 1, T_init=T_init)
    assert_ran_and_responded(res, n_expected=2 * 2 * 50)
    # the parsed programme loads to 0.02 and unloads to 0, twice
    np.testing.assert_allclose(res["Strain"][0].max(), 0.02, rtol=1e-8)
    np.testing.assert_allclose(res["Strain"][0, -1], 0.0, atol=1e-12)


def test_convert_to_json(tmp_path):
    """A legacy data directory converted once: the JSON says what the text said."""
    flags = ["E"] + ["S"] * 5
    steps = [_step_text(flags, [0.02, 0, 0, 0, 0, 0], ninc=50),
             _step_text(flags, [0.0, 0, 0, 0, 0, 0], ninc=50)]
    data = write_path_file(tmp_path, _path_text(steps, ncycle=2))
    (Path(data) / "material.dat").write_text(
        "Material\nName\tELISO\nNumber_of_material_parameters\t3\n"
        "Number_of_internal_variables\t1\n\n#Orientation\npsi\t0\ntheta\t0\nphi\t0\n\n"
        "#Mechanical\nE\t70000\nnu\t0.3\nalpha\t1e-5\n")
    (Path(data) / "solver_essentials.inp").write_text(
        "Solver_type_0_Newton_tangent_1_RNL\n0\nRate_type\n2\n")

    written = legacy.convert_to_json(data)
    assert sorted(Path(w).name for w in written) == ["material.json", "path.json"]

    from_json = slv.load_simulation_json(Path(data) / "material.json",
                                                Path(data) / "path.json")
    blocks, T_init = legacy.from_file(data, "path.txt")
    assert from_json["corate"] == "logarithmic"          # Rate_type 2 of the .inp
    assert from_json["T_init"] == T_init
    as_text = lambda bs: json.dumps([b.to_dict(T_init) for b in bs],
                                    default=lambda o: np.asarray(o).tolist())
    assert as_text(from_json["blocks"]) == as_text(blocks)
    np.testing.assert_array_equal(from_json["props"], ELISO_PROPS)
    res = solve(**from_json)
    assert_ran_and_responded(res, n_expected=2 * 2 * 50)


def test_from_file_thermomechanical(tmp_path):
    flags = ["S"] * 6
    steps = [_thermo_step_text(flags, [0.0] * 6, "T 340")]
    data = write_path_file(tmp_path, _path_text(steps, loading_type=2))

    blocks, T_init = legacy.from_file(data, "path.txt")
    s = blocks[0].steps[0]
    assert isinstance(s, StepThermomeca) and s.thermal_control == "temperature"
    assert s.T_final == 340.0

    res = solve(blocks, "ELISO", ELISO_T_PROPS, 1, T_init=T_init)
    assert res.status == 0
    # the parsed ramp lands on its target, stress-free, with free thermal expansion
    np.testing.assert_allclose(res["Temp"][-1], 340.0, rtol=1e-10)
    np.testing.assert_allclose(res["Strain"][0, -1], 1.0e-5 * 50.0, rtol=1e-8)
    np.testing.assert_allclose(res["Stress"][:, -1], 0.0, atol=1e-6)


def test_from_file_tabular(tmp_path):
    t = np.linspace(0.02, 1.0, 50)
    e11 = 0.015 * np.sin(np.pi * t)
    tab_lines = "\n".join(f"{i+1} {t[i]:.16g} {e11[i]:.16g}" for i in range(len(t)))
    step3 = """#Mode
3
#File
tab.txt
#Dn_init 1.
#Dn_mini 1.
#prescribed_mechanical_state
E
0 0
0 0 0
#T_is_set
0
"""
    data = write_path_file(tmp_path, _path_text([step3]),
                           extra_files={"tab.txt": tab_lines})

    blocks, T_init = legacy.from_file(data, "path.txt")
    s = blocks[0].steps[0]
    assert s.mode == 3 and not s.tabular_T
    assert s.control[0] == "strain" and s.control[1:] == ["zero"] * 5
    np.testing.assert_allclose(s.tabular[:, 0], t, atol=1e-14)

    res = solve(blocks, "ELISO", ELISO_PROPS, 1, T_init=T_init)
    assert_ran_and_responded(res, n_expected=len(t))
    # the table read from disk drives the run, in its own absolute time
    np.testing.assert_allclose(res["Time"], t, atol=1e-12)
    np.testing.assert_allclose(res["Strain"][0], e11, atol=1e-10)


def test_material_from_file(tmp_path):
    (tmp_path / "material.dat").write_text("""Material
Name\tELISO
Number_of_material_parameters\t3
Number_of_internal_variables\t1

#Orientation
psi\t0.1
theta\t0.2
phi\t0.3

#Mechanical
E 70000.
nu 0.3
alpha 1.E-5
""")
    kw = legacy.material_from_file(str(tmp_path), "material.dat")
    assert kw["umat_name"] == "ELISO"
    assert kw["nstatev"] == 1
    np.testing.assert_allclose(kw["props"], [70000.0, 0.3, 1.0e-5])
    assert kw["orientation"] == (0.1, 0.2, 0.3)
    res = solve(StepMeca(control=_UNIAXIAL, value=[0.01, 0, 0, 0, 0, 0], ninc=10),
                T_init=290.0, **kw)
    assert res.status == 0


# No other reader of the legacy grammar exists to compare with, so the parse itself is
# pinned: every field below is a place where a silent drift would produce a different
# loading programme that "status == 0 and the stress is finite" would not notice.
PATH_FIXTURE = """#Initial_temperature
300.0
#Number_of_blocks
2

#Block
1
#Loading_type
1
#Control_type(NLGEOM)
1
#Repeat
3
#Steps
2

#Mode
1
#Dn_init 1.
#Dn_mini 0.1
#Dn_inc 0.02
#time
7.5
#Consigne
E 0.11
S 0.12 E 0.22
S 0.13 S 0.23 E 0.33
#Consigne_T
T 305.

#Mode
2
#Dn_init 0.5
#Dn_mini 0.01
#Dn_inc 0.25
#time
2.
#Consigne
S 0. E 0. S 0. E 0. S 0. E 0.
#Consigne_T
T 310.

#Block
2
#Loading_type
2
#Control_type(NLGEOM)
1
#Repeat
1
#Steps
1

#Mode
1
#Dn_init 1.
#Dn_mini 0.05
#Dn_inc 0.1
#time
4.
#Consigne
E 0.01
S 0. S 0.
S 0. S 0. S 0.
#Consigne_T
Q 1500.
"""


def test_from_file_parses_the_documented_grammar(tmp_path):
    """Field-by-field pin of the legacy path-file grammar as from_file reads it."""
    (tmp_path / "path.txt").write_text(PATH_FIXTURE)
    blocks, T_init = legacy.from_file(str(tmp_path), "path.txt")

    assert T_init == 300.0
    assert len(blocks) == 2

    mech = blocks[0]
    assert mech.ncycle == 3 and mech.control_type == 1
    assert len(mech.steps) == 2

    first = mech.steps[0]
    assert first.mode in (1, "linear")
    assert first.time == 7.5
    assert first.ninc == 50                       # round(1 / Dn_inc)
    assert first.Dn_init == 1.0 and first.Dn_mini == 0.1
    assert first.T_final == 305.0
    # the file lists 11, 12, 22, 13, 23, 33 (lower triangle, row-wise); Voigt is
    # 11, 22, 33, 12, 13, 23. Swapping the two would scramble the shear components.
    assert list(first.control) == ["strain", "strain", "strain",
                                   "stress", "stress", "stress"]
    np.testing.assert_allclose(np.asarray(first.value, dtype=float),
                               [0.11, 0.22, 0.33, 0.12, 0.13, 0.23])

    second = mech.steps[1]
    assert second.mode in (2, "sinusoidal")
    assert second.ninc == 4 and second.time == 2.0
    # S E S E S E read as 11, 12, 22, 13, 23, 33 lands in Voigt as
    # 11=S, 22=S, 33=E, 12=E, 13=E, 23=S — the interleaving is the point
    assert list(second.control) == ["stress", "stress", "strain",
                                    "strain", "strain", "stress"]

    thermo = blocks[1]
    assert thermo.ncycle == 1
    heat = thermo.steps[0]
    assert heat.ninc == 10 and heat.time == 4.0
    # a 'Q' thermal condition is a heat flux, not a temperature target
    assert getattr(heat, "thermal_control", None) in ("heat_flux", 1)
    assert float(getattr(heat, "Q", 0.0)) == 1500.0


# The legacy .dat readers arrived with the JSON-only migration, when the parsing
# of those files moved out of the C++ src/Simulation/Phase/read.cpp.

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
        phases = legacy.load_phases_dat(_write_dat(tmp_path, "Nphases0.dat", PHASES_DAT))
        assert len(phases) == 2
        assert phases[0].umat_name == "ELISO"
        assert phases[0].save == 1
        assert phases[0].concentration == 0.8
        assert phases[0].nstatev == 1
        np.testing.assert_allclose(phases[0].props, [70000, 0.4, 0])
        np.testing.assert_allclose(phases[1].props, [3000, 0.4, 0])

    def test_layers_carry_geometry_orientation(self, tmp_path):
        layers = legacy.load_layers_dat(_write_dat(tmp_path, "Nlayers0.dat", LAYERS_DAT))
        assert len(layers) == 2
        assert euler_angles(layers[0].geometry_orientation)["theta"] == pytest.approx(90)
        assert euler_angles(layers[0].geometry_orientation)["phi"] == pytest.approx(-90)
        assert layers[0].material_orientation.is_identity()

    def test_ellipsoids_semi_axes_and_shape(self, tmp_path):
        ells = legacy.load_ellipsoids_dat(_write_dat(tmp_path, "Nellipsoids0.dat", ELLIPSOIDS_DAT))
        assert [e.number for e in ells] == [0, 1]
        assert ells[0].shape_type == "sphere"
        assert (ells[1].a1, ells[1].a2, ells[1].a3) == (50, 1, 1)
        assert ells[1].shape_type == "prolate_spheroid"
        assert ells[1].coatingof == 0

    def test_cylinders_length_and_radius(self, tmp_path):
        cyls = legacy.load_cylinders_dat(_write_dat(tmp_path, "Ncylinders0.dat", CYLINDERS_DAT))
        assert (cyls[0].L, cyls[0].R) == (1, 1)
        assert cyls[1].aspect_ratio == 50

    def test_sections_keep_their_name(self, tmp_path):
        secs = legacy.load_sections_dat(_write_dat(tmp_path, "Nsections0.dat", SECTIONS_DAT))
        assert [s.name for s in secs] == ["Yarn0", "Yarn1"]
        assert secs[0].umat_name == "ELISO"
        np.testing.assert_allclose(secs[0].props, [70000, 0.4, 0])

    def test_trailing_blank_lines_are_ignored(self, tmp_path):
        # PHASES_DAT ends with three empty lines, as the shipped fixture does.
        assert len(legacy.load_phases_dat(_write_dat(tmp_path, "Nphases0.dat", PHASES_DAT))) == 2

    def test_row_contradicting_its_nprops_is_rejected(self, tmp_path):
        truncated = PHASES_DAT.replace("\t3\t1\t70000\t0.4\t0", "\t3\t1\t70000\t0.4")
        with pytest.raises(ValueError, match="announces"):
            legacy.load_phases_dat(_write_dat(tmp_path, "Nphases0.dat", truncated))

    def test_non_integer_nprops_is_rejected(self, tmp_path):
        broken = PHASES_DAT.replace("\t3\t1\t70000", "\tthree\t1\t70000")
        with pytest.raises(ValueError, match="nprops"):
            legacy.load_phases_dat(_write_dat(tmp_path, "Nphases0.dat", broken))

    def test_empty_file_is_rejected(self, tmp_path):
        with pytest.raises(ValueError, match="header"):
            legacy.load_phases_dat(_write_dat(tmp_path, "Nphases0.dat", ""))

    def test_missing_file(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            legacy.load_phases_dat(tmp_path / "absent.dat")

    def test_abaqus_deck_is_named_for_what_it_is(self, tmp_path):
        # testBin ships Nsections1.dat in this form; it is not tabular data.
        deck = "** ==================\n*Material, name=ELISO-0\n*Depvar\n     5\n"
        with pytest.raises(ValueError, match="Abaqus"):
            legacy.load_sections_dat(_write_dat(tmp_path, "Nsections1.dat", deck))

    def test_kind_inferred_from_name(self):
        assert legacy.kind_from_dat_name("Nellipsoids0.dat") == "ellipsoids"
        assert legacy.kind_from_dat_name("Nphases12.dat") == "phases"
        with pytest.raises(ValueError):
            legacy.kind_from_dat_name("whatever.dat")


class TestDatToJsonConversion:
    """The one-way door: a legacy file in, the JSON we now write out."""

    def test_ellipsoids_round_trip_through_json(self, tmp_path):
        dat = _write_dat(tmp_path, "Nellipsoids0.dat", ELLIPSOIDS_DAT)
        out = legacy.convert_dat_to_json(dat)
        assert out.name == "ellipsoids0.json"

        from_dat = legacy.load_ellipsoids_dat(dat)
        from_json = load_ellipsoids_json(out)
        assert len(from_json) == len(from_dat)
        for a, b in zip(from_dat, from_json):
            assert (a.number, a.umat_name, a.coatingof) == (b.number, b.umat_name, b.coatingof)
            assert (a.a1, a.a2, a.a3) == (b.a1, b.a2, b.a3)
            assert a.concentration == b.concentration
            assert a.geometry_orientation.equals(b.geometry_orientation)
            np.testing.assert_allclose(a.props, b.props)

    def test_layers_conversion_to_an_explicit_path(self, tmp_path):
        dat = _write_dat(tmp_path, "Nlayers0.dat", LAYERS_DAT)
        out = legacy.convert_dat_to_json(dat, tmp_path / "chosen.json")
        assert out.name == "chosen.json"
        layers = load_layers_json(out)
        assert [euler_angles(lay.geometry_orientation)["phi"] for lay in layers] == pytest.approx([-90, -90])

    def test_unknown_kind_is_refused(self, tmp_path):
        dat = _write_dat(tmp_path, "mystery.dat", PHASES_DAT)
        with pytest.raises(ValueError):
            legacy.convert_dat_to_json(dat)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])


# ---------------------------------------------------------------------------
# mean-field models: the legacy props open with [nphases, file number]
# ---------------------------------------------------------------------------

_MATERIAL = ("Material\nName\t{name}\nNumber_of_material_parameters\t{n}\n"
             "Number_of_internal_variables\t1\n\n#Orientation\npsi\t0\ntheta\t0\nphi\t0\n\n#Mechanical\n{props}\n")

_ELLIPSOID_ROW = "{i}\t0\t{name}\t1\t{c}\t0\t0\t0\t{a1}\t1\t1\t0\t0\t0\t{n}\t1\t{props}\n"


def _ellipsoids(rows):
    return "header\n" + "".join(_ELLIPSOID_ROW.format(i=i, name=name, c=c, a1=a1, n=len(props),
                                                     props="\t".join(str(p) for p in props))
                                for i, (name, c, a1, props) in enumerate(rows))


@pytest.mark.parametrize("name, legacy_props, expected", [
    ("MIMTN", [2, 0, 20, 20, 0], [20.0, 20.0, 0.0]),
    ("MISCN", [2, 0, 20, 20, 0, 1], [20.0, 20.0, 0.0, 1.0]),
    ("MIPLN", [2, 0], []),
    ("ELISO", [70000.0, 0.3, 1e-5], [70000.0, 0.3, 1e-5]),
])
def test_material_props_lose_the_legacy_mean_field_slots(tmp_path, name, legacy_props, expected):
    body = "\n".join(f"p{i} {v}" for i, v in enumerate(legacy_props))
    (tmp_path / "material.dat").write_text(_MATERIAL.format(name=name, n=len(legacy_props), props=body))
    kw = legacy.material_from_file(str(tmp_path), "material.dat")
    np.testing.assert_array_equal(kw["props"], expected)


def test_nested_mean_field_rows_are_stripped_and_nested(tmp_path):
    """A MIMTN row names its own sub-phase file in props[1]: the converted phase carries
    them, with the 2.0 props, and the result solves."""
    (tmp_path / "Nellipsoids0.dat").write_text(_ellipsoids([
        ("MIMTN", 0.8, 1, [2, 1, 20, 20, 0]), ("ELISO", 0.2, 50, [50000.0, 0.3, 0.0])]))
    (tmp_path / "Nellipsoids1.dat").write_text(_ellipsoids([
        ("ELISO", 0.7, 1, [5000.0, 0.3, 0.0]), ("ELISO", 0.3, 1, [20000.0, 0.3, 0.0])]))
    out = legacy.convert_dat_to_json(tmp_path / "Nellipsoids0.dat", tmp_path / "ellipsoids0.json")
    phases = load_ellipsoids_json(out)
    np.testing.assert_array_equal(phases[0].props, [20.0, 20.0, 0.0])
    assert [q.concentration for q in phases[0].phases] == [0.7, 0.3]
    L = sim.L_eff("MIMTN", [20.0, 20.0, 0.0], 1, phases=phases)
    assert np.all(np.linalg.eigvalsh(0.5 * (L + L.T)) > 0.0)


def test_missing_nested_file_warns(tmp_path):
    (tmp_path / "Nellipsoids0.dat").write_text(_ellipsoids([
        ("MIMTN", 0.8, 1, [2, 7, 20, 20, 0]), ("ELISO", 0.2, 50, [50000.0, 0.3, 0.0])]))
    with pytest.warns(UserWarning, match="Nellipsoids7.dat"):
        phases = legacy.load_ellipsoids_dat(tmp_path / "Nellipsoids0.dat")
    assert phases[0].phases == []


def test_cycled_tabular_block_is_refused(tmp_path):
    step3 = "#Mode\n3\n#File\ntab.txt\n#Dn_init 1.\n#Dn_mini 1.\n#prescribed_mechanical_state\nE\n0 0\n0 0 0\n#T_is_set\n0\n"
    data = write_path_file(tmp_path, _path_text([step3], ncycle=2),
                           extra_files={"tab.txt": "1 0.5 0.001\n2 1.0 0.002\n"})
    with pytest.raises(ValueError, match="Unroll the cycles"):
        legacy.from_file(data, "path.txt")


def test_legacy_solver_settings_are_reported(tmp_path):
    (tmp_path / "solver_control.inp").write_text(
        "div_tnew_dt_solver\n0.5\nmul_tnew_dt_solver\n2\nminiter_solver\n10\nmaxiter_solver\n100\n"
        "inforce_solver\n1\nprecision_solver\n1.E-5\nlambda_solver\n10000.\n")
    (tmp_path / "solver_essentials.inp").write_text("Solver_type_0_Newton_tangent_1_RNL\n1\nCorate_type\n2\n")
    assert legacy.solver_kwargs_of(str(tmp_path)) == {"precision": 1e-5, "solver_type": 1}


def test_a_latin1_label_does_not_hide_the_path(tmp_path):
    flags = ["E"] + ["S"] * 5
    text = _path_text([_step_text(flags, [0.01, 0, 0, 0, 0, 0], ninc=10)]).replace("#time", "#dur\u00e9e")
    data = tmp_path / "data"
    data.mkdir()
    (data / "path.txt").write_bytes(text.encode("latin-1"))
    assert [Path(w).name for w in legacy.convert_to_json(str(data))] == ["path.json"]
