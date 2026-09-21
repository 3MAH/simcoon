"""sim.L_eff with its sub-phases given in memory.

The reference values pin a Mori-Tanaka composite (ELISO matrix at 80 %, ELISO fibre at
20 % with a1 = 50 and a 45 deg geometry angle); they were produced by the 1.x file-based
binding, which agreed with this one to the last bit.
"""

import numpy as np
import pytest

import simcoon as sim
from simcoon.solver.micromechanics import (
    Ellipsoid,
    to_phase_dicts,
)

# props of the RVE: [mp, np, index of the matrix phase] (the phase count is the list's)
MIMTN_PROPS = np.array([20.0, 20.0, 0.0])
NSTATEV = 10000

# Values of the file-based binding on Nellipsoids1.dat.
L00 = 10745.836046
L01 = 5402.394823
L33 = 4344.046376


def two_phase_composite():
    matrix = Ellipsoid(number=0, umat_name="ELISO", save=1, concentration=0.8, nstatev=1,
                       props=np.array([5000.0, 0.3, 0.0]), a1=1.0, a2=1.0, a3=1.0)
    fibre = Ellipsoid(number=1, umat_name="ELISO", save=1, concentration=0.2, nstatev=1,
                      props=np.array([50000.0, 0.3, 0.0]), a1=50.0, a2=1.0, a3=1.0,
                      geometry_orientation=(45.0, 0.0, 0.0))
    return [matrix, fibre]


class TestLeffInMemoryPhases:

    def test_matches_the_file_based_reference(self):
        L = np.asarray(sim._core.L_eff("MIMTN", MIMTN_PROPS, NSTATEV,
                                       phases=to_phase_dicts(two_phase_composite())))
        assert L.shape == (6, 6)
        assert L[0, 0] == pytest.approx(L00, rel=1e-9)
        assert L[0, 1] == pytest.approx(L01, rel=1e-9)
        assert L[3, 3] == pytest.approx(L33, rel=1e-9)

    def test_stiffness_is_symmetric(self):
        L = np.asarray(sim._core.L_eff("MIMTN", MIMTN_PROPS, NSTATEV,
                                       phases=to_phase_dicts(two_phase_composite())))
        np.testing.assert_allclose(L, L.T, rtol=1e-12, atol=1e-9)

    def test_geometry_angle_is_read_in_degrees(self):
        # The files stored degrees and the C++ side radians; a fibre at 45 deg must not
        # give the same stiffness as one left at 0.
        aligned = two_phase_composite()
        aligned[1].geometry_orientation = None   # coerced to the identity
        L_45 = np.asarray(sim._core.L_eff("MIMTN", MIMTN_PROPS, NSTATEV,
                                          phases=to_phase_dicts(two_phase_composite())))
        L_0 = np.asarray(sim._core.L_eff("MIMTN", MIMTN_PROPS, NSTATEV,
                                         phases=to_phase_dicts(aligned)))
        assert np.abs(L_45 - L_0).max() > 1.0

    def test_mean_field_model_without_phases_is_refused(self):
        with pytest.raises(Exception, match="phases"):
            sim._core.L_eff("MIMTN", MIMTN_PROPS, NSTATEV)

    def test_matrix_index_must_be_a_phase(self):
        # one phase at 100 %: the concentration check passes, the matrix index must exist
        single = two_phase_composite()[:1]
        single[0].concentration = 1.0
        L = np.asarray(sim._core.L_eff("MIMTN", MIMTN_PROPS, NSTATEV,
                                       phases=to_phase_dicts(single)))
        L_iso = np.asarray(sim.L_iso([5000.0, 0.3], "Enu"))
        np.testing.assert_allclose(L, L_iso, rtol=1e-9)   # a one-phase composite is that phase

    @pytest.mark.parametrize("props, message", [
        ([2.0, 0.0, 20.0, 20.0, 0.0], "got 5 values"),      # the pre-2.0 layout
        ([20.0, 20.0], "got 2 values"),
        ([20.0, 20.0, 2.0], "n_matrix = 2.* is not one of the 2 phases"),
        ([20.0, 20.0, -1.0], "n_matrix = -1"),
        ([20.0, 20.0, 1e300], "n_matrix"),
        ([20.0, 20.0, float("nan")], "n_matrix"),
        ([0.0, 20.0, 0.0], "mp and np"),
        ([1e12, 20.0, 0.0], "mp and np"),
        ([float("nan"), 20.0, 0.0], "mp and np"),
    ])
    def test_props_are_validated(self, props, message):
        with pytest.raises(ValueError, match=message):
            sim.L_eff("MIMTN", props, NSTATEV, phases=two_phase_composite())

    def test_self_consistent_start_option(self):
        phases = two_phase_composite()
        L_default = sim.L_eff("MISCN", [20.0, 20.0, 0.0], NSTATEV, phases=phases)
        L_mt = sim.L_eff("MISCN", [20.0, 20.0, 0.0, 1.0], NSTATEV, phases=phases)
        np.testing.assert_allclose(L_default, L_mt, rtol=1e-12)
        with pytest.raises(ValueError, match="start = 2"):
            sim.L_eff("MISCN", [20.0, 20.0, 0.0, 2.0], NSTATEV, phases=phases)
        with pytest.raises(ValueError, match="takes n_matrix < 0"):
            sim.L_eff("MISCN", [20.0, 20.0, 0.0, 0.0], NSTATEV, phases=phases)

    # Each case below used to crash the interpreter or return a silently wrong stiffness.
    def test_self_referencing_phases_are_refused(self):
        p = to_phase_dicts(two_phase_composite())
        p[1]["umat_name"], p[1]["props"], p[1]["phases"] = "MIMTN", MIMTN_PROPS, p
        with pytest.raises(ValueError, match="16 levels"):
            sim._core.L_eff("MIMTN", MIMTN_PROPS, NSTATEV, None, p)

    @pytest.mark.parametrize("mutate, message", [
        (lambda p: p[1].update(semi_axis=p[1].pop("semi_axes")), "unknown entry 'semi_axis'"),
        (lambda p: p[1].update(geometry_orientaton=p[1].pop("geometry_orientation")), "unknown entry"),
        (lambda p: p[1]["geometry_orientation"].update(psy=45.0), "unknown entry 'psy'"),
        (lambda p: p[1]["semi_axes"].update(a4=1.0), "unknown entry 'a4'"),
        (lambda p: p[0].update(concentration=float("nan")), "concentration"),
        (lambda p: (p[0].update(concentration=1.5), p[1].update(concentration=-0.5)), "concentration"),
        (lambda p: p[1].update(nstatev=-1), "nstatev"),
        (lambda p: p[1].update(props=[[1.0, 2.0], [3.0]]), "not a sequence of numbers"),
    ])
    def test_malformed_phase_dicts_are_refused(self, mutate, message):
        p = to_phase_dicts(two_phase_composite())
        mutate(p)
        with pytest.raises(ValueError, match=message):
            sim._core.L_eff("MIMTN", MIMTN_PROPS, NSTATEV, None, p)

    def test_orientation_keys_and_nstatev_are_checked(self):
        phases = to_phase_dicts(two_phase_composite())
        with pytest.raises(ValueError, match="unknown entry 'psy'"):
            sim._core.L_eff("MIMTN", MIMTN_PROPS, NSTATEV, {"psy": 45.0}, phases)
        with pytest.raises(ValueError, match="nstatev"):
            sim._core.L_eff("MIMTN", MIMTN_PROPS, -1, None, phases)

    def test_non_elastic_phase_is_refused(self):
        phases = two_phase_composite()
        phases[1].umat_name, phases[1].props, phases[1].nstatev = "EPICP", np.array([50000., .3, 0., 300., 1000., .3]), 8
        with pytest.raises(ValueError, match="not a linear elastic model"):
            sim.L_eff("MIMTN", MIMTN_PROPS, NSTATEV, phases=phases)
        with pytest.raises(ValueError, match="not a linear elastic model"):
            sim.L_eff("XXXXX", [1.0], 1)

    def test_homogeneous_model_needs_no_phases(self):
        L = np.asarray(sim._core.L_eff("ELISO", np.array([70000.0, 0.3, 1.0e-5]), 1))
        np.testing.assert_allclose(L, np.asarray(sim.L_iso([70000.0, 0.3], "Enu")),
                                   rtol=1e-12, atol=1e-9)
