"""sim.L_eff with its sub-phases given in memory, instead of read from Nellipsoids<N>.dat.

The reference values were produced by the previous, file-based binding on the shipped
fixture testBin/Umats/MIMTN/data/Nellipsoids1.dat (Mori-Tanaka, ELISO matrix at 80 %,
ELISO fibre at 20 % with a1 = 50 and a 45 deg geometry angle). Both paths agreed to the
last bit, so these numbers pin the migration.
"""

import numpy as np
import pytest

import simcoon as sim
from simcoon.solver.micromechanics import (
    Ellipsoid,
    GeometryOrientation,
    to_phase_dicts,
)

# inspect.signature() raises on a pybind11 builtin; the signature lives in the docstring.
pytestmark = pytest.mark.skipif(
    "phases" not in (sim._core.L_eff.__doc__ or ""),
    reason="the built _core predates the in-memory phases argument",
)

# props of the RVE: [nphases, file number (ignored now), mp, np, index of the matrix phase]
MIMTN_PROPS = np.array([2.0, 1.0, 20.0, 20.0, 0.0])
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
                      geometry_orientation=GeometryOrientation(psi=45.0, theta=0.0, phi=0.0))
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
        aligned[1].geometry_orientation = GeometryOrientation()
        L_45 = np.asarray(sim._core.L_eff("MIMTN", MIMTN_PROPS, NSTATEV,
                                          phases=to_phase_dicts(two_phase_composite())))
        L_0 = np.asarray(sim._core.L_eff("MIMTN", MIMTN_PROPS, NSTATEV,
                                         phases=to_phase_dicts(aligned)))
        assert np.abs(L_45 - L_0).max() > 1.0

    def test_mean_field_model_without_phases_is_refused(self):
        with pytest.raises(Exception, match="phases"):
            sim._core.L_eff("MIMTN", MIMTN_PROPS, NSTATEV)

    def test_phase_count_must_match_props(self):
        with pytest.raises(Exception, match="announces"):
            sim._core.L_eff("MIMTN", MIMTN_PROPS, NSTATEV,
                            phases=to_phase_dicts(two_phase_composite()[:1]))

    def test_homogeneous_model_needs_no_phases(self):
        L = np.asarray(sim._core.L_eff("ELISO", np.array([70000.0, 0.3, 1.0e-5]), 1))
        np.testing.assert_allclose(L, np.asarray(sim.L_iso([70000.0, 0.3], "Enu")),
                                   rtol=1e-12, atol=1e-9)
