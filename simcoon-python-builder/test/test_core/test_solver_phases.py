"""The solver driving a mean-field model whose sub-phases are given in memory.

umat_multi used to read Nellipsoids<N>.dat at the first increment, from a "data" directory
relative to the working directory. The phases now travel with the call.

The cross-check that matters: a uniaxial tension with stress-free lateral faces must return
the effective Young's modulus that L_eff computes on the same microstructure, although the
two go through entirely different code (incremental solver vs direct homogenization).
"""

from pathlib import Path

import numpy as np
import pytest

import simcoon as sim
from simcoon import solver as slv
from simcoon.solver.micromechanics import (Cylinder, Ellipsoid, Layer, Phase,
                                           load_ellipsoids_json, to_phase_dicts)

#: the historical reference case of the C++ test TMIMTN: a two-level composite
MIMTN_CASE = Path(__file__).resolve().parents[3] / "testBin" / "Umats" / "MIMTN"

# [nphases, unused (was the Nellipsoids file number), mp, np, index of the matrix phase]
MIMTN_PROPS = np.array([2.0, 0.0, 20.0, 20.0, 0.0])
NSTATEV = 10000

# Glass/epoxy of the shipped examples: 80 % matrix, 20 % spherical reinforcement.
E_EFF = 3293.160186


def two_phase_composite():
    matrix = Ellipsoid(number=0, umat_name="ELISO", save=1, concentration=0.8, nstatev=1,
                       props=np.array([2250.0, 0.19, 8.8e-5]))
    reinforcement = Ellipsoid(number=1, umat_name="ELISO", save=1, concentration=0.2, nstatev=1,
                              props=np.array([73000.0, 0.19, 0.5e-6]))
    return [matrix, reinforcement]


def uniaxial_step(eps=0.01, ninc=5):
    return slv.StepMeca(control=["strain"] + ["stress"] * 5,
                        value=[eps, 0.0, 0.0, 0.0, 0.0, 0.0],
                        time=1.0, ninc=ninc)


class TestSolverWithInMemoryPhases:

    def test_apparent_modulus_matches_L_eff(self):
        phases = two_phase_composite()
        res = slv.solve(uniaxial_step(), "MIMTN", MIMTN_PROPS, NSTATEV, phases=phases)
        stress, strain = res["Stress"], res["Strain"]
        apparent = stress[0, -1] / strain[0, -1]
        assert apparent == pytest.approx(E_EFF, rel=1e-6)

        direct = sim.L_iso_props(
            sim.L_eff("MIMTN", MIMTN_PROPS, NSTATEV, phases=phases)
        ).flatten()[0]
        assert apparent == pytest.approx(direct, rel=1e-9)

    def test_lateral_faces_stay_stress_free(self):
        res = slv.solve(uniaxial_step(), "MIMTN", MIMTN_PROPS, NSTATEV,
                        phases=two_phase_composite())
        assert np.abs(res["Stress"][1:, -1]).max() < 1e-9

    def test_dicts_are_accepted_as_well_as_dataclasses(self):
        phases = two_phase_composite()
        from_objects = slv.solve(uniaxial_step(), "MIMTN", MIMTN_PROPS, NSTATEV, phases=phases)
        from_dicts = slv.solve(uniaxial_step(), "MIMTN", MIMTN_PROPS, NSTATEV,
                               phases=to_phase_dicts(phases))
        np.testing.assert_allclose(from_objects["Stress"], from_dicts["Stress"], rtol=1e-12)

    def test_mean_field_model_without_phases_is_refused(self):
        with pytest.raises(Exception, match="phases"):
            slv.solve(uniaxial_step(), "MIMTN", MIMTN_PROPS, NSTATEV)

    def test_phase_count_must_match_props(self):
        #one phase at 100 %: the concentration check passes, the count check must not
        single = two_phase_composite()[:1]
        single[0].concentration = 1.0
        with pytest.raises(Exception, match="sub-phases props\\[0\\] announces"):
            slv.solve(uniaxial_step(), "MIMTN", MIMTN_PROPS, NSTATEV,
                      phases=single)

    def test_single_phase_model_is_untouched(self):
        # ELISO takes no sub-phases: the new argument must not disturb it.
        res = slv.solve(uniaxial_step(), "ELISO", np.array([70000.0, 0.3, 1.0e-5]), 1)
        apparent = res["Stress"][0, -1] / res["Strain"][0, -1]
        assert apparent == pytest.approx(70000.0, rel=1e-9)

    def test_phases_given_to_a_single_phase_model_are_refused(self):
        with pytest.raises(Exception, match="homogeneous"):
            slv.solve(uniaxial_step(), "ELISO", np.array([70000.0, 0.3, 1.0e-5]), 1,
                      phases=two_phase_composite())


class TestPhaseChecks:
    """The binding refuses what the schemes would otherwise homogenise silently."""

    def test_geometry_kind_must_match_the_model(self):
        for wrong in (Layer, Cylinder, Phase):
            phases = [wrong(number=0, concentration=0.8, props=[2250.0, 0.19, 0.0]),
                      wrong(number=1, concentration=0.2, props=[73000.0, 0.19, 0.0])]
            with pytest.raises(Exception, match="builds ellipsoids"):
                slv.solve(uniaxial_step(), "MIMTN", MIMTN_PROPS, NSTATEV, phases=phases)

    def test_concentration_is_required(self):
        dicts = to_phase_dicts(two_phase_composite())
        del dicts[1]["concentration"]
        with pytest.raises(Exception, match="concentration"):
            slv.solve(uniaxial_step(), "MIMTN", MIMTN_PROPS, NSTATEV, phases=dicts)

    def test_concentrations_must_sum_to_one(self):
        phases = two_phase_composite()
        phases[1].concentration = 0.1
        with pytest.raises(Exception, match="sum to"):
            slv.solve(uniaxial_step(), "MIMTN", MIMTN_PROPS, NSTATEV, phases=phases)


class TestNestedComposite:
    """A sub-phase that is itself a mean-field model, as Nellipsoids0.dat -> Nellipsoids1.dat
    chained through props[1] in the file era."""

    def test_reference_case_of_the_cpp_test(self):
        """TMIMTN's fixture, driven from Python through its JSON files (converted from
        the .dat / .txt the C++ test still reads), against its committed reference."""
        data = MIMTN_CASE / "data"
        outer = load_ellipsoids_json(data / "ellipsoids0.json")
        outer[0].phases = load_ellipsoids_json(data / "ellipsoids1.json")
        assert outer[0].umat_name == "MIMTN"

        kwargs = slv.load_simulation_json(data / "material.json", data / "path.json")
        res = slv.solve(phases=outer, **kwargs)

        #the reference carries 6 significant digits
        ref = np.loadtxt(MIMTN_CASE / "comparison" / "results_job_global-0.txt")
        assert len(res) == ref.shape[0]
        np.testing.assert_allclose(res["Strain"].T, ref[:, 8:14], rtol=1e-5, atol=1e-9)
        np.testing.assert_allclose(res["Stress"].T, ref[:, 14:20], rtol=1e-5, atol=1e-6)

    def test_inner_phases_are_required(self):
        outer = load_ellipsoids_json(MIMTN_CASE / "data" / "ellipsoids0.json")
        with pytest.raises(Exception, match="MIMTN is a mean-field model"):
            slv.solve(uniaxial_step(), "MIMTN", MIMTN_PROPS, NSTATEV, phases=outer)

    def test_homogeneous_phase_takes_no_inner_phases(self):
        phases = two_phase_composite()
        phases[1].phases = two_phase_composite()
        with pytest.raises(Exception, match="homogeneous"):
            slv.solve(uniaxial_step(), "MIMTN", MIMTN_PROPS, NSTATEV, phases=phases)
