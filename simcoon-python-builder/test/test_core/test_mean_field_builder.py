"""The mean-field schemes by name (MoriTanaka, SelfConsistent, Homogeneous, PeriodicLayers): the
builder only names the props slots the C++ reads (multiphase.hpp), so every result is the one of
the hand-written list, bit for bit."""

import re

import numpy as np
import pytest

import simcoon as sim
from simcoon import solver as slv
from simcoon.solver.micromechanics import (
    Ellipsoid,
    Homogeneous,
    Layer,
    MeanFieldScheme,
    MoriTanaka,
    PeriodicLayers,
    SelfConsistent,
    save_ellipsoids_json,
    load_ellipsoids_json,
    to_phase_dict,
)


def eliso(n, c, E, nu, **kw):
    return Ellipsoid(number=n, umat_name="ELISO", save=1, concentration=c, nstatev=1,
                     props=np.array([E, nu, 1e-5]), **kw)


def composite(c=0.3, a1=1.0):
    return [eliso(0, 1 - c, 2250.0, 0.19), eliso(1, c, 73000.0, 0.19, a1=a1)]


# (scheme, the list a script writes for it)
PAIRS = [
    (Homogeneous(), [20, 20]),
    (Homogeneous(mp=30, np=15), [30, 15]),
    (MoriTanaka(50, 50, 0), [50, 50, 0]),
    (MoriTanaka(20, 20, 0, closed_form=True), [20, 20, 0, 1]),
    (MoriTanaka(50, 50, 0, rule="half_sphere"), [50, 50, 0, 0, 1]),
    (MoriTanaka(50, 50, 1, closed_form=True, rule=1), [50, 50, 1, 1, 1]),
    (SelfConsistent(20, 20, 0), [20, 20, 0]),
    (SelfConsistent(20, 20, -1, start="homogeneous"), [20, 20, -1, 0]),
    (SelfConsistent(20, 20, 0, closed_form=True), [20, 20, 0, 1, 1]),
    (SelfConsistent(20, 20, 0, closed_form=True, anderson=2, rtol=1e-12), [20, 20, 0, 1, 1, 0, 2, 1e-12]),
    (SelfConsistent(20, 20, 0, rtol=1e-9), [20, 20, 0, 1, 0, 0, 0, 1e-9]),
    (PeriodicLayers(), []),
]


@pytest.mark.parametrize("scheme, slots", PAIRS, ids=[repr(p[0]) for p in PAIRS])
def test_props_are_the_hand_written_list(scheme, slots):
    """The trailing slots at their C++ default are omitted, as a script would."""
    np.testing.assert_array_equal(scheme.props, np.asarray(slots, dtype=float))
    assert scheme.props.dtype == float
    assert MeanFieldScheme.from_props(scheme.umat_name, slots) == scheme


def test_schemes_are_immutable_hashable_values():
    a, b = MoriTanaka(50, 50, 0, rule=1), MoriTanaka(50, 50, 0, rule="half_sphere")
    assert a == b and hash(a) == hash(b) and a != MoriTanaka(50, 50, 0)
    assert a.rule == "half_sphere"                       # the code is stored by name
    with pytest.raises(Exception):
        a.mp = 10
    assert len({a, b, SelfConsistent(), PeriodicLayers()}) == 3


@pytest.mark.parametrize("build, message", [
    (lambda: MoriTanaka(0, 20, 0), "mp = 0"),
    (lambda: MoriTanaka(20, 20.5, 0), "np = 20.5"),
    (lambda: MoriTanaka(20, 20, -1), "n_matrix = -1"),
    (lambda: MoriTanaka(20, 20, 0, rule="simpson"), "rule = 'simpson'"),
    (lambda: MoriTanaka(20, 20, 0, closed_form=2), "closed_form = 2"),
    (lambda: SelfConsistent(20, 20, 0, start="homogeneous"), "n_matrix = -1"),
    (lambda: SelfConsistent(20, 20, -1), "n_matrix >= 0"),
    (lambda: SelfConsistent(20, 20, 0, anderson=11), "anderson = 11"),
    (lambda: SelfConsistent(20, 20, 0, rtol=1.0), "rtol = 1.0"),
    (lambda: SelfConsistent(20, 20, 0, start=2), "start = 2"),
    (lambda: MeanFieldScheme.from_props("MIMTN", [20, 20]), "MIMTN takes props"),
    (lambda: MeanFieldScheme.from_props("MIPLN", [1]), "MIPLN takes no props"),
    (lambda: MeanFieldScheme.from_props("ELISO", [1, 2]), "not a mean-field model"),
])
def test_options_are_checked_by_name(build, message):
    with pytest.raises(ValueError, match=re.escape(message)):
        build()


def test_L_eff_takes_the_scheme_in_place_of_the_pair():
    ph = composite(a1=5.0)
    ref = sim.L_eff("MIMTN", [50, 50, 0, 1], 1, phases=ph)
    np.testing.assert_array_equal(sim.L_eff(MoriTanaka(50, 50, 0, closed_form=True), phases=ph), ref)
    np.testing.assert_array_equal(sim.L_eff("MIMTN", MoriTanaka(50, 50, 0, closed_form=True), 1, phases=ph), ref)
    sc = SelfConsistent(20, 20, 0, closed_form=True, anderson=2, rtol=1e-12)
    np.testing.assert_array_equal(sim.L_eff(sc, phases=composite()),
                                  sim.L_eff("MISCN", [20, 20, 0, 1, 1, 0, 2, 1e-12], 0, phases=composite()))
    with pytest.raises(TypeError, match="carries its own props"):
        sim.L_eff(MoriTanaka(), [20, 20, 0], phases=ph)
    with pytest.raises(ValueError, match="does not match"):
        sim.L_eff("MISCN", MoriTanaka(), 1, phases=ph)


def test_solve_takes_the_scheme():
    phases = [eliso(0, 0.7, 70000.0, 0.3), Ellipsoid(number=1, umat_name="EPICP", save=1, concentration=0.3, nstatev=8,
                                                    props=np.array([70000.0, 0.3, 1e-5, 200.0, 1000.0, 0.3]))]
    step = slv.StepMeca(control=["strain"] + ["stress"] * 5, value=[0.01, 0, 0, 0, 0, 0], ninc=5)
    ref = slv.solve(step, "MIMTN", np.array([20.0, 20.0, 0.0, 0.0, 1.0]), 1, phases=phases)
    for args in ((MoriTanaka(20, 20, 0, rule="half_sphere"), None), ("MIMTN", MoriTanaka(20, 20, 0, rule="half_sphere"))):
        res = slv.solve(step, args[0], args[1], 1, phases=phases)
        np.testing.assert_array_equal(res["Stress"], ref["Stress"])
        np.testing.assert_array_equal(res["Strain"], ref["Strain"])


def test_nested_scheme_names_the_sub_phase():
    """A sub-phase that is itself a mean-field model takes its scheme as props (and its name)."""
    inner = Ellipsoid(number=1, concentration=0.3, a1=3.0, props=MoriTanaka(20, 20, 0), phases=composite(a1=3.0))
    assert inner.umat_name == "MIMTN"
    np.testing.assert_array_equal(inner.props, [20.0, 20.0, 0.0])
    with pytest.raises(ValueError, match="does not match"):
        Ellipsoid(number=1, umat_name="MISCN", props=MoriTanaka(20, 20, 0))
    hand = Ellipsoid(number=1, umat_name="MIMTN", concentration=0.3, a1=3.0, props=[20, 20, 0], phases=composite(a1=3.0))
    assert inner == hand
    assert to_phase_dict(inner)["umat_name"] == to_phase_dict(hand)["umat_name"] == "MIMTN"
    two_level = [eliso(0, 0.7, 2250.0, 0.19), inner]
    np.testing.assert_array_equal(sim.L_eff(MoriTanaka(20, 20, 0), phases=two_level),
                                  sim.L_eff("MIMTN", [20, 20, 0], 0, phases=[eliso(0, 0.7, 2250.0, 0.19), hand]))


def test_nested_scheme_round_trips_through_json(tmp_path):
    inner = Ellipsoid(number=1, concentration=0.3, props=SelfConsistent(20, 20, 0, anderson=2), phases=composite())
    save_ellipsoids_json(tmp_path / "phases.json", [eliso(0, 0.7, 2250.0, 0.19), inner])
    back = load_ellipsoids_json(tmp_path / "phases.json")
    assert back[1].umat_name == "MISCN"
    assert MeanFieldScheme.from_props(back[1].umat_name, back[1].props) == SelfConsistent(20, 20, 0, anderson=2)


def test_layers_scheme_has_no_props():
    layers = [Layer(number=0, umat_name="ELISO", save=1, concentration=0.5, nstatev=1, props=[70000.0, 0.3, 0.0]),
              Layer(number=1, umat_name="ELISO", save=1, concentration=0.5, nstatev=1, props=[7000.0, 0.3, 0.0])]
    np.testing.assert_array_equal(sim.L_eff(PeriodicLayers(), phases=layers), sim.L_eff("MIPLN", [], 0, phases=layers))
