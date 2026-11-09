"""Numerical options of the ellipsoidal mean-field schemes (trailing props slots) and the
Fourier weights of discretize_odf. Every option defaults to the historical behaviour.

MIMTN props: [mp, np, n_matrix, closed_form, rule]
MISCN props: [mp, np, n_matrix, start, closed_form, rule, anderson, rtol]
"""

import warnings

import numpy as np
import pytest

import simcoon as sim
from simcoon import solver as slv
from simcoon.solver.micromechanics import Ellipsoid, Peak, discretize_odf


def eliso(n, c, E, nu, **kw):
    return Ellipsoid(number=n, umat_name="ELISO", save=1, concentration=c, nstatev=1,
                     props=np.array([E, nu, 1e-5]), **kw)


def composite(c=0.3, a1=1.0, E_m=2250.0):
    return [eliso(0, 1 - c, E_m, 0.19), eliso(1, c, E_m * 73000.0 / 2250.0, 0.19, a1=a1)]


def L_eff(scheme, props, phases):
    return np.asarray(sim.L_eff(scheme, np.asarray(props, dtype=float), 1, phases=phases))


def rel(a, b):
    return np.linalg.norm(a - b) / np.linalg.norm(b)


class TestDefaults:

    @pytest.mark.parametrize("scheme, base, zeros", [("MIMTN", [30, 30, 0], [0, 0]),
                                                    ("MISCN", [30, 30, 0], [1, 0, 0, 0, 0])])
    def test_explicit_zero_options_are_the_default(self, scheme, base, zeros):
        ph = composite(a1=4.0)
        np.testing.assert_array_equal(L_eff(scheme, base + zeros, ph), L_eff(scheme, base, ph))

    @pytest.mark.parametrize("props, message", [
        ([30, 30, 0, 2], "closed_form"), ([30, 30, 0, 0, 0.5], "rule"), ([30, 30, 0, 0, 0, 0], "MIMTN takes"),
    ])
    def test_mori_tanaka_options_are_checked(self, props, message):
        with pytest.raises(ValueError, match=message):
            L_eff("MIMTN", props, composite())

    @pytest.mark.parametrize("props, message", [
        ([30, 30, 0, 1, 0, 0, 11], "anderson"), ([30, 30, 0, 1, 0, 0, 1.5], "anderson"),
        ([30, 30, 0, 1, 0, 0, 2, -1e-8], "rtol"), ([30, 30, 0, 1, 0, 0, 2, 1e-8, 0], "MISCN takes"),
    ])
    def test_self_consistent_options_are_checked(self, props, message):
        with pytest.raises(ValueError, match=message):
            L_eff("MISCN", props, composite())


class TestClosedFormEshelby:

    @pytest.mark.parametrize("a1", [1.0, 3.0, 0.3])
    def test_mori_tanaka_matches_a_fine_quadrature(self, a1):
        ph = composite(a1=a1)
        closed = L_eff("MIMTN", [20, 20, 0, 1], ph)
        fine = L_eff("MIMTN", [200, 200, 0], ph)
        assert rel(closed, fine) < 1e-8

    def test_self_consistent_spheres(self):
        """Spheres in isotropic phases: the self-consistent medium stays isotropic, every
        Eshelby tensor of the fixed point is the closed form."""
        ph = composite(c=0.5)
        closed = L_eff("MISCN", [20, 20, 0, 1, 1], ph)
        assert rel(closed, L_eff("MISCN", [50, 50, 0], ph)) < 1e-10

    def test_anisotropic_medium_is_integrated(self):
        ph = composite(a1=4.0)
        ph[0] = Ellipsoid(number=0, umat_name="ELIST", save=1, concentration=0.7, nstatev=1,
                          props=np.array([3, 20000.0, 3000.0, 0.3, 0.35, 1500.0, 1e-5, 1e-5]))
        np.testing.assert_array_equal(L_eff("MIMTN", [30, 30, 0, 1], ph), L_eff("MIMTN", [30, 30, 0], ph))


class TestHalfSphereRule:

    def test_fibre_against_the_closed_form(self):
        """Isotropic matrix, fibre of aspect ratio 10: the closed form is exact. The half-sphere
        rule evaluates half the points of Gauss-Legendre for an error of the same order, and
        converges."""
        ph = composite(a1=10.0)
        exact = L_eff("MIMTN", [20, 20, 0, 1], ph)
        gauss = rel(L_eff("MIMTN", [50, 50, 0], ph), exact)
        half = rel(L_eff("MIMTN", [50, 50, 0, 0, 1], ph), exact)
        finer = rel(L_eff("MIMTN", [100, 100, 0, 0, 1], ph), exact)
        assert half < 1.5 * gauss
        assert finer < 0.2 * half

    def test_solver_takes_it(self):
        phases = [Ellipsoid(number=0, umat_name="EPICP", save=1, nstatev=8, concentration=0.7,
                            props=np.array([70000.0, 0.3, 1e-5, 200.0, 1000.0, 0.3])),
                  eliso(1, 0.3, 400000.0, 0.2)]
        step = slv.StepMeca(control=["strain"] + ["stress"] * 5, value=[0.02, 0, 0, 0, 0, 0],
                            time=1.0, ninc=20)
        ref = slv.solve(step, "MIMTN", np.array([50.0, 50.0, 0.0]), 1, phases=phases)["Stress"]
        half = slv.solve(step, "MIMTN", np.array([50.0, 50.0, 0.0, 0.0, 1.0]), 1, phases=phases)["Stress"]
        np.testing.assert_allclose(half[0], ref[0], rtol=1e-6)    # lateral stresses: solver residuals


class TestSelfConsistentFixedPoint:

    def test_anderson_reaches_the_same_fixed_point(self):
        ph = composite(c=0.5)
        converged = L_eff("MISCN", [20, 20, 0, 1, 1, 0, 0, 1e-15], ph)   # plain iterations, tight
        anderson = L_eff("MISCN", [20, 20, 0, 1, 1, 0, 2, 1e-12], ph)
        assert rel(anderson, converged) < 1e-11
        assert rel(L_eff("MISCN", [20, 20, 0, 1, 1], ph), converged) < 1e-9

    def test_relative_tolerance_does_not_depend_on_the_unit(self):
        mpa = L_eff("MISCN", [20, 20, 0, 1, 1, 0, 2, 1e-12], composite(c=0.5))
        pa = L_eff("MISCN", [20, 20, 0, 1, 1, 0, 2, 1e-12], composite(c=0.5, E_m=2250.0e6))
        assert rel(pa / 1e6, mpa) < 1e-11


class TestFourierOrientationWeights:
    """Mori-Tanaka, isotropic matrix, fibres swept about z: the average is exact with Fourier
    weights, whatever the number of phases."""

    matrix = Ellipsoid(number=0, umat_name="ELISO", save=1, nstatev=1, concentration=0.7,
                       props=np.array([3000.0, 0.35, 1e-5]))
    fibre = Ellipsoid(number=1, umat_name="ELISO", save=1, nstatev=1, concentration=0.3,
                      props=np.array([72000.0, 0.22, 1e-6]), a1=20.0)
    props = np.array([30.0, 30.0, 0.0])

    def average(self, n, weights, s_dev):
        peaks = [Peak(method=3, mean=47.0, s_dev=s_dev)]
        ph = discretize_odf([self.matrix, self.fibre], 1, peaks, n, axis=(0, 0, 1), weights=weights)
        return np.asarray(sim.L_eff("MIMTN", self.props, 0, phases=ph)), len(ph) - 1

    @pytest.mark.parametrize("s_dev", [40.0, 20.0])  # projected weights, then non-negative ones
    def test_exact_with_few_phases(self, s_dev):
        L6, _ = self.average(6, "fourier", s_dev)
        L36, kept = self.average(36, "fourier", s_dev)
        assert kept <= 36
        assert rel(L6, L36) < 1e-12
        simpson, _ = self.average(720, "simpson", s_dev)
        assert rel(L6, simpson) < 5e-6           # the error of 720 Simpson phases

    def test_too_peaked_density_falls_back_on_simpson(self):
        with pytest.warns(UserWarning, match="Simpson weights used instead"):
            fourier, n = self.average(6, "fourier", 2.0)
        simpson, _ = self.average(6, "simpson", 2.0)
        assert n == 6
        np.testing.assert_array_equal(fourier, simpson)

    def test_arguments_are_checked(self):
        peaks = [Peak(method=3, mean=47.0, s_dev=20.0)]
        with pytest.raises(ValueError, match="half turn"):
            discretize_odf([self.matrix, self.fibre], 1, peaks, 12, angle_range=(0, 360), weights="fourier")
        with pytest.raises(ValueError, match="nphases >= 5"):
            discretize_odf([self.matrix, self.fibre], 1, peaks, 4, weights="fourier")
        with pytest.raises(ValueError, match="'simpson' or 'fourier'"):
            discretize_odf([self.matrix, self.fibre], 1, peaks, 12, weights="gauss")
