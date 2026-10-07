/* This file is part of simcoon.

simcoon is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

simcoon is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with simcoon.  If not, see <http://www.gnu.org/licenses/>.

*/

///@file Ttutorial_umat.cpp
///@brief Validation of the tutorial J2 UMAT (examples/umat_tutorial/).
///
/// Three checks, mirroring the tutorial's claims:
///  1. Equivalence with the kept EPICP reference (m = 1 makes the legacy
///     power law R = k p^m exactly the tutorial's linear hardening H = k).
///  2. The hand-built analytical continuum tangent equals the shared
///     assemble_continuum_tangent operator AND its closed-form
///     specialization L - 4mu^2/(3mu+H) Lambda o Lambda.
///  3. The algorithmic tangent equals the finite-difference Jacobian of the
///     discrete stress update (exact for J2 radial return).

#include <gtest/gtest.h>
#include <armadillo>

#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Continuum_mechanics/Functions/criteria.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Continuum_mechanics/Umat/tangent_assembly.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Plasticity/plastic_isotropic_ccp.hpp>

#include "umat_tutorial_J2.hpp"

using namespace arma;
using namespace simcoon;

namespace {

const vec PROPS_TUTORIAL = {210000., 0.3, 1.2e-5, 300., 10000.};        // E nu alpha sY H
const vec PROPS_EPICP = {210000., 0.3, 1.2e-5, 300., 10000., 1.0};      // ... k m=1

struct Out {
    vec sigma = zeros(6);
    mat Lt = zeros(6, 6);
    vec statev;
};

// Drive a kernel over a strain increment sequence; return the final state.
template <typename Fn>
Out run(Fn&& umat, const vec& props, int nstatev,
        const std::vector<vec>& increments, int tangent_mode) {
    Out o;
    o.statev = zeros(nstatev);
    vec Etot = zeros(6);
    mat L(6, 6);
    const mat DR = eye(3, 3);
    double T = 293., DT = 0., Time = 0., DTime = 1.;
    double Wm = 0., Wm_r = 0., Wm_ir = 0., Wm_d = 0.;
    double tnew_dt = 1.;
    bool start = true;
    for (const vec& DEtot : increments) {
        umat("TUTO", Etot, DEtot, o.sigma, o.Lt, L, DR,
             static_cast<int>(props.n_elem), props, nstatev, o.statev,
             T, DT, Time, DTime, Wm, Wm_r, Wm_ir, Wm_d,
             3, 3, start, tnew_dt, tangent_mode);
        Etot += DEtot;
        Time += DTime;
        start = false;
    }
    return o;
}

std::vector<vec> cyclic_path(double amp, int n_quarter) {
    std::vector<vec> incs;
    const double de = amp / n_quarter;
    for (int i = 0; i < 4 * n_quarter; ++i) {
        vec d = zeros(6);
        d(0) = (i < n_quarter) ? de : (i < 3 * n_quarter) ? -de : de;
        d(3) = 0.3 * d(0);  // add shear so the flow direction is non-trivial
        incs.push_back(d);
    }
    return incs;
}

}  // namespace

TEST(Ttutorial_umat, matches_EPICP_reference) {
    const auto path = cyclic_path(0.02, 25);
    const Out tut = run(umat_tutorial_J2, PROPS_TUTORIAL, 8, path,
                        tangent_continuum);
    const Out ref = run(umat_plasticity_iso_CCP, PROPS_EPICP, 8, path,
                        tangent_continuum);
    const double peak = std::max(std::abs(ref.sigma.max()),
                                 std::abs(ref.sigma.min()));
    ASSERT_GT(peak, 300.);
    // EPICP iterates its CCP loop to 1e-9; the tutorial's radial return is
    // exact — they agree to the CCP tolerance.
    EXPECT_LT(norm(tut.sigma - ref.sigma, 2) / peak, 1e-8);
    EXPECT_NEAR(tut.statev(1), ref.statev(1), 1e-8);  // accumulated p
}

TEST(Ttutorial_umat, analytic_continuum_tangent_closed_form) {
    // Single plastic step; compare the UMAT's tangent to the closed form
    // L - 4 mu^2/(3 mu + H) Lambda o Lambda built independently here.
    const double E = PROPS_TUTORIAL(0), nu = PROPS_TUTORIAL(1);
    const double H = PROPS_TUTORIAL(4);
    const double mu = E / (2. * (1. + nu));

    vec DE = {0.004, -0.001, -0.001, 0.001, 0., 0.};
    const Out o = run(umat_tutorial_J2, PROPS_TUTORIAL, 8, {DE},
                      tangent_continuum);

    const vec Lam = eta_stress(o.sigma);  // radial return: same direction
    const mat L = L_iso(E, nu, "Enu");
    // CONVENTION LESSON: the closed form L - 4mu^2/(3mu+H) Lambda o Lambda
    // holds for the TENSORIAL Lambda. In Voigt components the outer product
    // must use the STRESS-Voigt normal (shear NOT doubled): kappa = L:Lambda
    // = 2 mu Lambda lands in stress Voigt, and eta_stress() returns the
    // ENGINEERING (doubled-shear) normal — using it directly overweights the
    // shear block of the rank-1 correction by up to 4x.
    vec Lam_sv = Lam;
    Lam_sv.rows(3, 5) *= 0.5;  // engineering -> stress Voigt
    const mat Lt_closed =
        L - (4. * mu * mu / (3. * mu + H)) * (Lam_sv * Lam_sv.t());

    EXPECT_LT(norm(o.Lt - Lt_closed, "fro") / norm(L, "fro"), 1e-9)
        << "hand-built tangent != closed form";

    // ... and to the shared continuum assembly used by production kernels.
    const vec kappa = L * Lam;
    const ContinuumTangent ct = assemble_continuum_tangent(
        dot(Lam, kappa) + H, kappa, Lam, 1.0, L);
    EXPECT_LT(norm(o.Lt - ct.Lt, "fro") / norm(L, "fro"), 1e-9)
        << "hand-built tangent != assemble_continuum_tangent";
}

TEST(Ttutorial_umat, algorithmic_tangent_matches_finite_difference) {
    // Build a plastic state, then FD the discrete map sigma(Etot + DEtot)
    // around a further plastic increment; the algorithmic tangent must be
    // its exact Jacobian (J2 radial return: CCP == CPP).
    vec DE1 = {0.003, -0.0009, -0.0009, 0.0005, 0., 0.};
    vec DE2 = {0.001, -0.0003, -0.0003, 0.0002, 0., 0.};

    const Out base = run(umat_tutorial_J2, PROPS_TUTORIAL, 8, {DE1, DE2},
                         tangent_algorithmic);

    const double h = 1e-7;
    mat J(6, 6);
    for (int j = 0; j < 6; ++j) {
        vec dp = DE2, dm = DE2;
        dp(j) += h;
        dm(j) -= h;
        const Out op = run(umat_tutorial_J2, PROPS_TUTORIAL, 8, {DE1, dp},
                           tangent_algorithmic);
        const Out om = run(umat_tutorial_J2, PROPS_TUTORIAL, 8, {DE1, dm},
                           tangent_algorithmic);
        J.col(j) = (op.sigma - om.sigma) / (2. * h);
    }
    EXPECT_LT(norm(base.Lt - J, "fro") / norm(J, "fro"), 1e-5)
        << "algorithmic tangent is not the discrete Jacobian";
}
