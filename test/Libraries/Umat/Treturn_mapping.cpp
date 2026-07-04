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

///@file Treturn_mapping.cpp
///@brief Unit tests for the closest-point projection helper (return_mapping.hpp):
///       closed-form J2 radial return, inactive-mechanism masking, Q-quadratic rate,
///       exact consistent tangent, elastic guard.
///@version 1.0

#include <gtest/gtest.h>
#include <armadillo>
#include <vector>
#include <cmath>

#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Continuum_mechanics/Functions/criteria.hpp>
#include <simcoon/Continuum_mechanics/Umat/return_mapping.hpp>

using namespace std;
using namespace arma;
using namespace simcoon;

namespace {

constexpr double E_mod = 70000.;
constexpr double nu_mod = 0.3;
constexpr double sigmaY = 200.;
constexpr double Hlin = 1000.;   // linear isotropic hardening

// Single-mechanism J2 + linear hardening: everything the helper needs.
struct J2Problem {
    mat L = L_iso(E_mod, nu_mod, "Enu");
    double p_n = 0.;
    double p = 0.;   // caller-owned state

    ReturnMechanism mech() {
        ReturnMechanism m;
        m.Phi = [this](const vec &sig) { return Mises_stress(sig) - sigmaY - Hlin * p; };
        m.dPhi_dsigma = [](const vec &sig) { return eta_stress(sig); };
        m.dLambda_dsigma = [](const vec &sig) { return deta_stress(sig); };
        return m;
    }
    std::function<bool(const vec &, double)> update_state() {
        return [this](const vec &, double Dl) { p = p_n + Dl; return true; };
    }
    std::function<double(const vec &, double)> K_scalar() {
        return [](const vec &, double) { return -Hlin; };
    }
};

} // namespace

// ---------------------------------------------------------------------
// Test 1: closed-form J2 radial return. For J2 + linear hardening the CPP
// solution is analytic: Dlambda = (Mises(sig_tr) - sigmaY - H p_n)/(3 mu + H),
// with the converged deviator colinear to the trial deviator.
// ---------------------------------------------------------------------
TEST(Treturn_mapping, j2_radial_return_closed_form)
{
    J2Problem pb;
    const double mu = E_mod / (2. * (1. + nu_mod));

    const vec Eps = {5.e-3, -1.e-3, 0.5e-3, 2.e-3, 0., 0.};
    const vec sigma_tr = pb.L * Eps;
    ASSERT_GT(Mises_stress(sigma_tr) - sigmaY, 0.) << "test set-up must yield";

    auto r = closest_point_return_mapping(sigma_tr, pb.L, pb.mech(),
                                          pb.update_state(), pb.K_scalar(), sigmaY);
    ASSERT_TRUE(r.converged);

    const double Dl_ref = (Mises_stress(sigma_tr) - sigmaY) / (3. * mu + Hlin);
    EXPECT_NEAR(r.Dlambda(0), Dl_ref, 1e-10 * Dl_ref);

    // Consistency: converged stress sits on the updated yield surface...
    EXPECT_NEAR(Mises_stress(r.sigma), sigmaY + Hlin * r.Dlambda(0), 1e-7);
    // ...and the return is radial: dev(sigma) colinear to dev(sigma_tr).
    const vec d1 = dev(r.sigma), d2 = dev(sigma_tr);
    const double cosang = sum(d1 % d2) / (norm(d1, 2) * norm(d2, 2));
    EXPECT_NEAR(cosang, 1., 1e-10);
    // Pressure untouched by the deviatoric flow.
    EXPECT_NEAR(tr(r.sigma), tr(sigma_tr), 1e-8 * fabs(tr(sigma_tr)));
}

// ---------------------------------------------------------------------
// Test 2: elastic guard — admissible trial returns bit-identical trial state.
// ---------------------------------------------------------------------
TEST(Treturn_mapping, elastic_guard_returns_trial)
{
    J2Problem pb;
    const vec Eps = {1.e-4, -0.3e-4, -0.3e-4, 0., 0., 0.};
    const vec sigma_tr = pb.L * Eps;
    ASSERT_LT(Mises_stress(sigma_tr) - sigmaY, 0.);

    auto r = closest_point_return_mapping(sigma_tr, pb.L, pb.mech(),
                                          pb.update_state(), pb.K_scalar(), sigmaY);
    ASSERT_TRUE(r.converged);
    EXPECT_EQ(r.niter, 0);
    EXPECT_LT(norm(r.sigma - sigma_tr, 2), 1e-14);
    EXPECT_EQ(r.Dlambda(0), 0.);
    // Elastic tangent from the returned pieces (Ds=0 masks all mechanisms -> L).
    ContinuumTangent ct = cpp_consistent_tangent(r, pb.L);
    EXPECT_LT(norm(ct.Lt - pb.L, "fro"), 1e-12);
}

// ---------------------------------------------------------------------
// Test 3: two mechanisms, one inactive. The second (very high yield) criterion
// must stay inactive: Dlambda_2 = 0 and the active solution matches Test 1.
// ---------------------------------------------------------------------
TEST(Treturn_mapping, inactive_mechanism_masked)
{
    J2Problem pb;
    const vec Eps = {5.e-3, -1.e-3, 0.5e-3, 2.e-3, 0., 0.};
    const vec sigma_tr = pb.L * Eps;

    std::vector<ReturnMechanism> mechs(2);
    mechs[0] = pb.mech();
    mechs[1].Phi = [](const vec &sig) { return Mises_stress(sig) - 1.e6; };  // never active
    mechs[1].dPhi_dsigma = [](const vec &sig) { return eta_stress(sig); };
    mechs[1].dLambda_dsigma = [](const vec &sig) { return deta_stress(sig); };

    ReturnStateHooks hooks;
    hooks.update_state = [&pb](const vec &, const vec &Dl) { pb.p = pb.p_n + Dl(0); return true; };
    hooks.K = [](const vec &, const vec &) {
        mat K = zeros(2, 2);
        K(0, 0) = -Hlin;
        return K;
    };
    vec Y_crit = {sigmaY, 1.e6};

    auto r = closest_point_return_mapping(sigma_tr, pb.L, mechs, hooks, Y_crit);
    ASSERT_TRUE(r.converged);
    EXPECT_EQ(r.Dlambda(1), 0.);

    const double mu = E_mod / (2. * (1. + nu_mod));
    const double Dl_ref = (Mises_stress(sigma_tr) - sigmaY) / (3. * mu + Hlin);
    EXPECT_NEAR(r.Dlambda(0), Dl_ref, 1e-9 * Dl_ref);
}

// ---------------------------------------------------------------------
// Test 4: Q-quadratic convergence on an ANISOTROPIC Hill mechanism (the flow
// rotates during the return — the case CCP cannot integrate consistently).
// ---------------------------------------------------------------------
TEST(Treturn_mapping, hill_anisotropic_q_quadratic)
{
    const mat L = L_iso(E_mod, nu_mod, "Enu");
    const vec hill_p = {0.5, 0.6, 0.7, 1.5, 1.4, 1.6};
    double p = 0.;

    ReturnMechanism m;
    m.Phi = [&](const vec &sig) { return Hill_stress(sig, hill_p) - sigmaY - Hlin * p; };
    m.dPhi_dsigma = [&](const vec &sig) { return dHill_stress(sig, hill_p); };
    m.dLambda_dsigma = [&](const vec &sig) { return ddHill_stress(sig, hill_p); };

    const vec Eps = {5.e-3, -1.e-3, 0.5e-3, 2.e-3, 0., 0.};
    const vec sigma_tr = L * Eps;
    ASSERT_GT(Hill_stress(sigma_tr, hill_p) - sigmaY, 0.);

    auto r = closest_point_return_mapping(
        sigma_tr, L, m,
        [&](const vec &, double Dl) { p = Dl; return true; },
        [](const vec &, double) { return -Hlin; }, sigmaY);
    ASSERT_TRUE(r.converged);
    EXPECT_LE(r.niter, 10);

    // Q-quadratic tail: once the combined error is below 1e-2, each iteration
    // must gain more than a superlinear factor.
    const auto &e = r.error_history;
    ASSERT_GE(e.size(), 2u);
    for (size_t k = 0; k + 1 < e.size(); k++) {
        if (e[k] < 1e-2 && e[k + 1] > 0.) {
            EXPECT_LT(e[k + 1], pow(e[k], 1.7))
                << "iteration " << k << ": " << e[k] << " -> " << e[k + 1];
        }
    }

    // Converged state satisfies the CPP system: Phi = 0 and R_sigma = 0.
    EXPECT_NEAR(Hill_stress(r.sigma, hill_p) - sigmaY - Hlin * r.Dlambda(0), 0., 1e-6);
    const vec R = r.sigma - sigma_tr + r.Dlambda(0) * (L * dHill_stress(r.sigma, hill_p));
    EXPECT_LT(norm(R, 2), 1e-6);
}

// ---------------------------------------------------------------------
// Test 5: the consistent tangent of the converged CPP map is the exact FD
// Jacobian d(sigma)/d(Eps) — for the ANISOTROPIC Hill mechanism, where the
// CCP-based algorithmic tangent could not be exact. THE PAYOFF.
// ---------------------------------------------------------------------
TEST(Treturn_mapping, hill_anisotropic_exact_consistent_tangent)
{
    const mat L = L_iso(E_mod, nu_mod, "Enu");
    const vec hill_p = {0.5, 0.6, 0.7, 1.5, 1.4, 1.6};
    double p = 0.;

    auto solve_cpp = [&](const vec &Eps) {
        ReturnMechanism m;
        m.Phi = [&](const vec &sig) { return Hill_stress(sig, hill_p) - sigmaY - Hlin * p; };
        m.dPhi_dsigma = [&](const vec &sig) { return dHill_stress(sig, hill_p); };
        m.dLambda_dsigma = [&](const vec &sig) { return ddHill_stress(sig, hill_p); };
        p = 0.;
        return closest_point_return_mapping(
            L * Eps, L, m,
            [&](const vec &, double Dl) { p = Dl; return true; },
            [](const vec &, double) { return -Hlin; }, sigmaY);
    };

    const vec Eps = {5.e-3, -1.e-3, 0.5e-3, 2.e-3, 0., 0.};
    auto r = solve_cpp(Eps);
    ASSERT_TRUE(r.converged);
    ContinuumTangent ct = cpp_consistent_tangent(r, L);

    const double h = 1.e-7;
    mat J_ref(6, 6);
    for (int j = 0; j < 6; j++) {
        vec Ep = Eps, Em = Eps;
        Ep(j) += h;
        Em(j) -= h;
        J_ref.col(j) = (solve_cpp(Ep).sigma - solve_cpp(Em).sigma) / (2. * h);
    }

    const double err = norm(ct.Lt - J_ref, "fro");
    EXPECT_LT(err, 1.) << "CPP consistent tangent is not the exact discrete Jacobian";
    // And the discrete map itself is symmetric (associated flow + CPP).
    EXPECT_LT(norm(J_ref - J_ref.t(), "fro"), 1e-4 * norm(J_ref, "fro"));
}
