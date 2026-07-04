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

///@file Ttangent_SMAUT.cpp
///@brief Tier-D acceptance: SMA unified_T (SMADI) under tangent_mode 2 — the CPP
///       integrator's consistent tangent is the exact FD Jacobian of the discrete map
///       (NON-associated mechanisms: no symmetry assertion), and the converged response
///       stays consistent with the legacy CCP integrator at moderate steps.
///@version 1.0

#include <gtest/gtest.h>
#include <armadillo>
#include <vector>
#include <cstdlib>

#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/SMA/unified_T.hpp>

using namespace std;
using namespace arma;
using namespace simcoon;

namespace {

// Canonical superelastic NiTi (SMADI layout: flagT, EA, EM, nuA, nuM, then common params).
// T = 320 K > Af0 = 280 K -> stress-induced, fully recoverable transformation.
struct Out { vec sigma; mat Lt; double xi; double tnew_dt; };

Out run_smadi(const vec &DEtot, int tangent_mode) {
    vec props = {0., 67538., 67538., 0.349, 0.349,
                 1.e-6, 1.e-6,             // alphaA, alphaM
                 0., 0.0418, 0.021, 0.,    // Hmin, Hmax, k1, sigmacrit
                 10., 10., 250., 230., 260., 280.,  // C_A, C_M, Ms0, Mf0, As0, Af0
                 0.2, 0.2, 0.2, 0.2,       // n1..n4
                 300., 1.4, 2.,            // sigmacaliber, prager_b, prager_n
                 1.e-6, 1.e-3, 1., 1.e8};  // lagrange
    const int nstatev = 17;

    vec Etot = zeros(6), statev = zeros(nstatev), sigma = zeros(6);
    mat Lt = zeros(6, 6), L = zeros(6, 6), DR = eye(3, 3);
    double T = 320., DT = 0., Time = 0., DTime = 1.;
    double Wm = 0., Wm_r = 0., Wm_ir = 0., Wm_d = 0., tnew_dt = 1.;

    umat_sma_unified_T("SMADI", Etot, DEtot, sigma, Lt, L, DR, int(props.n_elem), props,
                       nstatev, statev, T, DT, Time, DTime,
                       Wm, Wm_r, Wm_ir, Wm_d, 3, 3, true, tnew_dt, tangent_mode);
    return {sigma, Lt, statev(1), tnew_dt};
}

bool verbose() { return std::getenv("SIMCOON_TEST_VERBOSE") != nullptr; }

} // namespace

// ---------------------------------------------------------------------
// Test 1: elastic parity — a sub-transformation step returns the identical
// trial state across integrators.
// ---------------------------------------------------------------------
TEST(Ttangent_SMAUT, parity_elastic_step)
{
    const vec Deps = {5.e-4, -1.7e-4, -1.7e-4, 0., 0., 0.};
    Out c = run_smadi(Deps, 0);
    Out a = run_smadi(Deps, 2);
    ASSERT_LT(c.xi, 0.011);  // start value xi = limit; no transformation fired
    EXPECT_LT(norm(a.sigma - c.sigma, 2), 1e-10 * (norm(c.sigma, 2) + 1.));
}

// ---------------------------------------------------------------------
// Test 2: mode-2 consistent tangent == exact FD Jacobian of the CPP map
// (non-associated transformation mechanisms: NO symmetry expectation).
// ---------------------------------------------------------------------
TEST(Ttangent_SMAUT, mode2_exact_consistent_tangent)
{
    const vec Deps = {1.5e-2, 0., 0., 3.e-3, 0., 0.};  // fires forward transformation
    Out a = run_smadi(Deps, 2);
    ASSERT_EQ(a.tnew_dt, 1.);
    ASSERT_GT(a.xi, 0.02) << "transformation did not fire";

    const double h = 1.e-6;
    mat J_ref(6, 6);
    for (int j = 0; j < 6; ++j) {
        vec Dp = Deps, Dm = Deps;
        Dp(j) += h;
        Dm(j) -= h;
        J_ref.col(j) = (run_smadi(Dp, 2).sigma - run_smadi(Dm, 2).sigma) / (2. * h);
    }
    const double err = norm(a.Lt - J_ref, "fro");
    const double err_rel = err / norm(J_ref, "fro");
    if (verbose()) {
        std::cout << "\n[Ttangent_SMAUT mode2] ||Lt - J_FD||_F = " << err
                  << "  (rel " << err_rel << ")  xi = " << a.xi << "\n";
    }
    EXPECT_LT(err_rel, 1e-3) << "CPP consistent tangent is not the discrete Jacobian";
}

// ---------------------------------------------------------------------
// Test 3: CCP-vs-CPP consistency — converged stress and martensite fraction
// agree at a moderate transforming step (integrator difference is O(step^2)).
// ---------------------------------------------------------------------
TEST(Ttangent_SMAUT, ccp_cpp_consistency)
{
    const vec Deps = {1.5e-2, 0., 0., 3.e-3, 0., 0.};
    Out c = run_smadi(Deps, 0);
    Out a = run_smadi(Deps, 2);
    ASSERT_EQ(a.tnew_dt, 1.);
    if (verbose()) {
        std::cout << "\n[Ttangent_SMAUT ccp-cpp] |dsigma| = " << norm(a.sigma - c.sigma, 2)
                  << " (" << norm(a.sigma - c.sigma, 2) / norm(c.sigma, 2) * 100. << "%)"
                  << "  xi_ccp = " << c.xi << "  xi_cpp = " << a.xi << "\n";
    }
    EXPECT_LT(norm(a.sigma - c.sigma, 2), 0.05 * norm(c.sigma, 2));
    EXPECT_NEAR(a.xi, c.xi, 0.05 * std::max(c.xi, 1e-3));
}
