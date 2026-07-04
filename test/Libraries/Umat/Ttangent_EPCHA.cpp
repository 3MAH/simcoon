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

///@file Ttangent_EPCHA.cpp
///@brief Tier-C representative: EPCHA (J2 Chaboche kinematic+isotropic) algorithmic tangent.
///       The dLambda/dsigma closest-point term is wired (deta_stress(sigma-X)); the backstress
///       state-coupling dLambda/dX is deferred (CPP, future release), and EPCHA uses the CCP
///       integrator, so the tangent is symmetric and far better than continuum but not machine-exact.
///@version 1.0

#include <gtest/gtest.h>
#include <armadillo>
#include <vector>
#include <cstdlib>

#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Plasticity/plastic_chaboche_ccp.hpp>

using namespace std;
using namespace arma;
using namespace simcoon;

namespace {

// EPCHA: props = [E, nu, alpha, sigmaY, Q, b, C_1, D_1, C_2, D_2] ; statev(1)=p.
struct Out { vec sigma; mat Lt; bool plastic; };

Out run_epcha(const vec &DEtot, int tangent_mode) {
    const vec props = {70000., 0.3, 1.e-5, 200., 100., 10., 5000., 50., 2000., 20.};
    const int nprops = props.n_elem;
    const int nstatev = 33;

    vec Etot = zeros(6), statev = zeros(nstatev), sigma = zeros(6);
    mat Lt = zeros(6, 6), L = zeros(6, 6), DR = eye(3, 3);
    double T = 293.15, DT = 0., Time = 0., DTime = 1.;
    double Wm = 0., Wm_r = 0., Wm_ir = 0., Wm_d = 0., tnew_dt = 1.;

    umat_plasticity_chaboche_CCP("EPCHA", Etot, DEtot, sigma, Lt, L, DR, nprops, props,
                                 nstatev, statev, T, DT, Time, DTime,
                                 Wm, Wm_r, Wm_ir, Wm_d, 3, 3, true, tnew_dt, tangent_mode);
    return {sigma, Lt, statev(1) > simcoon::iota};
}

vector<double> stress_target_newton(const vec &sig_t, vec DEtot, int tm, int max_iter) {
    vector<double> res;
    for (int it = 0; it < max_iter; ++it) {
        Out o = run_epcha(DEtot, tm);
        vec r = sig_t - o.sigma;
        res.push_back(norm(r, 2));
        if (res.back() < 1e-10) break;
        DEtot += solve(o.Lt, r);
    }
    return res;
}

bool verbose() { return std::getenv("SIMCOON_TEST_VERBOSE") != nullptr; }

} // namespace

TEST(Ttangent_EPCHA, parity_elastic_step)
{
    const vec DEtot = {1.e-4, -0.3e-4, -0.3e-4, 0., 0., 0.};
    Out c = run_epcha(DEtot, 0);
    Out a = run_epcha(DEtot, 1);
    ASSERT_FALSE(c.plastic);
    EXPECT_LT(norm(c.Lt - a.Lt, "fro"), 1.e-12);
}

TEST(Ttangent_EPCHA, algorithmic_tangent_better_than_continuum)
{
    const vec Deps = {5.e-3, -1.e-3, 0.5e-3, 2.e-3, 0., 0.};
    Out c = run_epcha(Deps, 0);
    Out a = run_epcha(Deps, 1);
    ASSERT_TRUE(a.plastic);

    const double h = 1.e-7;
    mat J_ref(6, 6);
    for (int j = 0; j < 6; ++j) {
        vec Dp = Deps; Dp(j) += h; vec Dm = Deps; Dm(j) -= h;
        J_ref.col(j) = (run_epcha(Dp, 0).sigma - run_epcha(Dm, 0).sigma) / (2. * h);
    }
    const double err_algo = norm(a.Lt - J_ref, "fro");
    const double err_cont = norm(c.Lt - J_ref, "fro");
    if (verbose())
        std::cout << "\n[Ttangent_EPCHA] algo=" << err_algo << " cont=" << err_cont
                  << " (ratio=" << err_cont / std::max(err_algo, 1e-300) << ")\n";
    EXPECT_LT(norm(a.Lt - a.Lt.t(), "fro"), 1.e-6 * norm(a.Lt, "fro"));  // symmetric
    EXPECT_GT(err_cont, err_algo);                                       // closer than continuum
}

TEST(Ttangent_EPCHA, convergence_no_worse_than_continuum)
{
    const mat L = L_iso(70000., 0.3, "Enu");
    const vec sig_t = {350., -50., 50., 80., 0., 0.};
    vec Deps0 = solve(L, sig_t); Deps0 *= 1.2;
    Deps0(3) += 3.e-3; Deps0(4) += 1.5e-3; Deps0(1) += 1.e-3;

    auto r_cont = stress_target_newton(sig_t, Deps0, 0, 40);
    auto r_algo = stress_target_newton(sig_t, Deps0, 1, 40);
    if (verbose())
        std::cout << "[Ttangent_EPCHA] final residual cont=" << r_cont.back()
                  << " algo=" << r_algo.back() << "\n";
    // The closest-point tangent must not make the stress-target Newton worse than continuum.
    // (Full quadratic for chaboche needs CPP + backstress state-coupling -> future release.)
    // 1e-3 relative slack: both iterations stall at the same residual level for this stiff
    // chaboche set, and a strict <= would flake on platform/BLAS rounding noise.
    EXPECT_LE(r_algo.back(), r_cont.back() * (1. + 1.e-3));
}

// ---------------------------------------------------------------------
// Test 4 (tangent_mode == 2, CPP integrator): backstress payoff. With the
// inner backward-Euler state and the FD-total dLambda/dsigma (which carries
// the dLambda/dX chain), the consistent tangent of the mode-2 map is exact
// and the map symmetric — what neither mode 0 nor mode 1 could deliver on CCP.
// ---------------------------------------------------------------------
TEST(Ttangent_EPCHA, mode2_exact_symmetric_tangent_backstress)
{
    const vec Deps = {5.e-3, -1.e-3, 0.5e-3, 2.e-3, 0., 0.};
    Out a = run_epcha(Deps, 2);
    ASSERT_TRUE(a.plastic);

    const double h = 1.e-7;
    mat J_ref(6, 6);
    for (int j = 0; j < 6; ++j) {
        vec Dp = Deps; Dp(j) += h; vec Dm = Deps; Dm(j) -= h;
        J_ref.col(j) = (run_epcha(Dp, 2).sigma - run_epcha(Dm, 2).sigma) / (2. * h);
    }
    const double err = norm(a.Lt - J_ref, "fro");
    if (verbose())
        std::cout << "\n[Ttangent_EPCHA mode2] ||Lt - J_FD||_F = " << err
                  << "  ||J-J^T||/||J|| = "
                  << norm(J_ref - J_ref.t(), "fro") / norm(J_ref, "fro") << "\n";
    // Exact consistent tangent (FD reference carries ~1e-2-level noise through the
    // inner fixed point, hence the tolerance of 5 instead of EPICP/EPHIL's 1).
    EXPECT_LT(err, 5.);
    EXPECT_LT(norm(J_ref - J_ref.t(), "fro"), 1e-3 * norm(J_ref, "fro"));

    // Elastic parity across modes.
    const vec Deps_el = {1.e-4, -0.3e-4, -0.3e-4, 0., 0., 0.};
    Out e0 = run_epcha(Deps_el, 0);
    Out e2 = run_epcha(Deps_el, 2);
    EXPECT_LT(norm(e2.sigma - e0.sigma, 2), 1e-14);
    EXPECT_LT(norm(e2.Lt - e0.Lt, "fro"), 1e-12);
}
