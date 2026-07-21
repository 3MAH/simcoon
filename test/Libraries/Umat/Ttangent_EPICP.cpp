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

///@file Ttangent_EPICP.cpp
///@brief Tier-A acceptance test: the EPICP (J2 isotropic plasticity) algorithmic
///       tangent (algorithmic mode (tangent_algorithmic)) is the exact Jacobian of the discrete return map
///       and converges Q-quadratically, while tangent_continuum stays the continuum tangent.
///       Drives the real umat_plasticity_iso_CCP, toggling tangent_mode.
///@version 1.0

#include <gtest/gtest.h>
#include <armadillo>
#include <vector>
#include <cstdlib>

#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Plasticity/plastic_isotropic_ccp.hpp>

using namespace std;
using namespace arma;
using namespace simcoon;

namespace {

// EPICP material: props = [E, nu, alpha, sigmaY, k, m] ; nstatev = 8.
struct EpicpOut { vec sigma; mat Lt; bool plastic; };

EpicpOut run_epicp(const vec &DEtot, int tangent_mode) {
    const vec props = {70000., 0.3, 1.e-5, 200., 1000., 0.3};
    const int nprops = props.n_elem;
    const int nstatev = 8;

    vec Etot = zeros(6);
    vec statev = zeros(nstatev);
    vec sigma = zeros(6);
    mat Lt = zeros(6, 6), L = zeros(6, 6);
    mat DR = eye(3, 3);
    double T = 293.15, DT = 0.;
    double Time = 0., DTime = 1.;
    double Wm = 0., Wm_r = 0., Wm_ir = 0., Wm_d = 0.;
    double tnew_dt = 1.;

    umat_plasticity_iso_CCP("EPICP", Etot, DEtot, sigma, Lt, L, DR, nprops, props,
                            nstatev, statev, T, DT, Time, DTime,
                            Wm, Wm_r, Wm_ir, Wm_d, 3, 3, true, tnew_dt, tangent_mode);
    // EPICP statev layout: statev(0)=T_init, statev(1)=accumulated plastic strain p.
    return {sigma, Lt, statev(1) > simcoon::iota};
}

// Stress-target outer Newton using the UMAT's Lt as the Jacobian; returns the
// residual history ||sigma_target - sigma||.
vector<double> stress_target_newton(const vec &sigma_target, vec DEtot,
                                    int tangent_mode, int max_iter) {
    vector<double> res;
    for (int it = 0; it < max_iter; ++it) {
        EpicpOut o = run_epicp(DEtot, tangent_mode);
        vec r = sigma_target - o.sigma;
        res.push_back(norm(r, 2));
        if (res.back() < 1e-10) break;
        DEtot += solve(o.Lt, r);
    }
    return res;
}

bool verbose() { return std::getenv("SIMCOON_TEST_VERBOSE") != nullptr; }

} // namespace

// ---------------------------------------------------------------------
// Test 1: PARITY. An elastic (sub-yield) step has no active mechanism, so the
// algorithmic and continuum tangents must be bit-identical.
// ---------------------------------------------------------------------
TEST(Ttangent_EPICP, parity_elastic_step)
{
    const vec DEtot = {1.e-4, -0.3e-4, -0.3e-4, 0., 0., 0.}; // Mises ~ 7 MPa << sigmaY
    EpicpOut c = run_epicp(DEtot, simcoon::tangent_continuum);
    EpicpOut a = run_epicp(DEtot, simcoon::tangent_algorithmic);
    ASSERT_FALSE(c.plastic) << "Parity step unexpectedly yielded.";
    EXPECT_LT(norm(c.Lt - a.Lt, "fro"), 1.e-12);
}

// ---------------------------------------------------------------------
// Test 2: the algorithmic tangent is the exact Jacobian d(sigma)/d(DEtot)
// of the discrete return map (Simo-Hughes), far closer than the continuum one.
// ---------------------------------------------------------------------
TEST(Ttangent_EPICP, algorithmic_tangent_is_exact_jacobian)
{
    const vec Deps = {5.e-3, -1.e-3, 0.5e-3, 2.e-3, 0., 0.}; // non-radial -> yields
    EpicpOut c = run_epicp(Deps, simcoon::tangent_continuum);
    EpicpOut a = run_epicp(Deps, simcoon::tangent_algorithmic);
    ASSERT_TRUE(a.plastic) << "Test set-up did not yield - increase Deps.";

    const double h = 1.e-7;
    mat J_ref(6, 6);
    for (int j = 0; j < 6; ++j) {
        vec Dp = Deps; Dp(j) += h;
        vec Dm = Deps; Dm(j) -= h;
        J_ref.col(j) = (run_epicp(Dp, simcoon::tangent_continuum).sigma - run_epicp(Dm, simcoon::tangent_continuum).sigma) / (2. * h);
    }

    const double err_algo = norm(a.Lt - J_ref, "fro");
    const double err_cont = norm(c.Lt - J_ref, "fro");
    if (verbose()) {
        std::cout << "\n[Ttangent_EPICP] ||Lt - dsigma/dDeps||_F: algo=" << err_algo
                  << "  cont=" << err_cont
                  << "  (ratio=" << err_cont / std::max(err_algo, 1e-300) << ")\n";
    }
    EXPECT_LT(err_algo, 1.);
    EXPECT_GT(err_cont, 10. * err_algo);
}

// ---------------------------------------------------------------------
// Test 3: Newton convergence rate on a stress-target outer iteration.
// ---------------------------------------------------------------------
TEST(Ttangent_EPICP, convergence_rate_continuum_vs_algorithmic)
{
    const mat L = L_iso(70000., 0.3, "Enu");
    const vec sigma_target = {350., -50., 50., 80., 0., 0.};

    vec Deps0 = solve(L, sigma_target);
    Deps0 *= 1.2;
    Deps0(3) += 3.e-3;   // off-radial perturbation -> direction must be corrected
    Deps0(4) += 1.5e-3;
    Deps0(1) += 1.e-3;

    const int max_iter = 30;
    auto r_cont = stress_target_newton(sigma_target, Deps0, simcoon::tangent_continuum, max_iter);
    auto r_algo = stress_target_newton(sigma_target, Deps0, simcoon::tangent_algorithmic, max_iter);

    if (verbose()) {
        std::cout << "\n[Ttangent_EPICP] residual history (cont vs algo):\n";
        const size_t n = std::max(r_cont.size(), r_algo.size());
        for (size_t k = 0; k < n; ++k) {
            std::cout << "  " << k << "  "
                      << (k < r_cont.size() ? r_cont[k] : -1.) << "   "
                      << (k < r_algo.size() ? r_algo[k] : -1.) << "\n";
        }
    }

    ASSERT_GE(r_algo.size(), 2u);
    EXPECT_LT(r_algo.back(), 1e-8);
    EXPECT_LT(r_algo.size(), r_cont.size())
        << "algorithmic iters " << r_algo.size()
        << " not fewer than continuum " << r_cont.size();
}
