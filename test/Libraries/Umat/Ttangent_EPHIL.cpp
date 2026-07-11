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

///@file Ttangent_EPHIL.cpp
///@brief Tier-B acceptance test: the EPHIL (Hill anisotropic plasticity) algorithmic
///       tangent (tangent_mode==1) is symmetric, substantially closer to the discrete-map
///       Jacobian than the continuum tangent, and converges in fewer Newton iterations.
///       It is machine-exact for isotropic Hill (== J2); for ANISOTROPIC Hill the CCP
///       integrator differs from a closest-point projection, so exactness is deferred to
///       the CPP return-map rework (see the NOTE above the Jacobian test).
///       Drives the EPHIL name through the modular adapter path.
///@version 1.0

#include <gtest/gtest.h>
#include <armadillo>
#include <vector>
#include <cstdlib>

#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Continuum_mechanics/Umat/Modular/legacy_adapters.hpp>

using namespace std;
using namespace arma;
using namespace simcoon;

namespace {

// EPHIL: props = [E, nu, alpha, sigmaY, k, m, F, G, H, L, M, N] ; nstatev = 8, statev(1)=p.
struct Out { vec sigma; mat Lt; bool plastic; };

Out run_ephil(const vec &DEtot, int tangent_mode) {
    const vec props = {70000., 0.3, 1.e-5, 200., 1000., 0.3,
                       0.5, 0.6, 0.7, 1.5, 1.4, 1.6};
    const int nprops = props.n_elem;
    const int nstatev = 8;

    vec Etot = zeros(6);
    vec statev = zeros(nstatev);
    vec sigma = zeros(6);
    mat Lt = zeros(6, 6), L = zeros(6, 6);
    mat DR = eye(3, 3);
    double T = 293.15, DT = 0., Time = 0., DTime = 1.;
    double Wm = 0., Wm_r = 0., Wm_ir = 0., Wm_d = 0., tnew_dt = 1.;

    umat_legacy_modular("EPHIL", Etot, DEtot, sigma, Lt, L, DR, nprops, props,
                                  nstatev, statev, T, DT, Time, DTime,
                                  Wm, Wm_r, Wm_ir, Wm_d, 3, 3, true, tnew_dt, tangent_mode);
    return {sigma, Lt, statev(1) > simcoon::iota};
}

vector<double> stress_target_newton(const vec &sigma_target, vec DEtot, int tangent_mode, int max_iter) {
    vector<double> res;
    for (int it = 0; it < max_iter; ++it) {
        Out o = run_ephil(DEtot, tangent_mode);
        vec r = sigma_target - o.sigma;
        res.push_back(norm(r, 2));
        if (res.back() < 1e-10) break;
        DEtot += solve(o.Lt, r);
    }
    return res;
}

bool verbose() { return std::getenv("SIMCOON_TEST_VERBOSE") != nullptr; }

} // namespace

TEST(Ttangent_EPHIL, parity_elastic_step)
{
    const vec DEtot = {1.e-4, -0.3e-4, -0.3e-4, 0., 0., 0.};
    Out c = run_ephil(DEtot, 0);
    Out a = run_ephil(DEtot, 1);
    ASSERT_FALSE(c.plastic) << "Parity step unexpectedly yielded.";
    EXPECT_LT(norm(c.Lt - a.Lt, "fro"), 1.e-12);
}

// NOTE on anisotropic Hill: EPHIL uses the Convex Cutting Plane (CCP) integrator. For an
// ANISOTROPIC Hill criterion the deviatoric flow normal rotates during the return, so CCP differs
// from a closest-point projection: the (correct, symmetric) closest-point algorithmic tangent does
// NOT exactly match the CCP map's non-symmetric finite-difference Jacobian. It IS exact for
// J2/isotropic and remains a large improvement here. Machine-exact anisotropic tangents require
// switching EPHIL to a closest-point return map (future release; see project_future_cpp_return_map).
TEST(Ttangent_EPHIL, algorithmic_tangent_better_than_continuum)
{
    const vec Deps = {5.e-3, -1.e-3, 0.5e-3, 2.e-3, 0., 0.};
    Out c = run_ephil(Deps, 0);
    Out a = run_ephil(Deps, 1);
    ASSERT_TRUE(a.plastic) << "Test set-up did not yield - increase Deps.";

    const double h = 1.e-7;
    mat J_ref(6, 6);
    for (int j = 0; j < 6; ++j) {
        vec Dp = Deps; Dp(j) += h;
        vec Dm = Deps; Dm(j) -= h;
        J_ref.col(j) = (run_ephil(Dp, 0).sigma - run_ephil(Dm, 0).sigma) / (2. * h);
    }
    const double err_algo = norm(a.Lt - J_ref, "fro");
    const double err_cont = norm(c.Lt - J_ref, "fro");
    if (verbose()) {
        std::cout << "\n[Ttangent_EPHIL] ||Lt - dsigma/dDeps||_F: algo=" << err_algo
                  << "  cont=" << err_cont
                  << "  (ratio=" << err_cont / std::max(err_algo, 1e-300) << ")\n";
    }
    // The closest-point algorithmic tangent is symmetric (consistent tangent of associated flow)...
    EXPECT_LT(norm(a.Lt - a.Lt.t(), "fro"), 1.e-6 * norm(a.Lt, "fro"));
    // ...and substantially closer to the discrete-map Jacobian than the continuum tangent.
    // (Not machine-exact for anisotropic CCP -- see the note above; exact for J2/isotropic.)
    EXPECT_GT(err_cont, 3. * err_algo);
}

TEST(Ttangent_EPHIL, convergence_rate_continuum_vs_algorithmic)
{
    const mat L = L_iso(70000., 0.3, "Enu");
    const vec sigma_target = {350., -50., 50., 80., 0., 0.};

    vec Deps0 = solve(L, sigma_target);
    Deps0 *= 1.2;
    Deps0(3) += 3.e-3;
    Deps0(4) += 1.5e-3;
    Deps0(1) += 1.e-3;

    const int max_iter = 30;
    auto r_cont = stress_target_newton(sigma_target, Deps0, 0, max_iter);
    auto r_algo = stress_target_newton(sigma_target, Deps0, 1, max_iter);

    if (verbose()) {
        std::cout << "\n[Ttangent_EPHIL] residual history (cont vs algo):\n";
        const size_t n = std::max(r_cont.size(), r_algo.size());
        for (size_t k = 0; k < n; ++k)
            std::cout << "  " << k << "  " << (k < r_cont.size() ? r_cont[k] : -1.)
                      << "   " << (k < r_algo.size() ? r_algo[k] : -1.) << "\n";
    }
    ASSERT_GE(r_algo.size(), 2u);
    EXPECT_LT(r_algo.back(), 1e-8);
    EXPECT_LT(r_algo.size(), r_cont.size());
}
