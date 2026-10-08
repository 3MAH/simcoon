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

///@file Tcpp_ccp_robustness.cpp
///@brief Robustness comparison of the closest-point (CPP) helper against a faithful
///       replica of the legacy convex-cutting-plane (CCP) local loop, on identical
///       mechanism definitions, swept over step size / anisotropy / hardening.
///       Run with SIMCOON_TEST_VERBOSE=1 to print the full comparison table.
///       Asserts: CPP converges wherever CCP does (no robustness regression).
///@version 1.0

#include <gtest/gtest.h>
#include <armadillo>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <vector>

#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Continuum_mechanics/Functions/criteria.hpp>
#include <simcoon/Simulation/Maths/num_solve.hpp>
#include <simcoon/Continuum_mechanics/Umat/return_mapping.hpp>

using namespace std;
using namespace arma;
using namespace simcoon;

namespace {

constexpr double E_mod = 70000.;
constexpr double nu_mod = 0.3;
constexpr double sigY = 200.;

struct RunStats { bool ok = false; int iters = 0; double err = 0.; vec sigma = zeros(6); };

// ------------------------------------------------------------------ //
// Case definition: one mechanism family with optional backstress.
// kind: 0 = J2 iso-hardening, 1 = Hill(aniso) iso-hardening,
//       2 = J2 chaboche (1 backstress, linear recall) + linear iso.
// ------------------------------------------------------------------ //
struct CaseDef {
    int kind = 0;
    double H = 1000.;        // isotropic hardening slope
    double C_kin = 0.;       // backstress modulus (kind 2)
    double D_kin = 0.;       // backstress recall (kind 2)
    vec hill = {0.5, 0.6, 0.7, 1.5, 1.4, 1.6};

    double eqs(const vec &xi) const {
        return (kind == 1) ? Hill_stress(xi, hill) : Mises_stress(xi);
    }
    vec deqs(const vec &xi) const {
        return (kind == 1) ? dHill_stress(xi, hill) : eta_stress(xi);
    }
    mat ddeqs(const vec &xi) const {
        return (kind == 1) ? ddHill_stress(xi, hill) : deta_stress(xi);
    }
};

// ------------------------------------------------------------------ //
// Faithful CCP replica (mirrors the UMAT loops: accumulate EP with the
// current-iterate flow; Fischer-Burmeister on the 1x1 system).
// ------------------------------------------------------------------ //
RunStats run_ccp(const CaseDef &cd, const mat &L, const vec &Eps) {
    RunStats st;
    vec EP = zeros(6), a = zeros(6), X = zeros(6);
    double p = 0.;
    vec sigma = L * Eps;

    vec Phi = zeros(1), Y_crit = {sigY}, Ds_j = zeros(1), ds_j = zeros(1);
    mat B = zeros(1, 1);
    double error = 1.;
    int it = 0;
    for (it = 0; (it < simcoon::maxiter_umat) && (error > simcoon::precision_umat); it++) {
        const vec xi = sigma - X;
        const vec n = cd.deqs(xi);
        const vec Lambdap = n;
        const vec Lambdaa = (cd.kind == 2) ? vec(n - cd.D_kin * a) : n;
        const vec kappa = L * Lambdap;

        Phi(0) = cd.eqs(xi) - cd.H * p - sigY;
        double K = -cd.H;
        if (cd.kind == 2) {
            const vec dPhida = -1. * cd.C_kin * (n % Ir05());
            K += sum(dPhida % Lambdaa);
        }
        B(0, 0) = -1. * sum(n % kappa) + K;

        Fischer_Burmeister_m(Phi, Y_crit, B, Ds_j, ds_j, error);

        p += ds_j(0);
        EP += ds_j(0) * Lambdap;
        if (cd.kind == 2) {
            a += ds_j(0) * Lambdaa;
            X = cd.C_kin * (a % Ir05());
        }
        sigma = L * (Eps - EP);
        if (sigma.has_nan()) { st.ok = false; st.iters = it + 1; return st; }
    }
    st.ok = (error <= simcoon::precision_umat);
    st.iters = it;
    st.err = error;
    st.sigma = sigma;
    return st;
}

// ------------------------------------------------------------------ //
// CPP via the helper (same mechanism definitions; inner backward-Euler
// closed form for the chaboche backstress).
// ------------------------------------------------------------------ //
RunStats run_cpp(const CaseDef &cd, const mat &L, const vec &Eps) {
    RunStats st;
    double p = 0.;
    vec a = zeros(6), X = zeros(6);
    const vec sigma_tr = L * Eps;

    ReturnMechanism m;
    m.Phi = [&](const vec &sig) { return cd.eqs(sig - X) - cd.H * p - sigY; };
    m.dPhi_dsigma = [&](const vec &sig) { return cd.deqs(sig - X); };
    if (cd.kind == 2) {
        // Total derivative through the inner-consistent backstress (central FD of the map).
        m.dLambda_dsigma = [&, sigma_tr](const vec &sig) {
            // capture-by-ref of p, a, X; the FD probes re-run the inner solve at fixed Dlambda = p
            const double h = 1.e-5 * (norm(sig, 2) + 1.);
            mat D(6, 6);
            const double Dl = p;  // p = p_n(=0) + Dlambda at the current iterate
            auto solveX = [&](const vec &s) {
                // fixed point on X: a' = (a_n + Dl*n(s-X'))/(1+Dl*D_kin), X' = C*(a'%Ir05)
                vec Xl = X;
                for (int k = 0; k < 10; k++) {
                    const vec n = cd.deqs(s - Xl);
                    const vec al = (Dl * n) / (1. + Dl * cd.D_kin);   // a_n = 0
                    const vec Xn = cd.C_kin * (al % Ir05());
                    if (norm(Xn - Xl, 2) < 1e-12 * (norm(Xn, 2) + 1.)) { Xl = Xn; break; }
                    Xl = Xn;
                }
                return cd.deqs(s - Xl);
            };
            for (int c6 = 0; c6 < 6; c6++) {
                vec sp = sig, sm = sig;
                sp(c6) += h;
                sm(c6) -= h;
                D.col(c6) = (solveX(sp) - solveX(sm)) / (2. * h);
            }
            return D;
        };
    } else {
        m.dLambda_dsigma = [&](const vec &sig) { return cd.ddeqs(sig - X); };
    }

    auto update_state = [&](const vec &sig, double Dl) {
        p = Dl;  // p_n = 0
        if (cd.kind == 2) {
            // backward-Euler fixed point: a = (a_n + Dl*n(sig-X))/(1+Dl*D_kin)
            for (int k = 0; k < 20; k++) {
                const vec n = cd.deqs(sig - X);
                const vec an = (Dl * n) / (1. + Dl * cd.D_kin);
                const vec Xn = cd.C_kin * (an % Ir05());
                const double dX = norm(Xn - X, 2);
                a = an;
                X = Xn;
                if (dX < 1e-12 * (norm(Xn, 2) + 1.)) break;
            }
        }
        return true;
    };
    auto K_scalar = [&](const vec &sig, double) {
        double K = -cd.H;
        if (cd.kind == 2) {
            const vec n = cd.deqs(sig - X);
            const vec dPhida = -1. * cd.C_kin * (n % Ir05());
            const vec Lambdaa = n - cd.D_kin * a;
            K += sum(dPhida % Lambdaa);
        }
        return K;
    };

    auto r = closest_point_return_mapping(sigma_tr, L, m, update_state, K_scalar, sigY);
    st.ok = r.converged;
    st.iters = r.niter;
    st.err = r.error;
    st.sigma = r.sigma;
    return st;
}

bool verbose() { return std::getenv("SIMCOON_TEST_VERBOSE") != nullptr; }

} // namespace

// ---------------------------------------------------------------------
// The sweep: cases x step scales. Asserts CPP converges wherever CCP does.
// ---------------------------------------------------------------------
TEST(Tcpp_ccp_robustness, sweep_no_robustness_regression)
{
    const mat L = L_iso(E_mod, nu_mod, "Enu");
    const vec Eps0 = {5.e-3, -1.e-3, 0.5e-3, 2.e-3, 0., 0.};
    const std::vector<double> scales = {1., 2., 5., 10., 20., 50., 100.};

    struct Named { const char *name; CaseDef cd; };
    std::vector<Named> cases = {
        {"J2  H=1000        ", {0, 1000., 0., 0.}},
        {"J2  H=10 (plateau)", {0, 10., 0., 0.}},
        {"Hill-aniso H=1000 ", {1, 1000., 0., 0.}},
        {"Hill-aniso H=10   ", {1, 10., 0., 0.}},
        {"Chaboche C=5e3 D=50", {2, 100., 5000., 50.}},
        {"Chaboche C=2e4 D=200 (stiff)", {2, 10., 20000., 200.}},
    };

    if (verbose()) {
        std::cout << "\n[CPP vs CCP robustness] iterations to converge ('F' = failed within "
                  << simcoon::maxiter_umat << ")\n";
        std::cout << std::left << std::setw(30) << "case" << " | scale:";
        for (double s : scales) std::cout << std::setw(10) << s;
        std::cout << "\n";
    }

    int ccp_fail = 0, cpp_fail = 0, cpp_only_fail = 0;
    long long ccp_iters_tot = 0, cpp_iters_tot = 0;
    int both_ok = 0;

    for (auto &c : cases) {
        std::string line_ccp, line_cpp;
        for (double s : scales) {
            const vec Eps = s * Eps0;
            RunStats rc = run_ccp(c.cd, L, Eps);
            RunStats rp = run_cpp(c.cd, L, Eps);

            if (!rc.ok) ccp_fail++;
            if (!rp.ok) cpp_fail++;
            if (rc.ok && !rp.ok) cpp_only_fail++;
            if (rc.ok && rp.ok) {
                both_ok++;
                ccp_iters_tot += rc.iters;
                cpp_iters_tot += rp.iters;
                // The CCP/CPP difference is O(||Deps||^2): assert closeness only at
                // moderate steps. At 10-100x yield-scale increments the two (both
                // consistent) discretisations may differ arbitrarily -- that is the
                // expected integrator difference, recorded in the table, not a bug.
                if (s <= 2.) {
                    EXPECT_LT(norm(rp.sigma - rc.sigma, 2), 0.15 * norm(rc.sigma, 2) + 1.)
                        << c.name << " scale " << s;
                }
            }
            char buf[16];
            snprintf(buf, sizeof(buf), "%s%-9d", rc.ok ? "" : "F", rc.iters);
            line_ccp += buf;
            snprintf(buf, sizeof(buf), "%s%-9d", rp.ok ? "" : "F", rp.iters);
            line_cpp += buf;

            // THE assertion: no robustness regression.
            EXPECT_TRUE(rp.ok || !rc.ok)
                << "CPP failed where CCP converged: " << c.name << " scale " << s;
        }
        if (verbose()) {
            std::cout << std::left << std::setw(30) << c.name << " |  CCP: " << line_ccp << "\n";
            std::cout << std::left << std::setw(30) << " " << " |  CPP: " << line_cpp << "\n";
        }
    }

    if (verbose()) {
        std::cout << "\nsummary: CCP failures " << ccp_fail << ", CPP failures " << cpp_fail
                  << " (CPP-only failures " << cpp_only_fail << ")\n";
        if (both_ok > 0) {
            std::cout << "mean iterations where both converged: CCP "
                      << double(ccp_iters_tot) / both_ok << "  CPP "
                      << double(cpp_iters_tot) / both_ok << "\n";
        }
    }
}

// ---------------------------------------------------------------------
// Order of consistency: the CCP/CPP converged-stress difference shrinks as
// O(||Deps_plastic||^2). Strain = (exact-yield strain) + h*(plastic excess);
// halving h must divide the difference by ~4 (>= 3 with safety).
// ---------------------------------------------------------------------
TEST(Tcpp_ccp_robustness, ccp_cpp_second_order_consistency)
{
    const mat L = L_iso(E_mod, nu_mod, "Enu");
    const vec Eps0 = {5.e-3, -1.e-3, 0.5e-3, 2.e-3, 0., 0.};

    struct Named { const char *name; CaseDef cd; };
    std::vector<Named> cases = {
        {"Hill-aniso H=1000", {1, 1000., 0., 0.}},
        {"Chaboche stiff   ", {2, 10., 20000., 200.}},
    };

    for (auto &c : cases) {
        // Strain that lands exactly on the initial yield surface, plus scaled excess.
        const double f_y = sigY / c.cd.eqs(L * Eps0);
        const vec Eps_y = f_y * Eps0;
        auto d_of_h = [&](double h) {
            const vec Eps = Eps_y + h * Eps0;
            RunStats rc = run_ccp(c.cd, L, Eps);
            RunStats rp = run_cpp(c.cd, L, Eps);
            EXPECT_TRUE(rc.ok && rp.ok) << c.name << " h=" << h;
            return norm(rp.sigma - rc.sigma, 2);
        };
        const double d1 = d_of_h(1.), d05 = d_of_h(0.5), d025 = d_of_h(0.25),
                     d0125 = d_of_h(0.125);
        const double r1 = d1 / d05, r2 = d05 / d025, r3 = d025 / d0125;
        if (verbose()) {
            std::cout << "[order] " << c.name << ": d(1)=" << d1 << " d(1/2)=" << d05
                      << " d(1/4)=" << d025 << " d(1/8)=" << d0125 << "  ratios " << r1 << ", "
                      << r2 << ", " << r3 << "\n";
        }
        // Second-order consistency: halving ratios approach 4 from below as h -> 0.
        // Assert clear superlinearity at every level and near-asymptotic rate at the finest.
        EXPECT_GT(r1, 2.) << c.name;
        EXPECT_GT(r2, 2.5) << c.name;
        EXPECT_GT(r3, 3.) << c.name;
    }
}
