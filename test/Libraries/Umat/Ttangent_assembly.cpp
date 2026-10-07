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

///@file Ttangent_assembly.cpp
///@brief Tests for assemble_continuum_tangent / assemble_algorithmic_tangent.
///       The J2 closed-form return map lives only inside this test file —
///       the library itself stays general.
///@version 1.0

#include <gtest/gtest.h>
#include <armadillo>
#include <vector>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <cstdlib>

#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Continuum_mechanics/Umat/tangent_assembly.hpp>

using namespace std;
using namespace arma;
using namespace simcoon;

namespace {

// Diagnostic printing is off by default — ctest stdout stays clean. Set
// SIMCOON_TEST_VERBOSE=1 to see the FD-vs-Lt ratio and the residual history.
static bool verbose() {
    const char* v = std::getenv("SIMCOON_TEST_VERBOSE");
    return v && *v && *v != '0';
}

struct J2Result {
    vec sigma;
    double Dp;
    vec Lambda_eps;
    vec dPhidsigma;
    mat dLambda_dsigma;
    double K_lj;
    double Bhat;
    bool plastic;
};

// Closed-form J2 + linear iso hardening return map from σ_n = 0. The
// non-obvious bit is K_lj = -H' to match the UMAT convention
// Bhat = ∂Φ/∂σ · κ − K (where K = ∂Φ/∂p · ∂p/∂s).
static J2Result j2_radial_return(const mat& L, const vec& Deps,
                                 double sigmaY0, double Hp, double E, double nu)
{
    const vec sigma_trial = L * Deps;
    const double sigma_eq_trial = Mises_stress(sigma_trial);

    J2Result R;
    R.plastic = sigma_eq_trial > sigmaY0;

    if (!R.plastic) {
        R.sigma = sigma_trial;
        R.Dp = 0.;
        R.Lambda_eps = zeros(6);
        R.dPhidsigma = zeros(6);
        R.dLambda_dsigma = zeros(6, 6);
        R.K_lj = -Hp;
        R.Bhat = 1.;     // arbitrary > 0; inactive mechanism, masked out
        return R;
    }

    const double mu = E / (2. * (1. + nu));
    R.Dp = (sigma_eq_trial - sigmaY0) / (3. * mu + Hp);

    const vec s_trial_dev = dev(sigma_trial);
    R.sigma = sigma_trial - 2. * mu * R.Dp * (3./2.) * s_trial_dev / sigma_eq_trial;

    R.Lambda_eps = eta_stress(R.sigma);
    R.dPhidsigma = eta_stress(R.sigma);

    // ∂Λ_ε/∂σ via centred FD — test-only, the library stays general.
    R.dLambda_dsigma.set_size(6, 6);
    const double h = 1.0;
    for (int j = 0; j < 6; ++j) {
        vec sp = R.sigma; sp(j) += h;
        vec sm = R.sigma; sm(j) -= h;
        R.dLambda_dsigma.col(j) = (eta_stress(sp) - eta_stress(sm)) / (2. * h);
    }

    R.K_lj = -Hp;
    const vec kappa = L * R.Lambda_eps;
    R.Bhat = dot(R.dPhidsigma, kappa) - R.K_lj;
    return R;
}

// Drive one full local update from σ_n=0 with strain increment Δε,
// and assemble the requested tangent. Returns (σ_new, Lt).
struct UmatStep {
    vec sigma;
    mat Lt;
    bool plastic;
};

static UmatStep run_step(const mat& L, const vec& Deps,
                         double sigmaY0, double Hp, double E, double nu,
                         bool algorithmic)
{
    J2Result r = j2_radial_return(L, Deps, sigmaY0, Hp, E, nu);
    UmatStep s;
    s.sigma = r.sigma;
    s.plastic = r.plastic;

    if (!r.plastic) {
        s.Lt = L;
        return s;
    }

    const vec kappa = L * r.Lambda_eps;
    if (algorithmic) {
        ContinuumTangent ct = assemble_algorithmic_tangent(
            r.Bhat, kappa, r.dPhidsigma, r.Dp, L, r.dLambda_dsigma);
        s.Lt = ct.Lt;
    } else {
        ContinuumTangent ct = assemble_continuum_tangent(
            r.Bhat, kappa, r.dPhidsigma, r.Dp, L);
        s.Lt = ct.Lt;
    }
    return s;
}

} // anonymous namespace

// -------------------------------------------------------------------------
// Test 1: PARITY. When ∂Λ/∂σ = 0, algorithmic == continuum bit-for-bit.
// -------------------------------------------------------------------------
TEST(Ttangent_assembly, parity_zero_dLambda_singleMech)
{
    const mat L = L_iso(70000., 0.3, "Enu");
    const vec kappa = { 1., -0.5, -0.5, 0.3, 0.0, 0.1 };
    const vec dPhi  = { 1., -0.5, -0.5, 0.6, 0.0, 0.2 };
    const double Bhat = 1.2e3;
    const double Ds = 1e-3;

    ContinuumTangent ct = assemble_continuum_tangent(Bhat, kappa, dPhi, Ds, L);
    ContinuumTangent at = assemble_algorithmic_tangent(
        Bhat, kappa, dPhi, Ds, L, zeros(6, 6));

    EXPECT_LT(norm(ct.Lt - at.Lt, "fro"), 1e-12);
    EXPECT_LT(norm(ct.P_epsilon[0] - at.P_epsilon[0], 2), 1e-12);
    EXPECT_LT(std::abs(ct.invBhat(0,0) - at.invBhat(0,0)), 1e-12);
}

TEST(Ttangent_assembly, parity_zero_dLambda_multiMech)
{
    // Two-mechanism dummy problem; verifies the parity holds for Nmech > 1
    // and for an off-diagonal Bhat.
    const mat L = L_iso(70000., 0.3, "Enu");
    std::vector<vec> kappa = {
        L * vec({ 1., -0.5, -0.5, 0.2, 0.0, 0.0 }),
        L * vec({ 0., 0.0, 0.0, 0.5, 0.5, 0.5 })
    };
    std::vector<vec> dPhi = {
        { 1., -0.5, -0.5, 0.4, 0.0, 0.0 },
        { 0.0, 0.0, 0.0, 1.0, 1.0, 1.0 }
    };
    mat Bhat(2, 2);
    Bhat(0, 0) = dot(dPhi[0], kappa[0]) + 100.;
    Bhat(0, 1) = dot(dPhi[0], kappa[1]) - 10.;
    Bhat(1, 0) = dot(dPhi[1], kappa[0]) - 10.;
    Bhat(1, 1) = dot(dPhi[1], kappa[1]) + 50.;
    const vec Ds = { 1e-3, 5e-4 };

    std::vector<mat> dLambda_zero = { zeros(6, 6), zeros(6, 6) };

    auto ct = assemble_continuum_tangent(Bhat, kappa, dPhi, Ds, L);
    auto at = assemble_algorithmic_tangent(Bhat, kappa, dPhi, Ds, L, dLambda_zero);

    EXPECT_LT(norm(ct.Lt - at.Lt, "fro"), 1e-10);
    for (size_t l = 0; l < 2; ++l) {
        EXPECT_LT(norm(ct.P_epsilon[l] - at.P_epsilon[l], 2), 1e-10);
    }
    EXPECT_LT(norm(ct.invBhat - at.invBhat, "fro"), 1e-12);
}

TEST(Ttangent_assembly, parity_inactive_mechanism_masked_out)
{
    // Confirm the active-set mask still works in the algorithmic helper.
    const mat L = L_iso(70000., 0.3, "Enu");
    std::vector<vec> kappa = {
        L * vec({ 1., -0.5, -0.5, 0.0, 0.0, 0.0 }),
        L * vec({ 0., 0.0, 0.0, 1.0, 0.0, 0.0 })
    };
    std::vector<vec> dPhi = {
        { 1., -0.5, -0.5, 0.0, 0.0, 0.0 },
        { 0.0, 0.0, 0.0, 1.0, 0.0, 0.0 }
    };
    mat Bhat(2, 2);
    Bhat(0, 0) = 1e3;
    Bhat(1, 1) = 1e3;
    Bhat(0, 1) = -1.;
    Bhat(1, 0) = -1.;
    const vec Ds = { 1e-3, 0.0 }; // mechanism 2 inactive

    std::vector<mat> dL = { mat(6,6,fill::randu)*1e-4, zeros(6,6) };

    auto at = assemble_algorithmic_tangent(Bhat, kappa, dPhi, Ds, L, dL);

    // Inactive row/col of invBhat must be zero.
    EXPECT_LT(std::abs(at.invBhat(1, 0)), 1e-14);
    EXPECT_LT(std::abs(at.invBhat(0, 1)), 1e-14);
    EXPECT_LT(std::abs(at.invBhat(1, 1)), 1e-14);
    // Inactive P_eps must be zero.
    EXPECT_LT(norm(at.P_epsilon[1], 2), 1e-14);
}

// -------------------------------------------------------------------------
// Test 2: CONVERGENCE RATE. Drive a stress-target outer Newton iteration
// using each tangent as the iteration Jacobian, and verify the asymptotic
// order of convergence.
//
// Outer iteration:
//   F(Δε) := σ_target − σ_J2(Δε)
//   J(Δε) := dσ/dΔε  ≈  Lt(Δε)
//   Δε_{k+1} = Δε_k + J^{-1} F
//
// With the algorithmic tangent, Lt is the true derivative of the
// converged discrete update, so the iteration is the exact Newton
// scheme on F and converges quadratically (Simo–Hughes Thm 3.1).
// With the continuum tangent, Lt is off by Δp · ∂Λ/∂σ terms and the
// iteration converges linearly (Q-linear with a finite contraction
// factor that is not zero).
// -------------------------------------------------------------------------

static std::vector<double> stress_target_newton(const mat& L,
                                                const vec& sigma_target,
                                                const vec& Deps0,
                                                double sigmaY0, double Hp,
                                                double E, double nu,
                                                bool algorithmic,
                                                int max_iter)
{
    vec Deps = Deps0;
    std::vector<double> residuals;

    for (int k = 0; k < max_iter; ++k) {
        UmatStep s = run_step(L, Deps, sigmaY0, Hp, E, nu, algorithmic);
        vec F = sigma_target - s.sigma;
        residuals.push_back(norm(F, 2));
        // Converged once the residual reaches the test's convergence threshold (1e-8, as
        // asserted below). optimized BLAS reaches ~1e-13 but reference-netlib (conda Linux) floors above 1e-12.
        if (residuals.back() < 1e-8) break;
        // Newton step
        vec dDeps = solve(s.Lt, F);
        Deps += dDeps;
    }
    return residuals;
}

// ---------------------------------------------------------------------
// Test 2a: VERIFY THE EXACT DERIVATIVE PROPERTY
// The algorithmic tangent is the exact derivative of the discrete
// return-mapping update with respect to Δε (Simo–Hughes Thm 3.1).
// Compare both tangents to a centred-finite-difference Jacobian.
// ---------------------------------------------------------------------
TEST(Ttangent_assembly, algorithmic_tangent_is_exact_jacobian)
{
    const double E = 70000., nu = 0.3;
    const double sigmaY0 = 200., Hp = 1000.;
    const mat L = L_iso(E, nu, "Enu");
    // Non-aligned plastic strain increment: combination of normal and
    // shear components so ∂Λ/∂σ is non-trivial.
    const vec Deps = { 5e-3, -1e-3, 0.5e-3, 2e-3, 0.0, 0.0 };

    UmatStep s_cont = run_step(L, Deps, sigmaY0, Hp, E, nu, /*algorithmic=*/false);
    UmatStep s_algo = run_step(L, Deps, sigmaY0, Hp, E, nu, /*algorithmic=*/true);
    ASSERT_TRUE(s_algo.plastic) << "Test set-up did not yield — increase Δε.";

    // Centred-FD reference: J_ref(i,j) = (σ_i(Δε + h e_j) − σ_i(Δε − h e_j)) / (2h)
    const double h = 1e-7;
    mat J_ref(6, 6);
    for (int j = 0; j < 6; ++j) {
        vec Dp = Deps; Dp(j) += h;
        vec Dm = Deps; Dm(j) -= h;
        UmatStep sp = run_step(L, Dp, sigmaY0, Hp, E, nu, /*algorithmic=*/false);
        UmatStep sm = run_step(L, Dm, sigmaY0, Hp, E, nu, /*algorithmic=*/false);
        J_ref.col(j) = (sp.sigma - sm.sigma) / (2. * h);
    }

    const double err_algo = norm(s_algo.Lt - J_ref, "fro");
    const double err_cont = norm(s_cont.Lt - J_ref, "fro");

    if (verbose()) {
        std::cout << "\n[Ttangent_assembly] || Lt - dsigma/dDeps ||_F vs centred-FD:\n"
                  << "  algorithmic = " << err_algo
                  << " ;  continuum = " << err_cont
                  << "  (ratio cont/algo = " << err_cont / std::max(err_algo, 1e-300)
                  << ")\n";
    }

    EXPECT_LT(err_algo, 1.) << "Algorithmic tangent does not match the "
                               "discrete-update Jacobian.";
    EXPECT_GT(err_cont, 10. * err_algo)
        << "Continuum tangent should be visibly farther from the true "
           "Jacobian than the algorithmic one (non-radial Deps required).";
}

// ---------------------------------------------------------------------
// Test 2b: NEWTON CONVERGENCE RATE on a stress-target outer iteration.
// The initial guess is perturbed off the σ_target radial line so the
// continuum-vs-algorithmic difference shows up. (A purely radial scaling
// of the initial guess makes σ linear in the scale factor and any
// reasonable tangent converges in one Newton step.)
// ---------------------------------------------------------------------
TEST(Ttangent_assembly, convergence_rate_continuum_vs_algorithmic)
{
    const double E = 70000., nu = 0.3;
    const double sigmaY0 = 200., Hp = 1000.;
    const mat L = L_iso(E, nu, "Enu");
    const vec sigma_target = { 350., -50., 50., 80., 0., 0. };

    // Initial guess: the elastic predictor of σ_target shifted in a
    // non-radial direction (extra shear and off-axis normal terms) so
    // that the Newton outer loop must adjust both magnitude AND direction.
    vec Deps0 = solve(L, sigma_target);
    Deps0 *= 1.2;
    Deps0(3) += 3e-3;   // off-radial shear perturbation
    Deps0(4) += 1.5e-3;
    Deps0(1) += 1e-3;

    const int max_iter = 30;
    auto r_cont = stress_target_newton(L, sigma_target, Deps0, sigmaY0, Hp,
                                       E, nu, /*algorithmic=*/false, max_iter);
    auto r_algo = stress_target_newton(L, sigma_target, Deps0, sigmaY0, Hp,
                                       E, nu, /*algorithmic=*/true,  max_iter);

    if (verbose()) {
        std::cout << "\n[Ttangent_assembly] Residual history (||sigma_target - sigma||):\n";
        std::cout << "  iter | continuum                | algorithmic\n";
        std::cout << "  -----+--------------------------+--------------------------\n";
        const size_t n = std::max(r_cont.size(), r_algo.size());
        for (size_t k = 0; k < n; ++k) {
            std::cout << "  " << std::setw(4) << k << " | ";
            if (k < r_cont.size()) std::cout << std::setw(24) << r_cont[k];
            else                   std::cout << std::setw(24) << "-";
            std::cout << " | ";
            if (k < r_algo.size()) std::cout << std::setw(24) << r_algo[k];
            else                   std::cout << std::setw(24) << "-";
            std::cout << "\n";
        }
    }

    ASSERT_GE(r_algo.size(), 2u);
    EXPECT_LE(r_algo.back(), r_cont.back());

    // The algorithmic iteration reaches machine-precision residual.
    EXPECT_LT(r_algo.back(), 1e-8)
        << "Algorithmic Newton did not drive the residual to machine "
           "precision; check assemble_algorithmic_tangent.";

    // Quadratic asymptotic ratio for the algorithmic descent.
    if (r_algo.size() >= 3) {
        // pick the last two non-machine-precision points
        size_t k = 0;
        for (size_t i = 0; i + 1 < r_algo.size(); ++i) {
            if (r_algo[i] > 1e-6 && r_algo[i + 1] < r_algo[i]) k = i;
        }
        const double ratio_q = r_algo[k + 1] / (r_algo[k] * r_algo[k] + 1e-300);
        if (verbose()) {
            std::cout << "  algorithmic Q-quadratic ratio r[k+1]/r[k]^2 = "
                      << ratio_q << " (k=" << k << ")\n";
        }
        // Q-quadratic ratio O(1) is the smoking gun.
        EXPECT_LT(ratio_q, 1.)
            << "Algorithmic Newton failed Q-quadratic asymptotic test.";
    }

    // The algorithmic iteration converges to machine precision in
    // strictly fewer iterations than the continuum one — that's the
    // operational measure of the convergence-rate improvement.
    EXPECT_LT(r_algo.size(), r_cont.size())
        << "Algorithmic iteration count " << r_algo.size()
        << " is not strictly fewer than continuum " << r_cont.size()
        << "; the algorithmic tangent did not deliver its expected speed-up.";
}
