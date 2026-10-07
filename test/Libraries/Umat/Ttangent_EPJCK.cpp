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

///@file Ttangent_EPJCK.cpp
///@brief EPJCK (Johnson-Cook) kernel acceptance tests: rate-independent limit identical to
///       EPICP, rate and thermal sensitivities of the yield stress, tangent modes (none = L,
///       algorithmic = exact Jacobian of the discrete map at fixed DTime), and the
///       thermomechanical twin's mechanical block.
///@version 1.0

#include <gtest/gtest.h>
#include <armadillo>
#include <vector>
#include <cmath>

#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Plasticity/plastic_isotropic_ccp.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Plasticity/plastic_johnson_cook_ccp.hpp>
#include <simcoon/Continuum_mechanics/Umat/Thermomechanical/Plasticity/plastic_johnson_cook_ccp.hpp>

using namespace std;
using namespace arma;
using namespace simcoon;

namespace {

// AISI 4340 (Johnson & Cook, 1983): E, nu, alpha, A, B, n, C, edot0, m, T_ref, T_melt
const vec JC_PROPS = {200000., 0.33, 1.e-5, 792., 510., 0.26, 0.014, 1.0, 1.03, 293., 1793.};

struct Out { vec sigma; mat Lt; vec statev; double Wm_d; };

Out run_epjck(const vec &DEtot, double DTime, int tangent_mode, vec props = JC_PROPS,
              double T = 293.15, double DT = 0.) {
    const int nprops = props.n_elem;
    const int nstatev = 9;
    vec Etot = zeros(6);
    vec statev = zeros(nstatev);
    vec sigma = zeros(6);
    mat Lt = zeros(6, 6), L = zeros(6, 6);
    mat DR = eye(3, 3);
    double Time = 0.;
    double Wm = 0., Wm_r = 0., Wm_ir = 0., Wm_d = 0.;
    double tnew_dt = 1.;
    umat_plasticity_johnson_cook_CCP("EPJCK", Etot, DEtot, sigma, Lt, L, DR, nprops, props,
                                     nstatev, statev, T, DT, Time, DTime,
                                     Wm, Wm_r, Wm_ir, Wm_d, 3, 3, true, tnew_dt, tangent_mode);
    return {sigma, Lt, statev, Wm_d};
}

Out run_epicp(const vec &DEtot, int tangent_mode) {
    // EPICP props = [E, nu, alpha, sigmaY, k, m] : the rate-independent, isothermal JC law
    const vec props = {JC_PROPS(0), JC_PROPS(1), JC_PROPS(2), JC_PROPS(3), JC_PROPS(4), JC_PROPS(5)};
    const int nprops = props.n_elem;
    const int nstatev = 8;
    vec Etot = zeros(6);
    vec statev = zeros(nstatev);
    vec sigma = zeros(6);
    mat Lt = zeros(6, 6), L = zeros(6, 6);
    mat DR = eye(3, 3);
    double T = 293.15, DT = 0., Time = 0., DTime = 1.;
    double Wm = 0., Wm_r = 0., Wm_ir = 0., Wm_d = 0.;
    double tnew_dt = 1.;
    umat_plasticity_iso_CCP("EPICP", Etot, DEtot, sigma, Lt, L, DR, nprops, props,
                            nstatev, statev, T, DT, Time, DTime,
                            Wm, Wm_r, Wm_ir, Wm_d, 3, 3, true, tnew_dt, tangent_mode);
    return {sigma, Lt, statev, Wm_d};
}

const vec DEPS_YIELD = {8.e-3, -2.e-3, 1.e-3, 3.e-3, 0., 0.};   // non-radial, well past A

// Closed-form Johnson-Cook yield stress with the kernel's clamps (rate factor >= 1, T* in [0, 1))
double jc_yield(double p, double pdot, double T) {
    const double A = JC_PROPS(3), B = JC_PROPS(4), n = JC_PROPS(5), C = JC_PROPS(6),
                 edot0 = JC_PROPS(7), m = JC_PROPS(8), T_ref = JC_PROPS(9), T_melt = JC_PROPS(10);
    const double rate = 1. + C * std::max(std::log(pdot / edot0), 0.);
    const double Tstar = std::min(std::max((T - T_ref) / (T_melt - T_ref), 0.), 1.);
    return (A + B * std::pow(p, n)) * rate * (1. - std::pow(Tstar, m));
}

} // namespace

// ---------------------------------------------------------------------
// Rate-independent, isothermal limit: C = 0 and T = T_ref make the JC yield stress
// A + B p^n, i.e. EPICP's sigmaY + k p^m. Stress, tangent and plastic strain must agree
// to round-off in both tangent modes.
// ---------------------------------------------------------------------
TEST(Ttangent_EPJCK, reduces_to_EPICP_without_rate_and_thermal_effects)
{
    vec props = JC_PROPS;
    props(6) = 0.;   // C = 0
    for (int mode : {simcoon::tangent_continuum, simcoon::tangent_algorithmic}) {
        Out jc = run_epjck(DEPS_YIELD, 1.e-3, mode, props, 293.);   // T = T_ref: T* = 0
        Out ep = run_epicp(DEPS_YIELD, mode);
        ASSERT_GT(ep.statev(1), simcoon::iota) << "Test set-up did not yield";
        EXPECT_LT(norm(jc.sigma - ep.sigma, 2), 1.e-8 * norm(ep.sigma, 2));
        EXPECT_LT(norm(jc.Lt - ep.Lt, "fro"), 1.e-8 * norm(ep.Lt, "fro"));
        EXPECT_NEAR(jc.statev(1), ep.statev(1), 1.e-12);
        EXPECT_NEAR(jc.Wm_d, ep.Wm_d, 1.e-8 * std::abs(ep.Wm_d));
    }
}

// ---------------------------------------------------------------------
// Rate sensitivity: the same strain increment over a shorter DTime is a higher plastic
// strain rate and must give a higher Mises stress; below the reference rate the law is
// clamped at the quasi-static yield stress (rate factor = 1).
// ---------------------------------------------------------------------
TEST(Ttangent_EPJCK, yield_stress_increases_with_strain_rate_and_is_clamped_below_edot0)
{
    const double DTime_slow = 100.;     // pdot ~ 1e-4 /s << edot0 = 1 /s
    const double DTime_ref = 1.e-2;     // pdot ~ 1 /s   ~  edot0
    const double DTime_fast = 1.e-5;    // pdot ~ 1e3 /s >> edot0
    Out slow = run_epjck(DEPS_YIELD, DTime_slow, simcoon::tangent_algorithmic);
    Out ref = run_epjck(DEPS_YIELD, DTime_ref, simcoon::tangent_algorithmic);
    Out fast = run_epjck(DEPS_YIELD, DTime_fast, simcoon::tangent_algorithmic);
    Out stat = run_epjck(DEPS_YIELD, 0., simcoon::tangent_algorithmic);   // DTime = 0: rate-independent
    ASSERT_GT(slow.statev(1), simcoon::iota);

    // below edot0 the response is the quasi-static one
    EXPECT_LT(norm(slow.sigma - stat.sigma, 2), 1.e-10 * norm(stat.sigma, 2));
    EXPECT_GT(Mises_stress(fast.sigma), Mises_stress(ref.sigma));
    EXPECT_GT(Mises_stress(ref.sigma), Mises_stress(slow.sigma) - 1.e-6);

    // the fast step sits on the rate-dependent surface: Mises = (A + B p^n)(1 + C ln(pdot/edot0))
    const double p = fast.statev(1);
    const double pdot = fast.statev(8);
    EXPECT_NEAR(pdot, p / DTime_fast, 1.e-10 * pdot);
    const double sigmaY = jc_yield(p, pdot, 293.15);
    EXPECT_NEAR(Mises_stress(fast.sigma), sigmaY, 1.e-8 * sigmaY);
}

// ---------------------------------------------------------------------
// Thermal softening: a hotter step yields at a lower stress, consistent with (1 - T*^m).
// ---------------------------------------------------------------------
TEST(Ttangent_EPJCK, yield_stress_decreases_with_temperature)
{
    vec props = JC_PROPS;
    props(2) = 0.;   // no thermal expansion: isolate the softening
    const double DTime = 0.;   // rate-independent: Mises = (A + B p^n)(1 - T*^m)
    Out cold = run_epjck(DEPS_YIELD, DTime, simcoon::tangent_algorithmic, props, 293.);
    Out hot = run_epjck(DEPS_YIELD, DTime, simcoon::tangent_algorithmic, props, 793.);
    EXPECT_GT(Mises_stress(cold.sigma), Mises_stress(hot.sigma));
    const double sigmaY = jc_yield(hot.statev(1), 0., 793.);
    EXPECT_NEAR(Mises_stress(hot.sigma), sigmaY, 1.e-8 * sigmaY);
}

// ---------------------------------------------------------------------
// tangent_none: the kernel returns the elastic operator (explicit integration), while
// the stress update is unchanged.
// ---------------------------------------------------------------------
TEST(Ttangent_EPJCK, tangent_none_returns_elastic_operator)
{
    Out n = run_epjck(DEPS_YIELD, 1.e-3, simcoon::tangent_none);
    Out a = run_epjck(DEPS_YIELD, 1.e-3, simcoon::tangent_algorithmic);
    ASSERT_GT(n.statev(1), simcoon::iota);
    const mat L = L_iso(JC_PROPS(0), JC_PROPS(1), "Enu");
    EXPECT_LT(norm(n.Lt - L, "fro"), 1.e-12 * norm(L, "fro"));
    EXPECT_LT(norm(n.sigma - a.sigma, 2), 1.e-12 * norm(a.sigma, 2));
    EXPECT_GT(norm(a.Lt - L, "fro"), 1.e-3 * norm(L, "fro"));   // the plastic tangent differs
}

// ---------------------------------------------------------------------
// Algorithmic tangent: with the rate term in Bhat, it is the exact Jacobian of the
// discrete map d(sigma)/d(DEtot) at fixed DTime, far closer than the continuum one.
// ---------------------------------------------------------------------
TEST(Ttangent_EPJCK, algorithmic_tangent_is_exact_jacobian_at_fixed_DTime)
{
    const double DTime = 1.e-4;   // pdot ~ 100 /s: the rate term is active
    Out c = run_epjck(DEPS_YIELD, DTime, simcoon::tangent_continuum);
    Out a = run_epjck(DEPS_YIELD, DTime, simcoon::tangent_algorithmic);
    ASSERT_GT(a.statev(1), simcoon::iota);

    const double h = 1.e-7;
    mat J_ref(6, 6);
    for (int j = 0; j < 6; ++j) {
        vec Dp = DEPS_YIELD; Dp(j) += h;
        vec Dm = DEPS_YIELD; Dm(j) -= h;
        J_ref.col(j) = (run_epjck(Dp, DTime, simcoon::tangent_continuum).sigma
                      - run_epjck(Dm, DTime, simcoon::tangent_continuum).sigma) / (2. * h);
    }
    const double err_algo = norm(a.Lt - J_ref, "fro");
    const double err_cont = norm(c.Lt - J_ref, "fro");
    EXPECT_LT(err_algo, 1.e-4 * norm(J_ref, "fro"));
    EXPECT_GT(err_cont, 10. * err_algo);
}

// ---------------------------------------------------------------------
// Thermomechanical twin: same stress and mechanical tangent as the mechanical kernel at
// DT = 0, a positive heat source under plastic flow (dissipation; the increment is isochoric
// so the thermoelastic term -theta L:alpha:D vanishes), and dSdT below the pure thermoelastic
// -L:alpha (softening adds to the stress drop with temperature).
// ---------------------------------------------------------------------
TEST(Ttangent_EPJCK, thermomechanical_twin_matches_mechanical_block)
{
    const double rho = 7.85e-9, c_p = 4.75e8;   // t/mm^3, mJ/(t K)
    vec props_T(13);
    props_T(0) = rho; props_T(1) = c_p;
    props_T.subvec(2, 12) = JC_PROPS;
    const double DTime = 1.e-3;
    const double T = 500.;   // T* > 0 so that dPhi/dT is active
    const vec Deps = {8.e-3, -4.e-3, -4.e-3, 3.e-3, 0., 0.};   // isochoric

    vec Etot = zeros(6), statev = zeros(9), sigma = zeros(6);
    double r = 0.;
    mat dSdE = zeros(6, 6), dSdT = zeros(6, 1), drdE = zeros(6, 1), drdT = zeros(1, 1);
    mat DR = eye(3, 3);
    double Wm = 0., Wm_r = 0., Wm_ir = 0., Wm_d = 0., Wt = 0., Wt_r = 0., Wt_ir = 0.;
    double tnew_dt = 1.;
    umat_plasticity_johnson_cook_CCP_T(Etot, Deps, sigma, r, dSdE, dSdT, drdE, drdT, DR,
                                       13, props_T, 9, statev, T, 0., 0., DTime,
                                       Wm, Wm_r, Wm_ir, Wm_d, Wt, Wt_r, Wt_ir, 3, 3, true, tnew_dt,
                                       simcoon::tangent_algorithmic);
    Out m = run_epjck(Deps, DTime, simcoon::tangent_algorithmic, JC_PROPS, T);
    ASSERT_GT(statev(1), simcoon::iota);
    EXPECT_LT(norm(sigma - m.sigma, 2), 1.e-10 * norm(m.sigma, 2));
    EXPECT_LT(norm(dSdE - m.Lt, "fro"), 1.e-10 * norm(m.Lt, "fro"));
    EXPECT_NEAR(Wm_d, m.Wm_d, 1.e-10 * std::abs(m.Wm_d));
    EXPECT_GT(r, 0.) << "plastic dissipation must heat the material";

    const mat L = L_iso(JC_PROPS(0), JC_PROPS(1), "Enu");
    const vec alpha = JC_PROPS(2) * Ith();
    const vec dSdT_el = -L * alpha;
    // flow direction ~ sigma', so the softening term kappa P_theta reduces sigma_11 further
    EXPECT_LT(dSdT(0, 0), dSdT_el(0));
}
