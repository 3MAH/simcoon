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

///@file Tunified_TR.cpp
///@brief Smoke tests for the unified SMA transformation+reorientation UMAT
///       (umat_sma_unified_TR, variant SMRDI). Drives a single material point
///       directly through the constitutive routine with adaptive sub-stepping,
///       checking the qualitative physics rather than tabulated values:
///         - superelastic recovery (closed loop: ξ and σ return to ~0 at ε=0),
///         - dissipation Wm_d is positive and 2nd-law-monotone over the cycle,
///         - the energy partition closes (Wm = Wm_r + Wm_d, Wm_ir = 0),
///         - the finite-rotation (rotate_strain) path runs without blowing up.
///       Parameters follow the Chemisky-thesis / Grabe-Bruhns SMRDI set.
///@version 1.0

#include <gtest/gtest.h>
#include <armadillo>
#include <cmath>

#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/SMA/unified_TR.hpp>

using namespace std;
using namespace arma;
using namespace simcoon;

namespace {

// SMRDI (isotropic elasticity, Drucker criteria) — 35 properties.
// Layout per unified_TR.cpp: offset 5 (iso), aniso_offset 28, reorientation 28..34.
// Values: Chemisky-thesis-faithful set shifted to a fully-pseudoelastic margin.
vec smrdi_props() {
    vec p = zeros(35);
    p(0)  = 0.;                                  // flagT (smooth transition functions)
    p(1)  = 70000.;  p(2) = 70000.;              // E_A, E_M
    p(3)  = 0.30;    p(4) = 0.30;                // nu_A, nu_M
    p(5)  = 8.e-6;   p(6) = 8.e-6;              // alphaA, alphaM
    p(7)  = 0.0;     p(8) = 0.050;              // Hmin, Hmax
    p(9)  = 0.012;   p(10) = 0.0;               // k1, sigmacrit
    p(11) = 6.0;     p(12) = 5.0;               // C_A, C_M
    p(13) = 220.65;  p(14) = 200.65;            // Ms0, Mf0
    p(15) = 240.65;  p(16) = 260.65;            // As0, Af0   (Af0 = 40 K below T below)
    p(17) = 0.3; p(18) = 0.3; p(19) = 0.3; p(20) = 0.3;   // n1..n4
    p(21) = 400.0;                               // sigmacaliber
    p(22) = 0.79;    p(23) = 2.0;               // prager_b, prager_n (Sittner asymmetry)
    p(24) = 1.e-6; p(25) = 1.e-3; p(26) = 1.0; p(27) = 1.e8;   // c/p0/n/alpha lambda (transf.)
    p(28) = 100.0;   p(29) = 1500.0; p(30) = 0.10;            // YReo, HReo, ETRmax
    p(31) = 1.e-6; p(32) = 1.e-3; p(33) = 1.0; p(34) = 1.e8;   // c/p0/n/alpha lambdaReo
    return p;
}

// Single material-point driver with solver-style adaptive sub-stepping: if the
// UMAT requests a cut (tnew_dt < 1) the state is rolled back and the macro-step
// is halved. Mirrors what the global solver does, so the smoke test never
// commits a non-converged increment.
struct Driver {
    vec    props   = smrdi_props();
    int    nprops  = 35;
    int    nstatev = 30;
    vec    statev  = zeros(30);
    vec    stress  = zeros(6);
    vec    Etot    = zeros(6);
    mat    Lt = zeros(6, 6), L = zeros(6, 6);
    double T = 300.65, DT = 0.0;            // isothermal, ~40 K above Af0 -> pseudoelastic
    int    ndi = 3, nshr = 3;
    bool   start = true;
    double Wm = 0., Wm_r = 0., Wm_ir = 0., Wm_d = 0.;
    double Time = 0., DTime = 1.;

    // Advance by a prescribed strain increment vector (Voigt), sub-stepping as needed.
    // Returns false if the sub-step collapsed (model failed to converge).
    bool step(const vec& DEtot_macro) {
        double done = 0., frac = 1.;
        int guard = 0;
        while (done < 1.0 - 1e-9) {
            if (++guard > 2000) return false;
            const double s = std::min(frac, 1.0 - done);
            // snapshot
            const vec    ss = stress, sv = statev, et = Etot;
            const double w = Wm, wr = Wm_r, wir = Wm_ir, wd = Wm_d;
            const bool   st = start;
            vec    DE   = s * DEtot_macro;
            double tnew = 1.0;
            const mat DR = eye(3, 3);
            umat_sma_unified_TR("SMRDI", Etot, DE, stress, Lt, L, DR, nprops, props,
                                nstatev, statev, T, DT, Time, DTime, Wm, Wm_r, Wm_ir, Wm_d,
                                ndi, nshr, start, tnew, 0);
            if (!stress.is_finite() || !statev.is_finite()) return false;
            if (tnew < 1.0) {                     // rejected: roll back, halve
                stress = ss; statev = sv; Etot = et;
                Wm = w; Wm_r = wr; Wm_ir = wir; Wm_d = wd; start = st;
                frac *= 0.5;
                if (frac < 1e-5) return false;
            } else {                              // accepted
                Etot += DE; done += s; start = false;
                frac = std::min(1.0, frac * 1.5);
            }
        }
        return true;
    }

    double xi() const { return statev(1); }
};

} // namespace

// ---------------------------------------------------------------------------
// Superelastic uniaxial cycle: load E11 0 -> e_max -> 0 at constant T well above
// Af. The loop must close (transformation forward then fully reverse), the
// dissipation must be positive and monotone, and the energy partition must
// balance.
// ---------------------------------------------------------------------------
TEST(Tunified_TR, smrdi_superelastic_cycle)
{
    Driver d;
    const double e_max = 0.06;
    const int    N     = 120;
    const double de    = e_max / N;

    vec dir = zeros(6); dir(0) = 1.0;     // uniaxial axial strain control

    double xi_peak = 0.0, Wm_d_prev = 0.0;

    // Loading
    for (int k = 0; k < N; ++k) {
        ASSERT_TRUE(d.step(de * dir)) << "loading failed to converge at step " << k;
        xi_peak = std::max(xi_peak, d.xi());
        ASSERT_GE(d.Wm_d, Wm_d_prev - 1e-3) << "dissipation decreased (2nd law) on load step " << k;
        Wm_d_prev = d.Wm_d;
    }
    // Unloading
    for (int k = 0; k < N; ++k) {
        ASSERT_TRUE(d.step(-de * dir)) << "unloading failed to converge at step " << k;
        ASSERT_GE(d.Wm_d, Wm_d_prev - 1e-3) << "dissipation decreased (2nd law) on unload step " << k;
        Wm_d_prev = d.Wm_d;
    }

    // The model actually transformed (not a trivial elastic pass).
    EXPECT_GT(xi_peak, 0.3) << "forward transformation barely occurred (xi_peak)";
    // Superelastic recovery: martensite reverts and stress returns to ~0 at ε=0.
    EXPECT_LT(d.xi(), 0.05)            << "martensite did not revert (xi at end)";
    EXPECT_LT(Mises_stress(d.stress), 25.) << "stress did not recover at zero strain";
    // A hysteresis loop was traversed -> strictly positive dissipation.
    EXPECT_GT(d.Wm_d, 0.0) << "no dissipation: hysteresis missing";
    // Energy partition closes; SMA has no irrecoverable channel.
    EXPECT_NEAR(d.Wm, d.Wm_r + d.Wm_d, 1e-6 * std::max(1.0, std::abs(d.Wm)));
    EXPECT_NEAR(d.Wm_ir, 0.0, 1e-12);
}

// ---------------------------------------------------------------------------
// Finite-rotation path: form martensite, then apply a step carrying a genuine
// rigid rotation increment DR (exercising rotate_strain on ET/EReo/areo). The
// constitutive update must stay finite and stable.
// ---------------------------------------------------------------------------
TEST(Tunified_TR, smrdi_rotation_runs_clean)
{
    Driver d;
    vec dir = zeros(6); dir(0) = 1.0;

    // Partially transform.
    for (int k = 0; k < 40; ++k) {
        ASSERT_TRUE(d.step(0.001 * dir)) << "pre-rotation loading failed at step " << k;
    }
    ASSERT_GT(d.xi(), 0.05) << "no martensite formed before rotation step";

    // One increment carrying a 0.1 rad rotation about z, with a small strain bump.
    const double th = 0.1;
    mat DR = eye(3, 3);
    DR(0, 0) = std::cos(th); DR(0, 1) = -std::sin(th);
    DR(1, 0) = std::sin(th); DR(1, 1) =  std::cos(th);

    vec DE = zeros(6); DE(0) = 0.0005;
    double tnew = 1.0;
    bool start = false;
    ASSERT_NO_THROW(
        umat_sma_unified_TR("SMRDI", d.Etot, DE, d.stress, d.Lt, d.L, DR, d.nprops, d.props,
                            d.nstatev, d.statev, d.T, d.DT, d.Time, d.DTime,
                            d.Wm, d.Wm_r, d.Wm_ir, d.Wm_d, d.ndi, d.nshr, start, tnew, 0)
    );
    EXPECT_TRUE(d.stress.is_finite()) << "stress non-finite after rotation step";
    EXPECT_TRUE(d.statev.is_finite()) << "state non-finite after rotation step";
    EXPECT_TRUE(d.Lt.is_finite())     << "tangent non-finite after rotation step";
}
