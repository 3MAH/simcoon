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

///@file Treference_umats.cpp
///@brief Independent-physics validation of the modular adapters.
///
/// The legacy UMAT kernels removed in the modular rationalization are kept
/// VERBATIM under reference_kernels/ and compiled ONLY into this test. Each
/// adapter-served name is driven through a cyclic strain path side by side
/// with its reference kernel; the stress histories must match within the
/// tolerances established by the Phase-1 equivalence campaign (bit-identical
/// for the power-law/elastic families; Voce incremental-vs-closed-form
/// integration differences for the Chaboche family).

#include <gtest/gtest.h>
#include <armadillo>
#include <functional>
#include <string>

#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Umat/Modular/legacy_adapters.hpp>

#include "reference_kernels/elastic_isotropic.hpp"
#include "reference_kernels/elastic_transverse_isotropic.hpp"
#include "reference_kernels/elastic_orthotropic.hpp"
#include "reference_kernels/plastic_kin_iso_ccp.hpp"
#include "reference_kernels/Hill_isoh.hpp"
#include "reference_kernels/Hill_chaboche_ccp.hpp"
#include "reference_kernels/Ani_chaboche_ccp.hpp"
#include "reference_kernels/DFA_chaboche_ccp.hpp"
#include "reference_kernels/Generic_chaboche_ccp.hpp"
#include "reference_kernels/Hill_isoh_Nfast.hpp"

using namespace arma;
using namespace simcoon;

namespace {

// Standard-signature UMAT entry (reference kernel or adapter).
using UmatFn = std::function<void(
    const std::string&, const vec&, const vec&, vec&, mat&, mat&, const mat&,
    const int&, const vec&, const int&, vec&, const double&, const double&,
    const double&, const double&, double&, double&, double&, double&,
    const int&, const int&, const bool&, double&, const int&)>;

// Drive one kernel through a cyclic confined-uniaxial strain path
// (E11: 0 -> +amp -> -amp -> 0, other components held at zero) and return the
// sigma11 history. Continuum tangent on both sides — NOTE the reference
// kernels are verbatim PRE-2.0 copies, so their continuum mode is the OLD
// literal 0, not tangent_continuum (=1, which they would read as algorithmic).
// The mode only selects the Lt output in these kernels (stress is unaffected),
// but the intent is a continuum-vs-continuum comparison.
vec drive(const UmatFn& umat, const std::string& name, const vec& props,
          int nstatev, double amp, int n_quarter, int tangent_mode) {
    vec Etot = zeros(6);
    vec sigma = zeros(6);
    vec statev = zeros(nstatev);
    mat Lt(6, 6), L(6, 6);
    const mat DR = eye(3, 3);
    double T = 293., DT = 0., Time = 0., DTime = 1.;
    double Wm = 0., Wm_r = 0., Wm_ir = 0., Wm_d = 0.;
    double tnew_dt = 1.;

    const int n_total = 4 * n_quarter;
    vec s11(n_total);
    const double de = amp / n_quarter;
    bool start = true;
    for (int i = 0; i < n_total; ++i) {
        // quarter 1: load +, quarters 2-3: unload to -amp, quarter 4: back to 0
        const double d = (i < n_quarter) ? de
                       : (i < 3 * n_quarter) ? -de
                                             : de;
        vec DEtot = zeros(6);
        DEtot(0) = d;
        umat(name, Etot, DEtot, sigma, Lt, L, DR,
             static_cast<int>(props.n_elem), props, nstatev, statev,
             T, DT, Time, DTime, Wm, Wm_r, Wm_ir, Wm_d,
             3, 3, start, tnew_dt, tangent_mode);
        Etot += DEtot;
        Time += DTime;
        start = false;
        s11(i) = sigma(0);
    }
    return s11;
}

void expect_equiv(const UmatFn& reference, const std::string& name,
                  const vec& props, int nstatev, double rel_tol,
                  double amp = 0.02, int n_quarter = 25) {
    const vec ref = drive(reference, name, props, nstatev, amp, n_quarter,
                          /*pre-2.0 continuum literal*/ 0);
    const vec mod = drive(umat_legacy_modular, name, props, nstatev, amp,
                          n_quarter, tangent_continuum);
    const double peak = std::max(std::abs(ref.max()), std::abs(ref.min()));
    ASSERT_GT(peak, 0.);
    const double dev = abs(ref - mod).max() / peak;
    EXPECT_LT(dev, rel_tol) << name << ": adapter deviates from the reference "
                            << "kernel by " << dev;
}

}  // namespace

TEST(Treference_umats, ELISO_matches_reference) {
    expect_equiv(umat_elasticity_iso, "ELISO", {210000., 0.3, 1.2e-5}, 1, 1e-9);
}

TEST(Treference_umats, ELIST_matches_reference) {
    expect_equiv(umat_elasticity_trans_iso, "ELIST",
                 {3, 230000., 15000., 0.02, 0.4, 50000., 0., 0.}, 1, 1e-9);
}

TEST(Treference_umats, ELORT_matches_reference) {
    expect_equiv(umat_elasticity_ortho, "ELORT",
                 {70000., 30000., 15000., 0.3, 0.3, 0.3,
                  8000., 6000., 5000., 1e-5, 2e-5, 3e-5}, 1, 1e-9);
}

TEST(Treference_umats, EPKCP_matches_reference) {
    // m = 1 avoids the power-law p=0 tangent singularity (proven bit-identical)
    expect_equiv(umat_plasticity_kin_iso_CCP, "EPKCP",
                 {210000., 0.3, 0., 300., 1000., 1.0, 20000.}, 33, 1e-9);
}

TEST(Treference_umats, EPHIL_matches_reference) {
    expect_equiv(umat_plasticity_hill_isoh_CCP, "EPHIL",
                 {210000., 0.3, 0., 300., 5000., 1.0,
                  0.5, 0.4, 0.6, 1.5, 1.5, 1.5}, 33, 1e-9);
}

TEST(Treference_umats, EPHAC_matches_reference) {
    // Voce: legacy incremental Hp vs modular closed form -> small drift
    expect_equiv(umat_hill_chaboche_CCP, "EPHAC",
                 {210000., 0.3, 85000., 0., 300., 200., 20.,
                  30000., 172., 19500., 301.,
                  0.5, 0.4, 0.6, 1.5, 1.5, 1.5}, 33, 2e-3);
}

TEST(Treference_umats, EPANI_matches_reference) {
    expect_equiv(umat_ani_chaboche_CCP, "EPANI",
                 {210000., 0.3, 85000., 0., 300., 200., 20.,
                  30000., 172., 19500., 301.,
                  1.2, 1.1, 1.1, -0.6, -0.6, -0.5, 1.6, 1.5, 1.4}, 33, 2e-3);
}

TEST(Treference_umats, EPDFA_matches_reference) {
    expect_equiv(umat_dfa_chaboche_CCP, "EPDFA",
                 {210000., 0.3, 85000., 0., 300., 200., 20.,
                  30000., 172., 19500., 301.,
                  0.5, 0.4, 0.6, 1.5, 1.5, 1.5, 0.1}, 33, 2e-3);
}

TEST(Treference_umats, EPCHG_matches_reference) {
    expect_equiv(umat_generic_chaboche_CCP, "EPCHG",
                 {210000., 0.3, 85000., 0., 300., 2, 2, 0,
                  150., 15., 50., 40., 30000., 172., 19500., 301.}, 33, 2e-3);
}

TEST(Treference_umats, EPHIN_matches_reference) {
    // N = 1 only: the reference kernel is defective for N >= 2 (returns NaN
    // even for identical or never-active second surfaces) — the modular
    // engine is the correct implementation there.
    expect_equiv(umat_plasticity_hill_isoh_CCP_N, "EPHIN",
                 {210000., 0.3, 0., 1,
                  300., 3000., 1.0, 0.5, 0.4, 0.6, 1.5, 1.5, 1.5}, 33, 1e-9);
}
