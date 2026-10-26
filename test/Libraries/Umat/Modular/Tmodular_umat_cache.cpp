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

///@file Tmodular_umat_cache.cpp
///@brief umat_modular keeps configured ModularUMAT instances per thread: every call must
///       return, bit for bit, what a fresh instance returns (the behaviour before the cache)

#include <gtest/gtest.h>
#include <armadillo>
#include <atomic>
#include <cstring>
#include <stdexcept>
#include <thread>
#include <vector>

#include <simcoon/Continuum_mechanics/Umat/Modular/modular_umat.hpp>
#include <simcoon/Continuum_mechanics/Umat/Modular/damage_mechanism.hpp>
#include <simcoon/Continuum_mechanics/Umat/Modular/elasticity_module.hpp>
#include <simcoon/Continuum_mechanics/Umat/Modular/hardening.hpp>
#include <simcoon/Continuum_mechanics/Umat/Modular/strain_mechanism.hpp>
#include <simcoon/Continuum_mechanics/Umat/Modular/yield_criterion.hpp>
#include <simcoon/parameter.hpp>

using namespace std;
using namespace arma;
using namespace simcoon;

namespace {

double code(ElasticityType t) { return static_cast<double>(t); }
double code(MechanismType t) { return static_cast<double>(t); }
double code(YieldType t) { return static_cast<double>(t); }
double code(IsoHardType t) { return static_cast<double>(t); }
double code(KinHardType t) { return static_cast<double>(t); }
double code(DamageType t) { return static_cast<double>(t); }

struct material {
    vec props;
    int nstatev;
    int tangent_mode;
};

// Elastoplastic (von Mises, Voce) + linear damage: the damage compliance cache M_cached_.
material plastic_damage(double sigma_Y) {
    return {vec{code(ElasticityType::ISOTROPIC), 0.0, 210000.0, 0.3, 1e-5,
                2.0,
                code(MechanismType::PLASTICITY), code(YieldType::VON_MISES), code(IsoHardType::VOCE),
                code(KinHardType::NONE), 1.0, 0.0, sigma_Y, 200.0, 10.0,
                code(MechanismType::DAMAGE), code(DamageType::LINEAR), 0.05, 5.0},
            10, tangent_algorithmic};
}

// Linear isotropic + Armstrong-Frederick kinematic hardening, closest-point integrator.
material plastic_af() {
    return {vec{code(ElasticityType::ISOTROPIC), 0.0, 70000.0, 0.33, 0.0,
                1.0,
                code(MechanismType::PLASTICITY), code(YieldType::VON_MISES), code(IsoHardType::LINEAR),
                code(KinHardType::ARMSTRONG_FREDERICK), 1.0, 1.0, 150.0, 500.0, 20000.0, 100.0},
            14, tangent_closest_point};
}

// One Prony branch + power-law plasticity: the viscoelastic dEV/dE map of predict().
material visco_plastic() {
    return {vec{code(ElasticityType::ISOTROPIC), 0.0, 3000.0, 0.4, 0.0,
                2.0,
                code(MechanismType::VISCOELASTICITY), 1.0, 2000.0, 0.4, 1.0e4, 1.0e4,
                code(MechanismType::PLASTICITY), code(YieldType::VON_MISES), code(IsoHardType::POWER_LAW),
                code(KinHardType::NONE), 1.0, 0.0, 40.0, 50.0, 0.4},
            15, tangent_continuum};
}

struct point {
    vec Etot = zeros(6), sigma = zeros(6), statev;
    mat Lt = zeros(6, 6), L = zeros(6, 6);
    double Wm = 0., Wm_r = 0., Wm_ir = 0., Wm_d = 0., tnew_dt = 1.;
    explicit point(int nstatev) : statev(zeros(nstatev)) {}
};

// One increment through umat_modular (the cached instances) ...
void call_cached(const material& m, point& p, const vec& DEtot, double Time, bool start) {
    p.tnew_dt = 1.;
    umat_modular("MODUL", p.Etot, DEtot, p.sigma, p.Lt, p.L, eye(3, 3), static_cast<int>(m.props.n_elem),
                 m.props, m.nstatev, p.statev, 293., 0., Time, 0.1, p.Wm, p.Wm_r, p.Wm_ir, p.Wm_d,
                 3, 3, start, p.tnew_dt, m.tangent_mode);
    p.Etot += DEtot;
}

// ... and through a fresh instance, what umat_modular did before it kept them.
void call_fresh(const material& m, point& p, const vec& DEtot, double Time, bool start) {
    p.tnew_dt = 1.;
    ModularUMAT fresh;
    fresh.set_kinematic_form(KinematicForm::rate);
    fresh.run("MODUL", p.Etot, DEtot, p.sigma, p.Lt, p.L, eye(3, 3), static_cast<int>(m.props.n_elem),
              m.props, m.nstatev, p.statev, 293., 0., Time, 0.1, p.Wm, p.Wm_r, p.Wm_ir, p.Wm_d,
              3, 3, start, p.tnew_dt, m.tangent_mode);
    p.Etot += DEtot;
}

template <typename T>
bool same_bits(const T& a, const T& b) {
    return a.n_rows == b.n_rows && a.n_cols == b.n_cols
        && std::memcmp(a.memptr(), b.memptr(), a.n_elem * sizeof(double)) == 0;
}

bool same_bits(const point& a, const point& b) {
    const double wa[5] = {a.Wm, a.Wm_r, a.Wm_ir, a.Wm_d, a.tnew_dt};
    const double wb[5] = {b.Wm, b.Wm_r, b.Wm_ir, b.Wm_d, b.tnew_dt};
    return same_bits(a.sigma, b.sigma) && same_bits(a.statev, b.statev) && same_bits(a.Lt, b.Lt)
        && same_bits(a.L, b.L) && std::memcmp(wa, wb, sizeof(wa)) == 0;
}

// Strain increment of a load / unload history, scaled per material.
vec increment(int inc, double scale) {
    const vec dir = {1.0, -0.35, -0.25, 0.4, -0.2, 0.3};
    return (inc < 8 ? 1.0 : -1.0) * scale * dir;
}

// Interleave the materials call by call (as the phases of a multiphase material do), each with
// two material points; the second point starts mid-history (start = true on an instance that
// already holds the first point's plastic state). Returns the number of mismatching calls.
int interleaved_mismatches(const std::vector<material>& mats, const std::vector<double>& scales) {
    std::vector<point> cached, fresh;
    for (const auto& m : mats) {
        for (int k = 0; k < 2; ++k) {
            cached.emplace_back(m.nstatev);
            fresh.emplace_back(m.nstatev);
        }
    }
    int mismatches = 0;
    for (int inc = 0; inc < 12; ++inc) {
        for (size_t i = 0; i < mats.size(); ++i) {
            for (int k = 0; k < 2; ++k) {
                const int first = (k == 0) ? 0 : 5;
                if (inc < first) continue;
                const size_t ip = 2 * i + k;
                const vec DEtot = increment(inc, scales[i]);
                call_cached(mats[i], cached[ip], DEtot, 0.1 * inc, inc == first);
                call_fresh(mats[i], fresh[ip], DEtot, 0.1 * inc, inc == first);
                if (!same_bits(cached[ip], fresh[ip])) ++mismatches;
            }
        }
    }
    return mismatches;
}

}  // namespace

// Three configurations alternating at every call, with damage, viscoelasticity, the cutting-plane
// and closest-point integrators: the cached instances return what fresh ones do, bit for bit.
TEST(ModularUMATCache, AlternatingPropsMatchFreshInstances) {
    const std::vector<material> mats = {plastic_damage(300.0), plastic_af(), visco_plastic()};
    const std::vector<double> scales = {1.2e-3, 1.5e-3, 4.0e-3};

    // the histories are not trivial: plasticity, damage and viscous flow all activate
    std::vector<point> pts;
    for (size_t i = 0; i < mats.size(); ++i) {
        pts.emplace_back(mats[i].nstatev);
        for (int inc = 0; inc < 8; ++inc) {
            call_fresh(mats[i], pts[i], increment(inc, scales[i]), 0.1 * inc, inc == 0);
        }
    }
    EXPECT_GT(pts[0].statev(1), 0.0) << "plasticity (damage material) never activated";
    EXPECT_GT(pts[0].statev(8), 0.0) << "damage never activated";
    EXPECT_GT(pts[1].statev(1), 0.0) << "plasticity (AF material) never activated";
    EXPECT_GT(pts[2].statev(8), 0.0) << "plasticity (viscoplastic material) never activated";

    EXPECT_EQ(interleaved_mismatches(mats, scales), 0);
}

// More configurations than cache entries: every call evicts the least recently used one.
TEST(ModularUMATCache, MoreConfigurationsThanEntries) {
    std::vector<material> mats;
    std::vector<double> scales;
    for (int i = 0; i < 11; ++i) {
        mats.push_back(plastic_damage(250.0 + 15.0 * i));
        scales.push_back(1.2e-3);
    }
    EXPECT_EQ(interleaved_mismatches(mats, scales), 0);
}

// The instances are per thread: concurrent interleaved histories stay exact.
TEST(ModularUMATCache, ConcurrentThreads) {
    const std::vector<material> mats = {plastic_damage(300.0), plastic_af(), visco_plastic()};
    const std::vector<double> scales = {1.2e-3, 1.5e-3, 4.0e-3};
    std::atomic<int> mismatches{0};
    std::vector<std::thread> threads;
    for (int t = 0; t < 4; ++t) {
        threads.emplace_back([&]() { mismatches += interleaved_mismatches(mats, scales); });
    }
    for (auto& th : threads) th.join();
    EXPECT_EQ(mismatches.load(), 0);
}

// A call that throws must not leave its instance in the cache: here initialize() fails half-way
// (variables registered, statev too short to unpack), and reusing that instance would register
// the variables a second time.
TEST(ModularUMATCache, FailedCallIsNotReused) {
    const material m = plastic_damage(333.0);   // props used by no other test: a new entry
    point bad(1);   // nstatev says 10, statev holds 1 value
    EXPECT_THROW(call_cached(m, bad, increment(0, 1.2e-3), 0.0, true), std::runtime_error);

    point cached(m.nstatev), fresh(m.nstatev);
    for (int inc = 0; inc < 10; ++inc) {
        ASSERT_NO_THROW(call_cached(m, cached, increment(inc, 1.2e-3), 0.1 * inc, inc == 0));
        call_fresh(m, fresh, increment(inc, 1.2e-3), 0.1 * inc, inc == 0);
        EXPECT_TRUE(same_bits(cached, fresh)) << "increment " << inc;
    }

    // the same on an instance already in use: the failing call drops it, the next one rebuilds it
    EXPECT_THROW(call_cached(m, bad, increment(3, 1.2e-3), 0.3, false), std::runtime_error);
    for (int inc = 10; inc < 12; ++inc) {
        ASSERT_NO_THROW(call_cached(m, cached, increment(inc, 1.2e-3), 0.1 * inc, false));
        call_fresh(m, fresh, increment(inc, 1.2e-3), 0.1 * inc, false);
        EXPECT_TRUE(same_bits(cached, fresh)) << "increment " << inc;
    }
}
