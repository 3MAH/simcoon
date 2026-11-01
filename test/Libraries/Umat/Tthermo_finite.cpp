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

///@file Tthermo_finite.cpp
///@brief The thermomechanical dispatch under finite strain (run_umat_T, control types 2-6).
///@version 1.0

#include <gtest/gtest.h>
#include <map>
#include <string>
#include <armadillo>

#include <simcoon/parameter.hpp>
#include <simcoon/exception.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Simulation/Phase/phase_characteristics.hpp>
#include <simcoon/Simulation/Phase/state_variables_T.hpp>
#include <simcoon/Simulation/Maths/rotation.hpp>
#include <simcoon/Continuum_mechanics/Umat/umat_smart.hpp>

using namespace std;
using namespace arma;
using namespace simcoon;

namespace {

const double rho = 7.8e-9, c_p = 4.6e8, E_ = 210000., nu_ = 0.3, alpha = 1.2e-5;

struct law { vec props; int nstatev; double T0; };

const std::map<string, law> &laws() {
    static const std::map<string, law> l = {
        {"ELISO", {{rho, c_p, E_, nu_, alpha}, 1, 293.15}},
        {"ELIST", {{rho, c_p, 1, 140000., 9000., 0.3, 0.4, 4500., 1e-6, 3e-5}, 1, 293.15}},
        {"ELORT", {{rho, c_p, 140000., 9000., 9500., 0.3, 0.32, 0.4, 4500., 4600., 3000., 1e-6, 3e-5, 2.5e-5}, 1, 293.15}},
        {"EPICP", {{rho, c_p, E_, nu_, alpha, 300., 500., 0.4}, 8, 293.15}},
        {"EPKCP", {{rho, c_p, E_, nu_, alpha, 300., 500., 0.4, 20000.}, 14, 293.15}},
        {"EPJCK", {{7.85e-9, 4.75e8, 200000.0, 0.33, 1.0e-5, 792.0, 510.0, 0.26, 0.014, 1.0, 1.03, 293.0, 1793.0}, 9, 293.15}},
        {"ZENER", {{4.4, 0.656, 3000., 0.35, 1e-4, 1500., 0.35, 3000., 1200.}, 8, 293.15}},
        {"PRONK", {{4.4, 0.656, 3000., 0.35, 1e-4, 1., 1500., 0.35, 3000., 1200.}, 14, 293.15}},
        {"SMADI", {{5.5, 0.5, 0.5, 0, 67538., 67538., 0.349, 0.349, 1e-6, 1e-6, 0., 0.0418, 0.021, 0., 10., 10.,
                    300., 290., 295., 305., 0.2, 0.2, 0.2, 0.2, 300., 0.2, 2., 1e-6, 1e-3, 1., 1e8}, 17, 320.}},
    };
    return l;
}

// One increment of the law from the virgin state: strain increment De (both Green-Lagrange and
// box slots), DT, under control type ct and corate; the start frame R and F0 = F1 = F
shared_ptr<state_variables_T> one_increment(const string &name, const vec &De, const double &DT, const unsigned int &ct,
                                            const int &corate, const mat &R = eye(3,3), const mat &F = eye(3,3)) {
    const law &L = laws().at(name);
    phase_characteristics rve;
    rve.construct(0, 2);
    rve.sptr_matprops->update(0, name, 1., 0., 0., 0., L.props.n_elem, L.props);
    auto sv = std::dynamic_pointer_cast<state_variables_T>(rve.sptr_sv_global);
    sv->update(zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6),
               F, F, eye(3,3), eye(3,3), R, eye(3,3), L.T0, 0., L.nstatev, zeros(L.nstatev), zeros(L.nstatev));
    sv->dSdE = zeros(6,6);
    sv->dSdT = zeros(6,1);
    sv->tangent_mode = 2;
    sv->DEtot = De;
    sv->Detot = De;
    sv->DT = DT;
    bool start = true;
    double tnew_dt = 1.;
    run_umat_T(rve, eye(3,3), 0., 0.1, 3, 3, start, 0, ct, corate, tnew_dt);
    EXPECT_GE(tnew_dt, 1.);
    return sv;
}

void expect_same(const shared_ptr<state_variables_T> &a, const shared_ptr<state_variables_T> &b, const string &name) {
    const double s = std::max(norm(a->sigma, "inf"), 1.);
    EXPECT_LT(norm(a->tau - b->tau, "inf"), 1e-12*s) << name;
    EXPECT_LT(norm(a->dSdE - b->dSdE, "inf"), 1e-12*std::max(norm(a->dSdE, "inf"), 1.)) << name;
    EXPECT_LT(norm(vectorise(a->dSdT) - vectorise(b->dSdT), "inf"), 1e-12*std::max(norm(vectorise(a->dSdT), "inf"), 1e-12)) << name;
    EXPECT_LT(fabs(a->r - b->r), 1e-12*std::max(fabs(a->r), 1e-12)) << name;
    // Wm(0), Wm(1): under log_R the work is replaced by the stress power of F0 -> F1, here none
    EXPECT_LT(norm(a->Wm.tail(2) - b->Wm.tail(2), "inf"), 1e-12*std::max(norm(a->Wm, "inf"), 1e-12)) << name;
    EXPECT_LT(norm(a->Wt - b->Wt, "inf"), 1e-12*std::max(norm(a->Wt, "inf"), 1e-12)) << name;
}

const vec De_probe = {4e-3, -1e-3, 5e-4, 2e-3, -1e-3, 5e-4};

}  // namespace

// F = I: the finite dispatch has nothing to transport or rotate, so every law answers as under
// small strain, in both forms (the kernel strain is the one handed over, F carries none of it)
TEST(ThermoFinite, AtTheIdentityEveryLawIsItsSmallStrainSelf) {
    for (const auto &kv : laws()) {
        const auto small = one_increment(kv.first, De_probe, 2., 1, 0);
        expect_same(small, one_increment(kv.first, De_probe, 2., 2, 3), kv.first + " rate");
        expect_same(small, one_increment(kv.first, De_probe, 2., 2, corate_type_total), kv.first + " total");
    }
}

// The kernels run in the material frame R_hat = DR R: a law whose frame has turned by Q answers
// with the stress turned by Q, to round-off, the anisotropic ones included; heat and works unchanged
TEST(ThermoFinite, TheKernelsRunInTheMaterialFrame) {
    Rotation rot = Rotation::from_euler(0.7, -0.4, 1.1, "zxz");
    const mat Q = rot.as_matrix();
    for (const auto &kv : laws()) {
        const auto ref = one_increment(kv.first, De_probe, 2., 2, 3);
        const auto turned = one_increment(kv.first, rotate_strain(De_probe, Q), 2., 2, 3, Q);
        const double s = std::max(norm(ref->tau, "inf"), 1.);
        EXPECT_LT(norm(turned->tau - rotate_stress(ref->tau, Q), "inf"), 1e-10*s) << kv.first;
        EXPECT_LT(fabs(turned->r - ref->r), 1e-10*std::max(fabs(ref->r), 1e-12)) << kv.first;
        EXPECT_LT(norm(turned->dSdE - rotate_stiffness(ref->dSdE, Q), "inf"), 1e-9*norm(ref->dSdE, "inf")) << kv.first;
    }
}

// log_R and the total form only: the other rates would feed the kernels a strain their entropy
// and energies were not written for
TEST(ThermoFinite, OtherRatesAreRefused) {
    for (const int corate : {0, 1, 2, 4, 5})
        EXPECT_THROW(one_increment("ELISO", De_probe, 0., 2, corate), simcoon::exception_solver) << corate;
}
