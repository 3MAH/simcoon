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

///@file Tdamage_LLD_0.cpp
///@brief The Ladeveze-Le Dantec damage UMAT (LLDM0): tangent against finite differences, energy
///       split, irreversibility, transverse failure, fibre axis
///@version 1.0

#include <gtest/gtest.h>
#include <string>
#include <armadillo>

#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Simulation/Maths/rotation.hpp>
#include <simcoon/Simulation/Maths/lagrange.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Damage/damage_LLD_0.hpp>

using namespace std;
using namespace arma;
using namespace simcoon;

namespace {

// Carbon/epoxy: nuTL is the minor Poisson ratio (nu_LT E_T / E_L)
vec props_LLD(const double &axis = 1, const double &Y_22_u = 1.2) {
    return {axis, 138000., 11000., 0.0223, 0.40, 5500., -0.3e-6, 30e-6,
            0.10, 0.30, 0.15, 0.40, Y_22_u, 2.5, 2.0, 50., 0.5, 300.};
}

struct point {
    vec props;
    vec e = zeros(6), s = zeros(6), sv = zeros(10);
    double Wm[4] = {0., 0., 0., 0.};
    bool started = false;
    explicit point(const vec &p) : props(p) {}

    // one increment; commit = false leaves the point unchanged
    vec step(const vec &de, mat &Lt, mat &L, const bool &commit = true) {
        vec s1 = s, sv1 = sv;
        double W[4] = {Wm[0], Wm[1], Wm[2], Wm[3]};
        double tnew_dt = 1.;
        umat_damage_LLD_0("LLDM0", e, de, s1, Lt, L, eye(3, 3), props.n_elem, props, 10, sv1, 293.15, 0., 1., 0.01,
                          W[0], W[1], W[2], W[3], 3, 3, !started, tnew_dt);
        if (commit) {
            e += de; s = s1; sv = sv1; started = true;
            for (int i = 0; i < 4; i++) Wm[i] = W[i];
        }
        return s1;
    }
    void path(const vec &dir, const int &n, const double &h) {
        mat Lt, L;
        for (int k = 0; k < n; k++) step(dir*h, Lt, L);
    }
    // relative error of the tangent of the next increment de against central differences
    double tangent_error(const vec &de) {
        mat Lt, L, X, Y;
        step(de, Lt, L, false);
        mat fd(6, 6);
        const double h = 1.e-9;
        for (int j = 0; j < 6; j++) {
            vec dp = de, dm = de;
            dp(j) += h;
            dm(j) -= h;
            fd.col(j) = (step(dp, X, Y, false) - step(dm, X, Y, false))/(2.*h);
        }
        return norm(Lt - fd)/norm(fd);
    }
};

const vec transverse = {0., 1., 0., 0., 0., 0.};
const vec shear = {0., 0., 0., 1., 0., 0.};
const vec mixed = {0., 0.6, 0., 1., 0., 0.};
const vec mixed_neg = {0., 0.6, 0., -1., 0., 0.};

}  // namespace

TEST(Tdamage_LLD_0, tangent_is_exact_while_only_damage_evolves)
{
    // d22 and d12 coupled through Y_ts and the stress, the plasticity not yet active: the
    // tangent of the 2x2 damage system is exact
    for (const vec &dir : {shear, mixed, mixed_neg}) {
        point pt(props_LLD());
        pt.path(dir, 10, 2.5e-4);
        ASSERT_GT(pt.sv(2), 0.05);
        ASSERT_LT(pt.sv(3), 1.e-6);
        EXPECT_LT(pt.tangent_error(dir*2.5e-4), 1.e-6);
    }
}

TEST(Tdamage_LLD_0, tangent_with_plasticity_converges_with_the_increment)
{
    // With plasticity the direction is held fixed over the increment (continuum tangent): the
    // error is first order in the increment
    for (const vec &dir : {transverse, shear, mixed}) {
        point pt(props_LLD());
        pt.path(dir, 24, 2.5e-4);
        ASSERT_GT(pt.sv(3), 1.e-4);
        const double e4 = pt.tangent_error(dir*1.e-4), e5 = pt.tangent_error(dir*1.e-5);
        EXPECT_LT(e5, 1.e-3);
        EXPECT_LT(e5, 0.2*e4);
    }
}

TEST(Tdamage_LLD_0, unloading_is_elastic_on_the_damaged_stiffness)
{
    point pt(props_LLD());
    pt.path(mixed, 16, 2.5e-4);
    const vec sv_loaded = pt.sv;
    mat Lt, L;
    pt.step(-mixed*2.5e-4, Lt, L);
    EXPECT_LT(norm(pt.sv.subvec(1, 9) - sv_loaded.subvec(1, 9)), 1.e-12);   // nothing evolves
    EXPECT_LT(norm(Lt - L), 1.e-10*norm(L));
    EXPECT_LT(pt.tangent_error(-mixed*2.5e-4), 1.e-6);
}

TEST(Tdamage_LLD_0, work_split_and_irreversibility)
{
    // Wm = Wm_r + Wm_ir + Wm_d at every increment, the dissipation and the damage never decrease,
    // over loading, unloading and reloading
    point pt(props_LLD());
    double Wd = 0., d22 = 0., d12 = 0.;
    mat Lt, L;
    for (int k = 0; k < 60; k++) {
        const double sgn = (k < 25) ? 1. : (k < 40 ? -1. : 1.);
        pt.step(mixed*sgn*2.5e-4, Lt, L);
        EXPECT_NEAR(pt.Wm[0], pt.Wm[1] + pt.Wm[2] + pt.Wm[3], 1.e-12*fabs(pt.Wm[0]) + 1.e-14);
        EXPECT_GE(pt.Wm[3], Wd - 1.e-12);
        EXPECT_GE(pt.sv(1), d22);
        EXPECT_GE(pt.sv(2), d12);
        Wd = pt.Wm[3];
        d22 = pt.sv(1);
        d12 = pt.sv(2);
    }
    EXPECT_GT(pt.Wm[3], 0.);
    EXPECT_GT(pt.Wm[2], 0.);   // stored hardening
}

TEST(Tdamage_LLD_0, transverse_failure_breaks_d22_only)
{
    // Past Y_22_u the transverse damage jumps to 0.99 and stays; d12 keeps its own law,
    // d12 + lambda(d12) = <Y_ts - Y_12_0> / Y_12_c, while it grows
    const vec p = props_LLD(1, 0.08);
    point pt(p);
    pt.path(transverse, 20, 2.5e-4);
    ASSERT_NEAR(pt.sv(1), 0.99, 1.e-12);
    const double d22 = pt.sv(1), d12 = pt.sv(2);
    ASSERT_GT(d12, 0.1);
    ASSERT_LT(d12, 0.99);
    const double s22 = Macaulay_p(pt.s(1)), s33 = Macaulay_p(pt.s(2));
    const double Y22 = 0.5*(s22*s22 + s33*s33)/(p(2)*(1. - d22)*(1. - d22));
    const double Y12 = 0.5*(pt.s(3)*pt.s(3) + pt.s(4)*pt.s(4))/(p(5)*(1. - d12)*(1. - d12));
    const double Y_ts = sqrt(Y12 + p(13)*Y22);
    EXPECT_NEAR(d12 + lagrange_pow_1(d12, 1.e-2, 1.e-1, 1., 0.), Macaulay_p(Y_ts - p(8))/p(9), 1.e-6);
    pt.path(shear, 10, 2.5e-3);
    EXPECT_NEAR(pt.sv(1), 0.99, 1.e-12);
}

TEST(Tdamage_LLD_0, plastic_flow_is_associated)
{
    // The plastic strain increment is the gradient of the criterion in effective stress,
    // Theta eta(Theta sigma~), sigma~ = L0 : (eps - eps^p); the cutting-plane iterations turn the
    // direction over the increment, hence a small one
    point pt(props_LLD());
    pt.path(mixed, 24, 2.5e-4);
    const vec EP0 = pt.sv.subvec(4, 9);
    mat Lt, L;
    pt.step(mixed*1.e-8, Lt, L);
    const vec DEP = pt.sv.subvec(4, 9) - EP0;
    const vec p = props_LLD();
    const mat L0 = L_isotrans(p(1), p(2), p(3), p(4), p(5), 1);
    mat Theta = diagmat(vec({0., p(14), p(14), 1., 1., 0.}));
    const vec n = Theta*eta_stress(vec(Theta*(L0*(pt.e - pt.sv.subvec(4, 9)))));
    ASSERT_GT(norm(DEP), 0.);
    EXPECT_LT(norm(DEP/norm(DEP) - n/norm(n)), 1.e-5);
}

TEST(Tdamage_LLD_0, fibre_axis_rotates_the_material)
{
    // The same loading in the frames with the fibres along 2 and 3: the response is the rotated one
    const mat R2 = {{0., -1., 0.}, {1., 0., 0.}, {0., 0., 1.}};
    const mat R3 = {{0., 0., -1.}, {0., 1., 0.}, {1., 0., 0.}};
    for (const auto &[axis, R] : {std::pair<double, mat>{2., R2}, std::pair<double, mat>{3., R3}}) {
        point a(props_LLD(1)), b(props_LLD(axis));
        mat Lt_a, L_a, Lt_b, L_b;
        for (int k = 0; k < 24; k++) {
            a.step(mixed*2.5e-4, Lt_a, L_a);
            b.step(rotate_strain(vec(mixed*2.5e-4), R), Lt_b, L_b);
        }
        EXPECT_LT(norm(b.s - rotate_stress(a.s, R)), 1.e-9*norm(a.s));
        EXPECT_LT(norm(Lt_b - rotate_stiffness(Lt_a, R)), 1.e-9*norm(Lt_a));
        EXPECT_LT(norm(b.sv.subvec(1, 3) - a.sv.subvec(1, 3)), 1.e-12);
        EXPECT_LT(norm(b.sv.subvec(4, 9) - rotate_strain(vec(a.sv.subvec(4, 9)), R)), 1.e-12);
    }
    point bad(props_LLD(4));
    mat Lt, L;
    EXPECT_THROW(bad.step(mixed*2.5e-4, Lt, L), std::invalid_argument);
}

