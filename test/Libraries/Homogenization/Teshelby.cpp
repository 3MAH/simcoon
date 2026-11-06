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

///@file Teshelby.cpp
///@brief Test for Eshelby and Interaction tensors
///@version 1.0

#include <gtest/gtest.h>
#include <armadillo>

#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Homogenization/eshelby.hpp>

using namespace std;
using namespace arma;
using namespace simcoon;

TEST(Teshelby, S_TII)
{
    
    double a1 = 1.;
    double a2 = 1.;
    double a3 = 1.;
    int mp = 300;
    int np = 300;
    
    double E = 70000.;
    double nu = 0.3;
    
    double mu = E/(2.*(1+nu));
    double lambda = E*nu/((1.+nu)*(1.-2.*nu));
    
    mat S_anal = zeros(6,6);
    mat S_num = zeros(6,6);
    mat T_II_num = zeros(6,6);
    
    mat Lt = zeros(6,6);
    Lt(0,0) = 2.*mu+lambda;
    Lt(0,1) = lambda;
    Lt(0,2) = lambda;
    Lt(1,0) = lambda;
    Lt(1,1) = 2.*mu+lambda;
    Lt(1,2) = lambda;
    Lt(2,0) = lambda;
    Lt(2,1) = lambda;
    Lt(2,2) = 2.*mu+lambda;
    Lt(3,3) = mu;
    Lt(4,4) = mu;
    Lt(5,5) = mu;
    
    vec x = zeros(mp);
    vec wx = zeros(mp);
    vec y = zeros(np);
    vec wy = zeros(np);
    points(x, wx, y, wy, mp, np);
    
    S_anal = Eshelby_sphere(nu);
    S_num = Eshelby(Lt, a1, a2, a3, x, wx, y, wy, mp, np);
    EXPECT_LT(norm(S_num-S_anal,2),1.E-9);

    cout << S_num-S_anal << "\n";
    
    T_II_num = T_II(Lt, a1, a2, a3, x, wx, y, wy, mp, np);
    EXPECT_LT(norm(T_II_num*Lt-S_anal,2),1.E-9);
    
    a1 = 1.E16;
    S_anal = Eshelby_cylinder(nu);
    S_num = Eshelby(Lt, a1, a2, a3, x, wx, y, wy, mp, np);
    EXPECT_LT(norm(S_num-S_anal,2),1.E-4);
    
    cout << S_num-S_anal << "\n";
    
    T_II_num = T_II(Lt, a1, a2, a3, x, wx, y, wy, mp, np);
    EXPECT_LT(norm(T_II_num*Lt-S_anal,2),1.E-4);
    
}

#include <simcoon/Continuum_mechanics/Homogenization/ellipsoid_multi.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>

// eshelby_memo: the key is the exact (medium stiffness, semi-axes); a hit returns the stored
// tensor, any change of medium or geometry misses.
TEST(Teshelby, memo_keys_on_exact_medium_and_geometry)
{
    eshelby_memo memo;
    mat L = L_iso(70000., 0.3, "Enu");
    EXPECT_EQ(memo.find(L, 1., 1., 1.), nullptr);
    mat S = eye(6, 6) * 0.5;
    memo.store(L, 1., 1., 1., S);
    ASSERT_NE(memo.find(L, 1., 1., 1.), nullptr);
    EXPECT_LT(norm(*memo.find(L, 1., 1., 1.) - S, 2), simcoon::iota);
    EXPECT_EQ(memo.find(L, 2., 1., 1.), nullptr);                 // other geometry
    mat L2 = L; L2(0, 0) += 1.e-9;
    EXPECT_EQ(memo.find(L2, 1., 1., 1.), nullptr);                 // other medium, however close
    EXPECT_EQ(memo.size(), 1u);
    memo.clear();
    EXPECT_EQ(memo.size(), 0u);
}

#include <thread>
#include <vector>

// eshelby_quadrature_points: the tables of points(), one shared immutable instance per (mp, np),
// whatever the thread that asks for it.
TEST(Teshelby, quadrature_tables_are_shared_and_exact)
{
    const eshelby_quadrature &q = eshelby_quadrature_points(20, 30);
    vec x(20), wx(20), y(30), wy(30);
    points(x, wx, y, wy, 20, 30);
    EXPECT_EQ(q.mp, 20);
    EXPECT_EQ(q.np, 30);
    EXPECT_TRUE(approx_equal(q.x, x, "absdiff", 0.) && approx_equal(q.wx, wx, "absdiff", 0.));
    EXPECT_TRUE(approx_equal(q.y, y, "absdiff", 0.) && approx_equal(q.wy, wy, "absdiff", 0.));
    EXPECT_EQ(&eshelby_quadrature_points(20, 30), &q);
    EXPECT_NE(&eshelby_quadrature_points(30, 20), &q);
    EXPECT_EQ(&eshelby_quadrature_points(20, 30), &q);
    EXPECT_THROW(eshelby_quadrature_points(0, 10), std::invalid_argument);

    // threads alternating (mp, np) all get the tables of the size they asked for
    std::vector<int> wrong(8, 0);
    std::vector<std::thread> threads;
    for (int t = 0; t < 8; t++) {
        threads.emplace_back([t, &wrong]() {
            for (int it = 0; it < 200; it++) {
                const int m = 5 + (t + it) % 7;
                const eshelby_quadrature &qt = eshelby_quadrature_points(m, m + 1);
                if (qt.mp != m || qt.np != m + 1 || qt.x.n_elem != (uword)m || qt.wy.n_elem != (uword)(m + 1))
                    wrong[t]++;
            }
        });
    }
    for (auto &th : threads) th.join();
    for (int t = 0; t < 8; t++) EXPECT_EQ(wrong[t], 0);
}

// Gauss is the weighted sum of calG over the quadrature points (it inlines calG and hoists the
// azimuthal trigonometry out of the loops): same sums in the same order.
TEST(Teshelby, Gauss_is_the_weighted_sum_of_calG)
{
    Mat<int> Id = {{0, 3, 4}, {3, 1, 5}, {4, 5, 2}};
    mat L_ortho = zeros(6, 6);
    L_ortho.submat(0, 0, 2, 2) = mat({{120000., 30000., 20000.}, {30000., 80000., 25000.}, {20000., 25000., 40000.}});
    L_ortho(3, 3) = 30000.; L_ortho(4, 4) = 20000.; L_ortho(5, 5) = 15000.;
    const std::vector<mat> media = {L_iso(70000., 0.3, "Enu"), L_ortho};
    const std::vector<vec> shapes = {{1., 1., 1.}, {3., 1., 0.5}, {10., 1., 1.}};
    const int mp = 13, np = 7;
    vec x(mp), wx(mp), y(np), wy(np);
    points(x, wx, y, wy, mp, np);
    for (const mat &L : media) {
        for (const vec &a : shapes) {
            mat G = zeros(6, 6);
            Gauss(Id, L, G, a(0), a(1), a(2), x, wx, y, wy, mp, np);
            mat G_ref = zeros(6, 6);
            mat H = zeros(6, 6);
            for (int l = 0; l < mp; l++) {
                for (int i = 0; i < np; i++) {
                    calG(y(i), a(0), a(1), a(2), x(l), Id, L, H);
                    G_ref += wx(l) * (H * wy(i));
                }
            }
            EXPECT_LT(norm(G - G_ref, "fro"), 1.E-14 * norm(G_ref, "fro"));
        }
    }
}

// Eshelby_closed_form: the sphere / spheroid formulas for an isotropic medium, and nothing else.
TEST(Teshelby, closed_form_applies_to_isotropic_spheroids_only)
{
    const double nu = 0.3;
    mat L = L_iso(70000., nu, "Enu");
    mat S;
    ASSERT_TRUE(Eshelby_closed_form(L, 1., 1., 1., S));
    EXPECT_LT(norm(S - Eshelby_sphere(nu), "fro"), 1.E-14);
    ASSERT_TRUE(Eshelby_closed_form(L, 4., 2., 2., S));
    EXPECT_LT(norm(S - Eshelby_prolate(nu, 2.), "fro"), 1.E-14);
    ASSERT_TRUE(Eshelby_closed_form(L, 0.5, 2., 2., S));
    EXPECT_LT(norm(S - Eshelby_oblate(nu, 0.25), "fro"), 1.E-14);
    // and they are the limit of the integral
    const eshelby_quadrature &q = eshelby_quadrature_points(200, 200);
    for (double ar : {3., 0.3}) {
        ASSERT_TRUE(Eshelby_closed_form(L, ar, 1., 1., S));
        mat S_num = Eshelby(L, ar, 1., 1., q.x, q.wx, q.y, q.wy, q.mp, q.np);
        EXPECT_LT(norm(S - S_num, "fro"), 1.E-9 * norm(S, "fro")) << "aspect ratio " << ar;
    }
    mat S0 = eye(6, 6);
    S = S0;
    EXPECT_FALSE(Eshelby_closed_form(L, 1., 2., 3., S));            // general ellipsoid
    EXPECT_FALSE(Eshelby_closed_form(L, 2., 2., 1., S));            // spheroid of axis 3
    EXPECT_FALSE(Eshelby_closed_form(L, 1.005, 1., 1., S));         // too close to a sphere
    mat L_ti = L;
    L_ti(0, 0) *= 1.01;
    EXPECT_FALSE(Eshelby_closed_form(L_ti, 2., 1., 1., S));         // anisotropic medium
    EXPECT_EQ(norm(S - S0, "fro"), 0.);
}

// The penny-shaped crack (normal 1) is the ar -> 0 limit of the oblate tensor and of the integral.
TEST(Teshelby, penny_is_the_flat_oblate_limit)
{
    const eshelby_quadrature &q = eshelby_quadrature_points(200, 200);
    for (double nu : {0., 0.25, 0.3, 0.45}) {
        const mat S = Eshelby_penny(nu);
        EXPECT_LT(norm(S - Eshelby_oblate(nu, 1.E-6), "fro"), 1.E-4) << "nu " << nu;
        const mat L = L_iso(70000., nu, "Enu");
        const mat S_num = Eshelby(L, 1.E-5, 1., 1., q.x, q.wx, q.y, q.wy, q.mp, q.np);
        EXPECT_LT(norm(S - S_num, "fro"), 1.E-3) << "nu " << nu;
    }
}

// Rule 1 (half sphere, trapezoidal rule in phi): its own tables, the same integral.
TEST(Teshelby, half_sphere_rule)
{
    const eshelby_quadrature &g = eshelby_quadrature_points(40, 40);
    const eshelby_quadrature &b = eshelby_quadrature_points(40, 40, 1);
    const eshelby_quadrature &c = eshelby_quadrature_points(40, 40, 1, true);
    EXPECT_NE(&g, &b);
    EXPECT_NE(&b, &c);
    EXPECT_EQ(&eshelby_quadrature_points(40, 40, 1), &b);
    EXPECT_EQ(b.mp, 20);
    EXPECT_EQ(b.np, 40);
    EXPECT_NEAR(accu(b.wx), 2., 1.E-11);                            // twice the measure of [0, 1]
    EXPECT_NEAR(accu(b.wy), 2. * simcoon::pi, 1.E-13);
    EXPECT_GE(b.x.min(), 0.);
    EXPECT_FALSE(g.closed_form);
    EXPECT_TRUE(c.closed_form);
    EXPECT_THROW(eshelby_quadrature_points(40, 40, 2), std::invalid_argument);

    const double nu = 0.3;
    mat L = L_iso(70000., nu, "Enu");
    mat S_b = Eshelby(L, 1., 1., 1., b.x, b.wx, b.y, b.wy, b.mp, b.np);
    EXPECT_LT(norm(S_b - Eshelby_sphere(nu), "fro"), 1.E-12);
    // a prolate spheroid: the half-sphere rule with half the points is as accurate
    mat S_ref = Eshelby_prolate(nu, 10.);
    mat S_g = Eshelby(L, 10., 1., 1., g.x, g.wx, g.y, g.wy, g.mp, g.np);
    mat S_10 = Eshelby(L, 10., 1., 1., b.x, b.wx, b.y, b.wy, b.mp, b.np);
    EXPECT_LT(norm(S_10 - S_ref, "fro"), 1.5 * norm(S_g - S_ref, "fro"));
}
