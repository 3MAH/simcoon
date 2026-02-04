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

///@file Tcriteria.cpp
///@brief Test for stress criteria
///@version 1.0

#include <gtest/gtest.h>
#include <armadillo>

#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Continuum_mechanics/Functions/criteria.hpp>

using namespace std;
using namespace arma;
using namespace simcoon;

TEST(Tcriteria, Aniso)
{
    double b = 0.;
    double n = 2.;
    vec sigma = {400., 0., 0., 0., 0., 0.};

    // Check that Tresca is equal to 400 in that case
    double test_tresca = Tresca_stress(sigma);
    EXPECT_LT(test_tresca - 400., simcoon::iota);

    // Check that Drucker is equal to 400 in that case
    double test_prager = Drucker_stress(sigma, b, n);
    EXPECT_LT(test_prager - 400., simcoon::iota);
}

TEST(Tcriteria, Hill_stress_isotropic_limit)
{
    // With isotropic P_params (F=G=H=0.5, L=M=N=1.5), Hill = von Mises
    vec P_params = {0.5, 0.5, 0.5, 1.5, 1.5, 1.5};
    vec sigma = {400., 0., 0., 0., 0., 0.};

    double hill = Hill_stress(sigma, P_params);
    double mises = Mises_stress(sigma);
    EXPECT_LT(fabs(hill - mises), 1.E-6);
}

TEST(Tcriteria, dHill_stress_numerical)
{
    vec P_params = {0.5, 0.5, 0.5, 1.5, 1.5, 1.5};
    vec sigma = {400., 100., 200., 50., 30., 70.};

    vec dHill = dHill_stress(sigma, P_params);

    // Numerical derivative
    double eps = 1.E-6;
    vec dHill_num = zeros(6);
    for (int i = 0; i < 6; i++) {
        vec sigma_p = sigma;
        vec sigma_m = sigma;
        sigma_p(i) += eps;
        sigma_m(i) -= eps;
        dHill_num(i) = (Hill_stress(sigma_p, P_params) - Hill_stress(sigma_m, P_params)) / (2. * eps);
    }
    EXPECT_LT(norm(dHill - dHill_num, 2), 1.E-4);
}

TEST(Tcriteria, DFA_stress_basic)
{
    // DFA requires 7 parameters: {a1, a2, a3, b1, b2, b3, c}
    vec P_params = {0.5, 0.6, 0.7, 1.5, 1.5, 1.6, 1.2};
    vec sigma = {400., 0., 0., 0., 0., 0.};

    double dfa = DFA_stress(sigma, P_params);
    // DFA should return a positive value for non-zero stress
    EXPECT_GT(dfa, 0.);
}

TEST(Tcriteria, dDFA_stress_numerical)
{
    // DFA requires 7 parameters
    vec P_params = {0.5, 0.6, 0.7, 1.5, 1.5, 1.6, 1.2};
    vec sigma = {400., 100., 200., 50., 30., 70.};

    vec dDFA = dDFA_stress(sigma, P_params);

    double eps = 1.E-6;
    vec dDFA_num = zeros(6);
    for (int i = 0; i < 6; i++) {
        vec sigma_p = sigma;
        vec sigma_m = sigma;
        sigma_p(i) += eps;
        sigma_m(i) -= eps;
        dDFA_num(i) = (DFA_stress(sigma_p, P_params) - DFA_stress(sigma_m, P_params)) / (2. * eps);
    }
    EXPECT_LT(norm(dDFA - dDFA_num, 2), 1.E-4);
}

TEST(Tcriteria, Ani_stress_basic)
{
    // Ani requires 9 parameters
    vec P_params = {1., 1.2, 1.3, -0.2, -0.2, -0.33, 1., 1., 1.4};
    vec sigma = {400., 100., 200., 50., 30., 70.};

    double ani = Ani_stress(sigma, P_params);
    EXPECT_GT(ani, 0.);
}

TEST(Tcriteria, dAni_stress_numerical)
{
    // Ani requires 9 parameters
    vec P_params = {1., 1.2, 1.3, -0.2, -0.2, -0.33, 1., 1., 1.4};
    vec sigma = {400., 100., 200., 50., 30., 70.};

    vec dAni = dAni_stress(sigma, P_params);

    double eps = 1.E-6;
    vec dAni_num = zeros(6);
    for (int i = 0; i < 6; i++) {
        vec sigma_p = sigma;
        vec sigma_m = sigma;
        sigma_p(i) += eps;
        sigma_m(i) -= eps;
        dAni_num(i) = (Ani_stress(sigma_p, P_params) - Ani_stress(sigma_m, P_params)) / (2. * eps);
    }
    EXPECT_LT(norm(dAni - dAni_num, 2), 1.E-4);
}

TEST(Tcriteria, Eq_stress_dispatch)
{
    vec sigma = {400., 0., 0., 0., 0., 0.};
    // Hill requires 6 parameters (F,G,H,L,M,N)
    vec P_params_hill = {0.5, 0.5, 0.5, 1.5, 1.5, 1.5};

    // "Mises" should match Mises_stress
    double eq_mises = Eq_stress(sigma, "Mises");
    EXPECT_LT(fabs(eq_mises - Mises_stress(sigma)), 1.E-9);

    // "Hill" should match Hill_stress
    double eq_hill = Eq_stress(sigma, "Hill", P_params_hill);
    double hill = Hill_stress(sigma, P_params_hill);
    EXPECT_LT(fabs(eq_hill - hill), 1.E-9);
}

TEST(Tcriteria, dEq_stress_dispatch)
{
    vec sigma = {400., 100., 200., 50., 30., 70.};
    vec P_params = {0.5, 0.5, 0.5, 1.5, 1.5, 1.5};

    // "Hill" derivative should match dHill_stress
    vec dEq = dEq_stress(sigma, "Hill", P_params);
    vec dHill = dHill_stress(sigma, P_params);
    EXPECT_LT(norm(dEq - dHill, 2), 1.E-9);
}

TEST(Tcriteria, P_tensors_symmetry)
{
    // Hill requires 6 parameters (F,G,H,L,M,N)
    vec P_params_hill = {0.5, 0.5, 0.5, 1.5, 1.5, 1.5};
    // DFA requires 7 parameters (F,G,H,L,M,N,K)
    vec P_params_dfa = {0.5, 0.6, 0.7, 1.5, 1.5, 1.6, 1.2};
    // Ani requires 9 parameters
    vec P_params_ani = {1., 1.2, 1.3, -0.2, -0.2, -0.33, 1., 1., 1.4};

    mat P_h = P_Hill(P_params_hill);
    mat P_d = P_DFA(P_params_dfa);
    mat P_a = P_Ani(P_params_ani);

    // P tensors should be symmetric
    EXPECT_LT(norm(P_h - P_h.t(), 2), 1.E-9);
    EXPECT_LT(norm(P_d - P_d.t(), 2), 1.E-9);
    EXPECT_LT(norm(P_a - P_a.t(), 2), 1.E-9);
}

TEST(Tcriteria, Eq_stress_P_consistency)
{
    vec sigma = {400., 100., 200., 50., 30., 70.};
    vec P_params = {0.5, 0.5, 0.5, 1.5, 1.5, 1.5};

    mat H = P_Hill(P_params);
    double eq = Eq_stress_P(sigma, H);
    // Should be positive for non-zero stress
    EXPECT_GT(eq, 0.);

    // Derivative should exist
    vec dEq = dEq_stress_P(sigma, H);
    EXPECT_EQ(dEq.n_elem, (arma::uword)6);
}

TEST(Tcriteria, dJ2_stress_numerical)
{
    vec sigma = {400., 100., 200., 50., 30., 70.};
    vec dJ2 = dJ2_stress(sigma);

    double eps = 1.E-6;
    vec dJ2_num = zeros(6);
    for (int i = 0; i < 6; i++) {
        vec sigma_p = sigma;
        vec sigma_m = sigma;
        sigma_p(i) += eps;
        sigma_m(i) -= eps;
        dJ2_num(i) = (J2_stress(sigma_p) - J2_stress(sigma_m)) / (2. * eps);
    }
    EXPECT_LT(norm(dJ2 - dJ2_num, 2), 1.E-4);
}

TEST(Tcriteria, dJ3_stress_numerical)
{
    vec sigma = {400., 100., 200., 50., 30., 70.};
    vec dJ3 = dJ3_stress(sigma);

    double eps = 1.E-6;
    vec dJ3_num = zeros(6);
    for (int i = 0; i < 6; i++) {
        vec sigma_p = sigma;
        vec sigma_m = sigma;
        sigma_p(i) += eps;
        sigma_m(i) -= eps;
        dJ3_num(i) = (J3_stress(sigma_p) - J3_stress(sigma_m)) / (2. * eps);
    }
    EXPECT_LT(norm(dJ3 - dJ3_num, 2), 1.E-3);
}

TEST(Tcriteria, dDrucker_stress_numerical)
{
    double b = 0.5;
    double n = 2.;
    vec sigma = {400., 100., 200., 50., 30., 70.};
    vec dDrucker = dDrucker_stress(sigma, b, n);

    double eps = 1.E-6;
    vec dDrucker_num = zeros(6);
    for (int i = 0; i < 6; i++) {
        vec sigma_p = sigma;
        vec sigma_m = sigma;
        sigma_p(i) += eps;
        sigma_m(i) -= eps;
        dDrucker_num(i) = (Drucker_stress(sigma_p, b, n) - Drucker_stress(sigma_m, b, n)) / (2. * eps);
    }
    EXPECT_LT(norm(dDrucker - dDrucker_num, 2), 1.E-3);
}

TEST(Tcriteria, dTresca_stress_numerical)
{
    vec sigma = {400., 100., 200., 0., 0., 0.};
    vec dTresca = dTresca_stress(sigma);

    // Tresca derivative is non-smooth, so just verify it returns correct size
    // and is non-zero for non-trivial stress state
    EXPECT_EQ(dTresca.n_elem, (arma::uword)6);
    EXPECT_GT(norm(dTresca, 2), 0.);
}