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

///@file Tstress.cpp
///@brief Test for stress fonctions
///@version 1.0

#include <gtest/gtest.h>
#include <armadillo>

#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Functions/transfer.hpp>
#include <simcoon/Continuum_mechanics/Functions/kinematics.hpp>
#include <simcoon/Continuum_mechanics/Functions/stress.hpp>

using namespace std;
using namespace arma;
using namespace simcoon;

TEST(Tstress, round_trip_PKI)
{
    
    vec sigma_rand_vec = randu(6);
    mat F_rand = simcoon::v2t_strain(randu(6))+eye(3,3);

    mat PKI = Cauchy2PKI(v2t_stress(sigma_rand_vec), F_rand);
    mat sigma = PKI2Cauchy(PKI, F_rand);

    EXPECT_LT(norm(t2v_stress(sigma),2) - norm(sigma_rand_vec,2),1.E-9);
}

TEST(Tstress, round_trip_PKII)
{
    vec sigma_rand_vec = randu(6);
    mat F_rand = simcoon::v2t_strain(randu(6))+eye(3,3);

    mat PKII = Cauchy2PKII(v2t_stress(sigma_rand_vec), F_rand);
    mat sigma = PKII2Cauchy(PKII, F_rand);

    EXPECT_LT(norm(t2v_stress(sigma),2) - norm(sigma_rand_vec,2),1.E-9);    
}

TEST(Tstress, test_Biot)
{
    vec sigma_rand_vec = randu(6);
    mat F_rand = simcoon::v2t_strain(randu(6))+eye(3,3);

    mat R = zeros(3,3);
    mat U = zeros(3,3);
    RU_decomposition(R,U,F_rand);

    mat Biot = Cauchy2Biot(v2t_stress(sigma_rand_vec), F_rand);
    mat PKI = Cauchy2PKI(v2t_stress(sigma_rand_vec), F_rand);    
    mat Biot_test = 0.5*(R.t()*PKI + PKI.t()*R);

    EXPECT_LT(norm(Biot,2) - norm(Biot_test,2),1.E-9);
}

TEST(Tstress, Cauchy2Kirchoff_mat)
{
    mat sigma = {{100., 20., 10.}, {20., 200., 30.}, {10., 30., 300.}};
    mat F = eye(3, 3);
    F(0, 0) = 1.2;
    F(1, 1) = 0.9;
    F(2, 2) = 1.1;
    double J = det(F);

    // tau = J * sigma
    mat tau = Cauchy2Kirchoff(sigma, F);
    EXPECT_LT(norm(tau - J * sigma, 2), 1.E-9);
}

TEST(Tstress, Cauchy2Kirchoff_vec)
{
    vec sigma_v = {100., 200., 300., 20., 10., 30.};
    mat F = eye(3, 3);
    F(0, 0) = 1.2;
    F(1, 1) = 0.9;
    F(2, 2) = 1.1;
    double J = det(F);

    vec tau_v = Cauchy2Kirchoff(sigma_v, F);
    EXPECT_LT(norm(tau_v - J * sigma_v, 2), 1.E-9);
}

TEST(Tstress, Kirchoff2Cauchy_mat)
{
    mat sigma = {{100., 20., 10.}, {20., 200., 30.}, {10., 30., 300.}};
    mat F = eye(3, 3);
    F(0, 0) = 1.2;
    F(1, 1) = 0.9;
    F(2, 2) = 1.1;
    double J = det(F);
    mat tau = J * sigma;

    mat sigma_back = Kirchoff2Cauchy(tau, F);
    EXPECT_LT(norm(sigma_back - sigma, 2), 1.E-9);
}

TEST(Tstress, Kirchoff2Cauchy_vec)
{
    vec sigma_v = {100., 200., 300., 20., 10., 30.};
    mat F = eye(3, 3);
    F(0, 0) = 1.2;
    F(1, 1) = 0.9;
    F(2, 2) = 1.1;
    double J = det(F);
    vec tau_v = J * sigma_v;

    vec sigma_back = Kirchoff2Cauchy(tau_v, F);
    EXPECT_LT(norm(sigma_back - sigma_v, 2), 1.E-9);
}

TEST(Tstress, Cauchy_Kirchoff_roundtrip)
{
    vec sigma_v = randu(6);
    mat F = simcoon::v2t_strain(randu(6)) + eye(3, 3);
    double J = det(F);

    // mat version round-trip
    mat sigma_m = v2t_stress(sigma_v);
    mat tau_m = Cauchy2Kirchoff(sigma_m, F);
    mat sigma_m_back = Kirchoff2Cauchy(tau_m, F);
    EXPECT_LT(norm(sigma_m_back - sigma_m, 2), 1.E-9);

    // vec version round-trip
    vec tau_v = Cauchy2Kirchoff(sigma_v, F);
    vec sigma_v_back = Kirchoff2Cauchy(tau_v, F);
    EXPECT_LT(norm(sigma_v_back - sigma_v, 2), 1.E-9);
}

TEST(Tstress, Kirchoff2PKI_round_trip)
{
    vec sigma_v = randu(6);
    mat F = simcoon::v2t_strain(randu(6)) + eye(3, 3);
    double J = det(F);
    mat sigma_m = v2t_stress(sigma_v);
    mat tau = J * sigma_m;

    mat PKI = Kirchoff2PKI(tau, F);
    mat tau_back = PKI2Kirchoff(PKI, F);
    EXPECT_LT(norm(tau_back - tau, 2), 1.E-9);
}

TEST(Tstress, Kirchoff2PKII_mat_round_trip)
{
    vec sigma_v = randu(6);
    mat F = simcoon::v2t_strain(randu(6)) + eye(3, 3);
    double J = det(F);
    mat sigma_m = v2t_stress(sigma_v);
    mat tau = J * sigma_m;

    mat PKII = Kirchoff2PKII(tau, F);
    mat tau_back = PKII2Kirchoff(PKII, F);
    EXPECT_LT(norm(tau_back - tau, 2), 1.E-9);
}

TEST(Tstress, Kirchoff2PKII_vec)
{
    vec sigma_v = randu(6);
    mat F = simcoon::v2t_strain(randu(6)) + eye(3, 3);
    double J = det(F);
    vec tau_v = J * sigma_v;

    vec PKII_v = Kirchoff2PKII(tau_v, F);
    // Verify consistency with mat version
    mat tau_m = v2t_stress(tau_v);
    mat PKII_m = Kirchoff2PKII(tau_m, F);
    EXPECT_LT(norm(PKII_v - t2v_stress(PKII_m), 2), 1.E-9);
}

TEST(Tstress, PKII2Kirchoff)
{
    vec sigma_v = randu(6);
    mat F = simcoon::v2t_strain(randu(6)) + eye(3, 3);
    double J = det(F);
    mat sigma_m = v2t_stress(sigma_v);
    mat tau = J * sigma_m;

    // Forward and back: tau -> PKII -> tau
    mat S = Kirchoff2PKII(tau, F);
    mat tau_back = PKII2Kirchoff(S, F);
    EXPECT_LT(norm(tau_back - tau, 2), 1.E-9);
}