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