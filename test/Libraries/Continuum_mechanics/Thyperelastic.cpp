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

///@file Thyperelastic.cpp
///@brief Test for hyperelastic fonctions
///@version 1.0

#include <gtest/gtest.h>
#include <armadillo>

#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Functions/transfer.hpp>
#include <simcoon/Continuum_mechanics/Functions/hyperelastic.hpp>

using namespace std;
using namespace arma;
using namespace simcoon;

TEST(Thyperelastic, isochoric_invariants)
{
    mat b_test = eye(3,3);

    vec I_test = simcoon::isochoric_invariants(b_test);

    EXPECT_LT(I_test(0)-3.,1.E-9);
    EXPECT_LT(I_test(1)-3.,1.E-9);    
    EXPECT_LT(I_test(2)-1.,1.E-9);

    vec lambda = eig_sym(b_test);
    lambda.transform( [](double val) { return (sqrt(val)); } );

    vec I_test_lambda = simcoon::isochoric_invariants(lambda);

    EXPECT_LT(I_test_lambda(0)-3.,1.E-9);
    EXPECT_LT(I_test_lambda(1)-3.,1.E-9);    
    EXPECT_LT(I_test_lambda(2)-1.,1.E-9);

    mat b_rand = simcoon::v2t_strain(randu(6));

    double J = sqrt(det(b_rand));
    
    mat b_bar = pow(J,-2./3.)*b_rand;
    vec I = zeros(3);    
    I(0) = trace(b_bar);
    I(1) = 0.5*(pow(trace(b_bar),2.)-trace(powmat(b_bar,2)));
    I(2) = 1.;

    vec I_rand = simcoon::isochoric_invariants(b_rand);    

    cout << "b_rand = " << b_rand << endl;
    cout << "b_bar = " << b_bar << endl;
    cout << "J = " << J << endl;    

    cout << "I_rand = " << I_rand << endl;
    cout << "I = " << I << endl;    

    EXPECT_LT(I_rand(0)-I(0),1.E-9);
    EXPECT_LT(I_rand(1)-I(1),1.E-9);    
    EXPECT_LT(I_rand(2)-I(2),1.E-9);

    lambda = eig_sym(b_rand);
    lambda.transform( [](double val) { return (sqrt(val)); } );

    vec I_rand_lambda = simcoon::isochoric_invariants(lambda);

    EXPECT_LT(I_rand_lambda(0)-I(0),1.E-9);
    EXPECT_LT(I_rand_lambda(1)-I(1),1.E-9);    
    EXPECT_LT(I_rand_lambda(2)-I(2),1.E-9);

}

TEST(Thyperelastic, isochoric_pstretch)
{

    mat V_test = eye(3,3);
    mat b_test = powmat(V_test,2);    

    vec lambda_bar_test_from_V = simcoon::isochoric_pstretch_from_V(V_test);

    EXPECT_LT(lambda_bar_test_from_V(0)-3.,1.E-9);
    EXPECT_LT(lambda_bar_test_from_V(1)-3.,1.E-9);    
    EXPECT_LT(lambda_bar_test_from_V(2)-1.,1.E-9);

    vec lambda_bar_test_from_b = simcoon::isochoric_pstretch_from_b(b_test);

    EXPECT_LT(lambda_bar_test_from_b(0)-1.,1.E-9);
    EXPECT_LT(lambda_bar_test_from_b(1)-1.,1.E-9);    
    EXPECT_LT(lambda_bar_test_from_b(2)-1.,1.E-9);

    mat V_rand = simcoon::v2t_strain(randu(6)) + eye(3,3);
    mat b_rand = powmat(V_rand,2);    
    double J = det(V_rand);    

    vec lambda = eig_sym(V_rand);
    vec lambda_bar = pow(J,-1./3.)*lambda;

    lambda_bar_test_from_V = simcoon::isochoric_pstretch_from_V(V_rand);    

    cout << lambda_bar_test_from_V << endl;

    EXPECT_LT(lambda_bar_test_from_V(0)-lambda_bar(0),1.E-9);
    EXPECT_LT(lambda_bar_test_from_V(1)-lambda_bar(1),1.E-9);    
    EXPECT_LT(lambda_bar_test_from_V(2)-lambda_bar(2),1.E-9);

    lambda_bar_test_from_b = simcoon::isochoric_pstretch_from_b(b_rand);    

    EXPECT_LT(lambda_bar_test_from_b(0)-lambda_bar(0),1.E-9);
    EXPECT_LT(lambda_bar_test_from_b(1)-lambda_bar(1),1.E-9);    
    EXPECT_LT(lambda_bar_test_from_b(2)-lambda_bar(2),1.E-9);

}