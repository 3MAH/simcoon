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

///@file Tconstitutive.cpp
///@brief Test for Constitutive tensors in Voigt notation
///@version 1.0

#include <gtest/gtest.h>
#include <armadillo>

#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>

using namespace std;
using namespace arma;
using namespace simcoon;

TEST(Tconstitutive, L_iso_M_iso)
{
    double E = 70000.;
    double nu = 0.3;
    
    double mu = E/(2.*(1+nu));
    double lambda = E*nu/((1.+nu)*(1.-2.*nu));
    
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
    
    mat Mt = zeros(6,6);
    Mt(0,0) = 1./E;
    Mt(0,1) = -nu/E;
    Mt(0,2) = -nu/E;
    Mt(1,0) = -nu/E;
    Mt(1,1) = 1./E;
    Mt(1,2) = -nu/E;
    Mt(2,0) = -nu/E;
    Mt(2,1) = -nu/E;
    Mt(2,2) = 1./E;
    Mt(3,3) = 2.*(1.+nu)/E;
    Mt(4,4) = 2.*(1.+nu)/E;
    Mt(5,5) = 2.*(1.+nu)/E;
    
    //Test of L_iso function
    mat Ltest = L_iso(E, nu, "Enu");
    EXPECT_LT(norm(Ltest - Lt,2),1.E-9);
	Ltest = L_iso(lambda, mu, "lambdamu");
	EXPECT_LT(norm(Ltest - Lt,2),1.E-9);
    
	//Test of M_iso function
    mat Mtest = M_iso(E, nu, "Enu");
	EXPECT_LT(norm(Mtest - Mt,2),1.E-9);
	Mtest = M_iso(mu, lambda, "mulambda");
	EXPECT_LT(norm(Mtest - Mt,2),1.E-9);
}

TEST(Tconstitutive, L_cubic_M_cubic)
{
    double C11 = 1000;
    double C12 = 400;
    double C44 = 500;
    
    double nu = 1 / (1 + C11/C12);
    double E = C11*( 1 - 3*pow(nu,2) - 2*pow(nu,3) )/( 1-pow(nu,2));
    double G = C44;
    
    mat Lcub = zeros(6,6);
    Lcub(0,0) = C11;
    Lcub(0,1) = C12;
    Lcub(0,2) = C12;
    Lcub(1,0) = C12;
    Lcub(1,1) = C11;
    Lcub(1,2) = C12;
    Lcub(2,0) = C12;
    Lcub(2,1) = C12;
    Lcub(2,2) = C11;
    Lcub(3,3) = C44;
    Lcub(4,4) = C44;
    Lcub(5,5) = C44;
    
    //Test of L_cubic function
    mat Ltest = L_cubic(C11,C12,C44, "Cii");
    EXPECT_LT(norm(Ltest - Lcub,2),1.E-9);
    
    Ltest = L_cubic(E,nu,G, "EnuG");
    EXPECT_LT(norm(Ltest - Lcub,2),1.E-9);
    
    //Test of M_cubic function
    mat Mtest = M_cubic(C11,C12,C44, "Cii");
    EXPECT_LT(norm(Mtest - inv(Lcub),2),1.E-9);

    Mtest = M_cubic(E,nu,G, "EnuG");
    EXPECT_LT(norm(Mtest - inv(Lcub),2),1.E-9);
}

TEST(Tconstitutive, L_isotrans_M_isotrans)
{
    double EL = 10000.;
    double ET = 20000.;
    double nuTL = 0.3;
    double nuTT = 0.3;
    double GLT = 8000.;
    
    mat Mtrans = zeros(6,6);
    Mtrans(0,0) = 1./EL;
    Mtrans(0,1) = -nuTL/ET;
    Mtrans(0,2) = -nuTL/ET;
    Mtrans(1,0) = -nuTL/ET;
    Mtrans(1,1) = 1./ET;
    Mtrans(1,2) = -nuTT/ET;
    Mtrans(2,0) = -nuTL/ET;
    Mtrans(2,1) = -nuTT/ET;
    Mtrans(2,2) = 1./ET;
    Mtrans(3,3) = 1/GLT;
    Mtrans(4,4) = 1/GLT;
    Mtrans(5,5) = (2.*(1.+nuTT))/ET;
    
	//Test of L_isotrans function axis 1
    mat Ltest = L_isotrans(EL, ET, nuTL, nuTT, GLT, 1);
	EXPECT_LT(norm(Ltest - inv(Mtrans),2),1.E-9);
    //Test of M_isotrans function axis 1
    mat Mtest = M_isotrans(EL, ET, nuTL, nuTT, GLT, 1);
    EXPECT_LT(norm(Mtest - Mtrans,2),1.E-9);
    
}

TEST(Tconstitutive, L_ortho_M_ortho)
{
    double E1 = 10000.;
    double E2 = 20000.;
    double E3 = 30000.;
    double nu12 = 0.3;
    double nu13 = 0.3;
    double nu23 = 0.3;
    double G12 = 6000.;
    double G13 = 7000.;
    double G23 = 8000.;
    
    mat Mortho = zeros(6,6);
    Mortho(0,0) = 1/E1;
    Mortho(0,1) = -nu12/E1;
    Mortho(0,2) = -nu13/E1;
    Mortho(1,0) = -nu12/E1;
    Mortho(1,1) = 1/E2;
    Mortho(1,2) = -nu23/E2;
    Mortho(2,0) = -nu13/E1;
    Mortho(2,1) = -nu23/E2;
    Mortho(2,2) = 1/E3;
    Mortho(3,3) = 1/G12;
    Mortho(4,4) = 1/G13;
    Mortho(5,5) = 1/G23;
    
    //Test of L_ortho function
    mat Ltest = L_ortho(E1, E2, E3, nu12, nu13, nu23, G12, G13, G23, "EnuG");
    EXPECT_LT(norm(Ltest - inv(Mortho),2),1.E-9);
    //Test of M_ortho function
    mat Mtest = M_ortho(E1, E2, E3, nu12, nu13, nu23, G12, G13, G23, "EnuG");
    EXPECT_LT(norm(Mtest - Mortho,2),1.E-9);
    
    //Test of L_ortho function
    mat Ltest2 = L_ortho(Ltest(0,0), Ltest(0,1), Ltest(0,2), Ltest(1,1), Ltest(1,2), Ltest(2,2), Ltest(3,3), Ltest(4,4), Ltest(5,5), "Cii");
    EXPECT_LT(norm(Ltest2 - inv(Mortho),2),1.E-9);
    //Test of M_ortho function
    mat Mtest2 = M_ortho(Ltest(0,0), Ltest(0,1), Ltest(0,2), Ltest(1,1), Ltest(1,2), Ltest(2,2), Ltest(3,3), Ltest(4,4), Ltest(5,5), "Cii");
    EXPECT_LT(norm(Mtest2 - Mortho,2),1.E-9);

}
