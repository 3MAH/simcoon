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
#include <simcoon/Continuum_mechanics/Functions/recovery_props.hpp>

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

TEST(Tconstitutive, identity_tensors)
{
    mat Ir = Ireal();
    mat Iv = Ivol();
    mat Id = Idev();

    // Ireal: identity with 0.5 on shear diagonals
    mat Ir_ref = eye(6, 6);
    for (int i = 3; i < 6; i++) Ir_ref(i, i) = 0.5;
    EXPECT_LT(norm(Ir - Ir_ref, 2), sim_iota);

    // Ivol: 1/3 for first 3x3 block
    mat Iv_ref = zeros(6, 6);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            Iv_ref(i, j) = 1. / 3.;
    EXPECT_LT(norm(Iv - Iv_ref, 2), sim_iota);

    // Idev = Ireal - Ivol
    EXPECT_LT(norm(Id - (Ir - Iv), 2), sim_iota);
}

TEST(Tconstitutive, identity_tensors_alt)
{
    mat Ir2 = Ireal2();
    mat Id2 = Idev2();
    mat Iv = Ivol();

    // Ireal2: identity with 2.0 on shear diagonals
    mat Ir2_ref = eye(6, 6);
    for (int i = 3; i < 6; i++) Ir2_ref(i, i) = 2.;
    EXPECT_LT(norm(Ir2 - Ir2_ref, 2), sim_iota);

    // Idev2 = Ireal2 - Ivol
    EXPECT_LT(norm(Id2 - (Ir2 - Iv), 2), sim_iota);
}

TEST(Tconstitutive, Ith_Ir2_Ir05)
{
    vec th = Ith();
    vec r2 = Ir2();
    vec r05 = Ir05();

    // Ith = [1,1,1,0,0,0]
    vec th_ref = zeros(6);
    for (int i = 0; i < 3; i++) th_ref(i) = 1.;
    EXPECT_LT(norm(th - th_ref, 2), sim_iota);

    // Ir2 = [1,1,1,2,2,2]
    vec r2_ref = ones(6);
    for (int i = 3; i < 6; i++) r2_ref(i) = 2.;
    EXPECT_LT(norm(r2 - r2_ref, 2), sim_iota);

    // Ir05 = [1,1,1,0.5,0.5,0.5]
    vec r05_ref = ones(6);
    for (int i = 3; i < 6; i++) r05_ref(i) = 0.5;
    EXPECT_LT(norm(r05 - r05_ref, 2), sim_iota);

    // Ir2 % Ir05 = ones(6) (element-wise product)
    vec product = r2 % r05;
    EXPECT_LT(norm(product - ones(6), 2), sim_iota);
}

TEST(Tconstitutive, H_iso)
{
    double etaB = 1000.;
    double etaS = 500.;

    mat H = H_iso(etaB, etaS);

    // H should be symmetric
    EXPECT_LT(norm(H - H.t(), 2), 1.E-9);

    // H should be 6x6
    EXPECT_EQ(H.n_rows, (arma::uword)6);
    EXPECT_EQ(H.n_cols, (arma::uword)6);
}

TEST(Tconstitutive, el_pred_both_overloads)
{
    double E = 70000.;
    double nu = 0.3;
    mat L = L_iso(E, nu, "Enu");

    vec Eel = {0.001, -0.0003, -0.0003, 0., 0., 0.};
    vec sigma_start = zeros(6);

    // First overload: sigma = sigma_start + L * Eel
    vec sigma1 = el_pred(sigma_start, L, Eel);
    vec sigma_ref = sigma_start + L * Eel;
    EXPECT_LT(norm(sigma1 - sigma_ref, 2), 1.E-9);

    // Second overload: sigma = L * Eel (no sigma_start)
    vec sigma2 = el_pred(L, Eel);
    vec sigma_ref2 = L * Eel;
    EXPECT_LT(norm(sigma2 - sigma_ref2, 2), 1.E-9);

    // Both should give same result when sigma_start = 0
    EXPECT_LT(norm(sigma1 - sigma2, 2), 1.E-9);
}

TEST(Tconstitutive, Isotropize)
{
    // Build an orthotropic stiffness tensor
    double E1 = 10000.;
    double E2 = 20000.;
    double E3 = 30000.;
    double nu12 = 0.3;
    double nu13 = 0.3;
    double nu23 = 0.3;
    double G12 = 6000.;
    double G13 = 7000.;
    double G23 = 8000.;

    mat L = L_ortho(E1, E2, E3, nu12, nu13, nu23, G12, G13, G23, "EnuG");
    mat L_isotropized = Isotropize(L);

    // The isotropized tensor should be isotropic: check_symetries should report ELISO
    // Verify by extracting iso props
    vec props = L_iso_props(L_isotropized);
    // If L_iso_props succeeds and returns E > 0, nu in (-1, 0.5), it's isotropic
    EXPECT_GT(props(0), 0.);
    EXPECT_GT(props(1), -1.);
    EXPECT_LT(props(1), 0.5);
}
