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

#include <simcoon/FTensor.hpp>
#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Functions/transfer.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Continuum_mechanics/Functions/kinematics.hpp>
#include <simcoon/Continuum_mechanics/Functions/objective_rates.hpp>
#include <simcoon/Continuum_mechanics/Functions/stress.hpp>

using namespace std;
using namespace arma;
using namespace simcoon;
using namespace FTensor;

TEST(Tobjective_rates, get_B)
{

    mat F = zeros(3,3);
    
	F(0,0) = 4.;
	F(0,1) = 4.;
	F(0,2) = 1.5;
	F(1,0) = 1.75;
	F(1,1) = 2.;
	F(1,2) = 1.25;
	F(2,0) = 1.5;
	F(2,1) = 3.5;
	F(2,2) = 0.25;

    mat B = L_Cauchy_Green(F);
    vec bi = zeros(3);
    mat Bi;
    eig_sym(bi, Bi, B);

    mat Bij = zeros(3,3);
    Tensor2<double,3,3> Bij_ = mat_FTensor2(zeros(3,3));
    Tensor4<double,3,3,3,3> BBBB_ = mat_FTensor4(zeros(6,6));
    
    Index<'k', 3> k;
    Index<'l', 3> l;
    Index<'m', 3> m;
    Index<'n', 3> n;
    
    double f_z = 0.;
    for (unsigned int i=0; i<3; i++) {
        for (unsigned int j=0; j<3; j++) {
            if ((i!=j)&&(fabs(bi(i)-bi(j))>sim_iota)) {
                Bij = Bi.col(i)*(Bi.col(j)).t();
                Bij_ = mat_FTensor2(Bij);
                f_z = (1.+(bi(i)/bi(j)))/(1.-(bi(i)/bi(j)))+2./log(bi(i)/bi(j));
                //BBBB_(k,l,m,n) = b_[i](k)*b_[j](l)*b_[i](m)*b_[j](n);
                BBBB_(k,l,m,n) = BBBB_(k,l,m,n) + f_z*Bij_(k,l)*Bij_(m,n);
            }
        }
    }
    
    mat BBBB_test = FTensor4_mat(BBBB_);

    mat BBBB = get_BBBB(F);

    cout << "BBBB\n" << BBBB << endl;
    cout << "BBBB_test\n" << BBBB_test << endl;
    cout << "BBBB - BBBB_test \n" << BBBB - BBBB_test << endl;
    cout << "norm(BBBB - BBBB_test,2) = " << norm(BBBB - BBBB_test,2) << endl;
    
    //Test of BBBB function
    EXPECT_LT(norm(BBBB - BBBB_test,2),1.E-9);
}

TEST(Tobjective_rates, logarithmic_functions)
{
    mat F1 = zeros(3,3);
	F1(0,0) = 4.;
	F1(0,1) = 4.;
	F1(0,2) = 1.5;
	F1(1,0) = 1.75;
	F1(1,1) = 2.;
	F1(1,2) = 1.25;
	F1(2,0) = 1.5;
	F1(2,1) = 3.5;
	F1(2,2) = 0.25;    

    mat F0 = F1*(0.9);
	F0(0,1) = 3.;    
	F0(1,0) = 2.;

    double DTime = 1.E-3;

    //Compute log rate and increment of rotation
    mat I = eye(3,3);
    mat L = (1./DTime)*(F1-F0)*inv(F1);
    
    //decomposition of L
    mat D_test = 0.5*(L+L.t());
    mat W_test = 0.5*(L-L.t());

    mat B_test = L_Cauchy_Green(F1);
    
    vec bi = zeros(3);
    mat Bi;
    eig_sym(bi, Bi, B_test);

    mat Bij = zeros(3,3);
    Tensor2<double,3,3> Bij_ = mat_FTensor2(zeros(3,3));
    Tensor2<double,3,3> N_ = mat_FTensor2(zeros(3,3));
    Tensor2<double,3,3> D_= mat_FTensor2(zeros(3,3));
    Tensor4<double,3,3,3,3> BBBB_ = mat_FTensor4(zeros(6,6));
    
    D_ = mat_FTensor2(D_test);
    
    Index<'k', 3> k;
    Index<'l', 3> l;
    Index<'m', 3> m;
    Index<'n', 3> n;
    
    double f_z = 0.;
    BBBB_(k,l,m,n) = Bij_(k,l)*Bij_(m,n);
    for (unsigned int i=0; i<3; i++) {
        for (unsigned int j=0; j<3; j++) {
            if ((i!=j)&&(fabs(bi(i)-bi(j))>sim_iota)) {
                Bij = Bi.col(i)*(Bi.col(j)).t();
                Bij_ = mat_FTensor2(Bij);
                f_z = (1.+(bi(i)/bi(j)))/(1.-(bi(i)/bi(j)))+2./log(bi(i)/bi(j));
                //BBBB_(k,l,m,n) = b_[i](k)*b_[j](l)*b_[i](m)*b_[j](n);
                BBBB_(k,l,m,n) = BBBB_(k,l,m,n) + f_z*Bij_(k,l)*Bij_(m,n);
            }
        }
    }
    N_(k,l) = N_(k,l) + BBBB_(k,l,m,n)*D_(m,n);

    mat Omega_test = W_test + FTensor2_mat(N_);

    mat BBBB = get_BBBB(F1);
    mat Omega_test2 = W_test + v2t_skewsym(BBBB*t2v_strain(D_test));

    mat DR = zeros(3,3);
    mat D = zeros(3,3);
    mat Omega = zeros(3,3);
    
    logarithmic(DR, D, Omega, DTime, F0, F1);
    
    cout << "Omega\n" << BBBB << endl;
    cout << "Omega_test\n" << Omega_test << endl;
    cout << "Omega_test2\n" << Omega_test2 << endl;

    //Test of logarithmic function
    EXPECT_LT(norm(D - D_test,2),1.E-9);
    EXPECT_LT(norm(Omega - Omega_test,2),1.E-9);
    EXPECT_LT(norm(Omega - Omega_test2,2),1.E-9);
}

TEST(Tobjective_rates, Jaumann_rate)
{
    mat F0 = eye(3,3);
    mat F1 = eye(3,3);
    F1(0,0) = 1.01;
    F1(0,1) = 0.005;
    F1(1,0) = -0.005;
    F1(1,1) = 1.02;
    F1(2,2) = 0.98;
    double DTime = 0.001;

    mat DR = zeros(3,3);
    mat D = zeros(3,3);
    mat W = zeros(3,3);
    Jaumann(DR, D, W, DTime, F0, F1);

    // D should be symmetric
    EXPECT_LT(norm(D - D.t(), 2), 1.E-9);
    // W should be antisymmetric
    EXPECT_LT(norm(W + W.t(), 2), 1.E-9);
    // DR should be close to orthogonal for small increments
    EXPECT_LT(norm(DR.t() * DR - eye(3,3), 2), 1.E-3);
}

TEST(Tobjective_rates, Green_Naghdi_rate)
{
    mat F0 = eye(3,3);
    mat F1 = eye(3,3);
    F1(0,0) = 1.01;
    F1(0,1) = 0.005;
    F1(1,0) = -0.005;
    F1(1,1) = 1.02;
    F1(2,2) = 0.98;
    double DTime = 0.001;

    mat DR = zeros(3,3);
    mat D = zeros(3,3);
    mat Omega = zeros(3,3);
    Green_Naghdi(DR, D, Omega, DTime, F0, F1);

    // D should be symmetric
    EXPECT_LT(norm(D - D.t(), 2), 1.E-6);
    // Omega should be mostly antisymmetric (may have small numerical error)
    EXPECT_LT(norm(Omega + Omega.t(), 2), 0.1);
    // DR should be close to orthogonal
    EXPECT_LT(norm(DR.t() * DR - eye(3,3), 2), 1.E-3);
}

TEST(Tobjective_rates, logarithmic_R_rate)
{
    mat F0 = eye(3,3);
    mat F1 = eye(3,3);
    F1(0,0) = 1.01;
    F1(0,1) = 0.005;
    F1(1,0) = -0.005;
    F1(1,1) = 1.02;
    F1(2,2) = 0.98;
    double DTime = 0.001;

    mat DR = zeros(3,3);
    mat N_1 = zeros(3,3);
    mat N_2 = zeros(3,3);
    mat D = zeros(3,3);
    mat Omega = zeros(3,3);
    logarithmic_R(DR, N_1, N_2, D, Omega, DTime, F0, F1);

    // D should be symmetric
    EXPECT_LT(norm(D - D.t(), 2), 1.E-9);
    // DR should be close to orthogonal
    EXPECT_LT(norm(DR.t() * DR - eye(3,3), 2), 1.E-3);
}

TEST(Tobjective_rates, logarithmic_F_rate)
{
    mat F0 = eye(3,3);
    mat F1 = eye(3,3);
    F1(0,0) = 1.01;
    F1(0,1) = 0.005;
    F1(1,0) = -0.005;
    F1(1,1) = 1.02;
    F1(2,2) = 0.98;
    double DTime = 0.001;

    mat DF = zeros(3,3);
    mat N_1 = zeros(3,3);
    mat N_2 = zeros(3,3);
    mat D = zeros(3,3);
    mat L = zeros(3,3);
    logarithmic_F(DF, N_1, N_2, D, L, DTime, F0, F1);

    // D should be symmetric
    EXPECT_LT(norm(D - D.t(), 2), 1.E-9);
    // L should decompose into D + W
    mat W = 0.5 * (L - L.t());
    EXPECT_LT(norm(L - (D + W), 2), 1.E-9);
}

TEST(Tobjective_rates, Truesdell_rate)
{
    mat F0 = eye(3,3);
    mat F1 = eye(3,3);
    F1(0,0) = 1.01;
    F1(0,1) = 0.005;
    F1(1,0) = -0.005;
    F1(1,1) = 1.02;
    F1(2,2) = 0.98;
    double DTime = 0.001;

    mat DF = zeros(3,3);
    mat D = zeros(3,3);
    mat L = zeros(3,3);
    Truesdell(DF, D, L, DTime, F0, F1);

    // D should be symmetric part of L
    EXPECT_LT(norm(D - 0.5*(L + L.t()), 2), 1.E-9);
}

TEST(Tobjective_rates, Delta_log_strain_identity)
{
    mat D = zeros(3,3);
    mat Omega = zeros(3,3);
    double DTime = 0.001;

    // Zero D and Omega -> zero strain increment
    mat Deps = Delta_log_strain(D, Omega, DTime);
    EXPECT_LT(norm(Deps, 2), 1.E-9);
}

TEST(Tobjective_rates, get_BBBB_GN)
{
    mat F = zeros(3,3);
    F(0,0) = 4.;
    F(0,1) = 4.;
    F(0,2) = 1.5;
    F(1,0) = 1.75;
    F(1,1) = 2.;
    F(1,2) = 1.25;
    F(2,0) = 1.5;
    F(2,1) = 3.5;
    F(2,2) = 0.25;

    mat BBBB_GN = get_BBBB_GN(F);
    // Should be 6x6
    EXPECT_EQ(BBBB_GN.n_rows, (arma::uword)6);
    EXPECT_EQ(BBBB_GN.n_cols, (arma::uword)6);
}

TEST(Tobjective_rates, tangent_conversions_DtauDe_DSDE)
{
    // Use an isotropic material for testing
    double E = 70000.;
    double nu = 0.3;
    mat Lt = L_iso(E, nu, "Enu");

    mat F = eye(3,3);
    F(0,0) = 1.05;
    F(1,1) = 0.98;
    F(2,2) = 1.02;
    double J = det(F);

    mat B = get_BBBB(F);
    mat sigma = zeros(3,3);
    sigma(0,0) = 100.;
    sigma(1,1) = 50.;
    sigma(2,2) = 75.;
    mat tau = Cauchy2Kirchoff(sigma, F, J);

    // DtauDe -> DSDE should be 6x6
    mat DSDE = DtauDe_2_DSDE(Lt, B, F, tau);
    EXPECT_EQ(DSDE.n_rows, (arma::uword)6);
    EXPECT_EQ(DSDE.n_cols, (arma::uword)6);
}

TEST(Tobjective_rates, tangent_conversions_DsigmaDe_DSDE)
{
    double E = 70000.;
    double nu = 0.3;
    mat Lt = L_iso(E, nu, "Enu");

    mat F = eye(3,3);
    F(0,0) = 1.05;
    F(1,1) = 0.98;
    F(2,2) = 1.02;

    mat sigma = zeros(3,3);
    sigma(0,0) = 100.;
    sigma(1,1) = 50.;
    sigma(2,2) = 75.;

    mat DSDE = DsigmaDe_2_DSDE(Lt, F, sigma);
    EXPECT_EQ(DSDE.n_rows, (arma::uword)6);
    EXPECT_EQ(DSDE.n_cols, (arma::uword)6);
}

TEST(Tobjective_rates, DSDE_2_DtauDe_roundtrip)
{
    double E = 70000.;
    double nu = 0.3;
    mat Lt = L_iso(E, nu, "Enu");
    double J = 1.05;

    // DtauDe_2_DsigmaDe and DsigmaDe_2_DtauDe should be inverses
    mat Dsigma = DtauDe_2_DsigmaDe(Lt, J);
    mat Dtau = DsigmaDe_2_DtauDe(Dsigma, J);
    EXPECT_LT(norm(Dtau - Lt, 2), 1.E-9);
}

TEST(Tobjective_rates, Lie_Jaumann_conversion)
{
    double E = 70000.;
    double nu = 0.3;
    mat Lt = L_iso(E, nu, "Enu");

    mat tau = zeros(3,3);
    tau(0,0) = 100.;
    tau(1,1) = 50.;
    tau(2,2) = 75.;

    mat Jaumann_Lt = Dtau_LieDD_Dtau_JaumannDD(Lt, tau);
    EXPECT_EQ(Jaumann_Lt.n_rows, (arma::uword)6);
    EXPECT_EQ(Jaumann_Lt.n_cols, (arma::uword)6);
}

TEST(Tobjective_rates, DSDE_conversions_LieDD)
{
    double E = 70000.;
    double nu = 0.3;
    mat DSDE = L_iso(E, nu, "Enu");

    mat F = eye(3,3);
    F(0,0) = 1.05;
    F(1,1) = 0.98;
    F(2,2) = 1.02;

    mat Dtau_Lie = DSDE_2_Dtau_LieDD(DSDE, F);
    EXPECT_EQ(Dtau_Lie.n_rows, (arma::uword)6);
    EXPECT_EQ(Dtau_Lie.n_cols, (arma::uword)6);

    mat Dsigma_Lie = DSDE_2_Dsigma_LieDD(DSDE, F);
    EXPECT_EQ(Dsigma_Lie.n_rows, (arma::uword)6);
    EXPECT_EQ(Dsigma_Lie.n_cols, (arma::uword)6);

    // Dsigma_Lie should be (1/J) * Dtau_Lie
    double J = det(F);
    EXPECT_LT(norm(Dsigma_Lie - (1./J)*Dtau_Lie, 2), 1.E-6);
}
