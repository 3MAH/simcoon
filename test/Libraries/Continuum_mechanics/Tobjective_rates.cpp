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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "objective_rates"
#include <boost/test/unit_test.hpp>

#include <armadillo>
#include <FTensor.hpp>
#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Functions/transfer.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Continuum_mechanics/Functions/kinematics.hpp>
#include <simcoon/Continuum_mechanics/Functions/objective_rates.hpp>

using namespace std;
using namespace arma;
using namespace simcoon;
using namespace FTensor;

BOOST_AUTO_TEST_CASE( get_B )
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
	F(2,2) = 6.;

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
    BOOST_CHECK( norm(BBBB - BBBB_test,2) < 1.E-9 );
}

BOOST_AUTO_TEST_CASE( logarithmic_test )
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
	F1(2,2) = 6.;

    mat F0 = F1*(0.9);
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
    
    //Test of logarithmic function
    BOOST_CHECK( norm(D - D_test,2) < 1.E-9 );
    BOOST_CHECK( norm(Omega - Omega_test,2) < 1.E-9 );
    BOOST_CHECK( norm(Omega - Omega_test2,2) < 1.E-9 );
}
