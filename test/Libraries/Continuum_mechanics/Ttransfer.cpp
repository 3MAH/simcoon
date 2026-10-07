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

///@file Transfer.cpp
///@brief Test for transfer between tensors, vectors in armadillo (vec, mat) and Fastor
///@version 1.0

#include <gtest/gtest.h>
#include <armadillo>

#include <Fastor/Fastor.h>
#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Functions/transfer.hpp>
#include <simcoon/Continuum_mechanics/Functions/fastor_bridge.hpp>

using namespace std;
using namespace arma;
using namespace simcoon;

TEST(Ttransfer, v2t_func)
{
	
    vec test = zeros(6);
    test(0) = 4.;
    test(1) = 2.;
    test(2) = 6.;
    test(3) = 8.;
    test(4) = 3.;
    test(5) = 7.;
    
	mat strain = zeros(3,3);
	mat stress = zeros(3,3);
	mat teststrain = zeros(3,3);
	mat teststress = zeros(3,3);
	
	teststrain(0,0) = 4.;
	teststrain(0,1) = 4.;
	teststrain(0,2) = 1.5;
	teststrain(1,0) = 4.;
	teststrain(1,1) = 2.;
	teststrain(1,2) = 3.5;
	teststrain(2,0) = 1.5;
	teststrain(2,1) = 3.5;
	teststrain(2,2) = 6.;							

	teststress(0,0) = 4.;
	teststress(0,1) = 8.;
	teststress(0,2) = 3.;
	teststress(1,0) = 8.;
	teststress(1,1) = 2.;
	teststress(1,2) = 7.;
	teststress(2,0) = 3.;
	teststress(2,1) = 7.;
	teststress(2,2) = 6.;
	
	//Test of v2t_strain function
	strain = v2t_strain(zeros(6));
    EXPECT_LT(norm(strain,2),simcoon::iota);
	strain = v2t_strain(test);
    EXPECT_LT(norm(strain - teststrain,2),simcoon::iota);

    //Test of t2v_strain function
	vec dev1 = t2v_strain(zeros(3,3));
    vec dev2 = t2v_strain(teststrain);
    EXPECT_LT(norm(dev1,2),simcoon::iota);
    EXPECT_LT(norm(dev2 - test,2),simcoon::iota);

	//Test of v2t_stress function
	stress = v2t_stress(zeros(6));
    EXPECT_LT(norm(stress,2),simcoon::iota);
	stress = v2t_stress(test);
    EXPECT_LT(norm(stress - teststress,2),simcoon::iota);
    
	//Test of t2v_stress function
	dev1 = t2v_stress(zeros(3,3));
	dev2 = t2v_stress(teststress);
    EXPECT_LT(norm(dev1,2),simcoon::iota);
    EXPECT_LT(norm(dev2 - test,2),simcoon::iota);
    
}

TEST(Ttransfer, Fastor_transfer)
{
    // Test arma <-> Fastor rank-2 round-trip
    mat::fixed<3,3> testmat;
    testmat(0,0) = 4.; testmat(0,1) = 4.; testmat(0,2) = 1.5;
    testmat(1,0) = 4.; testmat(1,1) = 2.; testmat(1,2) = 3.5;
    testmat(2,0) = 1.5; testmat(2,1) = 3.5; testmat(2,2) = 6.;

    auto fastor2 = arma_to_fastor2(testmat);
    mat testmat2 = fastor2_to_arma(Fastor::Tensor<double,3,3>(fastor2));
    EXPECT_LT(norm(testmat2 - testmat, 2), simcoon::iota);

    // Test arma 6x6 <-> Fastor rank-4 round-trip (stiffness convention)
    mat L = zeros(6,6);
    L(0,0) = 16.; L(0,1) = 16.; L(0,2) = 2.25; L(0,3) = 16.; L(0,4) = 6.; L(0,5) = 6.;
    L(1,0) = 16.; L(1,1) = 4.; L(1,2) = 12.25; L(1,3) = 8.; L(1,4) = 14.; L(1,5) = 7.;
    L(2,0) = 2.25; L(2,1) = 12.25; L(2,2) = 36; L(2,3) = 5.25; L(2,4) = 9.; L(2,5) = 21.;
    L(3,0) = 16.; L(3,1) = 8.; L(3,2) = 5.25; L(3,3) = 12.; L(3,4) = 10.; L(3,5) = 8.5;
    L(4,0) = 6.; L(4,1) = 14.; L(4,2) = 9.; L(4,3) = 10.; L(4,4) = 13.125; L(4,5) = 14.625;
    L(5,0) = 6.; L(5,1) = 7.; L(5,2) = 21.; L(5,3) = 8.5; L(5,4) = 14.625; L(5,5) = 12.125;

    auto C = voigt_to_fastor4(mat::fixed<6,6>(L));
    mat L2 = fastor4_to_voigt(C);
    EXPECT_LT(norm(L2 - L, 2), simcoon::iota);
}

TEST(Ttransfer, t2v_sym_v2t_sym_roundtrip)
{
    // Symmetric tensor -> 6-component vector round-trip
    mat m = {{4., 2., 1.}, {2., 5., 3.}, {1., 3., 6.}};
    vec v = t2v_sym(m);
    EXPECT_EQ(v.n_elem, (arma::uword)6);

    mat m_back = v2t_sym(v);
    EXPECT_LT(norm(m_back - m, 2), simcoon::iota);
    // m_back should be symmetric
    EXPECT_LT(norm(m_back - m_back.t(), 2), simcoon::iota);
}

TEST(Ttransfer, v2t_skewsym)
{
    // v2t_skewsym requires 6 components: (11,22,33,12,13,23)
    // Returns a skew-symmetric matrix where off-diagonal terms have opposite signs
    vec v = {0., 0., 0., 1., 2., 3.};  // diag=0, off-diag=1,2,3
    mat m = v2t_skewsym(v);

    // Should be antisymmetric: m(i,j) = -m(j,i) for off-diagonal
    EXPECT_LT(fabs(m(0, 1) + m(1, 0)), simcoon::iota);
    EXPECT_LT(fabs(m(0, 2) + m(2, 0)), simcoon::iota);
    EXPECT_LT(fabs(m(1, 2) + m(2, 1)), simcoon::iota);
}

TEST(Ttransfer, v2t_9comp_roundtrip)
{
    // 9-component vector -> 3x3 matrix (general, not symmetric)
    // Mapping: v(i*3+j) -> m(i,j)
    vec v = {1., 2., 3., 4., 5., 6., 7., 8., 9.};
    mat m = v2t(v);
    EXPECT_EQ(m.n_rows, (arma::uword)3);
    EXPECT_EQ(m.n_cols, (arma::uword)3);

    // Verify entries: m(i,j) = v(i*3+j)
    EXPECT_LT(fabs(m(0, 0) - v(0)), simcoon::iota);  // v(0) -> m(0,0)
    EXPECT_LT(fabs(m(0, 1) - v(1)), simcoon::iota);  // v(1) -> m(0,1)
    EXPECT_LT(fabs(m(1, 1) - v(4)), simcoon::iota);  // v(4) -> m(1,1)
    EXPECT_LT(fabs(m(2, 2) - v(8)), simcoon::iota);  // v(8) -> m(2,2)
}

