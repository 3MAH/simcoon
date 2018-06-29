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

///@file Tcontimech.cpp
///@brief Test for continuum mechanics tensors in Voigt notation
///@version 1.0

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "contimech"
#include <boost/test/unit_test.hpp>

#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>

using namespace std;
using namespace arma;
using namespace simcoon;

BOOST_AUTO_TEST_CASE( tr_dev_Mises_eta )
{
    
    vec test = zeros(6);
    test(0) = 4.;
    test(1) = 2.;
    test(2) = 6.;
    test(3) = 8.;
    test(4) = 3.;
    test(5) = 7.;
    
    //Test of tr function
    double tr1 = tr(zeros(6));
    double tr2 = tr(test);
    BOOST_CHECK( fabs(tr1 - 0.) < iota );
    BOOST_CHECK( fabs(tr2 - 12.) < iota );

    //Test of dev function
    vec dev1 = dev(zeros(6));
    vec dev2 = dev(test);
    double trdev2 = tr(dev2);
	vec testdev = zeros(6);
    testdev(0) = 0.;
    testdev(1) = -2.;
    testdev(2) = 2.;
    testdev(3) = 8.;
    testdev(4) = 3.;
    testdev(5) = 7.;
    BOOST_CHECK( fabs(trdev2 - 0.) < iota );
    BOOST_CHECK( norm(dev2 - testdev,2) < iota );
    
    //Test of Mises_stress function
    double dstress = 3.*sqrt(42.);
    BOOST_CHECK( fabs(Mises_stress(test) - dstress) < iota );

    //Test of eta_stress function
    dev1 = eta_stress(zeros(6));
    dev2 = eta_stress(test);
	vec testeta = zeros(6);
    testeta(0) = 0.;
    testeta(1) = -2;
    testeta(2) = 2.;
    testeta(3) = 8.*2;
    testeta(4) = 3.*2;
    testeta(5) = 7.*2;
    BOOST_CHECK( norm(dev1,2) < iota );
    BOOST_CHECK( norm(dev2 - (3./2.)*testeta/dstress,2) < iota );
    
    //Test of Mises_strain function
    double dstrain = sqrt(46.);
    BOOST_CHECK( fabs(Mises_strain(test) - dstrain) < iota );

    //Test of eta_strain function
    dev1 = eta_strain(zeros(6));
    dev2 = eta_strain(test);
    BOOST_CHECK( norm(dev1,2) < iota );
    BOOST_CHECK( norm(dev2 - (2./3.)*testdev/dstrain,2) < iota );
    
}
    
BOOST_AUTO_TEST_CASE( J2_J3 )
{

    vec test = zeros(6);
    test(0) = 4.;
    test(1) = 2.;
    test(2) = 6.;
    test(3) = 8.;
    test(4) = 3.;
    test(5) = 7.;
    
	//Test of J2_stress function
    double J2_stress1 = J2_stress(zeros(6));
    double J2_stress2 = J2_stress(test);
    BOOST_CHECK( fabs(J2_stress1) < iota );
    BOOST_CHECK( fabs(J2_stress2 - 126.) < iota );
	
	//Test of J2_strain function
	double J2_strain1 = J2_strain(zeros(6));
	double J2_strain2 = J2_strain(test);
    BOOST_CHECK( fabs(J2_strain1) < iota );
    BOOST_CHECK( fabs(J2_strain2 - 34.5) < iota );

	//Test of J3_stress function
	double J3_stress1 = J3_stress(zeros(6));
	double J3_stress2 = J3_stress(test);
    BOOST_CHECK( fabs(J3_stress1) < iota );
    BOOST_CHECK( fabs(J3_stress2 - 226.) < iota );

	//Test of J3_stress function
	double J3_strain1 = J3_strain(zeros(6));
	double J3_strain2 = J3_strain(test);
    BOOST_CHECK( fabs(J3_strain1) < iota );
    BOOST_CHECK( fabs(J3_strain2 - 14.5) < iota );

}

BOOST_AUTO_TEST_CASE( Macaulay )
{
    BOOST_CHECK( Macaulay_p(-4.) == 0 );
    BOOST_CHECK( Macaulay_p(4.) - 4. == 0 );
    BOOST_CHECK( Macaulay_p(0.) == 0 );

    BOOST_CHECK( Macaulay_n(-4.) + 4 == 0 );
    BOOST_CHECK( Macaulay_n(4.) == 0 );
    BOOST_CHECK( Macaulay_n(0.) == 0 );

    BOOST_CHECK( sign(-4.) == -1. );
    BOOST_CHECK( sign(4.) == 1. );
    cout << "sign(0.)" << sign(0.) << "\n";
    BOOST_CHECK( sign(0.) == 0. );
}

BOOST_AUTO_TEST_CASE( ellipsoid )
{
    double a1 = 1.;
    double a2 = 10.;
    double a3 = 100.;

    double u = 0.;
    double v = 0.;
    vec test = { 0, 0, 1 };
    vec sigma_in = {4., 5., 6., 2., 1., 1.5};
    vec test_sig_int = {5.,2.5};
    vec normal = normal_ellipsoid(u,v,a1,a2,a3);
    BOOST_CHECK( norm(normal-test,2) < iota );
    
    u = 0.;
    v = pi/2.;
    test = { 1, 0, 0 };
    normal = normal_ellipsoid(u,v,a1,a2,a3);
    BOOST_CHECK( norm(normal-test,2) < iota );
    
    u = pi/2.;
    v = pi/2.;
    test = { 0, 1, 0 };
    normal = normal_ellipsoid(u,v,a1,a2,a3);
    BOOST_CHECK( norm(normal-test,2) < iota );
    vec sig_int = sigma_int(sigma_in,a1,a2,a3,u,v);
    BOOST_CHECK( norm(sig_int-test_sig_int,2) < iota );
}

BOOST_AUTO_TEST_CASE( P_ijkl )
{

	vec a = zeros(6);
	for (int i=0; i<3; i++)
		a(i) = 1.;
	
	vec b = ones(6);
    
	mat Ireal = eye(6,6);
	for (int i=3; i<6; i++)
		Ireal(i,i) = 0.5;
    
    vec test = zeros(6);
	test(0) = 4.;
	test(1) = 2.;
	test(2) = 6.;
	test(3) = 4.;
	test(4) = 1.5;
	test(5) = 3.5;
    
	//Test of p_ijkl function
	
	mat result = zeros(6,6);
	
	result(0,0) = 16.;
	result(0,1) = 16.;
	result(0,2) = 2.25;
	result(0,3) = 16.;
	result(0,4) = 6.;
	result(0,5) = 6.;
    
	result(1,0) = 16.;
	result(1,1) = 4.;
	result(1,2) = 12.25;
	result(1,3) = 8.;
	result(1,4) = 14.;
	result(1,5) = 7.;
    
	result(2,0) = 2.25;
	result(2,1) = 12.25;
	result(2,2) = 36;
	result(2,3) = 5.25;
	result(2,4) = 9.;
	result(2,5) = 21.;
	
	result(3,0) = 16.;
	result(3,1) = 8.;
	result(3,2) = 5.25;
	result(3,3) = 12.;
	result(3,4) = 10.;
	result(3,5) = 8.5;
	
	result(4,0) = 6.;
	result(4,1) = 14.;
	result(4,2) = 9.;
	result(4,3) = 10.;
	result(4,4) = 13.125;
	result(4,5) = 14.625;
    
	result(5,0) = 6.;
	result(5,1) = 7.;
	result(5,2) = 21.;
	result(5,3) = 8.5;
	result(5,4) = 14.625;
	result(5,5) = 12.125;
    
	//Test of p_ijkl function
	mat pikjlc = p_ikjl(test);
	mat pikjla = p_ikjl(a);
	mat pikjlb = p_ikjl(b);

    BOOST_CHECK( norm(pikjla - Ireal,2) < iota );
    BOOST_CHECK( norm(pikjlb - ones(6,6),2) < iota );
    BOOST_CHECK( norm(pikjlc - result,2) < iota );
    
}
