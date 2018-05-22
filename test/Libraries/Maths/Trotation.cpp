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
#define BOOST_TEST_MODULE "rotation"
#include <boost/test/unit_test.hpp>

#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Simulation/Maths/rotation.hpp>
#include <simcoon/Continuum_Mechanics/Homogenization/eshelby.hpp>

using namespace std;
using namespace arma;
using namespace simcoon;

BOOST_AUTO_TEST_CASE( basic_rotation )
{
    //Basic rotation
    vec u = {1,1,1};
    double psi = 23.;
    double theta = 42.;
    double phi = 165;
    
    psi *= (pi/180.);
    theta *= (pi/180.);
    phi *= (pi/180.);
    
    mat R1 = { {cos(psi),-sin(psi),0}, {sin(psi),cos(psi), 0}, {0,0,1}};
    mat R2 = { {1,0,0}, {0,cos(theta),-sin(theta)}, {0,sin(theta),cos(theta)} };
    mat R2b = { {cos(theta),0,sin(theta)}, {0,1,0}, {-sin(theta),0,cos(theta)} };
    mat R3 = { {cos(phi),-sin(phi),0}, {sin(phi),cos(phi), 0}, {0,0,1}};
    
    mat R_zxz = R3*R2*R1;
    mat R_zyz = R3*R2b*R1;
    
    mat R_zxz_a = fillR(psi, theta, phi, true, "zxz");
    mat R_zxz_p = fillR(psi, theta, phi, false, "zxz");
    mat R_zyz_a = fillR(psi, theta, phi, true, "zyz");
    mat R_zyz_p = fillR(psi, theta, phi, false, "zyz");
    
    cout << "R_zxz" << R_zxz << endl;
    cout << "R_zxz_a" << R_zxz_a << endl;
    cout << "R_zxz_p" << R_zxz_p << endl;
    
    cout << "R_zyz" << R_zyz << endl;
    cout << "R_zyz_a" << R_zyz_a << endl;
    cout << "R_zyz_p" << R_zyz_p << endl;
    
    vec upp = R_zxz*u;
    
    vec upp2 = R1*u;
    upp2 = R2*upp2;
    upp2 = R3*upp2;

    BOOST_CHECK( norm(R_zxz-R_zxz_a,2) < 1.E-9 );
    BOOST_CHECK( norm(trans(R_zxz)-R_zxz_p,2) < 1.E-9 );

    BOOST_CHECK( norm(R_zyz-R_zyz_a,2) < 1.E-9 );
    BOOST_CHECK( norm(trans(R_zyz)-R_zyz_p,2) < 1.E-9 );
    
    BOOST_CHECK( norm(upp-upp2,2) < 1.E-9 );
    
}

BOOST_AUTO_TEST_CASE( rotation )
{
    vec test = zeros(6);
    test(0) = 4.;
    test(1) = 2.;
    test(2) = 6.;
    test(3) = 8.;
    test(4) = 3.;
    test(5) = 7.;

    double psi = 12.5;
    double theta = 32.;
    double phi = -4.5;
    
    psi *= (pi/180.);
    theta *= (pi/180.);
    phi *= (pi/180.);
    
    //Rotation
    mat R1 = { {cos(psi),-sin(psi),0}, {sin(psi),cos(psi), 0}, {0,0,1}};
    mat R2 = { {1,0,0}, {0,cos(theta),-sin(theta)}, {0,sin(theta),cos(theta)} };
    mat R3 = { {cos(phi),-sin(phi),0}, {sin(phi),cos(phi), 0}, {0,0,1}};
    
    mat R = R3*R2*R1;
//    mat Rb = trans(R3)*trans(R2)*trans(R1);

    //test of rotate A
    mat S_c = Eshelby_cylinder(0.12);
    
    mat S_c1 = rotateA(S_c, psi, 3);
    mat S_c2 = rotateA(S_c, theta, 1);
    mat S_c3 = rotateA(S_c, phi, 3);
    
    mat S_c11 = rotateA(S_c,R1);
    mat S_c22 = rotateA(S_c,R2);
    mat S_c33 = rotateA(S_c,R3);

    mat S_cR = rotateA(S_c, R1);
    S_cR = rotateA(S_cR, R2);
    S_cR = rotateA(S_cR, R3);
    
    mat S_cRR = rotateA(S_c,R);
    
    BOOST_CHECK( norm(S_c33-S_c3,2) < 1.E-9 );
    BOOST_CHECK( norm(S_c22-S_c2,2) < 1.E-9 );
    BOOST_CHECK( norm(S_c11-S_c1,2) < 1.E-9 );
    
    BOOST_CHECK( norm(S_cRR-S_cR,2) < 1.E-9 );
    
/*    vec E_tot = {0.01,-0.033,-0.033,0.02,0.012,0.005}
    vec E_r = S_c*E_tot*/
    
    //test of rotate a matrix
    mat a = {{1.2,0.4,1.6}, {2.1, -0.9,1.3}, {0.8, -0.7, 0.02} };
    assert(axis_psi == 3);
    assert(axis_theta == 1);
    assert(axis_phi == 3);

    mat Rp2 = fillR(psi,theta,phi);
    mat Rp3 = fillR(psi,theta,phi, true, "zxz");
    
    mat a1 = rotate_mat(a, R);
    mat a2 = rotate_mat(a, Rp2);
    mat a3 = rotate_mat(a, Rp3);

    BOOST_CHECK( norm(a1-a2,2) < 1.E-9 );
    BOOST_CHECK( norm(a1-a3,2) < 1.E-9 );
    
}
