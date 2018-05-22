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
#define BOOST_TEST_MODULE "eshelby"
#include <boost/test/unit_test.hpp>

#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_Mechanics/Homogenization/eshelby.hpp>

using namespace std;
using namespace arma;
using namespace simcoon;

BOOST_AUTO_TEST_CASE( S_TII )
{
    
    double a1 = 1.;
    double a2 = 1.;
    double a3 = 1.;
    int mp = 300;
    int np = 300;
    
    double E = 70000.;
    double nu = 0.3;
    
    double mu = E/(2.*(1+nu));
    double lambda = E*nu/((1.+nu)*(1.-2.*nu));
    
    mat S_anal = zeros(6,6);
    mat S_num = zeros(6,6);
    mat T_II_num = zeros(6,6);
    
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
    
    vec x = zeros(mp);
    vec wx = zeros(mp);
    vec y = zeros(np);
    vec wy = zeros(np);
    points(x, wx, y, wy, mp, np);
    
    S_anal = Eshelby_sphere(nu);
    S_num = Eshelby(Lt, a1, a2, a3, x, wx, y, wy, mp, np);
    BOOST_CHECK( norm(S_num-S_anal,2) < 1.E-9 );

    cout << S_num-S_anal << "\n";
    
    T_II_num = T_II(Lt, a1, a2, a3, x, wx, y, wy, mp, np);
    BOOST_CHECK( norm(T_II_num*Lt-S_anal,2) < 1.E-9 );
    
    a1 = 1.E16;
    S_anal = Eshelby_cylinder(nu);
    S_num = Eshelby(Lt, a1, a2, a3, x, wx, y, wy, mp, np);
    BOOST_CHECK( norm(S_num-S_anal,2) < 1.E-4 );
    
    cout << S_num-S_anal << "\n";
    
    T_II_num = T_II(Lt, a1, a2, a3, x, wx, y, wy, mp, np);
    BOOST_CHECK( norm(T_II_num*Lt-S_anal,2) < 1.E-4 );
    
}
