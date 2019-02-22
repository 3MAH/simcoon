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

BOOST_AUTO_TEST_CASE( dev_sph )
{
    
    mat test = zeros(3,3);
    test(0,0) = 1.;
    
    mat testdev = zeros(3,3);
    testdev(0,0) = 2./3.;
    testdev(1,1) = -1./3.;
    testdev(1,1) = -1./3.;
    
    mat testsph = (1./3.)*eye(3,3);
 
    BOOST_CHECK( fnorm(test - testdev,2) < sim_iota );
    BOOST_CHECK( fnorm(test - testsph,2) < sim_iota );
}
