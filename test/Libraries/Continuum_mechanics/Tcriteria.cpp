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

///@file Tcriteria.cpp
///@brief Test for stress criteria
///@version 1.0

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "criteria"
#include <boost/test/unit_test.hpp>

#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Continuum_mechanics/Functions/criteria.hpp>

using namespace std;
using namespace arma;
using namespace simcoon;

BOOST_AUTO_TEST_CASE( Prager )
{
    double b = 0.;
    double n = 2.;
    vec sigma = {400.,0.,0.,0.,0.,0.};

    //Check that Prager is equal to 400 in that case
    double test_tresca = Tresca_stress(sigma);
    BOOST_CHECK( test_tresca - 400. < iota );
    
    //Check that Prager is equal to 400 in that case
    double test_prager = Prager_stress(sigma, b, n);
    BOOST_CHECK( test_prager - 400. < iota );
    
}