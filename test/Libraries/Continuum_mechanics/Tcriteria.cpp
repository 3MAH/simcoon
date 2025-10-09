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

#include <gtest/gtest.h>
#include <armadillo>

#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Continuum_mechanics/Functions/criteria.hpp>

using namespace std;
using namespace arma;
using namespace simcoon;

TEST(Tcriteria, Aniso)
{
    double b = 0.;
    double n = 2.;
    vec sigma = {400., 0., 0., 0., 0., 0.};

    // Check that Tresca is equal to 400 in that case
    double test_tresca = Tresca_stress(sigma);
    EXPECT_LT(test_tresca - 400., sim_iota);

    // Check that Drucker is equal to 400 in that case
    double test_prager = Drucker_stress(sigma, b, n);
    EXPECT_LT(test_prager - 400., sim_iota);
}