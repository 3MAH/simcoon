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

///@file Tnum_solve.cpp
///@brief Test for numerical solvers (Newton-Raphson, Fischer-Burmeister)
///@version 1.0

#include <gtest/gtest.h>
#include <armadillo>

#include <simcoon/parameter.hpp>
#include <simcoon/Simulation/Maths/num_solve.hpp>

using namespace std;
using namespace arma;
using namespace simcoon;

TEST(Tnum_solve, Newton_Raphon_simple)
{
    // Solve a simple 1D problem: Phi(x) = x^2 - 4 = 0, solution x=2
    // We simulate one Newton step from x=3
    // Phi = 3^2 - 4 = 5, dPhi/dx = 2*3 = 6
    vec Phi = {5.};
    vec Y_crit = {1.};
    mat denom = zeros(1, 1);
    denom(0, 0) = 6.;
    vec Dp = {3.};
    vec dp = zeros(1);
    double error = 0.;

    Newton_Raphon(Phi, Y_crit, denom, Dp, dp, error);

    // dp should be -Phi/denom = -5/6
    EXPECT_LT(fabs(dp(0) - (-5./6.)), 1.E-9);

    // Dp should be updated: 3 + (-5/6) = 13/6
    EXPECT_LT(fabs(Dp(0) - (3. - 5./6.)), 1.E-9);

    // error should be |Phi|/|Y_crit| = 5
    EXPECT_LT(fabs(error - 5.), 1.E-9);
}

TEST(Tnum_solve, Newton_Raphon_2D)
{
    // 2D system: Phi = [2, 3], Y_crit = [1, 1], denom = [[4, 0], [0, 4]]
    vec Phi = {2., 3.};
    vec Y_crit = {1., 1.};
    mat denom = {{4., 0.}, {0., 4.}};
    vec Dp = {1., 1.};
    vec dp = zeros(2);
    double error = 0.;

    Newton_Raphon(Phi, Y_crit, denom, Dp, dp, error);

    // dp = -solve(denom, Phi) = -[0.5, 0.75]
    EXPECT_LT(fabs(dp(0) + 0.5), 1.E-9);
    EXPECT_LT(fabs(dp(1) + 0.75), 1.E-9);

    // error = |2|/|1| + |3|/|1| = 5
    EXPECT_LT(fabs(error - 5.), 1.E-9);
}

TEST(Tnum_solve, Newton_Raphon_singular)
{
    // Singular denom should yield dp = 0
    vec Phi = {2.};
    vec Y_crit = {1.};
    mat denom = zeros(1, 1);
    vec Dp = {1.};
    vec dp = zeros(1);
    double error = 0.;

    Newton_Raphon(Phi, Y_crit, denom, Dp, dp, error);

    // dp should remain zero
    EXPECT_LT(fabs(dp(0)), 1.E-9);
    // Dp unchanged at 1
    EXPECT_LT(fabs(Dp(0) - 1.), 1.E-9);
}

TEST(Tnum_solve, Fischer_Burmeister_m_basic)
{
    // Test with a simple case: Phi > 0, Dp > 0
    vec Phi = {1.};
    vec Y_crit = {1.};
    mat denom = zeros(1, 1);
    denom(0, 0) = 2.;
    vec Dp = {1.};
    vec dp = zeros(1);
    double error = 0.;

    Fischer_Burmeister_m(Phi, Y_crit, denom, Dp, dp, error);

    // Should produce a correction step
    // FB = sqrt(Phi^2 + Dpstar^2) + Phi - Dpstar
    // Dpstar = Dp * |denom(0,0)| = 1 * 2 = 2
    double FB = sqrt(1. + 4.) + 1. - 2.;
    EXPECT_LT(fabs(error - fabs(FB)), 1.E-6);
}

TEST(Tnum_solve, denom_FB_m_basic)
{
    // Test denom_FB_m output
    vec Phi = {2.};
    mat denom = zeros(1, 1);
    denom(0, 0) = 3.;
    vec Dp = {1.};

    mat result = denom_FB_m(Phi, denom, Dp);
    EXPECT_EQ(result.n_rows, (arma::uword)1);
    EXPECT_EQ(result.n_cols, (arma::uword)1);

    // Dpstar = Dp * |denom(0,0)| = 1 * 3 = 3
    // dDpstar = 3
    // denomFB = (Phi/sqrt(Phi^2+Dpstar^2)+1)*denom + delta*dDpstar*(Dpstar/sqrt(Phi^2+Dpstar^2)-1)
    double Dpstar = 3.;
    double norm_val = sqrt(4. + 9.);
    double expected = (2./norm_val + 1.)*3. + 3.*(3./norm_val - 1.);
    EXPECT_LT(fabs(result(0,0) - expected), 1.E-9);
}

TEST(Tnum_solve, Fischer_Burmeister_zeros)
{
    // When both Phi and Dp are zero, FB = 0, error = 0
    vec Phi = {0.};
    vec Y_crit = {1.};
    mat denom = eye(1, 1);
    vec Dp = {0.};
    vec dp = zeros(1);
    double error = 0.;

    Fischer_Burmeister(Phi, Y_crit, denom, Dp, dp, error);

    EXPECT_LT(fabs(error), 1.E-9);
}
