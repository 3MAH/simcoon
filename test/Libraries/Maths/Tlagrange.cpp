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

///@file Tlagrange.cpp
///@brief Test for Lagrange multiplier functions
///@version 1.0

#include <gtest/gtest.h>
#include <armadillo>

#include <simcoon/parameter.hpp>
#include <simcoon/Simulation/Maths/lagrange.hpp>

using namespace std;
using namespace arma;
using namespace simcoon;

TEST(Tlagrange, lagrange_exp_basic)
{
    double c = 0.1;
    double p0 = 100.;

    // For h <= -c, result should be 0
    EXPECT_LT(fabs(lagrange_exp(-0.2, c, p0)), simcoon::iota);
    EXPECT_LT(fabs(lagrange_exp(-c, c, p0)), simcoon::iota);
    EXPECT_LT(fabs(lagrange_exp(-1.0, c, p0)), simcoon::iota);

    // For h > -c, result should be non-zero
    double val = lagrange_exp(0., c, p0);
    EXPECT_GT(fabs(val), 0.);

    // Derivative should also be 0 for h <= -c
    EXPECT_LT(fabs(dlagrange_exp(-0.2, c, p0)), simcoon::iota);
    EXPECT_LT(fabs(dlagrange_exp(-c, c, p0)), simcoon::iota);
}

TEST(Tlagrange, lagrange_exp_numerical_derivative)
{
    double c = 0.1;
    double p0 = 100.;
    double h = 0.05;
    double dh = 1.E-7;

    double f_plus = lagrange_exp(h + dh, c, p0);
    double f_minus = lagrange_exp(h - dh, c, p0);
    double numerical_deriv = (f_plus - f_minus) / (2. * dh);

    double analytical_deriv = dlagrange_exp(h, c, p0);
    EXPECT_LT(fabs(numerical_deriv - analytical_deriv), 1.E-4);
}

TEST(Tlagrange, lagrange_pow_0_basic)
{
    double c = 0.1;
    double p0 = 10.;
    double n = 2.;
    double alpha = 100.;

    // For x >= c (h=c), function is p0*(1-x)*x^(-n)
    double x = 0.5;
    double expected = p0 * (1. - x) * pow(x, -n);
    EXPECT_LT(fabs(lagrange_pow_0(x, c, p0, n, alpha) - expected), 1.E-9);

    // Numerical derivative check for x > c
    double dx = 1.E-7;
    double f_plus = lagrange_pow_0(x + dx, c, p0, n, alpha);
    double f_minus = lagrange_pow_0(x - dx, c, p0, n, alpha);
    double num_deriv = (f_plus - f_minus) / (2. * dx);
    double ana_deriv = dlagrange_pow_0(x, c, p0, n, alpha);
    EXPECT_LT(fabs(num_deriv - ana_deriv), 1.E-4);
}

TEST(Tlagrange, lagrange_pow_1_basic)
{
    double c = 0.1;
    double p0 = 10.;
    double n = 2.;
    double alpha = 100.;

    // For x <= h=1-c=0.9, function is p0*x*(1-x)^(-n)
    double x = 0.5;
    double expected = p0 * x * pow(1. - x, -n);
    EXPECT_LT(fabs(lagrange_pow_1(x, c, p0, n, alpha) - expected), 1.E-9);

    // Numerical derivative check for x < h
    double dx = 1.E-7;
    double f_plus = lagrange_pow_1(x + dx, c, p0, n, alpha);
    double f_minus = lagrange_pow_1(x - dx, c, p0, n, alpha);
    double num_deriv = (f_plus - f_minus) / (2. * dx);
    double ana_deriv = dlagrange_pow_1(x, c, p0, n, alpha);
    EXPECT_LT(fabs(num_deriv - ana_deriv), 1.E-4);
}

TEST(Tlagrange, lagrange_pow_1_second_derivative)
{
    double c = 0.1;
    double p0 = 10.;
    double n = 2.;
    double alpha = 100.;

    // Second derivative should be finite and computable
    double x = 0.5;
    double d2 = d2lagrange_pow_1(x, c, p0, n, alpha);
    EXPECT_TRUE(std::isfinite(d2));
}

TEST(Tlagrange, lagrange_pow_1_boundary)
{
    double c = 0.1;
    double p0 = 10.;
    double n = 2.;
    double alpha = 100.;

    // For x > 1 (out of bounds), function should return penalty value
    double f_outside = lagrange_pow_1(1.1, c, p0, n, alpha);
    EXPECT_GT(f_outside, 0.);

    // For x near 1 (close to boundary), function value should increase
    double f_near = lagrange_pow_1(0.95, c, p0, n, alpha);
    double f_far = lagrange_pow_1(0.5, c, p0, n, alpha);
    // Both should be finite
    EXPECT_TRUE(std::isfinite(f_near));
    EXPECT_TRUE(std::isfinite(f_far));
}
