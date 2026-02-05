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

///@file Tderivatives.cpp
///@brief Test for tensor derivative functions
///@version 1.0

#include <gtest/gtest.h>
#include <armadillo>

#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Functions/derivatives.hpp>
#include <simcoon/Continuum_mechanics/Functions/transfer.hpp>

using namespace std;
using namespace arma;
using namespace simcoon;

TEST(Tderivatives, dI1DS_identity)
{
    // dI1/dS = I (identity) for any input
    mat S = {{1., 2., 3.}, {4., 5., 6.}, {7., 8., 9.}};
    mat result = dI1DS(S);
    EXPECT_LT(norm(result - eye(3, 3), 2), simcoon::iota);

    // Also for identity input
    mat result2 = dI1DS(eye(3, 3));
    EXPECT_LT(norm(result2 - eye(3, 3), 2), simcoon::iota);

    // Also for zero input
    mat result3 = dI1DS(zeros(3, 3));
    EXPECT_LT(norm(result3 - eye(3, 3), 2), simcoon::iota);
}

TEST(Tderivatives, dI2DS_returns_S)
{
    // dI2/dS = S (the tensor itself)
    mat S = {{2., 1., 0.}, {1., 3., 1.}, {0., 1., 4.}};
    mat result = dI2DS(S);
    EXPECT_LT(norm(result - S, 2), simcoon::iota);

    // Identity case
    mat result2 = dI2DS(eye(3, 3));
    EXPECT_LT(norm(result2 - eye(3, 3), 2), simcoon::iota);
}

TEST(Tderivatives, dI3DS_returns_S_squared_transposed)
{
    // dI3/dS = (S*S)^T
    mat S = {{2., 1., 0.}, {1., 3., 1.}, {0., 1., 4.}};
    mat result = dI3DS(S);
    mat expected = (S * S).t();
    EXPECT_LT(norm(result - expected, 2), simcoon::iota);
}

TEST(Tderivatives, dtrSdS_identity)
{
    // dtr(S)/dS = I for any input (same as dI1DS)
    mat S = {{5., 2., 1.}, {2., 3., 0.}, {1., 0., 7.}};
    mat result = dtrSdS(S);
    EXPECT_LT(norm(result - eye(3, 3), 2), simcoon::iota);
}

TEST(Tderivatives, ddetSdS_cofactor)
{
    // ddet(S)/dS = det(S) * inv(S)^T = cofactor matrix
    mat S = {{2., 1., 0.}, {1., 3., 1.}, {0., 1., 4.}};
    mat result = ddetSdS(S);
    mat expected = det(S) * inv(S).t();
    EXPECT_LT(norm(result - expected, 2), 1.E-9);
}

TEST(Tderivatives, ddetSdS_identity)
{
    // For S = I: det(I)=1, inv(I)^T = I, so result = I
    mat result = ddetSdS(eye(3, 3));
    EXPECT_LT(norm(result - eye(3, 3), 2), 1.E-9);
}

TEST(Tderivatives, dinvSdSsym_basic)
{
    // Verify dinvSdSsym returns 6x6 matrix
    mat S = {{4., 1., 0.5}, {1., 3., 0.8}, {0.5, 0.8, 5.}};
    mat dinv = dinvSdSsym(S);

    // dinv is a 6x6 matrix (Voigt notation)
    EXPECT_EQ(dinv.n_rows, (arma::uword)6);
    EXPECT_EQ(dinv.n_cols, (arma::uword)6);

    // Should have major symmetry for elastic-type tensor
    // (Voigt notation 4th order tensor)
    EXPECT_LT(norm(dinv - dinv.t(), 2), 1.E-9);
}

TEST(Tderivatives, ddetSdS_numerical)
{
    // Verify ddetSdS with numerical differentiation
    mat S = {{4., 1., 0.5}, {1., 3., 0.8}, {0.5, 0.8, 5.}};
    mat analytical = ddetSdS(S);

    double eps = 1.E-7;
    mat numerical = zeros(3, 3);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            mat S_plus = S;
            mat S_minus = S;
            S_plus(i, j) += eps;
            S_minus(i, j) -= eps;
            numerical(i, j) = (det(S_plus) - det(S_minus)) / (2. * eps);
        }
    }
    EXPECT_LT(norm(analytical - numerical, 2), 1.E-4);
}
