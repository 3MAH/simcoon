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
    EXPECT_LT(norm(result - eye(3, 3), 2), sim_iota);

    // Also for identity input
    mat result2 = dI1DS(eye(3, 3));
    EXPECT_LT(norm(result2 - eye(3, 3), 2), sim_iota);

    // Also for zero input
    mat result3 = dI1DS(zeros(3, 3));
    EXPECT_LT(norm(result3 - eye(3, 3), 2), sim_iota);
}

TEST(Tderivatives, dI2DS_returns_S)
{
    // dI2/dS = S (the tensor itself)
    mat S = {{2., 1., 0.}, {1., 3., 1.}, {0., 1., 4.}};
    mat result = dI2DS(S);
    EXPECT_LT(norm(result - S, 2), sim_iota);

    // Identity case
    mat result2 = dI2DS(eye(3, 3));
    EXPECT_LT(norm(result2 - eye(3, 3), 2), sim_iota);
}

TEST(Tderivatives, dI3DS_returns_S_squared_transposed)
{
    // dI3/dS = (S*S)^T
    mat S = {{2., 1., 0.}, {1., 3., 1.}, {0., 1., 4.}};
    mat result = dI3DS(S);
    mat expected = (S * S).t();
    EXPECT_LT(norm(result - expected, 2), sim_iota);
}

TEST(Tderivatives, dtrSdS_identity)
{
    // dtr(S)/dS = I for any input (same as dI1DS)
    mat S = {{5., 2., 1.}, {2., 3., 0.}, {1., 0., 7.}};
    mat result = dtrSdS(S);
    EXPECT_LT(norm(result - eye(3, 3), 2), sim_iota);
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

TEST(Tderivatives, dinvSdSsym_finite_diff)
{
    // Verify dinvSdSsym using finite differences
    mat S = {{4., 1., 0.5}, {1., 3., 0.8}, {0.5, 0.8, 5.}};
    mat dinv = dinvSdSsym(S);

    // dinv is a 6x6 matrix (Voigt notation)
    // Verify by finite differences on the Voigt representation
    double eps = 1.E-6;
    mat dinv_num = zeros(6, 6);

    // Perturb each Voigt component and check
    for (int k = 0; k < 6; k++) {
        mat S_plus = S;
        mat S_minus = S;

        // Map Voigt index to (i,j)
        int row, col;
        if (k < 3) {
            row = k;
            col = k;
        }
        else if (k == 3) {
            row = 0;
            col = 1;
        }
        else if (k == 4) {
            row = 0;
            col = 2;
        }
        else {
            row = 1;
            col = 2;
        }

        S_plus(row, col) += eps;
        S_plus(col, row) += eps;
        S_minus(row, col) -= eps;
        S_minus(col, row) -= eps;

        // Factor 2 for off-diagonal in symmetric perturbation
        double factor = (k < 3) ? 1.0 : 2.0;
        S_plus(row, col) = S(row, col) + eps;
        S_plus(col, row) = S(col, row) + eps;
        S_minus(row, col) = S(row, col) - eps;
        S_minus(col, row) = S(col, row) - eps;

        mat inv_plus = inv(S_plus);
        mat inv_minus = inv(S_minus);
        mat diff = (inv_plus - inv_minus) / (2. * eps);

        // Extract Voigt vector from diff
        dinv_num(0, k) = diff(0, 0);
        dinv_num(1, k) = diff(1, 1);
        dinv_num(2, k) = diff(2, 2);
        dinv_num(3, k) = diff(0, 1);
        dinv_num(4, k) = diff(0, 2);
        dinv_num(5, k) = diff(1, 2);
    }

    EXPECT_LT(norm(dinv - dinv_num, 2), 1.E-4);
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
