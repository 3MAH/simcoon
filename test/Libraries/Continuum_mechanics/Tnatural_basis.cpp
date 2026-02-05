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

///@file Tnatural_basis.cpp
///@brief Test for natural_basis class
///@version 1.0

#include <gtest/gtest.h>
#include <armadillo>
#include <sstream>

#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Functions/natural_basis.hpp>

using namespace std;
using namespace arma;
using namespace simcoon;

TEST(Tnatural_basis, default_constructor)
{
    natural_basis nb;
    // Default: all vectors zero, metrics zero
    for (int i = 0; i < 3; i++) {
        EXPECT_LT(norm(nb.g_i[i], 2), simcoon::iota);
        EXPECT_LT(norm(nb.g0i[i], 2), simcoon::iota);
    }
    EXPECT_LT(norm(nb.g_ij, 2), simcoon::iota);
    EXPECT_LT(norm(nb.g0ij, 2), simcoon::iota);
}

TEST(Tnatural_basis, cartesian_constructor)
{
    // Cartesian basis: e1, e2, e3
    vector<vec> basis(3);
    basis[0] = {1., 0., 0.};
    basis[1] = {0., 1., 0.};
    basis[2] = {0., 0., 1.};

    natural_basis nb(basis);

    // Covariant metric should be identity
    EXPECT_LT(norm(nb.g_ij - eye(3, 3), 2), 1.E-9);

    // Contravariant metric should be identity
    EXPECT_LT(norm(nb.g0ij - eye(3, 3), 2), 1.E-9);

    // Covariant and contravariant vectors should match the basis
    for (int i = 0; i < 3; i++) {
        EXPECT_LT(norm(nb.g_i[i] - basis[i], 2), 1.E-9);
    }
}

TEST(Tnatural_basis, update)
{
    natural_basis nb;

    vector<vec> basis(3);
    basis[0] = {2., 0., 0.};
    basis[1] = {0., 3., 0.};
    basis[2] = {0., 0., 4.};

    nb.update(basis);

    // Covariant metric g_ij = diag(4, 9, 16)
    mat expected_g = diagmat(vec({4., 9., 16.}));
    EXPECT_LT(norm(nb.g_ij - expected_g, 2), 1.E-9);

    // Contravariant metric g0ij = diag(1/4, 1/9, 1/16)
    mat expected_g0 = diagmat(vec({0.25, 1. / 9., 1. / 16.}));
    EXPECT_LT(norm(nb.g0ij - expected_g0, 2), 1.E-9);
}

TEST(Tnatural_basis, from_F_identity)
{
    natural_basis nb;
    mat F = eye(3, 3);
    nb.from_F(F);

    // Identity F: basis vectors are Cartesian
    EXPECT_LT(norm(nb.g_ij - eye(3, 3), 2), 1.E-9);
    EXPECT_LT(norm(nb.g0ij - eye(3, 3), 2), 1.E-9);
}

TEST(Tnatural_basis, from_F_stretch)
{
    natural_basis nb;
    mat F = diagmat(vec({2., 3., 4.}));
    nb.from_F(F);

    // g_i[0] = F * e1 = (2,0,0), etc. -> columns of F
    vec expected0 = {2., 0., 0.};
    vec expected1 = {0., 3., 0.};
    vec expected2 = {0., 0., 4.};
    EXPECT_LT(norm(nb.g_i[0] - expected0, 2), 1.E-9);
    EXPECT_LT(norm(nb.g_i[1] - expected1, 2), 1.E-9);
    EXPECT_LT(norm(nb.g_i[2] - expected2, 2), 1.E-9);

    // g_ij = F^T * F = C (right Cauchy-Green)
    mat expected_g = F.t() * F;
    EXPECT_LT(norm(nb.g_ij - expected_g, 2), 1.E-9);
}

TEST(Tnatural_basis, copy_constructor)
{
    vector<vec> basis(3);
    basis[0] = {1., 0., 0.};
    basis[1] = {0., 2., 0.};
    basis[2] = {0., 0., 3.};

    natural_basis nb1(basis);
    natural_basis nb2(nb1);

    EXPECT_LT(norm(nb2.g_ij - nb1.g_ij, 2), simcoon::iota);
    EXPECT_LT(norm(nb2.g0ij - nb1.g0ij, 2), simcoon::iota);
    for (int i = 0; i < 3; i++) {
        EXPECT_LT(norm(nb2.g_i[i] - nb1.g_i[i], 2), simcoon::iota);
        EXPECT_LT(norm(nb2.g0i[i] - nb1.g0i[i], 2), simcoon::iota);
    }
}

TEST(Tnatural_basis, assignment_operator)
{
    vector<vec> basis(3);
    basis[0] = {1., 0., 0.};
    basis[1] = {0., 2., 0.};
    basis[2] = {0., 0., 3.};

    natural_basis nb1(basis);
    natural_basis nb2;
    nb2 = nb1;

    EXPECT_LT(norm(nb2.g_ij - nb1.g_ij, 2), simcoon::iota);
    EXPECT_LT(norm(nb2.g0ij - nb1.g0ij, 2), simcoon::iota);
}

TEST(Tnatural_basis, stream_operator)
{
    vector<vec> basis(3);
    basis[0] = {1., 0., 0.};
    basis[1] = {0., 1., 0.};
    basis[2] = {0., 0., 1.};

    natural_basis nb(basis);
    ostringstream oss;
    ASSERT_NO_THROW(oss << nb);
    // Stream should not be empty
    EXPECT_GT(oss.str().length(), (size_t)0);
}
