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

///@file Tkinematics.cpp
///@brief Test for finite strain kinematics functions
///@version 1.0

#include <gtest/gtest.h>
#include <armadillo>

#include <simcoon/parameter.hpp>
#include <simcoon/exception.hpp>
#include <simcoon/Simulation/Maths/eig_sym3.hpp>
#include <simcoon/Continuum_mechanics/Functions/kinematics.hpp>

using namespace std;
using namespace arma;
using namespace simcoon;

TEST(Tkinematics, ER_to_F_identity)
{
    // E = 0 (zero Green-Lagrange strain) and R = I -> F = I
    mat E = zeros(3, 3);
    mat R = eye(3, 3);
    mat F = ER_to_F(E, R);
    EXPECT_LT(norm(F - eye(3, 3), 2), 1.E-9);
}

TEST(Tkinematics, ER_to_F_pure_stretch)
{
    // Pure stretch: R = I, E = 0.5*(U^2 - I), so F = U
    mat U = diagmat(vec({1.5, 2.0, 0.8}));
    mat E = 0.5 * (U * U - eye(3, 3));
    mat R = eye(3, 3);
    mat F = ER_to_F(E, R);
    EXPECT_LT(norm(F - U, 2), 1.E-9);
}

TEST(Tkinematics, eR_to_F_identity)
{
    // e = 0 (zero log strain) and R = I -> F = I
    mat e = zeros(3, 3);
    mat R = eye(3, 3);
    mat F = eR_to_F(e, R);
    EXPECT_LT(norm(F - eye(3, 3), 2), 1.E-9);
}

TEST(Tkinematics, eR_to_F_pure_stretch)
{
    // Pure stretch: R = I, e = log(V), so F = V
    vec stretches = {1.5, 2.0, 0.8};
    mat V = diagmat(stretches);
    mat e = diagmat(log(stretches));
    mat R = eye(3, 3);
    mat F = eR_to_F(e, R);
    EXPECT_LT(norm(F - V, 2), 1.E-9);
}

TEST(Tkinematics, displacement_gradients)
{
    // G_UdX(F) = F - I, G_Udx(F) = I - inv(F)
    mat F = eye(3, 3);
    F(0, 0) = 1.1;
    F(0, 1) = 0.05;

    mat Lagrangian = G_UdX(F);
    EXPECT_LT(norm(Lagrangian - (F - eye(3, 3)), 2), simcoon::iota);

    mat Eulerian = G_Udx(F);
    EXPECT_LT(norm(Eulerian - (eye(3, 3) - inv(F)), 2), 1.E-9);

    // Identity F -> both gradients are zero
    mat Lag_id = G_UdX(eye(3, 3));
    mat Eul_id = G_Udx(eye(3, 3));
    EXPECT_LT(norm(Lag_id, 2), simcoon::iota);
    EXPECT_LT(norm(Eul_id, 2), 1.E-9);
}

TEST(Tkinematics, Cauchy_Green_tensors)
{
    mat F = eye(3, 3);
    F(0, 0) = 2.;
    F(0, 1) = 0.5;
    F(1, 1) = 1.5;
    F(2, 2) = 0.8;

    // Right Cauchy-Green: C = F^T * F (symmetric positive definite)
    mat C = R_Cauchy_Green(F);
    EXPECT_LT(norm(C - F.t() * F, 2), simcoon::iota);
    EXPECT_LT(norm(C - C.t(), 2), simcoon::iota); // symmetric

    // Left Cauchy-Green: B = F * F^T (symmetric positive definite)
    mat B = L_Cauchy_Green(F);
    EXPECT_LT(norm(B - F * F.t(), 2), simcoon::iota);
    EXPECT_LT(norm(B - B.t(), 2), simcoon::iota); // symmetric

    // Identity F -> C = B = I
    mat C_id = R_Cauchy_Green(eye(3, 3));
    mat B_id = L_Cauchy_Green(eye(3, 3));
    EXPECT_LT(norm(C_id - eye(3, 3), 2), simcoon::iota);
    EXPECT_LT(norm(B_id - eye(3, 3), 2), simcoon::iota);
}

TEST(Tkinematics, RU_VR_decomposition)
{
    mat F = zeros(3, 3);
    F(0, 0) = 2.;
    F(0, 1) = 0.5;
    F(1, 0) = 0.3;
    F(1, 1) = 1.5;
    F(2, 2) = 0.8;

    // RU decomposition: F = R * U
    mat R = zeros(3, 3);
    mat U = zeros(3, 3);
    RU_decomposition(R, U, F);
    EXPECT_LT(norm(F - R * U, 2), 1.E-9);
    // R should be orthogonal: R^T * R = I
    EXPECT_LT(norm(R.t() * R - eye(3, 3), 2), 1.E-9);
    // U should be symmetric positive definite
    EXPECT_LT(norm(U - U.t(), 2), 1.E-9);

    // VR decomposition: F = V * R
    mat V = zeros(3, 3);
    mat Rv = zeros(3, 3);
    VR_decomposition(V, Rv, F);
    EXPECT_LT(norm(F - V * Rv, 2), 1.E-9);
    // R should be orthogonal
    EXPECT_LT(norm(Rv.t() * Rv - eye(3, 3), 2), 1.E-9);
    // V should be symmetric
    EXPECT_LT(norm(V - V.t(), 2), 1.E-9);
}

TEST(Tkinematics, Inv_X)
{
    // Known eigenvalue case: diagonal matrix
    mat X = diagmat(vec({2., 3., 5.}));
    vec I = Inv_X(X);
    // I1 = trace = 10
    EXPECT_LT(fabs(I(0) - 10.), 1.E-9);
    // I2 = 0.5*(tr^2 - tr(X^2)) = 0.5*(100 - 38) = 31
    EXPECT_LT(fabs(I(1) - 31.), 1.E-9);
    // I3 = det = 30
    EXPECT_LT(fabs(I(2) - 30.), 1.E-9);

    // Identity: I1=3, I2=3, I3=1
    vec I_id = Inv_X(eye(3, 3));
    EXPECT_LT(fabs(I_id(0) - 3.), 1.E-9);
    EXPECT_LT(fabs(I_id(1) - 3.), 1.E-9);
    EXPECT_LT(fabs(I_id(2) - 1.), 1.E-9);
}

TEST(Tkinematics, Cauchy_tensor)
{
    // Cauchy(F) = inv(B) = inv(F*F^T)
    mat F = eye(3, 3);
    F(0, 0) = 2.;
    F(1, 1) = 1.5;
    F(2, 2) = 0.8;

    mat c = Cauchy(F);
    mat B = L_Cauchy_Green(F);
    EXPECT_LT(norm(c - inv(B), 2), 1.E-9);
}

TEST(Tkinematics, strain_tensors_identity)
{
    mat F = eye(3, 3);

    // Green-Lagrange: E = 0.5*(C-I) = 0 for F=I
    mat E = Green_Lagrange(F);
    EXPECT_LT(norm(E, 2), 1.E-9);

    // Euler-Almansi: e = 0.5*(I-c) = 0 for F=I
    mat e = Euler_Almansi(F);
    EXPECT_LT(norm(e, 2), 1.E-9);

    // Log strain: h = 0 for F=I
    mat h = Log_strain(F);
    EXPECT_LT(norm(h, 2), 1.E-9);
}

TEST(Tkinematics, strain_tensors_uniaxial)
{
    // Uniaxial stretch: F = diag(lambda, 1, 1)
    double lambda = 1.5;
    mat F = eye(3, 3);
    F(0, 0) = lambda;

    // Green-Lagrange: E_11 = 0.5*(lambda^2 - 1)
    mat E = Green_Lagrange(F);
    EXPECT_LT(fabs(E(0, 0) - 0.5 * (lambda * lambda - 1.)), 1.E-9);
    EXPECT_LT(fabs(E(1, 1)), 1.E-9);
    EXPECT_LT(fabs(E(2, 2)), 1.E-9);

    // Euler-Almansi: e_11 = 0.5*(1 - 1/lambda^2)
    mat e = Euler_Almansi(F);
    EXPECT_LT(fabs(e(0, 0) - 0.5 * (1. - 1. / (lambda * lambda))), 1.E-9);

    // Log strain: h_11 = log(lambda)
    mat h = Log_strain(F);
    EXPECT_LT(fabs(h(0, 0) - log(lambda)), 1.E-9);
}

TEST(Tkinematics, Green_Lagrange_definition)
{
    mat F = zeros(3, 3);
    F(0, 0) = 2.;
    F(0, 1) = 0.5;
    F(1, 0) = 0.3;
    F(1, 1) = 1.5;
    F(2, 2) = 0.8;

    mat E = Green_Lagrange(F);
    mat E_ref = 0.5 * (F.t() * F - eye(3, 3));
    EXPECT_LT(norm(E - E_ref, 2), 1.E-9);
}

TEST(Tkinematics, finite_velocity_gradient)
{
    mat F0 = eye(3, 3);
    mat F1 = eye(3, 3);
    F1(0, 0) = 1.01;
    F1(0, 1) = 0.005;
    double DTime = 0.01;

    // L = (1/DTime)*(F1-F0)*inv(F1)
    mat L = finite_L(F0, F1, DTime);
    mat L_ref = (1. / DTime) * (F1 - F0) * inv(F1);
    EXPECT_LT(norm(L - L_ref, 2), 1.E-9);

    // D should be symmetric part of L
    mat D = finite_D(F0, F1, DTime);
    EXPECT_LT(norm(D - 0.5 * (L + L.t()), 2), 1.E-9);
    EXPECT_LT(norm(D - D.t(), 2), 1.E-9);

    // W should be antisymmetric part of L
    mat W = finite_W(F0, F1, DTime);
    EXPECT_LT(norm(W - 0.5 * (L - L.t()), 2), 1.E-9);
    EXPECT_LT(norm(W + W.t(), 2), 1.E-9);

    // L = D + W
    EXPECT_LT(norm(L - (D + W), 2), 1.E-9);
}

TEST(Tkinematics, finite_Omega)
{
    mat F0 = eye(3, 3);
    mat F1 = eye(3, 3);
    F1(0, 0) = 1.01;
    F1(0, 1) = 0.005;
    F1(1, 0) = -0.005;
    double DTime = 0.01;

    mat Omega = finite_Omega(F0, F1, DTime);
    // Omega should be mostly antisymmetric (small numerical error possible)
    EXPECT_LT(norm(Omega + Omega.t(), 2), 0.01);
}

TEST(Tkinematics, finite_DQ)
{
    mat Omega = zeros(3, 3);
    Omega(0, 1) = 0.1;
    Omega(1, 0) = -0.1;
    double DTime = 0.01;

    // DQ = (I + 0.5*DTime*Omega0) * inv(I - 0.5*DTime*Omega1)
    mat DQ = finite_DQ(Omega, Omega, DTime);
    mat expected = (eye(3, 3) + 0.5 * DTime * Omega) * inv(eye(3, 3) - 0.5 * DTime * Omega);
    EXPECT_LT(norm(DQ - expected, 2), 1.E-9);

    // For zero Omega -> DQ = I
    mat DQ_zero = finite_DQ(zeros(3, 3), zeros(3, 3), DTime);
    EXPECT_LT(norm(DQ_zero - eye(3, 3), 2), 1.E-9);
}

TEST(Tkinematics, eig_sym_3x3_matches_eig_sym)
{
    // the fixed-size Jacobi solver against LAPACK: same ascending spectrum, orthonormal basis,
    // reconstruction to roundoff, on random, near-identity, repeated and stiff spectra
    arma_rng::set_seed(1);
    for (int t = 0; t < 400; t++) {
        mat F;
        mat Q0, R0;
        qr(Q0, R0, randn<mat>(3, 3));
        switch (t % 4) {
            case 0: F = eye(3, 3) + 0.3*randn<mat>(3, 3); break;
            case 1: F = eye(3, 3) + 1.e-7*randn<mat>(3, 3); break;
            case 2: F = Q0*diagmat(vec({1.2, 1.2, 0.9}))*Q0.t(); break;
            default: F = Q0*diagmat(vec({30., 1., 0.1}))*Q0.t(); break;
        }
        const mat A = F.t()*F;
        vec l_ref;
        mat Q_ref;
        ASSERT_TRUE(eig_sym(l_ref, Q_ref, A));
        vec::fixed<3> l;
        mat::fixed<3,3> Q;
        ASSERT_TRUE(eig_sym_3x3(l, Q, A));
        EXPECT_LT(norm(vec(l) - l_ref, "inf"), 1.E-12*norm(l_ref, "inf"));
        EXPECT_TRUE(l(0) <= l(1) && l(1) <= l(2));
        EXPECT_LT(abs(mat(Q.t()*Q) - eye(3, 3)).max(), 1.E-14);
        EXPECT_LT(abs(mat(spectral_compose_3x3(Q, l)) - A).max(), 1.E-13*abs(A).max());
    }
    // exact on a diagonal matrix
    vec::fixed<3> l;
    mat::fixed<3,3> Q;
    ASSERT_TRUE(eig_sym_3x3(l, Q, diagmat(vec({2., 0.5, 1.}))));
    EXPECT_EQ(l(0), 0.5);
    EXPECT_EQ(l(1), 1.);
    EXPECT_EQ(l(2), 2.);
    // failure on a non-finite input, as arma::eig_sym
    mat A = eye(3, 3);
    A(0, 1) = datum::nan;
    EXPECT_FALSE(eig_sym_3x3(l, Q, A));
}

TEST(Tkinematics, spectral_functions_3x3)
{
    arma_rng::set_seed(2);
    const mat F = eye(3, 3) + 0.2*randn<mat>(3, 3);
    const mat B = F*F.t();
    mat::fixed<3,3> X;
    ASSERT_TRUE(sqrtmat_sympd_3x3(X, B));
    EXPECT_LT(abs(mat(X) - sqrtmat_sympd(B)).max(), 1.E-13);
    ASSERT_TRUE(logmat_sympd_3x3(X, B));
    EXPECT_LT(abs(mat(X) - logmat_sympd(B)).max(), 1.E-13);
    ASSERT_TRUE(expmat_sym_3x3(X, 0.5*(F + F.t())));
    EXPECT_LT(abs(mat(X) - expmat_sym(0.5*(F + F.t()))).max(), 1.E-13);
    // not SPD: sqrt and log refuse, as their arma counterparts
    const mat N = diagmat(vec({1., -1.e-3, 2.}));
    EXPECT_FALSE(sqrtmat_sympd_3x3(X, N));
    EXPECT_FALSE(logmat_sympd_3x3(X, N));
}

TEST(Tkinematics, RU_spectral_basis)
{
    // one eigen-solve gives R, the eigenbasis of C and, through R N, the one of B = F F^T
    arma_rng::set_seed(3);
    const mat F = eye(3, 3) + 0.25*randn<mat>(3, 3);
    mat::fixed<3,3> R, N;
    vec::fixed<3> c;
    RU_spectral(R, c, N, F);
    mat R_ref, U_ref;
    RU_decomposition(R_ref, U_ref, F);
    EXPECT_LT(abs(mat(R) - R_ref).max(), 1.E-14);
    const mat Nb = R*N;
    EXPECT_LT(abs(mat(Nb*diagmat(vec(c))*Nb.t()) - F*F.t()).max(), 1.E-13*abs(F*F.t()).max());
    // a singular F is reported as RU_decomposition did
    mat F_sing = F;
    F_sing.col(2) = F_sing.col(0);
    EXPECT_THROW(RU_spectral(R, c, N, F_sing), simcoon::simcoon_error);
}
