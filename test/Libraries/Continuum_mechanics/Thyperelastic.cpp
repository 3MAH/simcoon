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

///@file Thyperelastic.cpp
///@brief Test for hyperelastic fonctions
///@version 1.0

#include <gtest/gtest.h>
#include <armadillo>

#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Functions/transfer.hpp>
#include <simcoon/Continuum_mechanics/Functions/hyperelastic.hpp>

using namespace std;
using namespace arma;
using namespace simcoon;

TEST(Thyperelastic, isochoric_invariants)
{
    mat b_test = eye(3,3);

    vec I_test = simcoon::isochoric_invariants(b_test);

    EXPECT_LT(I_test(0)-3.,1.E-9);
    EXPECT_LT(I_test(1)-3.,1.E-9);    
    EXPECT_LT(I_test(2)-1.,1.E-9);

    vec lambda = eig_sym(b_test);
    lambda.transform( [](double val) { return (sqrt(val)); } );

    vec I_test_lambda = simcoon::isochoric_invariants(lambda);

    EXPECT_LT(I_test_lambda(0)-3.,1.E-9);
    EXPECT_LT(I_test_lambda(1)-3.,1.E-9);    
    EXPECT_LT(I_test_lambda(2)-1.,1.E-9);

    mat b_rand = simcoon::v2t_strain(randu(6))+eye(3,3);

    double J = sqrt(det(b_rand));
    
    mat b_bar = pow(J,-2./3.)*b_rand;
    vec I = zeros(3);    
    I(0) = trace(b_bar);
    I(1) = 0.5*(pow(trace(b_bar),2.)-trace(powmat(b_bar,2)));
    I(2) = 1.;

    vec I_rand = simcoon::isochoric_invariants(b_rand);     

    EXPECT_LT(I_rand(0)-I(0),1.E-9);
    EXPECT_LT(I_rand(1)-I(1),1.E-9);    
    EXPECT_LT(I_rand(2)-I(2),1.E-9);

    lambda = eig_sym(b_rand);
    lambda.transform( [](double val) { return (sqrt(val)); } );

    vec I_rand_lambda = simcoon::isochoric_invariants(lambda);

    EXPECT_LT(I_rand_lambda(0)-I(0),1.E-9);
    EXPECT_LT(I_rand_lambda(1)-I(1),1.E-9);    
    EXPECT_LT(I_rand_lambda(2)-I(2),1.E-9);

}

TEST(Thyperelastic, isochoric_pstretch)
{

    mat V_test = eye(3,3);
    mat b_test = powmat(V_test,2);    

    vec lambda_bar_test_from_V = simcoon::isochoric_pstretch_from_V(V_test);

    EXPECT_LT(lambda_bar_test_from_V(0)-3.,1.E-9);
    EXPECT_LT(lambda_bar_test_from_V(1)-3.,1.E-9);    
    EXPECT_LT(lambda_bar_test_from_V(2)-1.,1.E-9);

    vec lambda_bar_test_from_b = simcoon::isochoric_pstretch_from_b(b_test);

    EXPECT_LT(lambda_bar_test_from_b(0)-1.,1.E-9);
    EXPECT_LT(lambda_bar_test_from_b(1)-1.,1.E-9);    
    EXPECT_LT(lambda_bar_test_from_b(2)-1.,1.E-9);

    mat V_rand = simcoon::v2t_strain(randu(6)) + eye(3,3);
    mat b_rand = powmat(V_rand,2);    
    double J = det(V_rand);    

    vec lambda = eig_sym(V_rand);
    vec lambda_bar = pow(J,-1./3.)*lambda;

    lambda_bar_test_from_V = simcoon::isochoric_pstretch_from_V(V_rand);    

    cout << lambda_bar_test_from_V << endl;

    EXPECT_LT(lambda_bar_test_from_V(0)-lambda_bar(0),1.E-9);
    EXPECT_LT(lambda_bar_test_from_V(1)-lambda_bar(1),1.E-9);    
    EXPECT_LT(lambda_bar_test_from_V(2)-lambda_bar(2),1.E-9);

    lambda_bar_test_from_b = simcoon::isochoric_pstretch_from_b(b_rand);

    EXPECT_LT(lambda_bar_test_from_b(0)-lambda_bar(0),1.E-9);
    EXPECT_LT(lambda_bar_test_from_b(1)-lambda_bar(1),1.E-9);
    EXPECT_LT(lambda_bar_test_from_b(2)-lambda_bar(2),1.E-9);

}

TEST(Thyperelastic, pstretch_from_b)
{
    // Create a known deformation: simple stretch
    mat F = diagmat(vec({1.5, 2.0, 0.8}));
    mat b = F * F.t();

    vec lambda;
    mat n_pvector;
    simcoon::pstretch(lambda, n_pvector, b, "b");

    // Eigenvalues should match the stretches (sorted ascending)
    vec expected = sort(vec({1.5, 2.0, 0.8}));
    EXPECT_LT(norm(sort(lambda) - expected, 2), 1.E-9);

    // Eigenvectors should be orthogonal
    EXPECT_LT(norm(n_pvector.t() * n_pvector - eye(3, 3), 2), 1.E-9);
}

TEST(Thyperelastic, pstretch_with_projectors)
{
    mat F = diagmat(vec({1.5, 2.0, 0.8}));
    mat b = F * F.t();

    vec lambda;
    mat n_pvector;
    vector<mat> N_projectors;
    simcoon::pstretch(lambda, n_pvector, N_projectors, b, "b");

    // Should have 3 projectors
    EXPECT_EQ(N_projectors.size(), (size_t)3);

    // Each projector N_i = n_i * n_i^T, and sum(N_i) = I
    mat sum_proj = zeros(3, 3);
    for (int i = 0; i < 3; i++) {
        sum_proj += N_projectors[i];
        // N_i should be symmetric
        EXPECT_LT(norm(N_projectors[i] - N_projectors[i].t(), 2), 1.E-9);
    }
    EXPECT_LT(norm(sum_proj - eye(3, 3), 2), 1.E-9);
}

TEST(Thyperelastic, isochoric_pstretch_with_directions)
{
    mat F = diagmat(vec({1.5, 2.0, 0.8}));
    mat b = F * F.t();

    vec lambda_bar;
    mat n_pvector;
    simcoon::isochoric_pstretch(lambda_bar, n_pvector, b, "b");

    // The product of isochoric stretches should be 1 (det=1 constraint)
    double J = det(F);
    vec expected_lambda = sort(vec({1.5, 2.0, 0.8})) / pow(J, 1. / 3.);
    double product = lambda_bar(0) * lambda_bar(1) * lambda_bar(2);
    EXPECT_LT(fabs(product - 1.0), 1.E-9);
}

TEST(Thyperelastic, beta_gamma_coefs)
{
    // Neo-Hookean: W = mu/2 * (I1_bar - 3), dW/dlambda_bar_i = mu * lambda_bar_i
    double mu = 1.0;
    vec lambda_bar = {1.0, 1.0, 1.0};
    vec dWdlambda_bar = mu * lambda_bar;

    vec beta = simcoon::beta_coefs(dWdlambda_bar, lambda_bar);
    EXPECT_EQ(beta.n_elem, (arma::uword)3);

    // gamma coefs
    mat dW2dlambda_bar2 = mu * eye(3, 3);
    mat gamma = simcoon::gamma_coefs(dWdlambda_bar, dW2dlambda_bar2, lambda_bar);
    EXPECT_EQ(gamma.n_rows, (arma::uword)3);
    EXPECT_EQ(gamma.n_cols, (arma::uword)3);
}

TEST(Thyperelastic, tau_iso_hyper_pstretch_identity)
{
    // For identity deformation (b = I), neo-Hookean stress should be zero
    double mu = 100.;
    mat b = eye(3, 3);
    vec dWdlambda_bar = mu * ones(3); // neo-Hookean

    mat tau = simcoon::tau_iso_hyper_pstretch(dWdlambda_bar, b);
    // For b=I, isochoric contribution vanishes
    // tau_iso should have zero trace or be deviatoric
    double trace_tau = tau(0, 0) + tau(1, 1) + tau(2, 2);
    EXPECT_LT(fabs(trace_tau), 1.E-6);
}

TEST(Thyperelastic, tau_vol_hyper)
{
    double dUdJ = 100.; // pressure-like term
    mat b = eye(3, 3);

    mat tau_vol = simcoon::tau_vol_hyper(dUdJ, b);

    // Volumetric Kirchoff stress should be hydrostatic: tau_vol = dU/dJ * J * I
    // For b = I, J = 1
    EXPECT_LT(fabs(tau_vol(0, 0) - dUdJ), 1.E-9);
    EXPECT_LT(fabs(tau_vol(1, 1) - dUdJ), 1.E-9);
    EXPECT_LT(fabs(tau_vol(2, 2) - dUdJ), 1.E-9);
    EXPECT_LT(fabs(tau_vol(0, 1)), 1.E-9);
}

TEST(Thyperelastic, sigma_vol_hyper)
{
    double dUdJ = 100.;
    mat b = eye(3, 3);

    mat sigma_vol = simcoon::sigma_vol_hyper(dUdJ, b);

    // For b=I (J=1), sigma_vol = tau_vol / J = dU/dJ * I
    EXPECT_LT(fabs(sigma_vol(0, 0) - dUdJ), 1.E-9);
    EXPECT_LT(fabs(sigma_vol(1, 1) - dUdJ), 1.E-9);
    EXPECT_LT(fabs(sigma_vol(2, 2) - dUdJ), 1.E-9);
}

TEST(Thyperelastic, L_vol_hyper_symmetry)
{
    double dUdJ = 100.;
    double dU2dJ2 = 50.;
    mat b = eye(3, 3);

    mat Lvol = simcoon::L_vol_hyper(dUdJ, dU2dJ2, b);

    // Tangent modulus should have major symmetry
    EXPECT_LT(norm(Lvol - Lvol.t(), 2), 1.E-9);
}

TEST(Thyperelastic, a_b_delta_coefs)
{
    // Test invariant-based coefficients
    double dWdI_1_bar = 1.0;
    double dWdI_2_bar = 0.0;
    vec I_bar = {3., 3., 1.};

    vec a = simcoon::a_coefs(dWdI_1_bar, dWdI_2_bar, I_bar);
    EXPECT_EQ(a.n_elem, (arma::uword)3);

    double dW2dI_11_bar = 0.;
    double dW2dI_12_bar = 0.;
    double dW2dI_22_bar = 0.;
    vec b = simcoon::b_coefs(dWdI_2_bar, dW2dI_11_bar, dW2dI_12_bar, dW2dI_22_bar, I_bar);
    EXPECT_EQ(b.n_elem, (arma::uword)6);

    mat b_tensor = eye(3, 3);
    vec delta = simcoon::delta_coefs(a, b, b_tensor);
    EXPECT_EQ(delta.n_elem, (arma::uword)6);
}

TEST(Thyperelastic, L_iso_hyper_pstretch_symmetry)
{
    // Build a non-trivial deformation
    mat F = diagmat(vec({1.5, 2.0, 0.8}));
    mat b = F * F.t();
    double J = sqrt(det(b));

    vec lambda;
    mat n_pvectors;
    simcoon::pstretch(lambda, n_pvectors, b, "b");

    vec lambda_bar = lambda / pow(J, 1. / 3.);
    double mu = 100.;
    vec dWdlambda_bar = mu * lambda_bar;
    mat dW2dlambda_bar2 = mu * eye(3, 3);

    mat L_iso = simcoon::L_iso_hyper_pstretch(dWdlambda_bar, dW2dlambda_bar2, b, J);

    // Tangent modulus should have major symmetry
    EXPECT_LT(norm(L_iso - L_iso.t(), 2), 1.E-6);
}