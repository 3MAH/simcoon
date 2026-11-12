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
#include <simcoon/exception.hpp>
#include <simcoon/Continuum_mechanics/Functions/transfer.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Continuum_mechanics/Functions/kinematics.hpp>
#include <simcoon/Continuum_mechanics/Functions/objective_rates.hpp>
#include <simcoon/Continuum_mechanics/Functions/hyperelastic.hpp>
#include <simcoon/Simulation/Maths/eig_sym3.hpp>
#include <simcoon/Continuum_mechanics/Umat/Finite/generic_hyper_invariants.hpp>

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

    double J;
    try {
        J = sqrt(det(b_rand));
    } catch (const std::runtime_error &e) {
        cerr << "Error in det or inv: throw inv exception " << e.what() << endl;
        throw simcoon::exception_det("Error in det function inside Thyperelastic.");
    }      
    
    mat b_bar = pow(J,-2./3.)*b_rand;
    vec I = zeros(3);    

    I(0) = trace(b_bar);
    mat b_bar2;
    try {
        I(1) = 0.5*(pow(trace(b_bar),2.)-trace(powmat(b_bar,2)));
    } catch (const std::runtime_error &e) {
        cerr << "Error in det: " << e.what() << endl;
        throw simcoon::exception_powmat("Error in powmat function inside Thyperelastic.");
    }            
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

    double J;
    try {
        J = det(V_rand);    
    } catch (const std::runtime_error &e) {
        cerr << "Error in det or inv: throw inv exception " << e.what() << endl;
        throw simcoon::exception_det("Error in det function inside Thyperelastic.");
    }      

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
    // N_projectors must be pre-sized to 3 elements
    vector<mat> N_projectors(3);
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

TEST(Thyperelastic, beta_gamma_coefs_values)
{
    // Test with a known stretch state
    vec lambda_bar = {1.2, 0.9, 0.926}; // product ~ 1
    double mu = 100.;
    vec dWdlambda_bar = mu * lambda_bar;

    vec beta = simcoon::beta_coefs(dWdlambda_bar, lambda_bar);
    // beta should have 3 components
    EXPECT_EQ(beta.n_elem, (arma::uword)3);

    mat dW2dlambda_bar2 = mu * eye(3, 3);
    mat gamma = simcoon::gamma_coefs(dWdlambda_bar, dW2dlambda_bar2, lambda_bar);
    // gamma should be 3x3
    EXPECT_EQ(gamma.n_rows, (arma::uword)3);
    EXPECT_EQ(gamma.n_cols, (arma::uword)3);
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
namespace {

// L_iso_hyper_invariants as the full expression of all five gamma terms (before the zero terms
// were skipped)
mat L_iso_hyper_invariants_full(const double &d1, const double &d2, const double &d11, const double &d12, const double &d22, const mat &b, const double &J) {
    const mat b_bar = pow(J,-2./3.)*b;
    const mat b_bar2 = powmat(b_bar,2);
    const mat dev_b_bar = dev(b_bar);
    const mat dev_b_bar2 = dev(b_bar2);
    const vec I_bar = isochoric_invariants(b,J);
    const mat Id = eye(3,3);
    const mat H_bar = auto_sym_dyadic_operator(b_bar);
    const mat devdevbb2 = I_bar(0)*dev(b_bar) - dev(b_bar2);
    const mat gamma_1 = (4./3.)*(I_bar(0)*Idev() - (sym_dyadic(dev_b_bar,Id)+sym_dyadic(Id,dev_b_bar)));
    const mat gamma_2 = (8./3.)*(I_bar(1)*(Ireal() - 2.*Ivol()) - I_bar(0)*(sym_dyadic(dev_b_bar,Id)+sym_dyadic(Id,dev_b_bar))
                    + (sym_dyadic(dev_b_bar2,Id)+sym_dyadic(Id,dev_b_bar2))) + 4*(auto_sym_dyadic(b_bar)-H_bar);
    const mat gamma_11 = 4.*auto_sym_dyadic(dev_b_bar);
    const mat gamma_22 = 4.*auto_sym_dyadic(devdevbb2);
    const mat gamma_12 = 4.*(sym_dyadic(dev_b_bar,devdevbb2)+sym_dyadic(devdevbb2,dev_b_bar));
    return (1./J)*(gamma_1*d1+gamma_2*d2+gamma_11*d11+gamma_12*d12+gamma_22*d22);
}

}  // namespace

TEST(Thyperelastic, L_iso_hyper_invariants_skips_only_zero_terms)
{
    arma_rng::set_seed(23);
    const mat F = eye(3,3) + 0.2*randn<mat>(3,3);
    const mat b = F*F.t();
    const double J = det(F);
    const std::vector<std::vector<double>> coefs = {
        {0.5, 0., 0., 0., 0.},           // Neo-Hookean: gamma_1 alone
        {0.26, -0.045, 0., 0., 0.},       // Mooney-Rivlin
        {0.4, 0., 0.03, 0., 0.},          // Yeoh-like: a second derivative in I_1
        {0.4, 0.1, 0.03, 0.02, 0.01}};    // every term
    for (const auto &c : coefs) {
        const mat L = L_iso_hyper_invariants(c[0], c[1], c[2], c[3], c[4], b, J);
        const mat L_ref = L_iso_hyper_invariants_full(c[0], c[1], c[2], c[3], c[4], b, J);
        EXPECT_LT(abs(L - L_ref).max(), 1.E-13*abs(L_ref).max());   // the reference keeps the difference form of the I_2 terms
    }
}

TEST(Thyperelastic, I2_terms_exact_at_large_stretch)
{
    // Uniaxial stretch lambda along e1 at constant volume: b_bar = diag(l^2, 1/l, 1/l), I_2 = 2 l + 1/l^2,
    // tr(b_bar) dev(b_bar) - dev(b_bar^2) = -dev(b_bar^-1) = -dev(diag(1/l^2, l, l)). The difference form
    // cancels two O(l^4) terms into O(l) and keeps 16 - 3 log10(l) digits; the cofactor form keeps 16.
    for (const double l : {10., 100., 1000.}) {
        const mat b = diagmat(vec({l*l, 1./l, 1./l}));
        const vec I = isochoric_invariants(b, 1.);
        EXPECT_LT(fabs(I(1) - (2.*l + 1./(l*l))), 1.E-15*I(1)) << l;
        const mat tau = tau_iso_hyper_invariants(0., 1., b, 1.);
        const mat inv_b = diagmat(vec({1./(l*l), l, l}));
        const mat tau_ref = -2.*dev(inv_b);
        EXPECT_LT(abs(tau - tau_ref).max(), 1.E-15*abs(tau_ref).max()) << l;
        // gamma_22 = 4 dd (x) dd with dd the same vector: its (1,1) entry is 4 (dd_11)^2, dd_11 = -(1/l^2 - I_2/3)
        const mat L = L_iso_hyper_invariants(0., 0., 0., 0., 1., b, 1.);
        const double dd_11 = -(1./(l*l) - I(1)/3.);
        EXPECT_LT(fabs(L(0,0) - 4.*dd_11*dd_11), 1.E-14*fabs(L(0,0))) << l;
    }
}

TEST(Thyperelastic, pure_stretch_box_tangent_closed_form)
{
    // Dtau_LieDD_2_DtauDe_stretch (Hencky block of the modular engine) against the pull-back +
    // exact spectral map of corate 3 at F = V, on an isotropic potential
    arma_rng::set_seed(29);
    const vec props = {0.5, 0.02, 0.003, 500.};   // YEOHH
    for (int t = 0; t < 10; t++) {
        mat e = 0.3*randn<mat>(3,3);
        e = 0.5*(e + e.t());
        if (t == 0) e = diagmat(vec({0.1, 0.1, -0.05}));   // repeated eigenvalue
        vec::fixed<3> ev;
        mat::fixed<3,3> N;
        ASSERT_TRUE(eig_sym_3x3(ev, N, e));
        const mat V = spectral_compose_3x3(N, exp(vec(ev)));
        const mat b = V*V;
        const double J = det(V);
        const hyper_invariants_dW dW = hyper_potential_derivatives(HyperPotential::YEOHH, props, isochoric_invariants(b, J), J);
        vec tau_v;
        mat Lt_ref, Lt_Lie;
        hyper_invariants_response(dW, b, J, V, 3, tau_v, Lt_ref, {});
        hyper_invariants_response(dW, b, J, V, 4, tau_v, Lt_Lie, {});
        const mat Lt = Dtau_LieDD_2_DtauDe_stretch(Lt_Lie, ev, N, v2t_stress(tau_v));
        EXPECT_LT(abs(Lt - Lt_ref).max(), 1.E-12*abs(Lt_ref).max());
    }
}

TEST(Thyperelastic, material_tangent_form)
{
    // corate_type_material_tangent: a kernel built on F returns dS/dE, the pull-back of its box
    // tangent in any corate, with the same stress and work
    arma_rng::set_seed(31);
    const mat F0 = eye(3,3);
    const mat F1 = eye(3,3) + 0.15*randn<mat>(3,3);
    const vec props = {0.26, -0.045, 10000.};   // MOORI
    const vec etot = zeros(6), Detot = t2v_strain(Log_strain(F1));
    for (int corate : {0, 1, 2, 3, 4, 5, corate_type_total}) {
        vec sig_box = zeros(6), sig_mat = zeros(6), sv_box = zeros(1), sv_mat = zeros(1);
        mat Lt_box(6,6), Lt_mat(6,6), L(6,6);
        double Wm_b = 0., Wr_b = 0., Wir_b = 0., Wd_b = 0., Wm_m = 0., Wr_m = 0., Wir_m = 0., Wd_m = 0., tnew = 1.;
        umat_generic_hyper_invariants("MOORI", etot, Detot, F0, F1, sig_box, Lt_box, L, eye(3,3), 3, props, 1, sv_box, 293., 0., 0., 1.,
                                      Wm_b, Wr_b, Wir_b, Wd_b, 3, 3, false, tnew, corate, tangent_default);
        umat_generic_hyper_invariants("MOORI", etot, Detot, F0, F1, sig_mat, Lt_mat, L, eye(3,3), 3, props, 1, sv_mat, 293., 0., 0., 1.,
                                      Wm_m, Wr_m, Wir_m, Wd_m, 3, 3, false, tnew, corate_type_material_tangent, tangent_default);
        EXPECT_TRUE(approx_equal(sig_box, sig_mat, "absdiff", 0.));
        EXPECT_EQ(Wm_b, Wm_m);
        const mat DSDE_ref = DtauDe_corate_2_DSDE(Lt_box, corate, F1, v2t_stress(sig_box));
        EXPECT_LT(abs(Lt_mat - DSDE_ref).max(), 1.E-12*abs(DSDE_ref).max());
        // and the dispatchers are the identity on that form
        EXPECT_TRUE(approx_equal(DtauDe_corate_2_DSDE(Lt_mat, corate_type_material_tangent, F1, v2t_stress(sig_mat)), Lt_mat, "absdiff", 0.));
        EXPECT_TRUE(approx_equal(DSDE_2_DtauDe_corate(Lt_mat, corate_type_material_tangent, F1, v2t_stress(sig_mat)), Lt_mat, "absdiff", 0.));
    }
}
