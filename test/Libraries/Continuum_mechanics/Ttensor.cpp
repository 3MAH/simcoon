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

///@file Ttensor.cpp
///@brief Tests for tensor2 and tensor4 classes
///@version 1.0

#include <gtest/gtest.h>
#include <armadillo>
#include <cmath>

#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Functions/tensor.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Continuum_mechanics/Functions/transfer.hpp>
#include <simcoon/Simulation/Maths/rotation.hpp>

using namespace std;
using namespace arma;
using namespace simcoon;

// ============================================================================
// tensor2 constructors
// ============================================================================

TEST(Ttensor2, DefaultConstructor)
{
    tensor2 t;
    EXPECT_TRUE(arma::approx_equal(t.mat(), mat::fixed<3,3>(fill::zeros), "absdiff", 1e-15));
    EXPECT_EQ(t.vtype(), VoigtType::stress);
}

TEST(Ttensor2, VoigtTypeConstructor)
{
    tensor2 t(VoigtType::strain);
    EXPECT_TRUE(arma::approx_equal(t.mat(), mat::fixed<3,3>(fill::zeros), "absdiff", 1e-15));
    EXPECT_EQ(t.vtype(), VoigtType::strain);
}

TEST(Ttensor2, MatrixConstructor)
{
    mat::fixed<3,3> m = {{1, 2, 3}, {2, 5, 6}, {3, 6, 9}};
    tensor2 t(m, VoigtType::stress);
    EXPECT_TRUE(arma::approx_equal(t.mat(), m, "absdiff", 1e-15));
}

TEST(Ttensor2, DynamicMatConstructor)
{
    mat m = {{1, 2, 3}, {2, 5, 6}, {3, 6, 9}};
    tensor2 t(m, VoigtType::stress);
    mat::fixed<3,3> expected = {{1, 2, 3}, {2, 5, 6}, {3, 6, 9}};
    EXPECT_TRUE(arma::approx_equal(t.mat(), expected, "absdiff", 1e-15));
}

TEST(Ttensor2, Factories)
{
    tensor2 z = tensor2::zeros(VoigtType::strain);
    EXPECT_TRUE(arma::approx_equal(z.mat(), mat::fixed<3,3>(fill::zeros), "absdiff", 1e-15));
    EXPECT_EQ(z.vtype(), VoigtType::strain);

    tensor2 id = tensor2::identity(VoigtType::stress);
    EXPECT_TRUE(arma::approx_equal(id.mat(), mat::fixed<3,3>(fill::eye), "absdiff", 1e-15));
}

// ============================================================================
// tensor2 Voigt conventions — cross-validate with transfer.hpp
// ============================================================================

TEST(Ttensor2, StressVoigtMatchesTransfer)
{
    // Create a symmetric stress tensor
    mat::fixed<3,3> sigma = {{100, 30, 20}, {30, 200, 40}, {20, 40, 300}};
    tensor2 t(sigma, VoigtType::stress);

    vec::fixed<6> v = t.voigt();
    vec v_ref = t2v_stress(mat(sigma));

    EXPECT_LT(norm(vec(v) - v_ref, 2), 1e-12);
}

TEST(Ttensor2, StrainVoigtMatchesTransfer)
{
    // Create a symmetric strain tensor
    mat::fixed<3,3> eps = {{0.01, 0.005, 0.003}, {0.005, 0.02, 0.004}, {0.003, 0.004, 0.03}};
    tensor2 t(eps, VoigtType::strain);

    vec::fixed<6> v = t.voigt();
    vec v_ref = t2v_strain(mat(eps));

    EXPECT_LT(norm(vec(v) - v_ref, 2), 1e-12);
}

TEST(Ttensor2, FromVoigtStress)
{
    vec::fixed<6> v = {100, 200, 300, 30, 20, 40};
    tensor2 t = tensor2::from_voigt(v, VoigtType::stress);

    // Should reconstruct the matrix
    mat m_ref = v2t_stress(vec(v));
    EXPECT_LT(norm(mat(t.mat()) - m_ref, "fro"), 1e-12);

    // And the voigt should roundtrip
    EXPECT_LT(norm(vec(t.voigt()) - vec(v), 2), 1e-12);
}

TEST(Ttensor2, FromVoigtStrain)
{
    vec::fixed<6> v = {0.01, 0.02, 0.03, 0.01, 0.006, 0.008};
    tensor2 t = tensor2::from_voigt(v, VoigtType::strain);

    // Should reconstruct the matrix (shear halved in matrix form)
    mat m_ref = v2t_strain(vec(v));
    EXPECT_LT(norm(mat(t.mat()) - m_ref, "fro"), 1e-12);

    // And the voigt should roundtrip
    EXPECT_LT(norm(vec(t.voigt()) - vec(v), 2), 1e-12);
}

TEST(Ttensor2, VoigtThrowsForNone)
{
    mat::fixed<3,3> F = {{1.1, 0.1, 0}, {0, 1, 0.05}, {0, 0, 0.95}};
    tensor2 t(F, VoigtType::none);
    EXPECT_THROW(t.voigt(), std::runtime_error);
}

// ============================================================================
// tensor2 set_voigt / set_mat
// ============================================================================

TEST(Ttensor2, SetVoigt)
{
    tensor2 t(VoigtType::stress);
    vec::fixed<6> v = {100, 200, 300, 30, 20, 40};
    t.set_voigt(v);
    EXPECT_LT(norm(vec(t.voigt()) - vec(v), 2), 1e-12);
}

TEST(Ttensor2, SetMat)
{
    tensor2 t(VoigtType::stress);
    mat::fixed<3,3> m = {{10, 3, 2}, {3, 20, 4}, {2, 4, 30}};
    t.set_mat(m);
    EXPECT_TRUE(arma::approx_equal(t.mat(), m, "absdiff", 1e-15));
}

// ============================================================================
// tensor2 symmetry check
// ============================================================================

TEST(Ttensor2, SymmetryCheck)
{
    mat::fixed<3,3> sym = {{1, 2, 3}, {2, 5, 6}, {3, 6, 9}};
    tensor2 ts(sym, VoigtType::stress);
    EXPECT_TRUE(ts.is_symmetric());

    mat::fixed<3,3> nonsym = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    tensor2 tn(nonsym, VoigtType::none);
    EXPECT_FALSE(tn.is_symmetric());
}

// ============================================================================
// tensor2 Fastor zero-copy map
// ============================================================================

TEST(Ttensor2, FastorMap)
{
    mat::fixed<3,3> m = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    tensor2 t(m, VoigtType::none);

    auto fmap = t.fastor();
    // Check a few elements (Fastor and armadillo are both column-major)
    EXPECT_DOUBLE_EQ(fmap(0,0), 1.0);
    EXPECT_DOUBLE_EQ(fmap(1,1), 5.0);
    EXPECT_DOUBLE_EQ(fmap(2,2), 9.0);
}

// ============================================================================
// tensor2 free functions
// ============================================================================

TEST(Ttensor2, StressStrainFreeFunctions)
{
    mat::fixed<3,3> m = {{100, 30, 20}, {30, 200, 40}, {20, 40, 300}};
    tensor2 s = stress(m);
    EXPECT_EQ(s.vtype(), VoigtType::stress);

    tensor2 e = strain(m);
    EXPECT_EQ(e.vtype(), VoigtType::strain);
}

TEST(Ttensor2, Trace)
{
    mat::fixed<3,3> m = {{100, 30, 20}, {30, 200, 40}, {20, 40, 300}};
    tensor2 t(m, VoigtType::stress);
    EXPECT_DOUBLE_EQ(trace(t), 600.0);
}

TEST(Ttensor2, DevStress)
{
    // Hydrostatic stress => deviatoric part is zero
    mat::fixed<3,3> hydro = {{100, 0, 0}, {0, 100, 0}, {0, 0, 100}};
    tensor2 t(hydro, VoigtType::stress);
    vec::fixed<6> d = dev(t);
    EXPECT_LT(norm(d, 2), 1e-12);
}

TEST(Ttensor2, Mises)
{
    // Pure tension sigma_11 = Y => von Mises = Y
    double Y = 250.0;
    mat::fixed<3,3> tension = fill::zeros;
    tension(0,0) = Y;
    tensor2 t(tension, VoigtType::stress);
    EXPECT_NEAR(Mises(t), Y, 1e-10);
}

// ============================================================================
// tensor2 arithmetic
// ============================================================================

TEST(Ttensor2, Addition)
{
    mat::fixed<3,3> a = {{1, 2, 3}, {2, 5, 6}, {3, 6, 9}};
    mat::fixed<3,3> b = {{10, 20, 30}, {20, 50, 60}, {30, 60, 90}};
    tensor2 ta(a, VoigtType::stress);
    tensor2 tb(b, VoigtType::stress);
    tensor2 tc = ta + tb;
    mat::fixed<3,3> expected;
    expected = a + b;
    EXPECT_TRUE(arma::approx_equal(tc.mat(), expected, "absdiff", 1e-14));
}

TEST(Ttensor2, ScalarMultiply)
{
    mat::fixed<3,3> a = {{1, 2, 3}, {2, 5, 6}, {3, 6, 9}};
    tensor2 ta(a, VoigtType::stress);
    tensor2 tb = ta * 3.0;
    tensor2 tc = 3.0 * ta;
    mat::fixed<3,3> expected;
    expected = a * 3.0;
    EXPECT_TRUE(arma::approx_equal(tb.mat(), expected, "absdiff", 1e-14));
    EXPECT_TRUE(arma::approx_equal(tc.mat(), expected, "absdiff", 1e-14));
}

TEST(Ttensor2, Equality)
{
    mat::fixed<3,3> a = {{1, 2, 3}, {2, 5, 6}, {3, 6, 9}};
    tensor2 ta(a, VoigtType::stress);
    tensor2 tb(a, VoigtType::stress);
    tensor2 tc(a, VoigtType::strain);
    EXPECT_TRUE(ta == tb);
    EXPECT_FALSE(ta == tc);
}

// ============================================================================
// tensor2 rotation
// ============================================================================

TEST(Ttensor2, RotationStressRoundtrip)
{
    // Rotate a stress tensor by some angle and back: should recover original
    mat::fixed<3,3> sigma = {{100, 30, 20}, {30, 200, 40}, {20, 40, 300}};
    tensor2 t(sigma, VoigtType::stress);

    Rotation R = Rotation::from_euler(0.3, 0.5, 0.7, "zxz");
    tensor2 t_rot = t.rotate(R, true);
    Rotation R_inv = Rotation::from_euler(-0.7, -0.5, -0.3, "zxz");
    tensor2 t_back = t_rot.rotate(R_inv, true);

    EXPECT_LT(norm(mat(t.mat()) - mat(t_back.mat()), "fro"), 1e-9);
}

TEST(Ttensor2, RotationStrainRoundtrip)
{
    mat::fixed<3,3> eps = {{0.01, 0.005, 0.003}, {0.005, 0.02, 0.004}, {0.003, 0.004, 0.03}};
    tensor2 t(eps, VoigtType::strain);

    Rotation R = Rotation::from_euler(0.3, 0.5, 0.7, "zxz");
    tensor2 t_rot = t.rotate(R, true);
    Rotation R_inv = Rotation::from_euler(-0.7, -0.5, -0.3, "zxz");
    tensor2 t_back = t_rot.rotate(R_inv, true);

    EXPECT_LT(norm(mat(t.mat()) - mat(t_back.mat()), "fro"), 1e-9);
}

// ============================================================================
// tensor2 push-forward / pull-back
// ============================================================================

TEST(Ttensor2, PushPullBackStress)
{
    // Push-forward then pull-back should be identity
    mat::fixed<3,3> sigma = {{100, 30, 20}, {30, 200, 40}, {20, 40, 300}};
    tensor2 t(sigma, VoigtType::stress);

    mat::fixed<3,3> F = {{1.1, 0.1, 0.05}, {0.02, 0.95, 0.03}, {0.01, 0.04, 1.05}};

    tensor2 t_push = t.push_forward(F);
    tensor2 t_back = t_push.pull_back(F);

    EXPECT_LT(norm(mat(t.mat()) - mat(t_back.mat()), "fro"), 1e-9);
}

TEST(Ttensor2, PushPullBackStrain)
{
    mat::fixed<3,3> eps = {{0.01, 0.005, 0.003}, {0.005, 0.02, 0.004}, {0.003, 0.004, 0.03}};
    tensor2 t(eps, VoigtType::strain);

    mat::fixed<3,3> F = {{1.1, 0.1, 0.05}, {0.02, 0.95, 0.03}, {0.01, 0.04, 1.05}};

    tensor2 t_push = t.push_forward(F);
    tensor2 t_back = t_push.pull_back(F);

    EXPECT_LT(norm(mat(t.mat()) - mat(t_back.mat()), "fro"), 1e-9);
}

// ============================================================================
// tensor4 constructors
// ============================================================================

TEST(Ttensor4, DefaultConstructor)
{
    tensor4 t;
    EXPECT_TRUE(arma::approx_equal(t.mat(), mat::fixed<6,6>(fill::zeros), "absdiff", 1e-15));
    EXPECT_EQ(t.type(), Tensor4Type::stiffness);
}

TEST(Ttensor4, MatrixConstructor)
{
    mat::fixed<6,6> L = L_iso(70000., 0.3, "Enu");
    tensor4 t(L, Tensor4Type::stiffness);
    EXPECT_TRUE(arma::approx_equal(t.mat(), L, "absdiff", 1e-15));
}

// ============================================================================
// tensor4 static factories cross-validate with constitutive.hpp
// ============================================================================

TEST(Ttensor4, IdentityMatchesIreal)
{
    tensor4 t = tensor4::identity();
    mat Iref = Ireal();
    EXPECT_LT(norm(mat(t.mat()) - Iref, "fro"), 1e-12);
}

TEST(Ttensor4, VolumetricMatchesIvol)
{
    tensor4 t = tensor4::volumetric();
    mat Iref = Ivol();
    EXPECT_LT(norm(mat(t.mat()) - Iref, "fro"), 1e-12);
}

TEST(Ttensor4, DeviatoricMatchesIdev)
{
    tensor4 t = tensor4::deviatoric();
    mat Iref = Idev();
    EXPECT_LT(norm(mat(t.mat()) - Iref, "fro"), 1e-12);
}

TEST(Ttensor4, Identity2MatchesIreal2)
{
    tensor4 t = tensor4::identity2();
    mat Iref = Ireal2();
    EXPECT_LT(norm(mat(t.mat()) - Iref, "fro"), 1e-12);
}

TEST(Ttensor4, Deviatoric2MatchesIdev2)
{
    tensor4 t = tensor4::deviatoric2();
    mat Iref = Idev2();
    EXPECT_LT(norm(mat(t.mat()) - Iref, "fro"), 1e-12);
}

// ============================================================================
// tensor4 contraction with type inference
// ============================================================================

TEST(Ttensor4, ContractStiffnessGivesStress)
{
    double E = 70000.0;
    double nu = 0.3;
    mat::fixed<6,6> L = L_iso(E, nu, "Enu");
    tensor4 stiffness(L, Tensor4Type::stiffness);

    // Uniaxial strain eps_11 = 0.01
    vec::fixed<6> eps_v = {0.01, 0.0, 0.0, 0.0, 0.0, 0.0};
    tensor2 eps = tensor2::from_voigt(eps_v, VoigtType::strain);

    tensor2 sigma = stiffness.contract(eps);
    EXPECT_EQ(sigma.vtype(), VoigtType::stress);

    // Cross-validate with raw matrix-vector product
    vec sig_ref = L * eps_v;
    EXPECT_LT(norm(vec(sigma.voigt()) - sig_ref, 2), 1e-10);
}

TEST(Ttensor4, ContractComplianceGivesStrain)
{
    double E = 70000.0;
    double nu = 0.3;
    mat::fixed<6,6> M = M_iso(E, nu, "Enu");
    tensor4 compliance(M, Tensor4Type::compliance);

    vec::fixed<6> sig_v = {100, 0, 0, 0, 0, 0};
    tensor2 sig = tensor2::from_voigt(sig_v, VoigtType::stress);

    tensor2 eps = compliance.contract(sig);
    EXPECT_EQ(eps.vtype(), VoigtType::strain);

    vec eps_ref = M * sig_v;
    EXPECT_LT(norm(vec(eps.voigt()) - eps_ref, 2), 1e-10);
}

TEST(Ttensor4, StiffnessComplianceInverse)
{
    double E = 70000.0;
    double nu = 0.3;
    mat::fixed<6,6> L = L_iso(E, nu, "Enu");
    mat::fixed<6,6> M = M_iso(E, nu, "Enu");

    tensor4 stiff(L, Tensor4Type::stiffness);
    tensor4 comp(M, Tensor4Type::compliance);

    // L * M should be the 6x6 identity
    mat LM = L * M;
    EXPECT_LT(norm(LM - eye(6,6), "fro"), 1e-9);
}

// ============================================================================
// tensor4 Fastor cache
// ============================================================================

TEST(Ttensor4, FastorCacheLazyCompute)
{
    mat::fixed<6,6> L = L_iso(70000., 0.3, "Enu");
    tensor4 t(L, Tensor4Type::stiffness);

    // First access computes the Fastor tensor
    const auto& f = t.fastor();
    // Check a diagonal element: C_0000 = L(0,0)
    EXPECT_NEAR(f(0,0,0,0), L(0,0), 1e-12);
    // Check a shear element: C_0101 = L(3,3) for stiffness convention
    EXPECT_NEAR(f(0,1,0,1), L(3,3), 1e-12);
}

// ============================================================================
// tensor4 arithmetic
// ============================================================================

TEST(Ttensor4, IsotropicFromProjectors)
{
    // L_iso = 3K * Ivol + 2mu * Idev
    double E = 70000.0;
    double nu = 0.3;
    double K = E / (3.0 * (1.0 - 2.0 * nu));
    double mu = E / (2.0 * (1.0 + nu));

    tensor4 L_from_proj = 3.0 * K * tensor4::volumetric() + 2.0 * mu * tensor4::deviatoric();

    mat::fixed<6,6> L_ref = L_iso(E, nu, "Enu");
    EXPECT_LT(norm(mat(L_from_proj.mat()) - mat(L_ref), "fro"), 1e-8);
}

TEST(Ttensor4, ScalarMultiply)
{
    mat::fixed<6,6> L = L_iso(70000., 0.3, "Enu");
    tensor4 t(L, Tensor4Type::stiffness);
    tensor4 t2 = t * 2.0;
    tensor4 t3 = 2.0 * t;
    mat::fixed<6,6> expected;
    expected = L * 2.0;
    EXPECT_TRUE(arma::approx_equal(t2.mat(), expected, "absdiff", 1e-12));
    EXPECT_TRUE(arma::approx_equal(t3.mat(), expected, "absdiff", 1e-12));
}

// ============================================================================
// tensor4 rotation
// ============================================================================

TEST(Ttensor4, RotationStiffnessRoundtrip)
{
    mat::fixed<6,6> L = L_iso(70000., 0.3, "Enu");
    tensor4 t(L, Tensor4Type::stiffness);

    // For isotropic material, rotation should give the same stiffness
    Rotation R = Rotation::from_euler(0.3, 0.5, 0.7, "zxz");
    tensor4 t_rot = t.rotate(R, true);

    EXPECT_LT(norm(mat(t.mat()) - mat(t_rot.mat()), "fro"), 1e-8);
}

TEST(Ttensor4, RotationComplianceRoundtrip)
{
    mat::fixed<6,6> M = M_iso(70000., 0.3, "Enu");
    tensor4 t(M, Tensor4Type::compliance);

    // Isotropic compliance should be invariant under rotation
    Rotation R = Rotation::from_euler(0.4, 0.6, 0.8, "zxz");
    tensor4 t_rot = t.rotate(R, true);

    EXPECT_LT(norm(mat(t.mat()) - mat(t_rot.mat()), "fro"), 1e-8);
}

// ============================================================================
// tensor4 push-forward / pull-back
// ============================================================================

TEST(Ttensor4, PushPullBackRoundtrip)
{
    mat::fixed<6,6> L = L_iso(70000., 0.3, "Enu");
    tensor4 t(L, Tensor4Type::stiffness);

    mat::fixed<3,3> F = {{1.1, 0.1, 0.05}, {0.02, 0.95, 0.03}, {0.01, 0.04, 1.05}};

    tensor4 t_push = t.push_forward(F);
    tensor4 t_back = t_push.pull_back(F);

    EXPECT_LT(norm(mat(t.mat()) - mat(t_back.mat()), "fro"), 1e-8);
}

// ============================================================================
// tensor4 dyadic product
// ============================================================================

TEST(Ttensor4, Dyadic)
{
    mat::fixed<3,3> a = {{1, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    mat::fixed<3,3> b = {{0, 0, 0}, {0, 1, 0}, {0, 0, 0}};
    tensor2 ta(a, VoigtType::stress);
    tensor2 tb(b, VoigtType::stress);

    tensor4 d = dyadic(ta, tb);
    // a has Voigt [1,0,0,0,0,0], b has Voigt [0,1,0,0,0,0]
    // dyadic -> 6x6 matrix with d(0,1) = 1, rest zero
    EXPECT_NEAR(d.mat()(0,1), 1.0, 1e-14);
    EXPECT_NEAR(d.mat()(0,0), 0.0, 1e-14);
    EXPECT_NEAR(d.mat()(1,1), 0.0, 1e-14);
}

TEST(Ttensor4, AutoDyadic)
{
    mat::fixed<3,3> hydro = fill::eye;
    tensor2 t(hydro, VoigtType::stress);

    tensor4 ad = auto_dyadic(t);
    // Identity tensor has Voigt [1,1,1,0,0,0]
    // auto_dyadic -> outer product with itself
    EXPECT_NEAR(ad.mat()(0,0), 1.0, 1e-14);
    EXPECT_NEAR(ad.mat()(0,1), 1.0, 1e-14);
    EXPECT_NEAR(ad.mat()(1,0), 1.0, 1e-14);
    EXPECT_NEAR(ad.mat()(3,3), 0.0, 1e-14);
}

// ============================================================================
// tensor4 equality
// ============================================================================

TEST(Ttensor4, Equality)
{
    mat::fixed<6,6> L = L_iso(70000., 0.3, "Enu");
    tensor4 a(L, Tensor4Type::stiffness);
    tensor4 b(L, Tensor4Type::stiffness);
    tensor4 c(L, Tensor4Type::compliance);
    EXPECT_TRUE(a == b);
    EXPECT_FALSE(a == c);
}

// ============================================================================
// Integration: L : epsilon = sigma workflow
// ============================================================================

TEST(Ttensor_integration, StiffnessContractStrain)
{
    double E = 70000.0;
    double nu = 0.3;

    // Build stiffness tensor
    tensor4 L(L_iso(E, nu, "Enu"), Tensor4Type::stiffness);

    // A general strain state
    vec::fixed<6> eps_v = {0.01, -0.003, -0.003, 0.005, 0.002, 0.001};
    tensor2 eps = tensor2::from_voigt(eps_v, VoigtType::strain);

    // Contract: sigma = L : epsilon
    tensor2 sigma = L.contract(eps);

    // Verify against raw arma computation
    vec sig_ref = L_iso(E, nu, "Enu") * eps_v;
    EXPECT_LT(norm(vec(sigma.voigt()) - sig_ref, 2), 1e-10);
    EXPECT_EQ(sigma.vtype(), VoigtType::stress);
}

TEST(Ttensor_integration, ComplianceContractStress)
{
    double E = 70000.0;
    double nu = 0.3;

    // Build compliance tensor
    tensor4 M(M_iso(E, nu, "Enu"), Tensor4Type::compliance);

    // A general stress state
    vec::fixed<6> sig_v = {100, 50, 75, 30, 20, 10};
    tensor2 sig = tensor2::from_voigt(sig_v, VoigtType::stress);

    // Contract: epsilon = M : sigma
    tensor2 eps = M.contract(sig);

    // Verify against raw arma computation
    vec eps_ref = M_iso(E, nu, "Enu") * sig_v;
    EXPECT_LT(norm(vec(eps.voigt()) - eps_ref, 2), 1e-10);
    EXPECT_EQ(eps.vtype(), VoigtType::strain);
}

TEST(Ttensor_integration, FullCycle_Stress_Rotate_Contract)
{
    double E = 70000.0;
    double nu = 0.3;

    // Build stiffness and strain
    tensor4 L(L_iso(E, nu, "Enu"), Tensor4Type::stiffness);
    vec::fixed<6> eps_v = {0.01, -0.003, -0.003, 0.005, 0.002, 0.001};
    tensor2 eps = tensor2::from_voigt(eps_v, VoigtType::strain);

    // Compute stress in material frame
    tensor2 sigma = L.contract(eps);

    // Rotate stiffness and strain to global frame
    Rotation R = Rotation::from_euler(0.3, 0.5, 0.7, "zxz");
    tensor4 L_rot = L.rotate(R, true);
    tensor2 eps_rot = eps.rotate(R, true);

    // Compute stress in global frame
    tensor2 sigma_rot_direct = L_rot.contract(eps_rot);

    // Rotate material-frame stress to global frame
    tensor2 sigma_rot_indirect = sigma.rotate(R, true);

    // Both should give the same result
    EXPECT_LT(norm(mat(sigma_rot_direct.mat()) - mat(sigma_rot_indirect.mat()), "fro"), 1e-8);
}
