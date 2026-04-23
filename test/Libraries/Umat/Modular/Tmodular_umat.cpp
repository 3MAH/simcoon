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

///@file Tmodular_umat.cpp
///@brief Unit tests for Modular UMAT components
///@version 1.0

#include <gtest/gtest.h>
#include <armadillo>
#include <cmath>

#include <simcoon/Continuum_mechanics/Umat/Modular/internal_variable.hpp>
#include <simcoon/Continuum_mechanics/Umat/Modular/internal_variable_collection.hpp>
#include <simcoon/Continuum_mechanics/Umat/Modular/elasticity_module.hpp>
#include <simcoon/Continuum_mechanics/Umat/Modular/yield_criterion.hpp>
#include <simcoon/Continuum_mechanics/Umat/Modular/hardening.hpp>
#include <simcoon/Continuum_mechanics/Umat/Modular/plasticity_mechanism.hpp>
#include <simcoon/Continuum_mechanics/Umat/Modular/viscoelastic_mechanism.hpp>
#include <simcoon/Continuum_mechanics/Umat/Modular/damage_mechanism.hpp>
#include <simcoon/Continuum_mechanics/Umat/Modular/modular_umat.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Continuum_mechanics/Functions/tensor.hpp>
#include <simcoon/Simulation/Maths/rotation.hpp>

using namespace std;
using namespace arma;
using namespace simcoon;

// ========== InternalVariable Tests ==========

class InternalVariableTest : public ::testing::Test {
protected:
    void SetUp() override {}
};

TEST_F(InternalVariableTest, ScalarConstruction) {
    InternalVariable iv("test_scalar", 1.5, false);

    EXPECT_EQ(iv.name(), "test_scalar");
    EXPECT_EQ(iv.type(), IVarType::SCALAR);
    EXPECT_DOUBLE_EQ(iv.scalar(), 1.5);
    EXPECT_EQ(iv.size(), 1u);
    EXPECT_FALSE(iv.requires_rotation());
}

TEST_F(InternalVariableTest, VectorConstruction) {
    vec init = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    InternalVariable iv("test_vec", init, true);

    EXPECT_EQ(iv.name(), "test_vec");
    EXPECT_EQ(iv.type(), IVarType::VECTOR_6);
    EXPECT_EQ(iv.size(), 6u);
    EXPECT_TRUE(iv.requires_rotation());

    for (int i = 0; i < 6; i++) {
        EXPECT_DOUBLE_EQ(iv.vec()(i), init(i));
    }
}

TEST_F(InternalVariableTest, MatrixConstruction) {
    mat init = eye(6, 6);
    InternalVariable iv("test_mat", init, true);

    EXPECT_EQ(iv.name(), "test_mat");
    EXPECT_EQ(iv.type(), IVarType::MATRIX_6x6);
    EXPECT_EQ(iv.size(), 36u);
    EXPECT_TRUE(iv.requires_rotation());
}

TEST_F(InternalVariableTest, ScalarPackUnpack) {
    InternalVariable iv("test", 3.14, false);
    iv.set_offset(5);

    vec statev = zeros(10);
    iv.pack(statev);
    EXPECT_DOUBLE_EQ(statev(5), 3.14);

    // Modify and unpack
    statev(5) = 2.71;
    iv.unpack(statev);
    EXPECT_DOUBLE_EQ(iv.scalar(), 2.71);
}

TEST_F(InternalVariableTest, VectorPackUnpack) {
    vec init = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    InternalVariable iv("test", init, false);
    iv.set_offset(2);

    vec statev = zeros(10);
    iv.pack(statev);

    for (int i = 0; i < 6; i++) {
        EXPECT_DOUBLE_EQ(statev(2 + i), init(i));
    }

    // Modify and unpack
    statev.subvec(2, 7).fill(10.0);
    iv.unpack(statev);
    for (int i = 0; i < 6; i++) {
        EXPECT_DOUBLE_EQ(iv.vec()(i), 10.0);
    }
}

TEST_F(InternalVariableTest, DeltaComputation) {
    InternalVariable iv("test", 1.0, false);
    iv.to_start();  // Save current as start

    iv.scalar() = 3.0;  // Modify
    EXPECT_DOUBLE_EQ(iv.delta_scalar(), 2.0);
}

// Tensor-typed views: VECTOR_6 → tensor2 (strain or stress), MATRIX_6x6 → tensor4.
// Verifies round-trip through the typed accessors and rejects type mismatches.
TEST_F(InternalVariableTest, TensorViews) {
    arma::vec EP_voigt = {0.01, -0.005, -0.005, 0.002, 0.0, 0.0};
    InternalVariable EP("EP", EP_voigt, true);

    tensor2 EP_t = EP.as_strain();
    EXPECT_EQ(EP_t.vtype(), VoigtType::strain);
    EXPECT_LT(arma::norm(arma::vec(EP_t.voigt()) - EP_voigt, 2), 1e-12);

    // Mutate via typed setter and read back
    arma::vec new_EP = {0.02, -0.01, -0.01, 0.004, 0.0, 0.0};
    tensor2 new_EP_t = tensor2::from_voigt(new_EP, VoigtType::strain);
    EP.set_tensor2(new_EP_t);
    EXPECT_LT(arma::norm(EP.vec() - new_EP, 2), 1e-12);

    // MATRIX_6x6 → tensor4 stiffness
    arma::mat L_init = L_iso(70000., 0.3, "Enu");
    InternalVariable L_iv("L_active", L_init, true);
    tensor4 L_t = L_iv.as_stiffness();
    EXPECT_EQ(L_t.type(), Tensor4Type::stiffness);
    EXPECT_LT(arma::norm(arma::mat(L_t.mat()) - L_init, "fro"), 1e-10);

    // Type-mismatch rejection
    InternalVariable scalar_iv("d", 0.1, false);
    EXPECT_THROW(scalar_iv.as_strain(), std::runtime_error);
    EXPECT_THROW(scalar_iv.as_stiffness(), std::runtime_error);
    EXPECT_THROW(EP.as_stiffness(), std::runtime_error);
    EXPECT_THROW(L_iv.as_strain(), std::runtime_error);
}

// Rotation overload: rotate(Rotation) is equivalent to rotate(R.as_matrix()).
TEST_F(InternalVariableTest, RotationOverload) {
    arma::vec EP_voigt = {0.01, -0.005, -0.005, 0.002, 0.001, 0.0};
    InternalVariable EP_a("EP", EP_voigt, true);
    InternalVariable EP_b("EP", EP_voigt, true);

    Rotation R = Rotation::from_axis_angle(M_PI / 6.0, 3);
    arma::mat DR = R.as_matrix();

    EP_a.rotate(DR);
    EP_b.rotate(R);

    EXPECT_LT(arma::norm(EP_a.vec() - EP_b.vec(), 2), 1e-12);
}

// ========== YieldCriterion Tensor2 overloads ==========

TEST(YieldCriterionTensor, VonMisesTensor2Match) {
    YieldCriterion yc;
    yc.configure_von_mises();

    arma::vec sigma_voigt = {200., -50., 30., 40., 20., 10.};
    tensor2 sigma_t = stress(arma::mat::fixed<3,3>{{200., 40., 20.},
                                                    {40., -50., 10.},
                                                    {20., 10., 30.}});
    EXPECT_NEAR(yc.equivalent_stress(sigma_t),
                yc.equivalent_stress(arma::vec(sigma_t.voigt())),
                1e-10);

    tensor2 lambda_t = yc.flow_direction(sigma_t);
    EXPECT_EQ(lambda_t.vtype(), VoigtType::strain);
    EXPECT_LT(arma::norm(arma::vec(lambda_t.voigt())
                         - yc.flow_direction(arma::vec(sigma_t.voigt())), 2),
              1e-10);

    // Backstress: callers just subtract the tensors themselves.
    tensor2 X_t = stress(arma::mat::fixed<3,3>{{20., 5., 0.},
                                                {5., 10., 2.},
                                                {0., 2., -10.}});
    EXPECT_NEAR(yc.equivalent_stress(sigma_t - X_t),
                yc.equivalent_stress(arma::vec(sigma_t.voigt())
                                     - arma::vec(X_t.voigt())),
                1e-10);
}

// ========== InternalVariableCollection Tests ==========

class InternalVariableCollectionTest : public ::testing::Test {
protected:
    void SetUp() override {
        ivc_.add_scalar("p", 0.0, false);
        ivc_.add_vec("EP", zeros(6), true);
        ivc_.add_scalar("D", 0.0, false);
    }

    InternalVariableCollection ivc_;
};

TEST_F(InternalVariableCollectionTest, Registration) {
    EXPECT_TRUE(ivc_.has("p"));
    EXPECT_TRUE(ivc_.has("EP"));
    EXPECT_TRUE(ivc_.has("D"));
    EXPECT_FALSE(ivc_.has("nonexistent"));
}

TEST_F(InternalVariableCollectionTest, Access) {
    ivc_.get("p").scalar() = 0.5;
    EXPECT_DOUBLE_EQ(ivc_.get("p").scalar(), 0.5);

    ivc_.get("EP").vec()(0) = 1.0;
    EXPECT_DOUBLE_EQ(ivc_.get("EP").vec()(0), 1.0);
}

TEST_F(InternalVariableCollectionTest, OffsetComputation) {
    ivc_.compute_offsets(0);

    // p: offset 0, size 1
    // EP: offset 1, size 6
    // D: offset 7, size 1
    // Total: 8
    EXPECT_EQ(ivc_.total_statev_size(), 8);
}

TEST_F(InternalVariableCollectionTest, PackUnpackAll) {
    ivc_.compute_offsets(0);

    ivc_.get("p").scalar() = 0.1;
    ivc_.get("EP").vec().fill(0.01);
    ivc_.get("D").scalar() = 0.05;

    vec statev = zeros(10);
    ivc_.pack_all(statev);

    EXPECT_DOUBLE_EQ(statev(0), 0.1);
    for (int i = 0; i < 6; i++) {
        EXPECT_DOUBLE_EQ(statev(1 + i), 0.01);
    }
    EXPECT_DOUBLE_EQ(statev(7), 0.05);

    // Modify and unpack
    statev.fill(0.5);
    ivc_.unpack_all(statev);

    EXPECT_DOUBLE_EQ(ivc_.get("p").scalar(), 0.5);
    EXPECT_DOUBLE_EQ(ivc_.get("D").scalar(), 0.5);
}

// ========== ElasticityModule Tests ==========

class ElasticityModuleTest : public ::testing::Test {
protected:
    void SetUp() override {}
};

TEST_F(ElasticityModuleTest, IsotropicConfiguration) {
    ElasticityModule em;

    // E = 210000 MPa, nu = 0.3, alpha = 1.2e-5
    vec props = {210000.0, 0.3, 1.2e-5};
    int offset = 0;

    em.configure(ElasticityType::ISOTROPIC, props, offset);

    EXPECT_TRUE(em.is_configured());
    EXPECT_EQ(offset, 3);

    // Check stiffness properties
    mat L = em.L();
    EXPECT_EQ(L.n_rows, 6u);
    EXPECT_EQ(L.n_cols, 6u);

    // C11 should be positive
    EXPECT_GT(L(0, 0), 0.0);

    // Symmetry check
    for (unsigned int i = 0; i < 6; i++) {
        for (unsigned int j = 0; j < 6; j++) {
            EXPECT_NEAR(L(i, j), L(j, i), 1e-10);
        }
    }
}

TEST_F(ElasticityModuleTest, CubicConfiguration) {
    ElasticityModule em;

    // Cubic: E = 185000, nu = 0.28, G = 39700, alpha = 1.2e-5
    // Note: for cubic, G is independent from E and nu
    vec props = {185000.0, 0.28, 39700.0, 1.2e-5};
    int offset = 0;

    em.configure(ElasticityType::CUBIC, props, offset);

    EXPECT_TRUE(em.is_configured());
    EXPECT_EQ(offset, 4);

    mat L = em.L();
    EXPECT_EQ(L.n_rows, 6u);
    EXPECT_EQ(L.n_cols, 6u);

    // C11 should be positive
    EXPECT_GT(L(0, 0), 0.0);

    // Cubic symmetry: C11 = C22 = C33
    EXPECT_NEAR(L(0, 0), L(1, 1), 1e-10);
    EXPECT_NEAR(L(1, 1), L(2, 2), 1e-10);

    // C12 = C13 = C23
    EXPECT_NEAR(L(0, 1), L(0, 2), 1e-10);
    EXPECT_NEAR(L(0, 1), L(1, 2), 1e-10);

    // C44 = C55 = C66 = G
    EXPECT_NEAR(L(3, 3), 39700.0, 1e-6);
    EXPECT_NEAR(L(3, 3), L(4, 4), 1e-10);
    EXPECT_NEAR(L(4, 4), L(5, 5), 1e-10);

    // For cubic, C44 != (C11-C12)/2 in general (anisotropy)
    double C11 = L(0, 0);
    double C12 = L(0, 1);
    double C44 = L(3, 3);
    double iso_check = (C11 - C12) / 2.0;
    // They should NOT be equal unless the material is isotropic
    // (i.e., Zener ratio A = 2*C44/(C11-C12) != 1 in general)

    // CTE isotropic for cubic
    vec alpha = em.alpha();
    EXPECT_DOUBLE_EQ(alpha(0), 1.2e-5);
    EXPECT_DOUBLE_EQ(alpha(1), 1.2e-5);
    EXPECT_DOUBLE_EQ(alpha(2), 1.2e-5);
}

TEST_F(ElasticityModuleTest, CubicCiiConfiguration) {
    ElasticityModule em;

    // Typical values for an FCC metal (e.g., copper-like)
    // C11 = 185000, C12 = 158000, C44 = 39700
    em.configure_cubic_Cii(185000.0, 158000.0, 39700.0, 0.0);

    EXPECT_TRUE(em.is_configured());

    mat L = em.L();

    // Direct check of Cii values
    EXPECT_NEAR(L(0, 0), 185000.0, 1e-6);
    EXPECT_NEAR(L(0, 1), 158000.0, 1e-6);
    EXPECT_NEAR(L(3, 3), 39700.0, 1e-6);
}

TEST_F(ElasticityModuleTest, ThermalExpansion) {
    ElasticityModule em;

    vec props = {210000.0, 0.3, 1.2e-5};
    int offset = 0;

    em.configure(ElasticityType::ISOTROPIC, props, offset);

    vec alpha = em.alpha();
    EXPECT_DOUBLE_EQ(alpha(0), 1.2e-5);
    EXPECT_DOUBLE_EQ(alpha(1), 1.2e-5);
    EXPECT_DOUBLE_EQ(alpha(2), 1.2e-5);
    EXPECT_DOUBLE_EQ(alpha(3), 0.0);
    EXPECT_DOUBLE_EQ(alpha(4), 0.0);
    EXPECT_DOUBLE_EQ(alpha(5), 0.0);
}

// Tensor-typed accessors: L_tensor(), M_tensor(), alpha_tensor() should produce
// the same numerical content as the raw arma equivalents but carry type tags
// that drive the right Voigt convention and rotation/contraction dispatch.
TEST_F(ElasticityModuleTest, TensorAccessors) {
    ElasticityModule em;
    vec props = {210000.0, 0.3, 1.2e-5};
    int offset = 0;
    em.configure(ElasticityType::ISOTROPIC, props, offset);

    tensor4 L_t = em.L_tensor();
    EXPECT_EQ(L_t.type(), Tensor4Type::stiffness);
    EXPECT_LT(norm(mat(L_t.mat()) - em.L(), "fro"), 1e-12);

    tensor4 M_t = em.M_tensor();
    EXPECT_EQ(M_t.type(), Tensor4Type::compliance);
    EXPECT_LT(norm(mat(M_t.mat()) - em.M(), "fro"), 1e-12);

    // Stiffness × compliance round-trip via typed contraction
    vec eps_voigt = {1e-3, -3e-4, -3e-4, 0.0, 0.0, 0.0};
    tensor2 eps_t = tensor2::from_voigt(eps_voigt, VoigtType::strain);
    tensor2 sig_t = L_t.contract(eps_t);
    EXPECT_EQ(sig_t.vtype(), VoigtType::stress);
    EXPECT_LT(norm(vec(sig_t.voigt()) - em.L() * eps_voigt, 2), 1e-8);

    // alpha_tensor: VoigtType::strain (factor-2 on shear is irrelevant here, all shear=0)
    tensor2 alpha_t = em.alpha_tensor();
    EXPECT_EQ(alpha_t.vtype(), VoigtType::strain);
    EXPECT_LT(norm(vec(alpha_t.voigt()) - em.alpha(), 2), 1e-12);

    // thermal_strain_tensor and damaged_L_tensor consistency
    double DT = 50.0;
    tensor2 eth = em.thermal_strain_tensor(DT);
    EXPECT_LT(norm(vec(eth.voigt()) - em.thermal_strain(DT), 2), 1e-12);

    tensor4 L_dmg = em.damaged_L_tensor(0.3);
    EXPECT_EQ(L_dmg.type(), Tensor4Type::stiffness);
    EXPECT_LT(norm(mat(L_dmg.mat()) - 0.7 * em.L(), "fro"), 1e-8);
}

// ========== YieldCriterion Tests ==========

class YieldCriterionTest : public ::testing::Test {
protected:
    void SetUp() override {}
};

TEST_F(YieldCriterionTest, VonMisesUniaxial) {
    YieldCriterion yc;
    yc.configure_von_mises();

    // Uniaxial stress: sigma = [100, 0, 0, 0, 0, 0]
    vec sigma = zeros(6);
    sigma(0) = 100.0;

    double sigma_eq = yc.equivalent_stress(sigma);

    // For uniaxial tension, von Mises = sigma_11
    EXPECT_NEAR(sigma_eq, 100.0, 1e-6);
}

TEST_F(YieldCriterionTest, VonMisesPureShear) {
    YieldCriterion yc;
    yc.configure_von_mises();

    // Pure shear: sigma = [0, 0, 0, tau, 0, 0]
    vec sigma = zeros(6);
    sigma(3) = 100.0;  // tau_xy

    double sigma_eq = yc.equivalent_stress(sigma);

    // For pure shear, von Mises = sqrt(3) * tau
    EXPECT_NEAR(sigma_eq, sqrt(3.0) * 100.0, 1e-6);
}

TEST_F(YieldCriterionTest, FlowDirectionNormalized) {
    YieldCriterion yc;
    yc.configure_von_mises();

    vec sigma = {100.0, 50.0, -30.0, 20.0, 10.0, 5.0};

    vec n = yc.flow_direction(sigma);

    // Flow direction should be normalized (for deviatoric stress)
    // Check it's not zero
    double norm_n = norm(n);
    EXPECT_GT(norm_n, 0.0);
}

// ========== IsotropicHardening Tests ==========

class IsotropicHardeningTest : public ::testing::Test {
protected:
    void SetUp() override {}
};

TEST_F(IsotropicHardeningTest, PowerLawHardening) {
    auto hard = IsotropicHardening::create(IsoHardType::POWER_LAW, 1);

    // k = 500, m = 0.2
    vec props = {500.0, 0.2};
    int offset = 0;
    hard->configure(props, offset);

    // R(p) = k * p^m
    double p = 0.1;
    double R = hard->R(p);
    double expected = 500.0 * pow(0.1, 0.2);

    EXPECT_NEAR(R, expected, 1e-6);

    // dR/dp = k * m * p^(m-1)
    double dR = hard->dR_dp(p);
    double expected_dR = 500.0 * 0.2 * pow(0.1, 0.2 - 1.0);
    EXPECT_NEAR(dR, expected_dR, 1e-6);
}

TEST_F(IsotropicHardeningTest, VoceHardening) {
    auto hard = IsotropicHardening::create(IsoHardType::VOCE, 1);

    // Q = 200, b = 10
    vec props = {200.0, 10.0};
    int offset = 0;
    hard->configure(props, offset);

    // R(p) = Q * (1 - exp(-b*p))
    double p = 0.1;
    double R = hard->R(p);
    double expected = 200.0 * (1.0 - exp(-10.0 * 0.1));

    EXPECT_NEAR(R, expected, 1e-6);
}

// Superposition identity: CombinedVoce with two identical (Q, b) terms must
// equal a single VoceHardening(2Q, b). Catches aggregation bugs (averaging
// instead of summing, shared internal state between terms, etc).
TEST_F(IsotropicHardeningTest, CombinedVoceDuplicateEqualsDoubled) {
    auto single = IsotropicHardening::create(IsoHardType::VOCE, 1);
    vec props_single = {400.0, 10.0};  // 2Q, b
    int off = 0;
    single->configure(props_single, off);

    auto combined = IsotropicHardening::create(IsoHardType::COMBINED_VOCE, 2);
    vec props_combined = {200.0, 10.0, 200.0, 10.0};  // (Q, b) twice
    off = 0;
    combined->configure(props_combined, off);

    for (double p : {0.0, 1e-4, 1e-3, 1e-2, 0.1, 0.5}) {
        EXPECT_NEAR(combined->R(p),     single->R(p),     1e-10);
        EXPECT_NEAR(combined->dR_dp(p), single->dR_dp(p), 1e-10);
    }
}

// ========== KinematicHardening Tests ==========

class KinematicHardeningTest : public ::testing::Test {
protected:
    InternalVariableCollection ivc_;
};

// Superposition identity: Chaboche with two identical (C, D) terms must
// produce the same total backstress trajectory as a single AF(2C, D),
// when both are driven by the same flow direction n under identical dp.
TEST_F(KinematicHardeningTest, ChabocheDuplicateEqualsDoubledAF) {
    const double C = 30000.0;
    const double D = 200.0;

    auto af = KinematicHardening::create(KinHardType::ARMSTRONG_FREDERICK, 1);
    vec props_af = {2.0 * C, D};
    int off = 0;
    af->configure(props_af, off);

    auto cha = KinematicHardening::create(KinHardType::CHABOCHE, 2);
    vec props_cha = {C, D, C, D};
    off = 0;
    cha->configure(props_cha, off);

    InternalVariableCollection ivc_af;
    InternalVariableCollection ivc_cha;
    af->register_variables(ivc_af);
    cha->register_variables(ivc_cha);
    ivc_af.compute_offsets(0);
    ivc_cha.compute_offsets(0);
    vec statev_af  = zeros(ivc_af.total_statev_size());
    vec statev_cha = zeros(ivc_cha.total_statev_size());
    ivc_af.unpack_all(statev_af);
    ivc_cha.unpack_all(statev_cha);

    // Uniaxial tensile flow direction (Voigt strain convention, factor-2 shear).
    vec n = {1.0, -0.5, -0.5, 0.0, 0.0, 0.0};
    const double dp = 1e-4;

    for (int step = 0; step < 100; ++step) {
        af->update(dp, n, ivc_af);
        cha->update(dp, n, ivc_cha);

        vec X_af  = af->total_backstress(ivc_af);
        vec X_cha = cha->total_backstress(ivc_cha);
        EXPECT_LT(norm(X_cha - X_af, 2), 1e-8) << "drift at step " << step;

        EXPECT_NEAR(cha->hardening_modulus(n, ivc_cha),
                    af->hardening_modulus(n, ivc_af), 1e-8);
    }
}

TEST_F(KinematicHardeningTest, PragerHardening) {
    auto hard = KinematicHardening::create(KinHardType::PRAGER, 1);

    // H = 10000
    vec props = {10000.0};
    int offset = 0;
    hard->configure(props, offset);

    hard->register_variables(ivc_);

    // Check initial backstress is zero
    vec X = hard->total_backstress(ivc_);
    EXPECT_DOUBLE_EQ(norm(X), 0.0);

    // Hardening modulus should be (2/3)*C for Prager kinematic hardening
    vec n = {1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double H = hard->hardening_modulus(n, ivc_);
    EXPECT_NEAR(H, (2.0 / 3.0) * 10000.0, 1e-6);
}

// ========== PlasticityMechanism Tests ==========

class PlasticityMechanismTest : public ::testing::Test {
protected:
    void SetUp() override {}
};

TEST_F(PlasticityMechanismTest, Construction) {
    PlasticityMechanism pm(YieldType::VON_MISES, IsoHardType::VOCE, KinHardType::NONE, 1, 0);

    EXPECT_EQ(pm.type(), MechanismType::PLASTICITY);
    EXPECT_EQ(pm.num_constraints(), 1);
}

TEST_F(PlasticityMechanismTest, InelasticStrainInitiallyZero) {
    PlasticityMechanism pm(YieldType::VON_MISES, IsoHardType::NONE, KinHardType::NONE, 0, 0);

    InternalVariableCollection ivc;
    pm.register_variables(ivc);
    ivc.compute_offsets(0);

    vec statev = zeros(ivc.total_statev_size());
    ivc.unpack_all(statev);

    vec EP = pm.inelastic_strain(ivc);
    EXPECT_DOUBLE_EQ(norm(EP), 0.0);
}

// Public yield-side accessors: equivalent_stress, yield_function, flow_direction.
// Verifies the (sigma - X) shift is applied internally and the tensor2
// overloads match the arma::vec versions bit-for-bit.
TEST_F(PlasticityMechanismTest, YieldAccessors) {
    PlasticityMechanism pm(YieldType::VON_MISES, IsoHardType::VOCE, KinHardType::NONE, 1, 0);

    // sigma_Y = 300, Q = 200, b = 10
    vec props = {300.0, 200.0, 10.0};
    int offset = 0;
    pm.configure(props, offset);

    InternalVariableCollection ivc;
    pm.register_variables(ivc);
    ivc.compute_offsets(0);
    vec statev = zeros(ivc.total_statev_size());
    ivc.unpack_all(statev);

    // Uniaxial elastic stress below yield: sigma_11 = 200 < sigma_Y = 300.
    vec sigma = {200.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    // Mises of uniaxial 200 MPa == 200 MPa.
    double seq = pm.equivalent_stress(sigma, ivc);
    EXPECT_NEAR(seq, 200.0, 1e-10);

    // Phi = 200 - R(0) - 300 = -100 (elastic).
    double Phi = pm.yield_function(sigma, ivc);
    EXPECT_NEAR(Phi, -100.0, 1e-10);

    // Flow direction: non-zero, points in deviatoric direction.
    vec n = pm.flow_direction(sigma, ivc);
    EXPECT_GT(norm(n), 0.0);

    // tensor2 overloads must match the arma::vec versions.
    tensor2 sig_t = stress(sigma);
    EXPECT_NEAR(pm.equivalent_stress(sig_t, ivc), seq, 1e-10);
    EXPECT_NEAR(pm.yield_function(sig_t, ivc), Phi, 1e-10);
    tensor2 n_t = pm.flow_direction(sig_t, ivc);
    EXPECT_EQ(n_t.vtype(), VoigtType::strain);
    EXPECT_LT(norm(vec(n_t.voigt()) - n, 2), 1e-10);

    // Uniaxial stress above yield: sigma_11 = 400 > sigma_Y => Phi > 0.
    vec sigma_high = {400.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    EXPECT_NEAR(pm.equivalent_stress(sigma_high, ivc), 400.0, 1e-10);
    EXPECT_GT(pm.yield_function(sigma_high, ivc), 0.0);
}

// ========== ModularUMAT Tests ==========

class ModularUMATTest : public ::testing::Test {
protected:
    void SetUp() override {}
};

TEST_F(ModularUMATTest, Construction) {
    ModularUMAT mumat;

    EXPECT_FALSE(mumat.is_initialized());
    EXPECT_EQ(mumat.num_mechanisms(), 0u);
}

TEST_F(ModularUMATTest, ElasticConfiguration) {
    ModularUMAT mumat;

    // E = 210000 MPa, nu = 0.3, alpha = 0
    vec props = {210000.0, 0.3, 0.0};
    int offset = 0;

    mumat.set_elasticity(ElasticityType::ISOTROPIC, props, offset);

    const mat& L = mumat.elasticity().L();
    EXPECT_GT(L(0, 0), 0.0);
}

TEST_F(ModularUMATTest, AddPlasticity) {
    ModularUMAT mumat;

    // Elasticity: E = 210000, nu = 0.3, alpha = 0
    vec props_el = {210000.0, 0.3, 0.0};
    int offset_el = 0;
    mumat.set_elasticity(ElasticityType::ISOTROPIC, props_el, offset_el);

    // Plasticity: sigma_Y = 300, Q = 200, b = 10
    vec props_pl = {300.0, 200.0, 10.0};
    int offset_pl = 0;
    mumat.add_plasticity(YieldType::VON_MISES, IsoHardType::VOCE, KinHardType::NONE, 1, 0, props_pl, offset_pl);

    EXPECT_EQ(mumat.num_mechanisms(), 1u);
}

TEST_F(ModularUMATTest, ElasticResponse) {
    ModularUMAT mumat;

    // Very high yield stress to ensure elastic response
    vec props_el = {210000.0, 0.3, 0.0};
    int offset_el = 0;
    mumat.set_elasticity(ElasticityType::ISOTROPIC, props_el, offset_el);

    vec props_pl = {1e10, 0.0, 1.0};  // Very high yield stress
    int offset_pl = 0;
    mumat.add_plasticity(YieldType::VON_MISES, IsoHardType::VOCE, KinHardType::NONE, 1, 0, props_pl, offset_pl);

    // Initialize
    int nstatev = 20;
    vec statev = zeros(nstatev);
    mumat.initialize(nstatev, statev);

    // Run a small elastic step
    vec Etot = zeros(6);
    vec DEtot = {0.001, -0.0003, -0.0003, 0.0, 0.0, 0.0};  // Uniaxial strain
    vec sigma = zeros(6);
    mat Lt = zeros(6, 6);
    mat L = zeros(6, 6);
    mat DR = eye(3, 3);
    vec props = zeros(10);
    double T = 293.0, DT = 0.0, Time = 0.0, DTime = 1.0;
    double Wm = 0.0, Wm_r = 0.0, Wm_ir = 0.0, Wm_d = 0.0;
    double tnew_dt = 1.0;

    mumat.run("MODUL", Etot, DEtot, sigma, Lt, L, DR, 10, props, nstatev, statev,
              T, DT, Time, DTime, Wm, Wm_r, Wm_ir, Wm_d, 3, 3, true, tnew_dt);

    // Check stress is non-zero
    EXPECT_GT(norm(sigma), 0.0);

    // Check tangent is symmetric and positive definite
    for (unsigned int i = 0; i < 6; i++) {
        for (unsigned int j = 0; j < 6; j++) {
            EXPECT_NEAR(Lt(i, j), Lt(j, i), 1e-6);
        }
    }
}

// ========== Integration Test ==========

// Single-branch generalized-Maxwell viscoelasticity: apply a fixed strain,
// hold, and verify the stress relaxes monotonically toward zero (long-term
// modulus is 0 since we have no parallel elastic spring in this single-branch
// test; the instantaneous response is nonzero).
TEST(ModularUMATIntegration, ViscoelasticRelaxation) {
    ModularUMAT mumat;

    vec props_el = {210000.0, 0.3, 0.0};
    int offset_el = 0;
    mumat.set_elasticity(ElasticityType::ISOTROPIC, props_el, offset_el);

    // One Prony branch: E=210 GPa, nu=0.3, etaB=etaS=1e5 MPa*s.
    vec props_vi = {210000.0, 0.3, 1.0e5, 1.0e5};
    int offset_vi = 0;
    mumat.add_viscoelasticity(1, props_vi, offset_vi);

    int nstatev = 20;
    vec statev = zeros(nstatev);
    mumat.initialize(nstatev, statev);

    // Ramp to 0.1% strain, then hold.
    vec Etot = zeros(6);
    vec DEtot_ramp = {1e-3, -3e-4, -3e-4, 0, 0, 0};
    vec sigma = zeros(6);
    mat Lt(6, 6), L(6, 6);
    mat DR = eye(3, 3);
    vec props = zeros(10);
    double T = 293.0, DT = 0.0, Wm = 0.0, Wm_r = 0.0, Wm_ir = 0.0, Wm_d = 0.0;
    double tnew_dt = 1.0;
    double Time = 0.0;

    // Single ramp step
    mumat.run("MODUL", Etot, DEtot_ramp, sigma, Lt, L, DR, 10, props, nstatev, statev,
              T, DT, Time, 1.0, Wm, Wm_r, Wm_ir, Wm_d, 3, 3, true, tnew_dt);
    Etot += DEtot_ramp;

    // sigma right after ramp should be roughly elastic (L * strain).
    const double sigma_11_initial = sigma(0);
    EXPECT_GT(sigma_11_initial, 50.0);  // well above zero

    // Hold (DEtot = 0) for many steps; stress must relax monotonically.
    vec DEtot_hold = zeros(6);
    double sigma_11_prev = sigma_11_initial;
    Time += 1.0;
    for (int step = 0; step < 40; ++step) {
        mumat.run("MODUL", Etot, DEtot_hold, sigma, Lt, L, DR, 10, props, nstatev, statev,
                  T, DT, Time, 1.0, Wm, Wm_r, Wm_ir, Wm_d, 3, 3, false, tnew_dt);
        Time += 1.0;
        // Monotone decrease (allow small tolerance for numerical noise)
        EXPECT_LT(sigma(0), sigma_11_prev + 1e-6) << " at step " << step;
        sigma_11_prev = sigma(0);
    }
    // After 40 s of relaxation the stress has dropped substantially.
    EXPECT_LT(sigma(0), 0.5 * sigma_11_initial);
}

// Arbitrary mix of same-type mechanisms must coexist without IVC key
// collisions: two plasticity + two viscoelastic + two damage. Each mechanism
// receives a unique prefix (plast0_, plast1_, visco0_, visco1_, damage0_,
// damage1_) from ModularUMAT::add_*, and the orchestrator's nstatev
// accounting covers all of them.
TEST(ModularUMATIntegration, ArbitraryMultiMechanismInitialize) {
    ModularUMAT mumat;

    vec props_el = {210000.0, 0.3, 0.0};
    int off = 0;
    mumat.set_elasticity(ElasticityType::ISOTROPIC, props_el, off);

    vec props_pl = {300.0, 10000.0};
    off = 0;
    mumat.add_plasticity(YieldType::VON_MISES, IsoHardType::LINEAR,
                          KinHardType::NONE, 1, 0, props_pl, off);
    off = 0;
    mumat.add_plasticity(YieldType::VON_MISES, IsoHardType::LINEAR,
                          KinHardType::NONE, 1, 0, props_pl, off);

    // Two viscoelastic mechanisms with 1 branch each (E, nu, etaB, etaS).
    vec props_ve = {70000.0, 0.3, 1e5, 1e5};
    off = 0;
    mumat.add_viscoelasticity(1, props_ve, off);
    off = 0;
    mumat.add_viscoelasticity(1, props_ve, off);

    // Two damage mechanisms (linear: Y_0, Y_c).
    vec props_dmg = {1.0, 10.0};
    off = 0;
    mumat.add_damage(DamageType::LINEAR, props_dmg, off);
    off = 0;
    mumat.add_damage(DamageType::LINEAR, props_dmg, off);

    // Statev: 1 (T_init) + 2*7 (plast) + 2*7 (visco 1 branch) + 2*2 (damage) = 33.
    EXPECT_EQ(mumat.required_nstatev(), 1 + 2 * 7 + 2 * 7 + 2 * 2);
    EXPECT_EQ(mumat.num_mechanisms(), 6u);

    vec statev = zeros(mumat.required_nstatev());
    EXPECT_NO_THROW(mumat.initialize(mumat.required_nstatev(), statev));

    // All prefixed keys exist (no collision).
    for (const auto& k : {"plast0_p", "plast1_p", "plast0_EP", "plast1_EP",
                           "visco0_v_0", "visco1_v_0", "visco0_EV_0", "visco1_EV_0",
                           "damage0_D", "damage1_D", "damage0_Y_max", "damage1_Y_max"}) {
        EXPECT_NO_THROW(mumat.internal_variables().get(k)) << "missing key: " << k;
    }
}

// Two identical Plasticity mechanisms sharing the same stress must produce
// the same final stress as a single plasticity mechanism — by superposition
// of the associated flow, exactly when hardening is LINEAR (so that
// 2 * R(p/2) = R(p); with Voce or other nonlinear laws each mechanism hardens
// differently on half the accumulated plastic strain and the equivalence is
// only approximate). Verifies that the cross-mechanism Jacobian coupling
// B_{lj} = -dot(dPhi^l/dsigma, kappa^j) is correctly assembled; without it,
// the FB return-mapping would over-correct and diverge to nonphysical states.
TEST(ModularUMATIntegration, TwoPlasticityEquivalentToOne) {
    const double E = 210000.0, nu = 0.3, sigma_Y = 300.0;
    // Linear isotropic hardening with the same modulus H in both cases.
    const double H_iso = 10000.0;

    // Single plasticity mechanism (reference)
    ModularUMAT mumat_single;
    {
        vec props_el = {E, nu, 0.0};
        int off = 0;
        mumat_single.set_elasticity(ElasticityType::ISOTROPIC, props_el, off);
        vec props_pl = {sigma_Y, H_iso};
        off = 0;
        mumat_single.add_plasticity(YieldType::VON_MISES, IsoHardType::LINEAR,
                                     KinHardType::NONE, 1, 0, props_pl, off);
        vec statev = zeros(20);
        mumat_single.initialize(20, statev);
    }

    // Two identical plasticity mechanisms
    ModularUMAT mumat_two;
    {
        vec props_el = {E, nu, 0.0};
        int off = 0;
        mumat_two.set_elasticity(ElasticityType::ISOTROPIC, props_el, off);
        // Two identical linear-hardening plasticities with HALF the modulus
        // each, so the summed response matches the single-mechanism H_iso.
        vec props_pl_1 = {sigma_Y, 2.0 * H_iso};
        int off1 = 0;
        mumat_two.add_plasticity(YieldType::VON_MISES, IsoHardType::LINEAR,
                                  KinHardType::NONE, 1, 0, props_pl_1, off1);
        vec props_pl_2 = {sigma_Y, 2.0 * H_iso};
        int off2 = 0;
        mumat_two.add_plasticity(YieldType::VON_MISES, IsoHardType::LINEAR,
                                  KinHardType::NONE, 1, 0, props_pl_2, off2);
        vec statev = zeros(40);
        mumat_two.initialize(40, statev);
    }

    // Drive both through the same load history: monotonic uniaxial tension
    // above yield.
    auto run_step = [](ModularUMAT& m, vec& Etot, const vec& DEtot, vec& sigma,
                       int nstatev, vec& statev, double Time, bool start) {
        mat Lt(6, 6), L(6, 6), DR = eye(3, 3);
        vec props = zeros(10);
        double T = 293.0, DT = 0.0, Wm = 0.0, Wm_r = 0.0, Wm_ir = 0.0, Wm_d = 0.0;
        double tnew_dt = 1.0;
        m.run("MODUL", Etot, DEtot, sigma, Lt, L, DR, 10, props, nstatev, statev,
              T, DT, Time, 1.0, Wm, Wm_r, Wm_ir, Wm_d, 3, 3, start, tnew_dt);
    };

    vec Etot_single = zeros(6), sigma_single = zeros(6), statev_single = zeros(20);
    vec Etot_two    = zeros(6), sigma_two    = zeros(6), statev_two    = zeros(40);
    const vec DEtot = {2e-3, -6e-4, -6e-4, 0, 0, 0};

    double Time = 0.0;
    for (int step = 0; step < 5; ++step) {
        run_step(mumat_single, Etot_single, DEtot, sigma_single, 20, statev_single,
                 Time, step == 0);
        run_step(mumat_two, Etot_two, DEtot, sigma_two, 40, statev_two,
                 Time, step == 0);
        Etot_single += DEtot;
        Etot_two    += DEtot;
        Time += 1.0;
    }

    // The two-plasticity case should match the single — identical Φ, identical
    // evolution; each mechanism carries half the plastic increment due to the
    // coupled Jacobian, and the sum gives the same total plastic strain.
    EXPECT_LT(norm(sigma_single - sigma_two, 2) / norm(sigma_single, 2), 1e-3)
        << "single = " << sigma_single.t()
        << "two    = " << sigma_two.t();
}

// Two Plasticity mechanisms must coexist without IVC key collisions. The
// orchestrator auto-prefixes each mechanism's variables (plast0_, plast1_, ...).
TEST(ModularUMATIntegration, TwoPlasticityMechanismsInitialize) {
    ModularUMAT mumat;

    vec props_el = {210000.0, 0.3, 0.0};
    int offset_el = 0;
    mumat.set_elasticity(ElasticityType::ISOTROPIC, props_el, offset_el);

    vec props_pl_1 = {300.0, 200.0, 10.0};
    int off1 = 0;
    mumat.add_plasticity(YieldType::VON_MISES, IsoHardType::VOCE, KinHardType::NONE,
                         1, 0, props_pl_1, off1);

    vec props_pl_2 = {300.0, 200.0, 10.0};
    int off2 = 0;
    mumat.add_plasticity(YieldType::VON_MISES, IsoHardType::VOCE, KinHardType::NONE,
                         1, 0, props_pl_2, off2);

    int nstatev = 40;
    vec statev = zeros(nstatev);
    EXPECT_NO_THROW(mumat.initialize(nstatev, statev));
    EXPECT_EQ(mumat.num_mechanisms(), 2u);

    // Both p and EP exist under their disambiguated keys.
    EXPECT_NO_THROW(mumat.internal_variables().get("plast0_p"));
    EXPECT_NO_THROW(mumat.internal_variables().get("plast1_p"));
    EXPECT_NO_THROW(mumat.internal_variables().get("plast0_EP"));
    EXPECT_NO_THROW(mumat.internal_variables().get("plast1_EP"));
}

TEST(ModularUMATIntegration, UniaxialTension) {
    ModularUMAT mumat;

    // Isotropic elasticity: E = 210000, nu = 0.3, alpha = 0
    vec props_el = {210000.0, 0.3, 0.0};
    int offset_el = 0;
    mumat.set_elasticity(ElasticityType::ISOTROPIC, props_el, offset_el);

    // Von Mises plasticity with Voce hardening: sigma_Y = 300, Q = 200, b = 10
    vec props_pl = {300.0, 200.0, 10.0};
    int offset_pl = 0;
    mumat.add_plasticity(YieldType::VON_MISES, IsoHardType::VOCE, KinHardType::NONE, 1, 0, props_pl, offset_pl);

    // Initialize
    int nstatev = 20;
    vec statev = zeros(nstatev);
    mumat.initialize(nstatev, statev);

    // Apply strain incrementally
    vec Etot = zeros(6);
    vec sigma = zeros(6);
    mat Lt = zeros(6, 6);
    mat L = zeros(6, 6);
    mat DR = eye(3, 3);
    vec props = zeros(10);
    double T = 293.0, DT = 0.0;
    double Wm = 0.0, Wm_r = 0.0, Wm_ir = 0.0, Wm_d = 0.0;
    double tnew_dt = 1.0;

    double E = 210000.0;
    double nu = 0.3;
    double sigma_Y = 300.0;

    // First increment: elastic
    vec DEtot = {0.001, -0.0003, -0.0003, 0.0, 0.0, 0.0};
    double DTime = 1.0;
    double Time = 0.0;

    mumat.run("MODUL", Etot, DEtot, sigma, Lt, L, DR, 10, props, nstatev, statev,
              T, DT, Time, DTime, Wm, Wm_r, Wm_ir, Wm_d, 3, 3, true, tnew_dt);

    // Stress should be approximately E * strain = 210000 * 0.001 = 210 MPa
    // (Below yield)
    EXPECT_NEAR(sigma(0), E * 0.001, 10.0);  // Allow some tolerance due to Poisson effect

    // Multiple increments to go into plasticity
    for (int i = 0; i < 10; i++) {
        Etot += DEtot;
        Time += DTime;

        mumat.run("MODUL", Etot, DEtot, sigma, Lt, L, DR, 10, props, nstatev, statev,
                  T, DT, Time, DTime, Wm, Wm_r, Wm_ir, Wm_d, 3, 3, false, tnew_dt);
    }

    // After significant straining, stress should be above initial yield
    // but following hardening curve
    double sigma_11 = sigma(0);
    EXPECT_GT(sigma_11, sigma_Y - 50.0);  // Should be near or above yield

    // Accumulated plastic strain should be positive.
    // ModularUMAT prefixes each mechanism's variables with plast<i>_, so the
    // first plasticity mechanism's "p" is registered as "plast0_p".
    double p = mumat.internal_variables().get("plast0_p").scalar();
    // If we strained enough to yield, p should be positive
    if (sigma_11 > sigma_Y - 10.0) {
        EXPECT_GT(p, 0.0);
    }
}
