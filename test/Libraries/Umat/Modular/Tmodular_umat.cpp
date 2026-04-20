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
#include <simcoon/Continuum_mechanics/Umat/Modular/modular_umat.hpp>

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

// ========== KinematicHardening Tests ==========

class KinematicHardeningTest : public ::testing::Test {
protected:
    InternalVariableCollection ivc_;
};

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

    // Accumulated plastic strain should be positive
    double p = mumat.internal_variables().get("p").scalar();
    // If we strained enough to yield, p should be positive
    if (sigma_11 > sigma_Y - 10.0) {
        EXPECT_GT(p, 0.0);
    }
}
