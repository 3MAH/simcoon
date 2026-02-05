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

///@file TRotationClass.cpp
///@brief Test for the Rotation class
///@version 1.0

#include <gtest/gtest.h>
#include <armadillo>

#include <simcoon/parameter.hpp>
#include <simcoon/Simulation/Maths/rotation_class.hpp>
#include <simcoon/Simulation/Maths/rotation.hpp>

using namespace std;
using namespace arma;
using namespace simcoon;

// =============================================================================
// Identity and Basic Tests
// =============================================================================

TEST(TRotationClass, identity)
{
    Rotation r = Rotation::identity();

    EXPECT_TRUE(r.is_identity());

    // Identity quaternion should be [0, 0, 0, 1]
    vec::fixed<4> q = r.as_quat();
    EXPECT_NEAR(q(0), 0.0, 1.E-12);
    EXPECT_NEAR(q(1), 0.0, 1.E-12);
    EXPECT_NEAR(q(2), 0.0, 1.E-12);
    EXPECT_NEAR(q(3), 1.0, 1.E-12);

    // Identity matrix should be eye(3,3)
    mat::fixed<3,3> R = r.as_matrix();
    EXPECT_LT(norm(R - eye(3,3), "fro"), 1.E-12);
}

TEST(TRotationClass, default_constructor)
{
    Rotation r;
    EXPECT_TRUE(r.is_identity());
}

// =============================================================================
// Factory Methods
// =============================================================================

TEST(TRotationClass, from_quat)
{
    // Test with a known quaternion (90 degrees around z-axis)
    double angle = simcoon::pi / 2.0;
    double s = sin(angle / 2.0);
    double c = cos(angle / 2.0);
    vec::fixed<4> q = {0, 0, s, c};

    Rotation r = Rotation::from_quat(q);
    vec::fixed<4> q_out = r.as_quat();

    // Should be normalized and equal
    EXPECT_NEAR(norm(q_out), 1.0, 1.E-12);
    EXPECT_LT(norm(q - q_out), 1.E-12);
}

TEST(TRotationClass, from_matrix)
{
    // Create a rotation matrix using existing fillR function
    double psi = 0.5;
    double theta = 0.3;
    double phi = 0.7;

    mat R_expected = fillR(psi, theta, phi, true, "zxz");

    mat::fixed<3,3> R_fixed;
    R_fixed = R_expected;

    Rotation r = Rotation::from_matrix(R_fixed);
    mat::fixed<3,3> R_actual = r.as_matrix();

    EXPECT_LT(norm(R_expected - R_actual, "fro"), 1.E-9);
}

TEST(TRotationClass, from_euler_zxz)
{
    double psi = 23.0 * (simcoon::pi / 180.0);
    double theta = 42.0 * (simcoon::pi / 180.0);
    double phi = 165.0 * (simcoon::pi / 180.0);

    // Use existing fillR to get expected result
    mat R_expected = fillR(psi, theta, phi, true, "zxz");

    // Create rotation using Rotation class
    Rotation r = Rotation::from_euler(psi, theta, phi, "zxz", true, false);
    mat::fixed<3,3> R_actual = r.as_matrix();

    EXPECT_LT(norm(R_expected - R_actual, "fro"), 1.E-9);
}

TEST(TRotationClass, from_euler_zyz)
{
    double psi = 23.0 * (simcoon::pi / 180.0);
    double theta = 42.0 * (simcoon::pi / 180.0);
    double phi = 165.0 * (simcoon::pi / 180.0);

    mat R_expected = fillR(psi, theta, phi, true, "zyz");

    Rotation r = Rotation::from_euler(psi, theta, phi, "zyz", true, false);
    mat::fixed<3,3> R_actual = r.as_matrix();

    EXPECT_LT(norm(R_expected - R_actual, "fro"), 1.E-9);
}

TEST(TRotationClass, from_euler_degrees)
{
    // Test with degrees flag
    Rotation r1 = Rotation::from_euler(90.0, 0.0, 0.0, "zxz", true, true);
    Rotation r2 = Rotation::from_euler(simcoon::pi/2.0, 0.0, 0.0, "zxz", true, false);

    EXPECT_TRUE(r1.equals(r2, 1.E-9));
}

TEST(TRotationClass, from_axis_angle)
{
    // Test rotation around z-axis
    double angle = simcoon::pi / 4.0;  // 45 degrees

    Rotation r = Rotation::from_axis_angle(angle, 3, false);
    mat::fixed<3,3> R_actual = r.as_matrix();

    mat R_expected = fillR(angle, 3, true);

    EXPECT_LT(norm(R_expected - R_actual, "fro"), 1.E-9);
}

TEST(TRotationClass, from_rotvec)
{
    // Rotation vector: 45 degrees around z-axis
    double angle = simcoon::pi / 4.0;
    vec::fixed<3> rotvec = {0, 0, angle};

    Rotation r = Rotation::from_rotvec(rotvec, false);
    mat::fixed<3,3> R_actual = r.as_matrix();

    mat R_expected = fillR(angle, 3, true);

    EXPECT_LT(norm(R_expected - R_actual, "fro"), 1.E-9);
}

// =============================================================================
// Conversion Methods
// =============================================================================

TEST(TRotationClass, as_euler_roundtrip)
{
    // Test round-trip: euler -> rotation -> euler
    double psi_in = 0.5;
    double theta_in = 0.3;
    double phi_in = 0.7;

    Rotation r = Rotation::from_euler(psi_in, theta_in, phi_in, "zxz", true, false);
    vec::fixed<3> euler_out = r.as_euler("zxz", true, false);

    // Reconstruct rotation from output euler angles
    Rotation r2 = Rotation::from_euler(euler_out(0), euler_out(1), euler_out(2), "zxz", true, false);

    // The two rotations should be equal
    EXPECT_TRUE(r.equals(r2, 1.E-9));
}

TEST(TRotationClass, as_rotvec_roundtrip)
{
    // Create a rotation
    Rotation r = Rotation::from_euler(0.5, 0.3, 0.7, "zxz");

    // Convert to rotvec and back
    vec::fixed<3> rotvec = r.as_rotvec();
    Rotation r2 = Rotation::from_rotvec(rotvec);

    EXPECT_TRUE(r.equals(r2, 1.E-9));
}

TEST(TRotationClass, as_QS_consistency)
{
    // Test that as_QS matches existing fillQS function
    Rotation r = Rotation::from_euler(0.5, 0.3, 0.7, "zxz");
    mat::fixed<3,3> R = r.as_matrix();

    mat::fixed<6,6> QS_class = r.as_QS(true);
    mat QS_func = fillQS(mat(R), true);

    EXPECT_LT(norm(mat(QS_class) - QS_func, "fro"), 1.E-9);
}

TEST(TRotationClass, as_QE_consistency)
{
    // Test that as_QE matches existing fillQE function
    Rotation r = Rotation::from_euler(0.5, 0.3, 0.7, "zxz");
    mat::fixed<3,3> R = r.as_matrix();

    mat::fixed<6,6> QE_class = r.as_QE(true);
    mat QE_func = fillQE(mat(R), true);

    EXPECT_LT(norm(mat(QE_class) - QE_func, "fro"), 1.E-9);
}

// =============================================================================
// Apply Methods
// =============================================================================

TEST(TRotationClass, apply_vector)
{
    // 90 degree rotation around z-axis
    Rotation r = Rotation::from_axis_angle(simcoon::pi / 2.0, 3);

    vec::fixed<3> v = {1, 0, 0};
    vec::fixed<3> v_rot = r.apply(v);

    // x-axis should rotate to y-axis
    EXPECT_NEAR(v_rot(0), 0.0, 1.E-9);
    EXPECT_NEAR(v_rot(1), 1.0, 1.E-9);
    EXPECT_NEAR(v_rot(2), 0.0, 1.E-9);
}

TEST(TRotationClass, apply_vector_inverse)
{
    Rotation r = Rotation::from_axis_angle(simcoon::pi / 2.0, 3);

    vec::fixed<3> v = {1, 0, 0};
    vec::fixed<3> v_rot = r.apply(v, false);
    vec::fixed<3> v_back = r.apply(v_rot, true);

    EXPECT_LT(norm(v - v_back), 1.E-9);
}

TEST(TRotationClass, apply_stress_consistency)
{
    // Test that apply_stress matches existing rotate_stress function
    Rotation r = Rotation::from_euler(0.5, 0.3, 0.7, "zxz");
    mat::fixed<3,3> R = r.as_matrix();

    vec::fixed<6> sigma = {100, 50, 0, 25, 0, 0};

    vec::fixed<6> sigma_class = r.apply_stress(sigma, true);
    vec sigma_func = rotate_stress(vec(sigma), mat(R), true);

    EXPECT_LT(norm(vec(sigma_class) - sigma_func), 1.E-9);
}

TEST(TRotationClass, apply_strain_consistency)
{
    Rotation r = Rotation::from_euler(0.5, 0.3, 0.7, "zxz");
    mat::fixed<3,3> R = r.as_matrix();

    vec::fixed<6> epsilon = {0.01, -0.005, -0.005, 0.002, 0.001, 0.003};

    vec::fixed<6> epsilon_class = r.apply_strain(epsilon, true);
    vec epsilon_func = rotate_strain(vec(epsilon), mat(R), true);

    EXPECT_LT(norm(vec(epsilon_class) - epsilon_func), 1.E-9);
}

TEST(TRotationClass, apply_stiffness_consistency)
{
    Rotation r = Rotation::from_euler(0.5, 0.3, 0.7, "zxz");
    mat::fixed<3,3> R = r.as_matrix();

    // Create a simple stiffness matrix (isotropic)
    mat L = eye(6, 6) * 100.0;
    L(0,1) = L(1,0) = 30.0;
    L(0,2) = L(2,0) = 30.0;
    L(1,2) = L(2,1) = 30.0;

    mat::fixed<6,6> L_fixed;
    L_fixed = L;

    mat::fixed<6,6> L_class = r.apply_stiffness(L_fixed, true);
    mat L_func = rotateL(L, mat(R), true);

    EXPECT_LT(norm(mat(L_class) - L_func, "fro"), 1.E-9);
}

// =============================================================================
// Operations
// =============================================================================

TEST(TRotationClass, inverse)
{
    Rotation r = Rotation::from_euler(0.5, 0.3, 0.7, "zxz");
    Rotation r_inv = r.inv();

    // r * r_inv should be identity
    Rotation product = r * r_inv;
    EXPECT_TRUE(product.is_identity(1.E-9));
}

TEST(TRotationClass, composition)
{
    // Test that (r1 * r2).apply(v) == r1.apply(r2.apply(v))
    Rotation r1 = Rotation::from_axis_angle(0.5, 1);
    Rotation r2 = Rotation::from_axis_angle(0.3, 2);

    vec::fixed<3> v = {1, 2, 3};

    // Method 1: compose then apply
    Rotation r_composed = r1 * r2;
    vec::fixed<3> v1 = r_composed.apply(v);

    // Method 2: apply sequentially (r2 first, then r1)
    vec::fixed<3> v2 = r1.apply(r2.apply(v));

    EXPECT_LT(norm(v1 - v2), 1.E-9);
}

TEST(TRotationClass, composition_inplace)
{
    Rotation r1 = Rotation::from_axis_angle(0.5, 1);
    Rotation r2 = Rotation::from_axis_angle(0.3, 2);

    Rotation r_composed = r1 * r2;

    Rotation r1_copy = r1;
    r1_copy *= r2;

    EXPECT_TRUE(r_composed.equals(r1_copy, 1.E-9));
}

TEST(TRotationClass, magnitude)
{
    // 45 degree rotation
    double angle = simcoon::pi / 4.0;
    Rotation r = Rotation::from_axis_angle(angle, 3);

    EXPECT_NEAR(r.magnitude(), angle, 1.E-9);
    EXPECT_NEAR(r.magnitude(true), 45.0, 1.E-9);  // degrees
}

TEST(TRotationClass, slerp)
{
    Rotation r1 = Rotation::identity();
    Rotation r2 = Rotation::from_axis_angle(simcoon::pi / 2.0, 3);

    // t=0 should give r1
    Rotation r_t0 = r1.slerp(r2, 0.0);
    EXPECT_TRUE(r_t0.equals(r1, 1.E-9));

    // t=1 should give r2
    Rotation r_t1 = r1.slerp(r2, 1.0);
    EXPECT_TRUE(r_t1.equals(r2, 1.E-9));

    // t=0.5 should give 45 degrees rotation
    Rotation r_t5 = r1.slerp(r2, 0.5);
    EXPECT_NEAR(r_t5.magnitude(), simcoon::pi / 4.0, 1.E-9);
}

TEST(TRotationClass, equals)
{
    Rotation r1 = Rotation::from_euler(0.5, 0.3, 0.7, "zxz");
    Rotation r2 = Rotation::from_euler(0.5, 0.3, 0.7, "zxz");
    Rotation r3 = Rotation::from_euler(0.6, 0.3, 0.7, "zxz");

    EXPECT_TRUE(r1.equals(r2, 1.E-9));
    EXPECT_FALSE(r1.equals(r3, 1.E-9));
}

// =============================================================================
// Edge Cases
// =============================================================================

TEST(TRotationClass, small_angle)
{
    // Very small rotation should not cause numerical issues
    double small_angle = 1.E-10;
    Rotation r = Rotation::from_axis_angle(small_angle, 3);

    vec::fixed<3> v = {1, 0, 0};
    vec::fixed<3> v_rot = r.apply(v);

    // Should be nearly unchanged
    EXPECT_LT(norm(v - v_rot), 1.E-8);
}

TEST(TRotationClass, rotation_180_degrees)
{
    // 180 degree rotation around z-axis
    Rotation r = Rotation::from_axis_angle(simcoon::pi, 3);

    vec::fixed<3> v = {1, 0, 0};
    vec::fixed<3> v_rot = r.apply(v);

    // x should become -x
    EXPECT_NEAR(v_rot(0), -1.0, 1.E-9);
    EXPECT_NEAR(v_rot(1), 0.0, 1.E-9);
    EXPECT_NEAR(v_rot(2), 0.0, 1.E-9);
}

TEST(TRotationClass, random_rotation_properties)
{
    // Random rotation should produce a valid rotation matrix
    Rotation r = Rotation::random();
    mat::fixed<3,3> R = r.as_matrix();

    // R * R^T should be identity
    EXPECT_LT(norm(R * R.t() - eye(3,3), "fro"), 1.E-9);

    // det(R) should be 1
    EXPECT_NEAR(det(R), 1.0, 1.E-9);
}

// =============================================================================
// Matrix Round-trip
// =============================================================================

TEST(TRotationClass, matrix_roundtrip)
{
    // Create rotation, convert to matrix, convert back
    Rotation r1 = Rotation::from_euler(0.5, 0.3, 0.7, "zxz");
    mat::fixed<3,3> R = r1.as_matrix();
    Rotation r2 = Rotation::from_matrix(R);

    EXPECT_TRUE(r1.equals(r2, 1.E-9));
}

TEST(TRotationClass, quat_matrix_quat_roundtrip)
{
    // quat -> matrix -> quat round-trip
    vec::fixed<4> q1 = {0.5, 0.5, 0.5, 0.5};  // This is already unit quaternion
    q1 /= norm(q1);  // Ensure normalization

    Rotation r1 = Rotation::from_quat(q1);
    mat::fixed<3,3> R = r1.as_matrix();
    Rotation r2 = Rotation::from_matrix(R);

    // Note: q and -q represent the same rotation
    EXPECT_TRUE(r1.equals(r2, 1.E-9));
}
