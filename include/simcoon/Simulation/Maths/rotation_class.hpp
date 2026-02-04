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

///@file rotation_class.hpp
///@brief Rotation class inspired by scipy.spatial.transform.Rotation
///@version 1.0

#pragma once

#include <armadillo>
#include <string>

namespace simcoon {

/**
 * @class Rotation
 * @brief A class representing 3D rotations using unit quaternions.
 *
 * Inspired by scipy.spatial.transform.Rotation, this class provides a unified
 * interface for working with rotations in various representations (quaternion,
 * rotation matrix, Euler angles, rotation vector) while maintaining high
 * numerical precision and performance.
 *
 * Internal representation uses unit quaternions in scalar-last convention [qx, qy, qz, qw].
 *
 * @details Example usage:
 * @code
 *      // Create rotation from Euler angles
 *      Rotation r = Rotation::from_euler(0.5, 0.3, 0.7, "zxz");
 *
 *      // Apply to a vector
 *      arma::vec::fixed<3> v = {1, 0, 0};
 *      arma::vec::fixed<3> v_rot = r.apply(v);
 *
 *      // Compose rotations
 *      Rotation r2 = Rotation::from_axis_angle(M_PI/4, 3);
 *      Rotation r_combined = r * r2;
 *
 *      // Apply to Voigt stress tensor
 *      arma::vec::fixed<6> sigma = {100, 50, 0, 25, 0, 0};
 *      arma::vec::fixed<6> sigma_rot = r.apply_stress(sigma);
 * @endcode
 */
class Rotation {
private:
    arma::vec::fixed<4> _quat;  // [qx, qy, qz, qw] - scalar-last convention, stack allocated

    // Private constructor from quaternion (assumes already normalized)
    explicit Rotation(const arma::vec::fixed<4>& quat);

public:
    /**
     * @brief Default constructor creating identity rotation.
     */
    Rotation();

    /**
     * @brief Copy constructor.
     */
    Rotation(const Rotation& other) = default;

    /**
     * @brief Move constructor.
     */
    Rotation(Rotation&& other) noexcept = default;

    /**
     * @brief Copy assignment operator.
     */
    Rotation& operator=(const Rotation& other) = default;

    /**
     * @brief Move assignment operator.
     */
    Rotation& operator=(Rotation&& other) noexcept = default;

    /**
     * @brief Destructor.
     */
    ~Rotation() = default;

    // =========================================================================
    // Factory Methods
    // =========================================================================

    /**
     * @brief Create an identity rotation.
     * @return Identity rotation (no rotation)
     */
    static Rotation identity();

    /**
     * @brief Create rotation from a unit quaternion.
     * @param quat Quaternion in [qx, qy, qz, qw] format (scalar-last)
     * @return Rotation object
     * @note The quaternion will be normalized if not already unit length
     */
    static Rotation from_quat(const arma::vec::fixed<4>& quat);

    /**
     * @brief Create rotation from a unit quaternion.
     * @param quat Quaternion in [qx, qy, qz, qw] format (scalar-last), dynamic size
     * @return Rotation object
     */
    static Rotation from_quat(const arma::vec& quat);

    /**
     * @brief Create rotation from a 3x3 rotation matrix.
     * @param R 3x3 rotation matrix (must be orthogonal with det=1)
     * @return Rotation object
     */
    static Rotation from_matrix(const arma::mat::fixed<3,3>& R);

    /**
     * @brief Create rotation from a 3x3 rotation matrix.
     * @param R 3x3 rotation matrix (dynamic size)
     * @return Rotation object
     */
    static Rotation from_matrix(const arma::mat& R);

    /**
     * @brief Create rotation from Euler angles.
     * @param psi First Euler angle in radians (or degrees if degrees=true)
     * @param theta Second Euler angle in radians (or degrees if degrees=true)
     * @param phi Third Euler angle in radians (or degrees if degrees=true)
     * @param conv Euler angle convention: "zxz", "zyz", "xyz", "xzy", "yxz", "yzx", "zxy", "zyx",
     *             "xyx", "xzx", "yxy", "yzy", or "user"
     * @param intrinsic If true (default), intrinsic rotations (rotating frame);
     *                  if false, extrinsic rotations (fixed frame)
     * @param degrees If true, angles are in degrees; if false (default), radians
     * @return Rotation object
     */
    static Rotation from_euler(double psi, double theta, double phi,
                               const std::string& conv = "zxz",
                               bool intrinsic = true, bool degrees = false);

    /**
     * @brief Create rotation from Euler angles (vector input).
     * @param angles 3-element vector of Euler angles
     * @param conv Euler angle convention
     * @param intrinsic If true, intrinsic rotations; if false, extrinsic rotations
     * @param degrees If true, angles are in degrees
     * @return Rotation object
     */
    static Rotation from_euler(const arma::vec::fixed<3>& angles,
                               const std::string& conv = "zxz",
                               bool intrinsic = true, bool degrees = false);

    /**
     * @brief Create rotation from a rotation vector (axis-angle representation).
     * @param rotvec Rotation vector where direction is axis and magnitude is angle
     * @param degrees If true, magnitude is in degrees; if false (default), radians
     * @return Rotation object
     */
    static Rotation from_rotvec(const arma::vec::fixed<3>& rotvec, bool degrees = false);

    /**
     * @brief Create rotation from a rotation vector (dynamic size).
     * @param rotvec Rotation vector (3 elements)
     * @param degrees If true, magnitude is in degrees
     * @return Rotation object
     */
    static Rotation from_rotvec(const arma::vec& rotvec, bool degrees = false);

    /**
     * @brief Create rotation around a single axis.
     * @param angle Rotation angle in radians (or degrees if degrees=true)
     * @param axis Axis of rotation: 1=x, 2=y, 3=z
     * @param degrees If true, angle is in degrees
     * @return Rotation object
     */
    static Rotation from_axis_angle(double angle, int axis, bool degrees = false);

    /**
     * @brief Create a random rotation (uniformly distributed).
     * @return Random rotation object
     */
    static Rotation random();

    // =========================================================================
    // Conversion Methods
    // =========================================================================

    /**
     * @brief Get rotation as a unit quaternion.
     * @return Quaternion in [qx, qy, qz, qw] format (scalar-last)
     */
    inline arma::vec::fixed<4> as_quat() const { return _quat; }

    /**
     * @brief Get rotation as a 3x3 rotation matrix.
     * @return 3x3 orthogonal rotation matrix
     */
    arma::mat::fixed<3,3> as_matrix() const;

    /**
     * @brief Get rotation as Euler angles.
     * @param conv Euler angle convention
     * @param intrinsic If true, return intrinsic angles; if false, extrinsic
     * @param degrees If true, return angles in degrees
     * @return 3-element vector of Euler angles [psi, theta, phi]
     */
    arma::vec::fixed<3> as_euler(const std::string& conv = "zxz",
                                  bool intrinsic = true, bool degrees = false) const;

    /**
     * @brief Get rotation as a rotation vector.
     * @param degrees If true, magnitude is in degrees
     * @return Rotation vector (axis * angle)
     */
    arma::vec::fixed<3> as_rotvec(bool degrees = false) const;

    /**
     * @brief Get 6x6 rotation matrix for stress tensors in Voigt notation.
     * @param active If true (default), active rotation; if false, passive
     * @return 6x6 stress rotation matrix (QS)
     */
    arma::mat::fixed<6,6> as_QS(bool active = true) const;

    /**
     * @brief Get 6x6 rotation matrix for strain tensors in Voigt notation.
     * @param active If true (default), active rotation; if false, passive
     * @return 6x6 strain rotation matrix (QE)
     */
    arma::mat::fixed<6,6> as_QE(bool active = true) const;

    // =========================================================================
    // Apply Methods (3D objects)
    // =========================================================================

    /**
     * @brief Apply rotation to a 3D vector.
     * @param v 3D vector to rotate
     * @param inverse If true, apply inverse rotation
     * @return Rotated vector
     * @details Uses direct quaternion-vector rotation (faster than matrix multiplication)
     */
    arma::vec::fixed<3> apply(const arma::vec::fixed<3>& v, bool inverse = false) const;

    /**
     * @brief Apply rotation to a 3D vector (dynamic size).
     * @param v 3D vector to rotate
     * @param inverse If true, apply inverse rotation
     * @return Rotated vector
     */
    arma::vec apply(const arma::vec& v, bool inverse = false) const;

    /**
     * @brief Apply rotation to a 3x3 tensor (matrix).
     * @param m 3x3 tensor to rotate
     * @param inverse If true, apply inverse rotation
     * @return Rotated tensor: R * m * R^T (or R^T * m * R for inverse)
     */
    arma::mat::fixed<3,3> apply_tensor(const arma::mat::fixed<3,3>& m, bool inverse = false) const;

    /**
     * @brief Apply rotation to a 3x3 tensor (dynamic size).
     * @param m 3x3 tensor to rotate
     * @param inverse If true, apply inverse rotation
     * @return Rotated tensor
     */
    arma::mat apply_tensor(const arma::mat& m, bool inverse = false) const;

    // =========================================================================
    // Apply Methods (Voigt notation)
    // =========================================================================

    /**
     * @brief Apply rotation to a stress vector in Voigt notation.
     * @param sigma 6-component stress vector [s11, s22, s33, s12, s13, s23]
     * @param active If true (default), active rotation
     * @return Rotated stress vector
     */
    arma::vec::fixed<6> apply_stress(const arma::vec::fixed<6>& sigma, bool active = true) const;

    /**
     * @brief Apply rotation to a stress vector (dynamic size).
     * @param sigma 6-component stress vector
     * @param active If true, active rotation
     * @return Rotated stress vector
     */
    arma::vec apply_stress(const arma::vec& sigma, bool active = true) const;

    /**
     * @brief Apply rotation to a strain vector in Voigt notation.
     * @param epsilon 6-component strain vector [e11, e22, e33, 2*e12, 2*e13, 2*e23]
     * @param active If true (default), active rotation
     * @return Rotated strain vector
     */
    arma::vec::fixed<6> apply_strain(const arma::vec::fixed<6>& epsilon, bool active = true) const;

    /**
     * @brief Apply rotation to a strain vector (dynamic size).
     * @param epsilon 6-component strain vector
     * @param active If true, active rotation
     * @return Rotated strain vector
     */
    arma::vec apply_strain(const arma::vec& epsilon, bool active = true) const;

    /**
     * @brief Apply rotation to a 6x6 stiffness matrix.
     * @param L 6x6 stiffness matrix in Voigt notation
     * @param active If true (default), active rotation
     * @return Rotated stiffness matrix: QS * L * QS^T
     */
    arma::mat::fixed<6,6> apply_stiffness(const arma::mat::fixed<6,6>& L, bool active = true) const;

    /**
     * @brief Apply rotation to a stiffness matrix (dynamic size).
     * @param L 6x6 stiffness matrix
     * @param active If true, active rotation
     * @return Rotated stiffness matrix
     */
    arma::mat apply_stiffness(const arma::mat& L, bool active = true) const;

    /**
     * @brief Apply rotation to a 6x6 compliance matrix.
     * @param M 6x6 compliance matrix in Voigt notation
     * @param active If true (default), active rotation
     * @return Rotated compliance matrix: QE * M * QE^T
     */
    arma::mat::fixed<6,6> apply_compliance(const arma::mat::fixed<6,6>& M, bool active = true) const;

    /**
     * @brief Apply rotation to a compliance matrix (dynamic size).
     * @param M 6x6 compliance matrix
     * @param active If true, active rotation
     * @return Rotated compliance matrix
     */
    arma::mat apply_compliance(const arma::mat& M, bool active = true) const;

    // =========================================================================
    // Operations
    // =========================================================================

    /**
     * @brief Get the inverse (conjugate) rotation.
     * @return Inverse rotation
     */
    Rotation inv() const;

    /**
     * @brief Get the magnitude (rotation angle) of this rotation.
     * @param degrees If true, return in degrees
     * @return Rotation angle in radians (or degrees)
     */
    double magnitude(bool degrees = false) const;

    /**
     * @brief Compose this rotation with another (this * other).
     * @param other Rotation to compose with
     * @return Composed rotation
     * @details The composition represents applying 'other' first, then 'this'
     */
    Rotation operator*(const Rotation& other) const;

    /**
     * @brief Compose this rotation with another in-place.
     * @param other Rotation to compose with
     * @return Reference to this rotation
     */
    Rotation& operator*=(const Rotation& other);

    /**
     * @brief Spherical linear interpolation between this rotation and another.
     * @param other Target rotation
     * @param t Interpolation parameter in [0, 1] (0 = this, 1 = other)
     * @return Interpolated rotation
     */
    Rotation slerp(const Rotation& other, double t) const;

    /**
     * @brief Check if this rotation equals another within tolerance.
     * @param other Rotation to compare with
     * @param tol Tolerance for comparison
     * @return True if rotations are equal (or antipodal quaternions)
     */
    inline bool equals(const Rotation& other, double tol = 1e-12) const {
        // Quaternions q and -q represent the same rotation
        double diff1 = arma::norm(_quat - other._quat);
        double diff2 = arma::norm(_quat + other._quat);
        return (diff1 < tol) || (diff2 < tol);
    }

    /**
     * @brief Check if this rotation is identity (no rotation).
     * @param tol Tolerance for comparison
     * @return True if this is the identity rotation
     */
    inline bool is_identity(double tol = 1e-12) const {
        // Identity quaternion is [0, 0, 0, 1] or [0, 0, 0, -1]
        return std::abs(std::abs(_quat(3)) - 1.0) < tol;
    }
};

} // namespace simcoon
