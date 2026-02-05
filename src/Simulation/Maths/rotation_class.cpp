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

///@file rotation_class.cpp
///@brief Implementation of the Rotation class
///@version 1.0

#include <cmath>
#include <stdexcept>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Simulation/Maths/rotation_class.hpp>
#include <simcoon/Simulation/Maths/rotation.hpp>

using namespace std;
using namespace arma;

namespace simcoon {

// =============================================================================
// Helper functions (internal)
// =============================================================================

namespace {

// Convert axis number (1,2,3) to axis index (0,1,2)
inline int axis_to_index(int axis) {
    if (axis < 1 || axis > 3) {
        throw invalid_argument("Axis must be 1 (x), 2 (y), or 3 (z)");
    }
    return axis - 1;
}

// Parse Euler convention string to get axis sequence
// Returns axes as 0=x, 1=y, 2=z
void parse_euler_convention(const string& conv, int& axis1, int& axis2, int& axis3) {
    if (conv.length() != 3) {
        throw invalid_argument("Euler convention must be 3 characters (e.g., 'zxz', 'xyz')");
    }

    auto char_to_axis = [](char c) -> int {
        switch (c) {
            case 'x': case 'X': return 0;
            case 'y': case 'Y': return 1;
            case 'z': case 'Z': return 2;
            default:
                throw invalid_argument("Invalid axis character in Euler convention");
        }
    };

    axis1 = char_to_axis(conv[0]);
    axis2 = char_to_axis(conv[1]);
    axis3 = char_to_axis(conv[2]);
}

// Normalize a quaternion to unit length
vec::fixed<4> normalize_quat(const vec::fixed<4>& q) {
    double n = norm(q);
    if (n < simcoon::iota) {
        // Return identity quaternion if input is near zero
        return {0.0, 0.0, 0.0, 1.0};
    }
    return q / n;
}

// Hamilton product of two quaternions (scalar-last convention)
// q = [qx, qy, qz, qw]
vec::fixed<4> quat_multiply(const vec::fixed<4>& q1, const vec::fixed<4>& q2) {
    double x1 = q1(0), y1 = q1(1), z1 = q1(2), w1 = q1(3);
    double x2 = q2(0), y2 = q2(1), z2 = q2(2), w2 = q2(3);

    return {
        w1*x2 + x1*w2 + y1*z2 - z1*y2,  // x
        w1*y2 - x1*z2 + y1*w2 + z1*x2,  // y
        w1*z2 + x1*y2 - y1*x2 + z1*w2,  // z
        w1*w2 - x1*x2 - y1*y2 - z1*z2   // w
    };
}

// Quaternion conjugate (inverse for unit quaternion)
vec::fixed<4> quat_conjugate(const vec::fixed<4>& q) {
    return {-q(0), -q(1), -q(2), q(3)};
}

} // anonymous namespace

// =============================================================================
// Constructors
// =============================================================================

Rotation::Rotation() : _quat({0.0, 0.0, 0.0, 1.0}) {
    // Identity rotation
}

Rotation::Rotation(const vec::fixed<4>& quat) : _quat(quat) {
    // Private constructor - assumes quaternion is already normalized
}

// =============================================================================
// Factory Methods
// =============================================================================

Rotation Rotation::identity() {
    return Rotation();
}

Rotation Rotation::from_quat(const vec::fixed<4>& quat) {
    return Rotation(normalize_quat(quat));
}

Rotation Rotation::from_quat(const vec& quat) {
    if (quat.n_elem != 4) {
        throw invalid_argument("Quaternion must have 4 elements");
    }
    vec::fixed<4> q_fixed = {quat(0), quat(1), quat(2), quat(3)};
    return from_quat(q_fixed);
}

Rotation Rotation::from_matrix(const mat::fixed<3,3>& R) {
    // Shepperd's method for numerical stability
    // Reference: S.W. Shepperd, "Quaternion from Rotation Matrix"
    // Journal of Guidance and Control, Vol. 1, No. 3, 1978

    double trace = R(0,0) + R(1,1) + R(2,2);
    vec::fixed<4> q;

    if (trace > 0) {
        double s = 0.5 / sqrt(trace + 1.0);
        q(3) = 0.25 / s;
        q(0) = (R(2,1) - R(1,2)) * s;
        q(1) = (R(0,2) - R(2,0)) * s;
        q(2) = (R(1,0) - R(0,1)) * s;
    } else if (R(0,0) > R(1,1) && R(0,0) > R(2,2)) {
        double s = 2.0 * sqrt(1.0 + R(0,0) - R(1,1) - R(2,2));
        q(3) = (R(2,1) - R(1,2)) / s;
        q(0) = 0.25 * s;
        q(1) = (R(0,1) + R(1,0)) / s;
        q(2) = (R(0,2) + R(2,0)) / s;
    } else if (R(1,1) > R(2,2)) {
        double s = 2.0 * sqrt(1.0 + R(1,1) - R(0,0) - R(2,2));
        q(3) = (R(0,2) - R(2,0)) / s;
        q(0) = (R(0,1) + R(1,0)) / s;
        q(1) = 0.25 * s;
        q(2) = (R(1,2) + R(2,1)) / s;
    } else {
        double s = 2.0 * sqrt(1.0 + R(2,2) - R(0,0) - R(1,1));
        q(3) = (R(1,0) - R(0,1)) / s;
        q(0) = (R(0,2) + R(2,0)) / s;
        q(1) = (R(1,2) + R(2,1)) / s;
        q(2) = 0.25 * s;
    }

    return Rotation(normalize_quat(q));
}

Rotation Rotation::from_matrix(const mat& R) {
    if (R.n_rows != 3 || R.n_cols != 3) {
        throw invalid_argument("Rotation matrix must be 3x3");
    }
    mat::fixed<3,3> R_fixed;
    R_fixed = R;
    return from_matrix(R_fixed);
}

Rotation Rotation::from_euler(double psi, double theta, double phi,
                              const string& conv, bool intrinsic, bool degrees) {
    // Convert to radians if needed
    if (degrees) {
        psi *= simcoon::pi / 180.0;
        theta *= simcoon::pi / 180.0;
        phi *= simcoon::pi / 180.0;
    }

    // Handle "user" convention using the existing fillR function
    if (conv == "user") {
        mat R = fillR(psi, theta, phi, true, "user");
        mat::fixed<3,3> R_fixed;
        R_fixed = R;
        return from_matrix(R_fixed);
    }

    // Parse the convention
    int axis1, axis2, axis3;
    parse_euler_convention(conv, axis1, axis2, axis3);

    // Build quaternions for each elementary rotation
    auto axis_quat = [](double angle, int axis) -> vec::fixed<4> {
        double half_angle = angle / 2.0;
        double s = sin(half_angle);
        double c = cos(half_angle);
        vec::fixed<4> q = {0, 0, 0, c};
        q(axis) = s;
        return q;
    };

    vec::fixed<4> q1 = axis_quat(psi, axis1);
    vec::fixed<4> q2 = axis_quat(theta, axis2);
    vec::fixed<4> q3 = axis_quat(phi, axis3);

    // Compose quaternions
    vec::fixed<4> result;
    if (intrinsic) {
        // Intrinsic: rotations about axes of the rotating frame
        // Order: q3 * q2 * q1 (applied right-to-left)
        result = quat_multiply(quat_multiply(q3, q2), q1);
    } else {
        // Extrinsic: rotations about fixed axes
        // Order: q1 * q2 * q3 (applied left-to-right)
        result = quat_multiply(quat_multiply(q1, q2), q3);
    }

    return Rotation(normalize_quat(result));
}

Rotation Rotation::from_euler(const vec::fixed<3>& angles, const string& conv,
                              bool intrinsic, bool degrees) {
    return from_euler(angles(0), angles(1), angles(2), conv, intrinsic, degrees);
}

Rotation Rotation::from_rotvec(const vec::fixed<3>& rotvec, bool degrees) {
    double angle = norm(rotvec);

    if (degrees) {
        angle *= simcoon::pi / 180.0;
    }

    if (angle < simcoon::iota) {
        return identity();
    }

    vec::fixed<3> axis = rotvec / (degrees ? norm(rotvec) * 180.0 / simcoon::pi : norm(rotvec));
    double half_angle = angle / 2.0;
    double s = sin(half_angle);
    double c = cos(half_angle);

    return Rotation({axis(0)*s, axis(1)*s, axis(2)*s, c});
}

Rotation Rotation::from_rotvec(const vec& rotvec, bool degrees) {
    if (rotvec.n_elem != 3) {
        throw invalid_argument("Rotation vector must have 3 elements");
    }
    vec::fixed<3> rv_fixed = {rotvec(0), rotvec(1), rotvec(2)};
    return from_rotvec(rv_fixed, degrees);
}

Rotation Rotation::from_axis_angle(double angle, int axis, bool degrees) {
    if (degrees) {
        angle *= simcoon::pi / 180.0;
    }

    int ax_idx = axis_to_index(axis);
    double half_angle = angle / 2.0;
    double s = sin(half_angle);
    double c = cos(half_angle);

    vec::fixed<4> q = {0, 0, 0, c};
    q(ax_idx) = s;

    return Rotation(q);
}

Rotation Rotation::random() {
    // Generate uniformly distributed random rotation using Shoemake's method
    // Reference: K. Shoemake, "Uniform random rotations", Graphics Gems III, 1992

    vec u = randu<vec>(3);

    double sqrt1_u0 = sqrt(1.0 - u(0));
    double sqrt_u0 = sqrt(u(0));
    double two_pi_u1 = 2.0 * simcoon::pi * u(1);
    double two_pi_u2 = 2.0 * simcoon::pi * u(2);

    vec::fixed<4> q = {
        sqrt1_u0 * sin(two_pi_u1),
        sqrt1_u0 * cos(two_pi_u1),
        sqrt_u0 * sin(two_pi_u2),
        sqrt_u0 * cos(two_pi_u2)
    };

    return Rotation(q);
}

// =============================================================================
// Conversion Methods
// =============================================================================

mat::fixed<3,3> Rotation::as_matrix() const {
    double qx = _quat(0), qy = _quat(1), qz = _quat(2), qw = _quat(3);

    // Precompute products
    double qx2 = qx * qx, qy2 = qy * qy, qz2 = qz * qz;
    double qxqy = qx * qy, qxqz = qx * qz, qxqw = qx * qw;
    double qyqz = qy * qz, qyqw = qy * qw, qzqw = qz * qw;

    mat::fixed<3,3> R;
    R(0,0) = 1.0 - 2.0*(qy2 + qz2);
    R(0,1) = 2.0*(qxqy - qzqw);
    R(0,2) = 2.0*(qxqz + qyqw);
    R(1,0) = 2.0*(qxqy + qzqw);
    R(1,1) = 1.0 - 2.0*(qx2 + qz2);
    R(1,2) = 2.0*(qyqz - qxqw);
    R(2,0) = 2.0*(qxqz - qyqw);
    R(2,1) = 2.0*(qyqz + qxqw);
    R(2,2) = 1.0 - 2.0*(qx2 + qy2);

    return R;
}

vec::fixed<3> Rotation::as_euler(const string& conv, bool intrinsic, bool degrees) const {
    mat::fixed<3,3> R = as_matrix();

    // Parse the convention
    int axis1, axis2, axis3;

    if (conv == "user") {
        // Use default axes from parameter.hpp
        axis1 = axis_psi - 1;   // Convert 1-indexed to 0-indexed
        axis2 = axis_theta - 1;
        axis3 = axis_phi - 1;
    } else {
        parse_euler_convention(conv, axis1, axis2, axis3);
    }

    // For extrinsic, we transpose R and swap axis order
    mat::fixed<3,3> R_work = intrinsic ? R : R.t();
    if (!intrinsic) {
        swap(axis1, axis3);
    }

    vec::fixed<3> angles;

    // Determine if proper Euler (axis1 == axis3) or Tait-Bryan (axis1 != axis3)
    bool proper_euler = (axis1 == axis3);

    // Find the "other" axis for proper Euler, or compute for Tait-Bryan
    int i = axis1;
    int j = axis2;
    int k = proper_euler ? (3 - i - j) : axis3;  // For proper Euler, k is the third axis

    // Sign based on cyclic permutation
    // Cyclic: (0,1,2), (1,2,0), (2,0,1) give +1
    // Anti-cyclic: (0,2,1), (2,1,0), (1,0,2) give -1
    double sign = ((j - i + 3) % 3 == 1) ? 1.0 : -1.0;

    if (proper_euler) {
        // Proper Euler angles (e.g., zxz, zyz)
        double sy = sqrt(R_work(i,j)*R_work(i,j) + R_work(i,k)*R_work(i,k));

        if (sy > simcoon::iota) {
            angles(0) = atan2(R_work(i,j), sign*R_work(i,k));
            angles(1) = atan2(sy, R_work(i,i));
            angles(2) = atan2(R_work(j,i), -sign*R_work(k,i));
        } else {
            // Gimbal lock
            angles(0) = atan2(-sign*R_work(j,k), R_work(j,j));
            angles(1) = atan2(sy, R_work(i,i));
            angles(2) = 0.0;
        }
    } else {
        // Tait-Bryan angles (e.g., xyz, zyx)
        double cy = sqrt(R_work(i,i)*R_work(i,i) + R_work(j,i)*R_work(j,i));

        if (cy > simcoon::iota) {
            angles(0) = atan2(sign*R_work(k,j), R_work(k,k));
            angles(1) = atan2(-sign*R_work(k,i), cy);
            angles(2) = atan2(sign*R_work(j,i), R_work(i,i));
        } else {
            // Gimbal lock
            angles(0) = atan2(-sign*R_work(j,k), R_work(j,j));
            angles(1) = atan2(-sign*R_work(k,i), cy);
            angles(2) = 0.0;
        }
    }

    // For extrinsic, swap back
    if (!intrinsic) {
        swap(angles(0), angles(2));
    }

    if (degrees) {
        angles *= 180.0 / simcoon::pi;
    }

    return angles;
}

vec::fixed<3> Rotation::as_rotvec(bool degrees) const {
    // Get rotation angle from quaternion
    double qw = _quat(3);
    vec::fixed<3> qv = {_quat(0), _quat(1), _quat(2)};

    double sin_half = norm(qv);

    if (sin_half < simcoon::iota) {
        // Near identity rotation
        vec::fixed<3> result = {0.0, 0.0, 0.0};
        return result;
    }

    // Angle = 2 * atan2(||qv||, qw)
    double angle = 2.0 * atan2(sin_half, qw);

    // Normalize to [-pi, pi]
    if (angle > simcoon::pi) {
        angle -= 2.0 * simcoon::pi;
    } else if (angle < -simcoon::pi) {
        angle += 2.0 * simcoon::pi;
    }

    // Axis is normalized qv
    vec::fixed<3> axis = qv / sin_half;

    if (degrees) {
        angle *= 180.0 / simcoon::pi;
    }

    return axis * angle;
}

mat::fixed<6,6> Rotation::as_QS(bool active) const {
    mat R = as_matrix();
    mat QS_dyn = fillQS(R, active);
    mat::fixed<6,6> QS;
    QS = QS_dyn;
    return QS;
}

mat::fixed<6,6> Rotation::as_QE(bool active) const {
    mat R = as_matrix();
    mat QE_dyn = fillQE(R, active);
    mat::fixed<6,6> QE;
    QE = QE_dyn;
    return QE;
}

// =============================================================================
// Apply Methods (3D objects)
// =============================================================================

vec::fixed<3> Rotation::apply(const vec::fixed<3>& v, bool inverse) const {
    // Direct quaternion-vector rotation using Rodrigues formula
    // v' = v + 2*qw*(qv × v) + 2*(qv × (qv × v))
    // This is equivalent to q * v * q^(-1) but more efficient

    double qx = _quat(0), qy = _quat(1), qz = _quat(2), qw = _quat(3);

    if (inverse) {
        // For inverse, negate the vector part of quaternion
        qx = -qx; qy = -qy; qz = -qz;
    }

    vec::fixed<3> qv = {qx, qy, qz};
    vec::fixed<3> uv = cross(qv, v);
    vec::fixed<3> uuv = cross(qv, uv);

    return v + 2.0 * (qw * uv + uuv);
}

vec Rotation::apply(const vec& v, bool inverse) const {
    if (v.n_elem != 3) {
        throw invalid_argument("Vector must have 3 elements");
    }
    vec::fixed<3> v_fixed = {v(0), v(1), v(2)};
    vec::fixed<3> result = apply(v_fixed, inverse);
    return vec(result);
}

mat::fixed<3,3> Rotation::apply_tensor(const mat::fixed<3,3>& m, bool inverse) const {
    mat::fixed<3,3> R = as_matrix();
    if (inverse) {
        return R.t() * m * R;
    } else {
        return R * m * R.t();
    }
}

mat Rotation::apply_tensor(const mat& m, bool inverse) const {
    if (m.n_rows != 3 || m.n_cols != 3) {
        throw invalid_argument("Tensor must be 3x3");
    }
    mat::fixed<3,3> m_fixed;
    m_fixed = m;
    mat::fixed<3,3> result = apply_tensor(m_fixed, inverse);
    return mat(result);
}

// =============================================================================
// Apply Methods (Voigt notation)
// =============================================================================

vec::fixed<6> Rotation::apply_stress(const vec::fixed<6>& sigma, bool active) const {
    mat::fixed<6,6> QS = as_QS(active);
    vec::fixed<6> result;
    result = QS * sigma;
    return result;
}

vec Rotation::apply_stress(const vec& sigma, bool active) const {
    if (sigma.n_elem != 6) {
        throw invalid_argument("Stress vector must have 6 elements");
    }
    mat R = as_matrix();
    return rotate_stress(sigma, R, active);
}

vec::fixed<6> Rotation::apply_strain(const vec::fixed<6>& epsilon, bool active) const {
    mat::fixed<6,6> QE = as_QE(active);
    vec::fixed<6> result;
    result = QE * epsilon;
    return result;
}

vec Rotation::apply_strain(const vec& epsilon, bool active) const {
    if (epsilon.n_elem != 6) {
        throw invalid_argument("Strain vector must have 6 elements");
    }
    mat R = as_matrix();
    return rotate_strain(epsilon, R, active);
}

mat::fixed<6,6> Rotation::apply_stiffness(const mat::fixed<6,6>& L, bool active) const {
    mat::fixed<6,6> QS = as_QS(active);
    mat::fixed<6,6> result;
    result = QS * L * QS.t();
    return result;
}

mat Rotation::apply_stiffness(const mat& L, bool active) const {
    if (L.n_rows != 6 || L.n_cols != 6) {
        throw invalid_argument("Stiffness matrix must be 6x6");
    }
    mat R = as_matrix();
    return rotateL(L, R, active);
}

mat::fixed<6,6> Rotation::apply_compliance(const mat::fixed<6,6>& M, bool active) const {
    mat::fixed<6,6> QE = as_QE(active);
    mat::fixed<6,6> result;
    result = QE * M * QE.t();
    return result;
}

mat Rotation::apply_compliance(const mat& M, bool active) const {
    if (M.n_rows != 6 || M.n_cols != 6) {
        throw invalid_argument("Compliance matrix must be 6x6");
    }
    mat R = as_matrix();
    return rotateM(M, R, active);
}

// =============================================================================
// Operations
// =============================================================================

Rotation Rotation::inv() const {
    return Rotation(quat_conjugate(_quat));
}

double Rotation::magnitude(bool degrees) const {
    // Rotation angle = 2 * acos(|qw|)
    double angle = 2.0 * acos(min(abs(_quat(3)), 1.0));

    if (degrees) {
        angle *= 180.0 / simcoon::pi;
    }

    return angle;
}

Rotation Rotation::operator*(const Rotation& other) const {
    return Rotation(normalize_quat(quat_multiply(_quat, other._quat)));
}

Rotation& Rotation::operator*=(const Rotation& other) {
    _quat = normalize_quat(quat_multiply(_quat, other._quat));
    return *this;
}

Rotation Rotation::slerp(const Rotation& other, double t) const {
    // Spherical linear interpolation
    // Reference: K. Shoemake, "Animating rotation with quaternion curves", SIGGRAPH 1985

    vec::fixed<4> q1 = _quat;
    vec::fixed<4> q2 = other._quat;

    // Compute dot product
    double dot = arma::dot(q1, q2);

    // If dot is negative, negate one quaternion to take shorter path
    if (dot < 0.0) {
        q2 = -q2;
        dot = -dot;
    }

    // Clamp dot to valid range
    dot = min(dot, 1.0);

    vec::fixed<4> result;

    if (dot > 0.9995) {
        // Quaternions are very close, use linear interpolation
        result = q1 + t * (q2 - q1);
    } else {
        // Standard SLERP
        double theta_0 = acos(dot);
        double theta = theta_0 * t;
        double sin_theta = sin(theta);
        double sin_theta_0 = sin(theta_0);

        double s0 = cos(theta) - dot * sin_theta / sin_theta_0;
        double s1 = sin_theta / sin_theta_0;

        result = s0 * q1 + s1 * q2;
    }

    return Rotation(normalize_quat(result));
}

} // namespace simcoon
