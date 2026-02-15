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

///@file rotation.cpp
///@brief Rotation class and Voigt tensor rotation utilities.
///@version 2.0

#include <cmath>
#include <cassert>
#include <stdexcept>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Simulation/Maths/rotation.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>

using namespace std;
using namespace arma;

namespace simcoon{

// =============================================================================
// Anonymous namespace: internal helpers
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
        return {0.0, 0.0, 0.0, 1.0};
    }
    return q / n;
}

// Hamilton product of two quaternions (scalar-last convention)
vec::fixed<4> quat_multiply(const vec::fixed<4>& q1, const vec::fixed<4>& q2) {
    double x1 = q1(0), y1 = q1(1), z1 = q1(2), w1 = q1(3);
    double x2 = q2(0), y2 = q2(1), z2 = q2(2), w2 = q2(3);

    return {
        w1*x2 + x1*w2 + y1*z2 - z1*y2,
        w1*y2 - x1*z2 + y1*w2 + z1*x2,
        w1*z2 + x1*y2 - y1*x2 + z1*w2,
        w1*w2 - x1*x2 - y1*y2 - z1*z2
    };
}

// Quaternion conjugate (inverse for unit quaternion)
vec::fixed<4> quat_conjugate(const vec::fixed<4>& q) {
    return {-q(0), -q(1), -q(2), q(3)};
}

// Internal: Generate a 3x3 rotation matrix for rotation around a single axis
mat fillR_internal(const double &alpha, const int &axis, const bool &active) {
    double act = active ? 1. : -1.;

    mat R = zeros(3,3);
    double c = cos(alpha);
    double s = act*sin(alpha);

    switch(axis) {
        case 1: {
            R = { {1,0,0}, {0,c,-s}, {0,s,c} };
            break;
        }
        case 2: {
            R = mat {{c, 0, s}, {0,1,0}, {-s,0,c}};
            break;
        }
        case 3: {
            R = {{c,-s,0}, {s,c, 0}, {0,0,1}};
            break;
        }
        default: {
            throw invalid_argument("Please choose the axis 1, 2 or 3");
        }
    }
    return R;
}

// Internal: Generate a 3x3 rotation matrix using Euler angles
mat fillR_euler_internal(const double &psi, const double &theta, const double &phi, const bool &active, const std::string &conv) {
    mat R = zeros(3,3);

    if(conv == "zxz") {
        double c1 = cos(psi);
        double s1 = sin(psi);
        double c2 = cos(theta);
        double s2 = sin(theta);
        double c3 = cos(phi);
        double s3 = sin(phi);

        if (active)
            R = { {c3*c1-c2*s3*s1,-c3*s1-c2*c1*s3,s3*s2}, {c1*s3+c3*c2*s1,c3*c2*c1-s3*s1,-c3*s2}, {s2*s1,c1*s2,c2} };
        else
            R = { {c3*c1-c2*s3*s1,c1*s3+c3*c2*s1,s2*s1}, {-c3*s1-c2*c1*s3,c3*c2*c1-s3*s1,c1*s2}, {s3*s2,-c3*s2,c2} };
    }
    else if(conv == "zyz") {
        double c1 = cos(psi);
        double s1 = sin(psi);
        double c2 = cos(theta);
        double s2 = sin(theta);
        double c3 = cos(phi);
        double s3 = sin(phi);

        if (active)
            R = { {c3*c2*c1-s3*s1,-c1*s3-c3*c2*s1,c3*s2}, {c3*s1+c2*c1*s3,c3*c1-c2*s3*s1,s3*s2}, {-c1*s2,s2*s1,c2} };
        else
            R = { {c3*c2*c1-s3*s1,c3*s1+c2*c1*s3,-c1*s2}, {-c1*s3-c3*c2*s1,c3*c1-c2*s3*s1,s2*s1}, {c3*s2,s3*s2,c2} };
    }
    else if(conv == "user") {
        mat R1 = fillR_internal(psi, simcoon::axis_psi, active);
        mat R2 = fillR_internal(theta, simcoon::axis_theta, active);
        mat R3 = fillR_internal(phi, simcoon::axis_phi, active);
        R = R3*R2*R1;
    }
    else {
        throw std::invalid_argument("error in Simulation/Maths/rotation.cpp: invalid Euler convention '" + conv + "'. Supported conventions: zxz, zyz, xyz, xzy, yxz, yzx, zxy, zyx, xyx, xzx, yxy, yzy, user");
    }

    return R;
}

// Internal: Generate 6x6 stress Voigt rotation matrix from angle/axis
mat fill_voigt_stress(const double &alpha, const int &axis, const bool &active) {
    double act = active ? 1. : -1.;

    mat voigt_stress = zeros(6,6);
    double c = cos(alpha);
    double s = act*sin(alpha);

    switch(axis) {
        case 1: {
            voigt_stress(0,0) = 1.;
            voigt_stress(1,1) = c*c;
            voigt_stress(1,2) = s*s;
            voigt_stress(1,5) = -2.*s*c;
            voigt_stress(2,1) = s*s;
            voigt_stress(2,2) = c*c;
            voigt_stress(2,5) = 2.*s*c;
            voigt_stress(3,3) = c;
            voigt_stress(3,4) = -s;
            voigt_stress(4,3) = s;
            voigt_stress(4,4) = c;
            voigt_stress(5,1) = s*c;
            voigt_stress(5,2) = -s*c;
            voigt_stress(5,5) = c*c-s*s;
            break;
        }
        case 2: {
            voigt_stress(0,0) = c*c;
            voigt_stress(0,2) = s*s;
            voigt_stress(0,4) = 2*c*s;
            voigt_stress(1,1) = 1.;
            voigt_stress(2,0) = s*s;
            voigt_stress(2,2) = c*c;
            voigt_stress(2,4) = -2*c*s;
            voigt_stress(3,3) = c;
            voigt_stress(3,5) = s;
            voigt_stress(4,0) = -c*s;
            voigt_stress(4,2) = c*s;
            voigt_stress(4,4) = c*c-s*s;
            voigt_stress(5,3) = -s;
            voigt_stress(5,5) = c;
            break;
        }
        case 3: {
            voigt_stress(0,0) = c*c;
            voigt_stress(0,1) = s*s;
            voigt_stress(0,3) = -2*s*c;
            voigt_stress(1,0) = s*s;
            voigt_stress(1,1) = c*c;
            voigt_stress(1,3) = 2*s*c;
            voigt_stress(2,2) = 1.;
            voigt_stress(3,0) = s*c;
            voigt_stress(3,1) = -s*c;
            voigt_stress(3,3) = c*c-s*s;
            voigt_stress(4,4) = c;
            voigt_stress(4,5) = -s;
            voigt_stress(5,4) = s;
            voigt_stress(5,5) = c;
            break;
        }
        default: {
            throw invalid_argument("Please choose the axis 1, 2 or 3");
        }
    }
    return voigt_stress;
}

// Internal: Generate 6x6 stress Voigt rotation matrix from rotation matrix
mat fill_voigt_stress(const mat &DR, const bool &active) {
    double a, b, c, d, e, f, g, h, i;

    if(active) {
        a = DR(0,0); b = DR(0,1); c = DR(0,2);
        d = DR(1,0); e = DR(1,1); f = DR(1,2);
        g = DR(2,0); h = DR(2,1); i = DR(2,2);
    } else {
        a = DR(0,0); d = DR(0,1); g = DR(0,2);
        b = DR(1,0); e = DR(1,1); h = DR(1,2);
        c = DR(2,0); f = DR(2,1); i = DR(2,2);
    }

    mat voigt_stress = zeros(6,6);
    voigt_stress(0,0) = a*a;    voigt_stress(0,1) = b*b;    voigt_stress(0,2) = c*c;
    voigt_stress(0,3) = 2.*a*b; voigt_stress(0,4) = 2.*a*c; voigt_stress(0,5) = 2.*b*c;
    voigt_stress(1,0) = d*d;    voigt_stress(1,1) = e*e;    voigt_stress(1,2) = f*f;
    voigt_stress(1,3) = 2.*d*e; voigt_stress(1,4) = 2.*d*f; voigt_stress(1,5) = 2.*e*f;
    voigt_stress(2,0) = g*g;    voigt_stress(2,1) = h*h;    voigt_stress(2,2) = i*i;
    voigt_stress(2,3) = 2.*g*h; voigt_stress(2,4) = 2.*g*i; voigt_stress(2,5) = 2.*h*i;
    voigt_stress(3,0) = a*d;    voigt_stress(3,1) = b*e;    voigt_stress(3,2) = c*f;
    voigt_stress(3,3) = d*b+a*e; voigt_stress(3,4) = d*c+a*f; voigt_stress(3,5) = e*c+b*f;
    voigt_stress(4,0) = a*g;    voigt_stress(4,1) = b*h;    voigt_stress(4,2) = c*i;
    voigt_stress(4,3) = g*b+a*h; voigt_stress(4,4) = g*c+a*i; voigt_stress(4,5) = h*c+b*i;
    voigt_stress(5,0) = d*g;    voigt_stress(5,1) = e*h;    voigt_stress(5,2) = f*i;
    voigt_stress(5,3) = g*e+d*h; voigt_stress(5,4) = g*f+d*i; voigt_stress(5,5) = h*f+e*i;

    return voigt_stress;
}

// Internal: Generate 6x6 strain Voigt rotation matrix from angle/axis
mat fill_voigt_strain(const double &alpha, const int &axis, const bool &active) {
    double act = active ? 1. : -1.;

    mat voigt_strain = zeros(6,6);
    double c = cos(alpha);
    double s = act*sin(alpha);

    switch(axis) {
        case 1: {
            voigt_strain(0,0) = 1.;
            voigt_strain(1,1) = c*c;
            voigt_strain(1,2) = s*s;
            voigt_strain(1,5) = -s*c;
            voigt_strain(2,1) = s*s;
            voigt_strain(2,2) = c*c;
            voigt_strain(2,5) = s*c;
            voigt_strain(3,3) = c;
            voigt_strain(3,4) = -s;
            voigt_strain(4,3) = s;
            voigt_strain(4,4) = c;
            voigt_strain(5,1) = 2.*s*c;
            voigt_strain(5,2) = -2.*s*c;
            voigt_strain(5,5) = c*c-s*s;
            break;
        }
        case 2: {
            voigt_strain(0,0) = c*c;
            voigt_strain(0,2) = s*s;
            voigt_strain(0,4) = c*s;
            voigt_strain(1,1) = 1.;
            voigt_strain(2,0) = s*s;
            voigt_strain(2,2) = c*c;
            voigt_strain(2,4) = -c*s;
            voigt_strain(3,3) = c;
            voigt_strain(3,5) = s;
            voigt_strain(4,0) = -2.*c*s;
            voigt_strain(4,2) = 2.*c*s;
            voigt_strain(4,4) = c*c-s*s;
            voigt_strain(5,3) = -s;
            voigt_strain(5,5) = c;
            break;
        }
        case 3: {
            voigt_strain(0,0) = c*c;
            voigt_strain(0,1) = s*s;
            voigt_strain(0,3) = -s*c;
            voigt_strain(1,0) = s*s;
            voigt_strain(1,1) = c*c;
            voigt_strain(1,3) = s*c;
            voigt_strain(2,2) = 1.;
            voigt_strain(3,0) = 2.*s*c;
            voigt_strain(3,1) = -2.*s*c;
            voigt_strain(3,3) = c*c-s*s;
            voigt_strain(4,4) = c;
            voigt_strain(4,5) = -s;
            voigt_strain(5,4) = s;
            voigt_strain(5,5) = c;
            break;
        }
        default: {
            throw invalid_argument("Please choose the axis 1, 2 or 3");
        }
    }
    return voigt_strain;
}

// Internal: Generate 6x6 strain Voigt rotation matrix from rotation matrix
mat fill_voigt_strain(const mat &DR, const bool &active) {
    double a, b, c, d, e, f, g, h, i;

    if(active) {
        a = DR(0,0); b = DR(0,1); c = DR(0,2);
        d = DR(1,0); e = DR(1,1); f = DR(1,2);
        g = DR(2,0); h = DR(2,1); i = DR(2,2);
    } else {
        a = DR(0,0); d = DR(0,1); g = DR(0,2);
        b = DR(1,0); e = DR(1,1); h = DR(1,2);
        c = DR(2,0); f = DR(2,1); i = DR(2,2);
    }

    mat voigt_strain = zeros(6,6);
    voigt_strain(0,0) = a*a;      voigt_strain(0,1) = b*b;      voigt_strain(0,2) = c*c;
    voigt_strain(0,3) = a*b;      voigt_strain(0,4) = a*c;      voigt_strain(0,5) = b*c;
    voigt_strain(1,0) = d*d;      voigt_strain(1,1) = e*e;      voigt_strain(1,2) = f*f;
    voigt_strain(1,3) = d*e;      voigt_strain(1,4) = d*f;      voigt_strain(1,5) = e*f;
    voigt_strain(2,0) = g*g;      voigt_strain(2,1) = h*h;      voigt_strain(2,2) = i*i;
    voigt_strain(2,3) = g*h;      voigt_strain(2,4) = g*i;      voigt_strain(2,5) = h*i;
    voigt_strain(3,0) = 2.*a*d;   voigt_strain(3,1) = 2.*b*e;   voigt_strain(3,2) = 2.*c*f;
    voigt_strain(3,3) = d*b+a*e;  voigt_strain(3,4) = d*c+a*f;  voigt_strain(3,5) = e*c+b*f;
    voigt_strain(4,0) = 2.*a*g;   voigt_strain(4,1) = 2.*b*h;   voigt_strain(4,2) = 2.*c*i;
    voigt_strain(4,3) = g*b+a*h;  voigt_strain(4,4) = g*c+a*i;  voigt_strain(4,5) = h*c+b*i;
    voigt_strain(5,0) = 2.*d*g;   voigt_strain(5,1) = 2.*e*h;   voigt_strain(5,2) = 2.*f*i;
    voigt_strain(5,3) = g*e+d*h;  voigt_strain(5,4) = g*f+d*i;  voigt_strain(5,5) = h*f+e*i;

    return voigt_strain;
}

} // anonymous namespace

// =============================================================================
// Rotation Class: Constructors
// =============================================================================

Rotation::Rotation() : _quat({0.0, 0.0, 0.0, 1.0}) {
}

Rotation::Rotation(const vec::fixed<4>& quat) : _quat(quat) {
}

// =============================================================================
// Rotation Class: Factory Methods
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
    if (degrees) {
        psi *= simcoon::pi / 180.0;
        theta *= simcoon::pi / 180.0;
        phi *= simcoon::pi / 180.0;
    }

    if (conv == "user") {
        mat R = fillR_euler_internal(psi, theta, phi, true, "user");
        mat::fixed<3,3> R_fixed;
        R_fixed = R;
        return from_matrix(R_fixed);
    }

    int axis1, axis2, axis3;
    parse_euler_convention(conv, axis1, axis2, axis3);

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

    vec::fixed<4> result;
    if (intrinsic) {
        result = quat_multiply(quat_multiply(q3, q2), q1);
    } else {
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
// Rotation Class: Conversion Methods
// =============================================================================

mat::fixed<3,3> Rotation::as_matrix() const {
    double qx = _quat(0), qy = _quat(1), qz = _quat(2), qw = _quat(3);

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

    int axis1, axis2, axis3;

    if (conv == "user") {
        axis1 = axis_psi - 1;
        axis2 = axis_theta - 1;
        axis3 = axis_phi - 1;
    } else {
        parse_euler_convention(conv, axis1, axis2, axis3);
    }

    mat::fixed<3,3> R_work = intrinsic ? R : R.t();
    if (!intrinsic) {
        swap(axis1, axis3);
    }

    vec::fixed<3> angles;

    bool proper_euler = (axis1 == axis3);

    int i = axis1;
    int j = axis2;
    int k = proper_euler ? (3 - i - j) : axis3;

    double sign = ((j - i + 3) % 3 == 1) ? 1.0 : -1.0;

    if (proper_euler) {
        double sy = sqrt(R_work(i,j)*R_work(i,j) + R_work(i,k)*R_work(i,k));

        if (sy > simcoon::iota) {
            angles(0) = atan2(R_work(i,j), sign*R_work(i,k));
            angles(1) = atan2(sy, R_work(i,i));
            angles(2) = atan2(R_work(j,i), -sign*R_work(k,i));
        } else {
            angles(0) = atan2(-sign*R_work(j,k), R_work(j,j));
            angles(1) = atan2(sy, R_work(i,i));
            angles(2) = 0.0;
        }
    } else {
        double cy = sqrt(R_work(i,i)*R_work(i,i) + R_work(j,i)*R_work(j,i));

        if (cy > simcoon::iota) {
            angles(0) = atan2(sign*R_work(k,j), R_work(k,k));
            angles(1) = atan2(-sign*R_work(k,i), cy);
            angles(2) = atan2(sign*R_work(j,i), R_work(i,i));
        } else {
            angles(0) = atan2(-sign*R_work(j,k), R_work(j,j));
            angles(1) = atan2(-sign*R_work(k,i), cy);
            angles(2) = 0.0;
        }
    }

    if (!intrinsic) {
        swap(angles(0), angles(2));
    }

    if (degrees) {
        angles *= 180.0 / simcoon::pi;
    }

    return angles;
}

vec::fixed<3> Rotation::as_rotvec(bool degrees) const {
    double qw = _quat(3);
    vec::fixed<3> qv = {_quat(0), _quat(1), _quat(2)};

    double sin_half = norm(qv);

    if (sin_half < simcoon::iota) {
        vec::fixed<3> result = {0.0, 0.0, 0.0};
        return result;
    }

    double angle = 2.0 * atan2(sin_half, qw);

    if (angle > simcoon::pi) {
        angle -= 2.0 * simcoon::pi;
    } else if (angle < -simcoon::pi) {
        angle += 2.0 * simcoon::pi;
    }

    vec::fixed<3> axis = qv / sin_half;

    if (degrees) {
        angle *= 180.0 / simcoon::pi;
    }

    return axis * angle;
}

mat::fixed<6,6> Rotation::as_voigt_stress_rotation(bool active) const {
    mat R = as_matrix();
    mat vs_dyn = fill_voigt_stress(R, active);
    mat::fixed<6,6> vs;
    vs = vs_dyn;
    return vs;
}

mat::fixed<6,6> Rotation::as_voigt_strain_rotation(bool active) const {
    mat R = as_matrix();
    mat ve_dyn = fill_voigt_strain(R, active);
    mat::fixed<6,6> ve;
    ve = ve_dyn;
    return ve;
}

// =============================================================================
// Rotation Class: Apply Methods (3D objects)
// =============================================================================

vec::fixed<3> Rotation::apply(const vec::fixed<3>& v, bool inverse) const {
    double qx = _quat(0), qy = _quat(1), qz = _quat(2), qw = _quat(3);

    if (inverse) {
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
// Rotation Class: Apply Methods (Voigt notation)
// =============================================================================

vec::fixed<6> Rotation::apply_stress(const vec::fixed<6>& sigma, bool active) const {
    mat::fixed<6,6> vs = as_voigt_stress_rotation(active);
    vec::fixed<6> result;
    result = vs * sigma;
    return result;
}

vec Rotation::apply_stress(const vec& sigma, bool active) const {
    if (sigma.n_elem != 6) {
        throw invalid_argument("Stress vector must have 6 elements");
    }
    mat R = as_matrix();
    mat vs = fill_voigt_stress(R, active);
    return vs * sigma;
}

vec::fixed<6> Rotation::apply_strain(const vec::fixed<6>& epsilon, bool active) const {
    mat::fixed<6,6> ve = as_voigt_strain_rotation(active);
    vec::fixed<6> result;
    result = ve * epsilon;
    return result;
}

vec Rotation::apply_strain(const vec& epsilon, bool active) const {
    if (epsilon.n_elem != 6) {
        throw invalid_argument("Strain vector must have 6 elements");
    }
    mat R = as_matrix();
    mat ve = fill_voigt_strain(R, active);
    return ve * epsilon;
}

mat::fixed<6,6> Rotation::apply_stiffness(const mat::fixed<6,6>& L, bool active) const {
    mat::fixed<6,6> vs = as_voigt_stress_rotation(active);
    mat::fixed<6,6> result;
    result = vs * L * vs.t();
    return result;
}

mat Rotation::apply_stiffness(const mat& L, bool active) const {
    if (L.n_rows != 6 || L.n_cols != 6) {
        throw invalid_argument("Stiffness matrix must be 6x6");
    }
    mat R = as_matrix();
    mat vs = fill_voigt_stress(R, active);
    return vs * (L * trans(vs));
}

mat::fixed<6,6> Rotation::apply_compliance(const mat::fixed<6,6>& M, bool active) const {
    mat::fixed<6,6> ve = as_voigt_strain_rotation(active);
    mat::fixed<6,6> result;
    result = ve * M * ve.t();
    return result;
}

mat Rotation::apply_compliance(const mat& M, bool active) const {
    if (M.n_rows != 6 || M.n_cols != 6) {
        throw invalid_argument("Compliance matrix must be 6x6");
    }
    mat R = as_matrix();
    mat ve = fill_voigt_strain(R, active);
    return ve * (M * trans(ve));
}

mat::fixed<6,6> Rotation::apply_strain_concentration(const mat::fixed<6,6>& A, bool active) const {
    mat::fixed<6,6> ve = as_voigt_strain_rotation(active);
    mat::fixed<6,6> vs = as_voigt_stress_rotation(active);
    mat::fixed<6,6> result;
    result = ve * mat(A) * vs.t();
    return result;
}

mat Rotation::apply_strain_concentration(const mat& A, bool active) const {
    if (A.n_rows != 6 || A.n_cols != 6) {
        throw invalid_argument("Strain concentration tensor must be 6x6");
    }
    mat R = as_matrix();
    mat ve = fill_voigt_strain(R, active);
    mat vs = fill_voigt_stress(R, active);
    return ve * (A * trans(vs));
}

mat::fixed<6,6> Rotation::apply_stress_concentration(const mat::fixed<6,6>& B, bool active) const {
    mat::fixed<6,6> vs = as_voigt_stress_rotation(active);
    mat::fixed<6,6> ve = as_voigt_strain_rotation(active);
    mat::fixed<6,6> result;
    result = vs * mat(B) * ve.t();
    return result;
}

mat Rotation::apply_stress_concentration(const mat& B, bool active) const {
    if (B.n_rows != 6 || B.n_cols != 6) {
        throw invalid_argument("Stress concentration tensor must be 6x6");
    }
    mat R = as_matrix();
    mat vs = fill_voigt_stress(R, active);
    mat ve = fill_voigt_strain(R, active);
    return vs * (B * trans(ve));
}

// =============================================================================
// Rotation Class: Operations
// =============================================================================

Rotation Rotation::inv() const {
    return Rotation(quat_conjugate(_quat));
}

double Rotation::magnitude(bool degrees) const {
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
    vec::fixed<4> q1 = _quat;
    vec::fixed<4> q2 = other._quat;

    double dot = arma::dot(q1, q2);

    if (dot < 0.0) {
        q2 = -q2;
        dot = -dot;
    }

    dot = min(dot, 1.0);

    vec::fixed<4> result;

    if (dot > 0.9995) {
        result = q1 + t * (q2 - q1);
    } else {
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

// =============================================================================
// Convenience Free Functions
// =============================================================================

vec rotate_vec(const vec &v, const mat &DR) {
    return DR*v;
}

vec rotate_vec(const vec &v, const double &alpha, const int &axis) {
    assert(axis>0);
    Rotation rot = Rotation::from_axis_angle(alpha, axis);
    return rot.apply(v);
}

mat rotate_mat(const mat &m, const mat &DR) {
    return trans(DR)*m*DR;
}

mat rotate_mat(const mat &m, const double &alpha, const int &axis) {
    Rotation rot = Rotation::from_axis_angle(alpha, axis);
    return rot.apply_tensor(m);
}

vec rotate_stress(const vec &V, const double &alpha, const int &axis, const bool &active) {
    mat vs = fill_voigt_stress(alpha, axis, active);
    return vs*V;
}

vec rotate_stress(const vec &V, const mat &DR, const bool &active) {
    mat vs = fill_voigt_stress(DR, active);
    return vs*V;
}

vec rotate_strain(const vec &V, const double &alpha, const int &axis, const bool &active) {
    mat ve = fill_voigt_strain(alpha, axis, active);
    return ve*V;
}

vec rotate_strain(const vec &V, const mat &DR, const bool &active) {
    mat ve = fill_voigt_strain(DR, active);
    return ve*V;
}

mat rotate_stiffness(const mat &L, const double &alpha, const int &axis, const bool &active) {
    mat vs = fill_voigt_stress(alpha, axis, active);
    return vs*(L*trans(vs));
}

mat rotate_stiffness(const mat &L, const mat &DR, const bool &active) {
    mat vs = fill_voigt_stress(DR, active);
    return vs*(L*trans(vs));
}

mat rotate_compliance(const mat &M, const double &alpha, const int &axis, const bool &active) {
    mat ve = fill_voigt_strain(alpha, axis, active);
    return ve*(M*trans(ve));
}

mat rotate_compliance(const mat &M, const mat &DR, const bool &active) {
    mat ve = fill_voigt_strain(DR, active);
    return ve*(M*trans(ve));
}

mat rotate_strain_concentration(const mat &A, const double &alpha, const int &axis, const bool &active) {
    mat ve = fill_voigt_strain(alpha, axis, active);
    mat vs = fill_voigt_stress(alpha, axis, active);
    return ve*(A*trans(vs));
}

mat rotate_strain_concentration(const mat &A, const mat &DR, const bool &active) {
    mat ve = fill_voigt_strain(DR, active);
    mat vs = fill_voigt_stress(DR, active);
    return ve*(A*trans(vs));
}

mat rotate_stress_concentration(const mat &B, const double &alpha, const int &axis, const bool &active) {
    mat ve = fill_voigt_strain(alpha, axis, active);
    mat vs = fill_voigt_stress(alpha, axis, active);
    return vs*(B*trans(ve));
}

mat rotate_stress_concentration(const mat &B, const mat &DR, const bool &active) {
    mat ve = fill_voigt_strain(DR, active);
    mat vs = fill_voigt_stress(DR, active);
    return vs*(B*trans(ve));
}

} //namespace simcoon
