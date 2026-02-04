#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <string>
#include <carma>
#include <armadillo>

#include <simcoon/Simulation/Maths/rotation_class.hpp>
#include <simcoon/python_wrappers/Libraries/Maths/rotation_class.hpp>

using namespace std;
using namespace arma;
namespace py = pybind11;

namespace simpy {

void register_rotation_class(py::module_& m) {
    py::class_<simcoon::Rotation>(m, "Rotation",
        R"doc(
        A class representing 3D rotations using unit quaternions.

        Inspired by scipy.spatial.transform.Rotation, this class provides a unified
        interface for working with rotations in various representations (quaternion,
        rotation matrix, Euler angles, rotation vector) while maintaining high
        numerical precision and performance.

        The internal representation uses unit quaternions in scalar-last convention [qx, qy, qz, qw].

        Example
        -------
        >>> import simcoon as smc
        >>> import numpy as np
        >>>
        >>> # Create rotation from Euler angles
        >>> r = smc.Rotation.from_euler(0.5, 0.3, 0.7, "zxz")
        >>>
        >>> # Apply to a vector
        >>> v = np.array([1.0, 0.0, 0.0])
        >>> v_rot = r.apply(v)
        >>>
        >>> # Compose rotations
        >>> r2 = smc.Rotation.from_axis_angle(np.pi/4, 3)
        >>> r_combined = r * r2
        >>>
        >>> # Apply to Voigt stress tensor
        >>> sigma = np.array([100.0, 50.0, 0.0, 25.0, 0.0, 0.0])
        >>> sigma_rot = r.apply_stress(sigma)
        )doc")

        // Default constructor
        .def(py::init<>(), "Create an identity rotation")

        // Static factory methods
        .def_static("identity", &simcoon::Rotation::identity,
            "Create an identity rotation (no rotation)")

        .def_static("from_quat",
            [](py::array_t<double> quat) {
                vec q = carma::arr_to_col(quat);
                return simcoon::Rotation::from_quat(q);
            },
            py::arg("quat"),
            R"doc(
            Create rotation from a unit quaternion.

            Parameters
            ----------
            quat : array_like
                Quaternion in [qx, qy, qz, qw] format (scalar-last).
                Will be normalized if not already unit length.

            Returns
            -------
            Rotation
                Rotation object
            )doc")

        .def_static("from_matrix",
            [](py::array_t<double> R) {
                mat R_mat = carma::arr_to_mat(R);
                return simcoon::Rotation::from_matrix(R_mat);
            },
            py::arg("R"),
            R"doc(
            Create rotation from a 3x3 rotation matrix.

            Parameters
            ----------
            R : array_like
                3x3 rotation matrix (must be orthogonal with det=1)

            Returns
            -------
            Rotation
                Rotation object
            )doc")

        .def_static("from_euler",
            py::overload_cast<double, double, double, const string&, bool, bool>(
                &simcoon::Rotation::from_euler),
            py::arg("psi"), py::arg("theta"), py::arg("phi"),
            py::arg("conv") = "zxz", py::arg("intrinsic") = true, py::arg("degrees") = false,
            R"doc(
            Create rotation from Euler angles.

            Parameters
            ----------
            psi : float
                First Euler angle
            theta : float
                Second Euler angle
            phi : float
                Third Euler angle
            conv : str, optional
                Euler angle convention: "zxz", "zyz", "xyz", "xzy", "yxz", "yzx", "zxy", "zyx",
                "xyx", "xzx", "yxy", "yzy", or "user". Default is "zxz".
            intrinsic : bool, optional
                If True (default), intrinsic rotations (rotating frame);
                if False, extrinsic rotations (fixed frame)
            degrees : bool, optional
                If True, angles are in degrees; if False (default), radians

            Returns
            -------
            Rotation
                Rotation object
            )doc")

        .def_static("from_rotvec",
            [](py::array_t<double> rotvec, bool degrees) {
                vec rv = carma::arr_to_col(rotvec);
                return simcoon::Rotation::from_rotvec(rv, degrees);
            },
            py::arg("rotvec"), py::arg("degrees") = false,
            R"doc(
            Create rotation from a rotation vector.

            Parameters
            ----------
            rotvec : array_like
                Rotation vector where direction is axis and magnitude is angle
            degrees : bool, optional
                If True, magnitude is in degrees; if False (default), radians

            Returns
            -------
            Rotation
                Rotation object
            )doc")

        .def_static("from_axis_angle", &simcoon::Rotation::from_axis_angle,
            py::arg("angle"), py::arg("axis"), py::arg("degrees") = false,
            R"doc(
            Create rotation around a single axis.

            Parameters
            ----------
            angle : float
                Rotation angle
            axis : int
                Axis of rotation: 1=x, 2=y, 3=z
            degrees : bool, optional
                If True, angle is in degrees; if False (default), radians

            Returns
            -------
            Rotation
                Rotation object
            )doc")

        .def_static("random", &simcoon::Rotation::random,
            "Create a uniformly distributed random rotation")

        // Conversion methods
        .def("as_quat",
            [](const simcoon::Rotation& self) {
                return carma::col_to_arr(vec(self.as_quat()));
            },
            R"doc(
            Get rotation as a unit quaternion.

            Returns
            -------
            ndarray
                Quaternion in [qx, qy, qz, qw] format (scalar-last)
            )doc")

        .def("as_matrix",
            [](const simcoon::Rotation& self) {
                return carma::mat_to_arr(mat(self.as_matrix()));
            },
            R"doc(
            Get rotation as a 3x3 rotation matrix.

            Returns
            -------
            ndarray
                3x3 orthogonal rotation matrix
            )doc")

        .def("as_euler",
            [](const simcoon::Rotation& self, const string& conv, bool intrinsic, bool degrees) {
                return carma::col_to_arr(vec(self.as_euler(conv, intrinsic, degrees)));
            },
            py::arg("conv") = "zxz", py::arg("intrinsic") = true, py::arg("degrees") = false,
            R"doc(
            Get rotation as Euler angles.

            Parameters
            ----------
            conv : str, optional
                Euler angle convention. Default is "zxz".
            intrinsic : bool, optional
                If True (default), return intrinsic angles; if False, extrinsic
            degrees : bool, optional
                If True, return angles in degrees; if False (default), radians

            Returns
            -------
            ndarray
                3-element array of Euler angles [psi, theta, phi]
            )doc")

        .def("as_rotvec",
            [](const simcoon::Rotation& self, bool degrees) {
                return carma::col_to_arr(vec(self.as_rotvec(degrees)));
            },
            py::arg("degrees") = false,
            R"doc(
            Get rotation as a rotation vector.

            Parameters
            ----------
            degrees : bool, optional
                If True, magnitude is in degrees; if False (default), radians

            Returns
            -------
            ndarray
                Rotation vector (axis * angle)
            )doc")

        .def("as_QS",
            [](const simcoon::Rotation& self, bool active) {
                return carma::mat_to_arr(mat(self.as_QS(active)));
            },
            py::arg("active") = true,
            R"doc(
            Get 6x6 rotation matrix for stress tensors in Voigt notation.

            Parameters
            ----------
            active : bool, optional
                If True (default), active rotation; if False, passive

            Returns
            -------
            ndarray
                6x6 stress rotation matrix (QS)
            )doc")

        .def("as_QE",
            [](const simcoon::Rotation& self, bool active) {
                return carma::mat_to_arr(mat(self.as_QE(active)));
            },
            py::arg("active") = true,
            R"doc(
            Get 6x6 rotation matrix for strain tensors in Voigt notation.

            Parameters
            ----------
            active : bool, optional
                If True (default), active rotation; if False, passive

            Returns
            -------
            ndarray
                6x6 strain rotation matrix (QE)
            )doc")

        // Apply methods
        .def("apply",
            [](const simcoon::Rotation& self, py::array_t<double> v, bool inverse) {
                vec v_cpp = carma::arr_to_col(v);
                vec result = self.apply(v_cpp, inverse);
                return carma::col_to_arr(result);
            },
            py::arg("v"), py::arg("inverse") = false,
            R"doc(
            Apply rotation to a 3D vector.

            Parameters
            ----------
            v : array_like
                3D vector to rotate
            inverse : bool, optional
                If True, apply inverse rotation. Default is False.

            Returns
            -------
            ndarray
                Rotated vector
            )doc")

        .def("apply_tensor",
            [](const simcoon::Rotation& self, py::array_t<double> m, bool inverse) {
                mat m_cpp = carma::arr_to_mat(m);
                mat result = self.apply_tensor(m_cpp, inverse);
                return carma::mat_to_arr(result);
            },
            py::arg("m"), py::arg("inverse") = false,
            R"doc(
            Apply rotation to a 3x3 tensor (matrix).

            Parameters
            ----------
            m : array_like
                3x3 tensor to rotate
            inverse : bool, optional
                If True, apply inverse rotation. Default is False.

            Returns
            -------
            ndarray
                Rotated tensor: R * m * R^T (or R^T * m * R for inverse)
            )doc")

        .def("apply_stress",
            [](const simcoon::Rotation& self, py::array_t<double> sigma, bool active) {
                vec sigma_cpp = carma::arr_to_col(sigma);
                vec result = self.apply_stress(sigma_cpp, active);
                return carma::col_to_arr(result);
            },
            py::arg("sigma"), py::arg("active") = true,
            R"doc(
            Apply rotation to a stress vector in Voigt notation.

            Parameters
            ----------
            sigma : array_like
                6-component stress vector [s11, s22, s33, s12, s13, s23]
            active : bool, optional
                If True (default), active rotation

            Returns
            -------
            ndarray
                Rotated stress vector
            )doc")

        .def("apply_strain",
            [](const simcoon::Rotation& self, py::array_t<double> epsilon, bool active) {
                vec epsilon_cpp = carma::arr_to_col(epsilon);
                vec result = self.apply_strain(epsilon_cpp, active);
                return carma::col_to_arr(result);
            },
            py::arg("epsilon"), py::arg("active") = true,
            R"doc(
            Apply rotation to a strain vector in Voigt notation.

            Parameters
            ----------
            epsilon : array_like
                6-component strain vector [e11, e22, e33, 2*e12, 2*e13, 2*e23]
            active : bool, optional
                If True (default), active rotation

            Returns
            -------
            ndarray
                Rotated strain vector
            )doc")

        .def("apply_stiffness",
            [](const simcoon::Rotation& self, py::array_t<double> L, bool active) {
                mat L_cpp = carma::arr_to_mat(L);
                mat result = self.apply_stiffness(L_cpp, active);
                return carma::mat_to_arr(result);
            },
            py::arg("L"), py::arg("active") = true,
            R"doc(
            Apply rotation to a 6x6 stiffness matrix.

            Parameters
            ----------
            L : array_like
                6x6 stiffness matrix in Voigt notation
            active : bool, optional
                If True (default), active rotation

            Returns
            -------
            ndarray
                Rotated stiffness matrix: QS * L * QS^T
            )doc")

        .def("apply_compliance",
            [](const simcoon::Rotation& self, py::array_t<double> M, bool active) {
                mat M_cpp = carma::arr_to_mat(M);
                mat result = self.apply_compliance(M_cpp, active);
                return carma::mat_to_arr(result);
            },
            py::arg("M"), py::arg("active") = true,
            R"doc(
            Apply rotation to a 6x6 compliance matrix.

            Parameters
            ----------
            M : array_like
                6x6 compliance matrix in Voigt notation
            active : bool, optional
                If True (default), active rotation

            Returns
            -------
            ndarray
                Rotated compliance matrix: QE * M * QE^T
            )doc")

        // Operations
        .def("inv", &simcoon::Rotation::inv,
            "Get the inverse (conjugate) rotation")

        .def("magnitude", &simcoon::Rotation::magnitude,
            py::arg("degrees") = false,
            R"doc(
            Get the magnitude (rotation angle) of this rotation.

            Parameters
            ----------
            degrees : bool, optional
                If True, return in degrees; if False (default), radians

            Returns
            -------
            float
                Rotation angle
            )doc")

        .def("__mul__", &simcoon::Rotation::operator*,
            py::arg("other"),
            "Compose this rotation with another (self * other)")

        .def("__imul__", &simcoon::Rotation::operator*=,
            py::arg("other"),
            "Compose this rotation with another in-place")

        .def("slerp", &simcoon::Rotation::slerp,
            py::arg("other"), py::arg("t"),
            R"doc(
            Spherical linear interpolation between this rotation and another.

            Parameters
            ----------
            other : Rotation
                Target rotation
            t : float
                Interpolation parameter in [0, 1] (0 = self, 1 = other)

            Returns
            -------
            Rotation
                Interpolated rotation
            )doc")

        .def("equals", &simcoon::Rotation::equals,
            py::arg("other"), py::arg("tol") = 1e-12,
            R"doc(
            Check if this rotation equals another within tolerance.

            Parameters
            ----------
            other : Rotation
                Rotation to compare with
            tol : float, optional
                Tolerance for comparison. Default is 1e-12.

            Returns
            -------
            bool
                True if rotations are equal (or antipodal quaternions)
            )doc")

        .def("is_identity", &simcoon::Rotation::is_identity,
            py::arg("tol") = 1e-12,
            R"doc(
            Check if this rotation is identity (no rotation).

            Parameters
            ----------
            tol : float, optional
                Tolerance for comparison. Default is 1e-12.

            Returns
            -------
            bool
                True if this is the identity rotation
            )doc")

        .def("__repr__",
            [](const simcoon::Rotation& self) {
                vec::fixed<4> q = self.as_quat();
                return "<Rotation: quat=[" +
                    to_string(q(0)) + ", " + to_string(q(1)) + ", " +
                    to_string(q(2)) + ", " + to_string(q(3)) + "]>";
            });
}

} // namespace simpy
