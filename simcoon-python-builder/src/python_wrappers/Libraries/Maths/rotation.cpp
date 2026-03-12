#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <string>
#include <carma>
#include <armadillo>

#include <simcoon/Simulation/Maths/rotation.hpp>
#include <simcoon/python_wrappers/Libraries/Maths/rotation.hpp>

using namespace std;
using namespace arma;
namespace py = pybind11;

namespace simpy {

namespace {
    void validate_vector_size(const py::array_t<double>& arr, size_t expected, const string& name) {
        if (arr.ndim() != 1) {
            throw invalid_argument(name + " must be a 1D array, got " + to_string(arr.ndim()) + "D");
        }
        if (static_cast<size_t>(arr.size()) != expected) {
            throw invalid_argument(name + " must have " + to_string(expected) +
                " elements, got " + to_string(arr.size()));
        }
    }

    void validate_matrix_size(const py::array_t<double>& arr, size_t rows, size_t cols, const string& name) {
        if (arr.ndim() != 2) {
            throw invalid_argument(name + " must be a 2D array, got " + to_string(arr.ndim()) + "D");
        }
        auto shape = arr.shape();
        if (static_cast<size_t>(shape[0]) != rows || static_cast<size_t>(shape[1]) != cols) {
            throw invalid_argument(name + " must have shape (" + to_string(rows) + ", " +
                to_string(cols) + "), got (" + to_string(shape[0]) + ", " + to_string(shape[1]) + ")");
        }
    }
} // anonymous namespace

void register_rotation(py::module_& m) {
    py::class_<simcoon::Rotation>(m, "_CppRotation",
        R"doc(
        Internal C++ rotation backend using unit quaternions (scalar-last).

        End users should use ``simcoon.Rotation`` instead, which inherits from
        ``scipy.spatial.transform.Rotation`` and delegates mechanics operations
        to this class.
        )doc")

        // The only factory method needed — Python Rotation._to_cpp() uses this
        .def_static("from_quat",
            [](py::array_t<double> quat) {
                validate_vector_size(quat, 4, "quat");
                vec q = carma::arr_to_col(quat);
                return simcoon::Rotation::from_quat(q);
            },
            py::arg("quat"),
            "Create rotation from quaternion [qx, qy, qz, qw] (scalar-last)")

        // Voigt rotation matrices
        .def("as_voigt_stress_rotation",
            [](const simcoon::Rotation& self, bool active) {
                return carma::mat_to_arr(mat(self.as_voigt_stress_rotation(active)));
            },
            py::arg("active") = true,
            "Get 6x6 rotation matrix for stress tensors in Voigt notation")

        .def("as_voigt_strain_rotation",
            [](const simcoon::Rotation& self, bool active) {
                return carma::mat_to_arr(mat(self.as_voigt_strain_rotation(active)));
            },
            py::arg("active") = true,
            "Get 6x6 rotation matrix for strain tensors in Voigt notation")

        // Apply methods — the core mechanics operations
        .def("apply_tensor",
            [](const simcoon::Rotation& self, py::array_t<double> m, bool inverse) {
                validate_matrix_size(m, 3, 3, "m");
                mat m_cpp = carma::arr_to_mat(m);
                mat result = self.apply_tensor(m_cpp, inverse);
                return carma::mat_to_arr(result);
            },
            py::arg("m"), py::arg("inverse") = false,
            "Apply rotation to a 3x3 tensor: R * m * R^T")

        .def("apply_stress",
            [](const simcoon::Rotation& self, py::array_t<double> sigma, bool active) {
                validate_vector_size(sigma, 6, "sigma");
                vec sigma_cpp = carma::arr_to_col(sigma);
                vec result = self.apply_stress(sigma_cpp, active);
                return carma::col_to_arr(result);
            },
            py::arg("sigma"), py::arg("active") = true,
            "Apply rotation to a 6-component stress vector in Voigt notation")

        .def("apply_strain",
            [](const simcoon::Rotation& self, py::array_t<double> epsilon, bool active) {
                validate_vector_size(epsilon, 6, "epsilon");
                vec epsilon_cpp = carma::arr_to_col(epsilon);
                vec result = self.apply_strain(epsilon_cpp, active);
                return carma::col_to_arr(result);
            },
            py::arg("epsilon"), py::arg("active") = true,
            "Apply rotation to a 6-component strain vector in Voigt notation")

        .def("apply_stiffness",
            [](const simcoon::Rotation& self, py::array_t<double> L, bool active) {
                validate_matrix_size(L, 6, 6, "L");
                mat L_cpp = carma::arr_to_mat(L);
                mat result = self.apply_stiffness(L_cpp, active);
                return carma::mat_to_arr(result);
            },
            py::arg("L"), py::arg("active") = true,
            "Apply rotation to a 6x6 stiffness matrix")

        .def("apply_compliance",
            [](const simcoon::Rotation& self, py::array_t<double> M, bool active) {
                validate_matrix_size(M, 6, 6, "M");
                mat M_cpp = carma::arr_to_mat(M);
                mat result = self.apply_compliance(M_cpp, active);
                return carma::mat_to_arr(result);
            },
            py::arg("M"), py::arg("active") = true,
            "Apply rotation to a 6x6 compliance matrix")

        .def("apply_strain_concentration",
            [](const simcoon::Rotation& self, py::array_t<double> A, bool active) {
                validate_matrix_size(A, 6, 6, "A");
                mat A_cpp = carma::arr_to_mat(A);
                mat result = self.apply_strain_concentration(A_cpp, active);
                return carma::mat_to_arr(result);
            },
            py::arg("A"), py::arg("active") = true,
            "Apply rotation to a 6x6 strain concentration tensor")

        .def("apply_stress_concentration",
            [](const simcoon::Rotation& self, py::array_t<double> B, bool active) {
                validate_matrix_size(B, 6, 6, "B");
                mat B_cpp = carma::arr_to_mat(B);
                mat result = self.apply_stress_concentration(B_cpp, active);
                return carma::mat_to_arr(result);
            },
            py::arg("B"), py::arg("active") = true,
            "Apply rotation to a 6x6 stress concentration tensor")

        ;
}

} // namespace simpy
