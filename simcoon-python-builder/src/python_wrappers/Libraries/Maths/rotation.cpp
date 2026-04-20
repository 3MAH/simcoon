#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <string>
#include <carma>
#include <armadillo>

#include <simcoon/Simulation/Maths/rotation.hpp>
#include <simcoon/python_wrappers/Libraries/Maths/rotation.hpp>
#include <simcoon/docs/Libraries/Maths/doc_rotation.hpp>

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
        simcoon_docs::CppRotation_class)

        // The only factory method needed — Python Rotation._to_cpp() uses this
        .def_static("from_quat",
            [](py::array_t<double> quat) {
                validate_vector_size(quat, 4, "quat");
                vec q = carma::arr_to_col(quat);
                return simcoon::Rotation::from_quat(q);
            },
            py::arg("quat"),
            simcoon_docs::CppRotation_from_quat)

        // Voigt rotation matrices
        .def("as_voigt_stress_rotation",
            [](const simcoon::Rotation& self, bool active) {
                return carma::mat_to_arr(mat(self.as_voigt_stress_rotation(active)));
            },
            py::arg("active") = true,
            simcoon_docs::as_voigt_stress_rotation)

        .def("as_voigt_strain_rotation",
            [](const simcoon::Rotation& self, bool active) {
                return carma::mat_to_arr(mat(self.as_voigt_strain_rotation(active)));
            },
            py::arg("active") = true,
            simcoon_docs::as_voigt_strain_rotation)

        // Apply methods — the core mechanics operations
        .def("apply_tensor",
            [](const simcoon::Rotation& self, py::array_t<double> m, bool inverse) {
                validate_matrix_size(m, 3, 3, "m");
                mat m_cpp = carma::arr_to_mat(m);
                mat result = self.apply_tensor(m_cpp, inverse);
                return carma::mat_to_arr(result);
            },
            py::arg("m"), py::arg("inverse") = false,
            simcoon_docs::apply_tensor)

        .def("apply_stress",
            [](const simcoon::Rotation& self, py::array_t<double> sigma, bool active) {
                validate_vector_size(sigma, 6, "sigma");
                vec sigma_cpp = carma::arr_to_col(sigma);
                vec result = self.apply_stress(sigma_cpp, active);
                return carma::col_to_arr(result);
            },
            py::arg("sigma"), py::arg("active") = true,
            simcoon_docs::apply_stress)

        .def("apply_strain",
            [](const simcoon::Rotation& self, py::array_t<double> epsilon, bool active) {
                validate_vector_size(epsilon, 6, "epsilon");
                vec epsilon_cpp = carma::arr_to_col(epsilon);
                vec result = self.apply_strain(epsilon_cpp, active);
                return carma::col_to_arr(result);
            },
            py::arg("epsilon"), py::arg("active") = true,
            simcoon_docs::apply_strain)

        .def("apply_stiffness",
            [](const simcoon::Rotation& self, py::array_t<double> L, bool active) {
                validate_matrix_size(L, 6, 6, "L");
                mat L_cpp = carma::arr_to_mat(L);
                mat result = self.apply_stiffness(L_cpp, active);
                return carma::mat_to_arr(result);
            },
            py::arg("L"), py::arg("active") = true,
            simcoon_docs::apply_stiffness)

        .def("apply_compliance",
            [](const simcoon::Rotation& self, py::array_t<double> M, bool active) {
                validate_matrix_size(M, 6, 6, "M");
                mat M_cpp = carma::arr_to_mat(M);
                mat result = self.apply_compliance(M_cpp, active);
                return carma::mat_to_arr(result);
            },
            py::arg("M"), py::arg("active") = true,
            simcoon_docs::apply_compliance)

        .def("apply_strain_concentration",
            [](const simcoon::Rotation& self, py::array_t<double> A, bool active) {
                validate_matrix_size(A, 6, 6, "A");
                mat A_cpp = carma::arr_to_mat(A);
                mat result = self.apply_strain_concentration(A_cpp, active);
                return carma::mat_to_arr(result);
            },
            py::arg("A"), py::arg("active") = true,
            simcoon_docs::apply_strain_concentration)

        .def("apply_stress_concentration",
            [](const simcoon::Rotation& self, py::array_t<double> B, bool active) {
                validate_matrix_size(B, 6, 6, "B");
                mat B_cpp = carma::arr_to_mat(B);
                mat result = self.apply_stress_concentration(B_cpp, active);
                return carma::mat_to_arr(result);
            },
            py::arg("B"), py::arg("active") = true,
            simcoon_docs::apply_stress_concentration)

        ;

    // Batch free functions — avoid Python loops for large N
    m.def("_batch_voigt_stress_rotation",
        [](py::array_t<double> quats, bool active) {
            // quats: (N, 4) numpy → (4, N) arma mat
            auto buf = quats.unchecked<2>();
            int N = buf.shape(0);
            mat q_cpp(4, N);
            for (int n = 0; n < N; n++) {
                for (int j = 0; j < 4; j++) {
                    q_cpp(j, n) = buf(n, j);
                }
            }
            cube result = simcoon::batch_voigt_stress_rotation(q_cpp, active);
            // (6,6,N) → (N,6,6) numpy
            py::array_t<double> out({N, 6, 6});
            auto out_buf = out.mutable_unchecked<3>();
            for (int n = 0; n < N; n++)
                for (int i = 0; i < 6; i++)
                    for (int j = 0; j < 6; j++)
                        out_buf(n, i, j) = result(i, j, n);
            return out;
        },
        py::arg("quats"), py::arg("active") = true,
        "Batch-build (N, 6, 6) stress Voigt rotation matrices from (N, 4) quaternions");

    m.def("_batch_voigt_strain_rotation",
        [](py::array_t<double> quats, bool active) {
            auto buf = quats.unchecked<2>();
            int N = buf.shape(0);
            mat q_cpp(4, N);
            for (int n = 0; n < N; n++) {
                for (int j = 0; j < 4; j++) {
                    q_cpp(j, n) = buf(n, j);
                }
            }
            cube result = simcoon::batch_voigt_strain_rotation(q_cpp, active);
            py::array_t<double> out({N, 6, 6});
            auto out_buf = out.mutable_unchecked<3>();
            for (int n = 0; n < N; n++)
                for (int i = 0; i < 6; i++)
                    for (int j = 0; j < 6; j++)
                        out_buf(n, i, j) = result(i, j, n);
            return out;
        },
        py::arg("quats"), py::arg("active") = true,
        "Batch-build (N, 6, 6) strain Voigt rotation matrices from (N, 4) quaternions");
}

} // namespace simpy
