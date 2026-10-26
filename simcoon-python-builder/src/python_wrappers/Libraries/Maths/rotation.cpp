#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <cstring>
#include <map>
#include <string>
#include <simcoon/python_wrappers/arma_to_numpy.hpp>
#include <simcoon/python_wrappers/numpy_to_arma.hpp>
#include <armadillo>

#include <simcoon/Simulation/Maths/rotation.hpp>
#include <simcoon/Continuum_mechanics/Functions/tensor.hpp>
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

    // A 2-D array as a cube of one slice (copy), a 3-D array as a cube (borrowed when
    // F-contiguous, copied otherwise). Returned as a prvalue so that the view stays a view.
    cube as_cube(const py::array_t<double>& arr) {
        if (arr.ndim() == 3) {
            return simpy::numpy_to_arma::arr_to_cube_view(arr);
        }
        const mat m = simpy::numpy_to_arma::arr_to_mat(arr);
        return cube(m.memptr(), m.n_rows, m.n_cols, 1);
    }

    void check_rotation_count(const cube& R, uword N) {
        if (R.n_slices != 1 && R.n_slices != N) {
            throw invalid_argument("rotate_voigt: DR holds " + to_string(R.n_slices)
                + " rotations, expected 1 or " + to_string(N) + " (one per tensor)");
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
                vec q = simpy::numpy_to_arma::arr_to_col(quat);
                return simcoon::Rotation::from_quat(q);
            },
            py::arg("quat"),
            simcoon_docs::CppRotation_from_quat)

        // Voigt rotation matrices
        .def("as_voigt_stress_rotation",
            [](const simcoon::Rotation& self, bool active) {
                return simpy::arma_to_numpy::mat_to_arr(mat(self.as_voigt_stress_rotation(active)));
            },
            py::arg("active") = true,
            simcoon_docs::as_voigt_stress_rotation)

        .def("as_voigt_strain_rotation",
            [](const simcoon::Rotation& self, bool active) {
                return simpy::arma_to_numpy::mat_to_arr(mat(self.as_voigt_strain_rotation(active)));
            },
            py::arg("active") = true,
            simcoon_docs::as_voigt_strain_rotation)

        // Apply methods — the core mechanics operations
        .def("apply_tensor",
            [](const simcoon::Rotation& self, py::array_t<double> m, bool inverse) {
                validate_matrix_size(m, 3, 3, "m");
                mat m_cpp = simpy::numpy_to_arma::arr_to_mat(m);
                mat result = self.apply_tensor(m_cpp, inverse);
                return simpy::arma_to_numpy::mat_to_arr(result);
            },
            py::arg("m"), py::arg("inverse") = false,
            simcoon_docs::apply_tensor)

        .def("apply_stress",
            [](const simcoon::Rotation& self, py::array_t<double> sigma, bool active) {
                validate_vector_size(sigma, 6, "sigma");
                vec sigma_cpp = simpy::numpy_to_arma::arr_to_col(sigma);
                vec result = self.apply_stress(sigma_cpp, active);
                return simpy::arma_to_numpy::col_to_arr(result);
            },
            py::arg("sigma"), py::arg("active") = true,
            simcoon_docs::apply_stress)

        .def("apply_strain",
            [](const simcoon::Rotation& self, py::array_t<double> epsilon, bool active) {
                validate_vector_size(epsilon, 6, "epsilon");
                vec epsilon_cpp = simpy::numpy_to_arma::arr_to_col(epsilon);
                vec result = self.apply_strain(epsilon_cpp, active);
                return simpy::arma_to_numpy::col_to_arr(result);
            },
            py::arg("epsilon"), py::arg("active") = true,
            simcoon_docs::apply_strain)

        .def("apply_stiffness",
            [](const simcoon::Rotation& self, py::array_t<double> L, bool active) {
                validate_matrix_size(L, 6, 6, "L");
                mat L_cpp = simpy::numpy_to_arma::arr_to_mat(L);
                mat result = self.apply_stiffness(L_cpp, active);
                return simpy::arma_to_numpy::mat_to_arr(result);
            },
            py::arg("L"), py::arg("active") = true,
            simcoon_docs::apply_stiffness)

        .def("apply_compliance",
            [](const simcoon::Rotation& self, py::array_t<double> M, bool active) {
                validate_matrix_size(M, 6, 6, "M");
                mat M_cpp = simpy::numpy_to_arma::arr_to_mat(M);
                mat result = self.apply_compliance(M_cpp, active);
                return simpy::arma_to_numpy::mat_to_arr(result);
            },
            py::arg("M"), py::arg("active") = true,
            simcoon_docs::apply_compliance)

        .def("apply_strain_concentration",
            [](const simcoon::Rotation& self, py::array_t<double> A, bool active) {
                validate_matrix_size(A, 6, 6, "A");
                mat A_cpp = simpy::numpy_to_arma::arr_to_mat(A);
                mat result = self.apply_strain_concentration(A_cpp, active);
                return simpy::arma_to_numpy::mat_to_arr(result);
            },
            py::arg("A"), py::arg("active") = true,
            simcoon_docs::apply_strain_concentration)

        .def("apply_stress_concentration",
            [](const simcoon::Rotation& self, py::array_t<double> B, bool active) {
                validate_matrix_size(B, 6, 6, "B");
                mat B_cpp = simpy::numpy_to_arma::arr_to_mat(B);
                mat result = self.apply_stress_concentration(B_cpp, active);
                return simpy::arma_to_numpy::mat_to_arr(result);
            },
            py::arg("B"), py::arg("active") = true,
            simcoon_docs::apply_stress_concentration)

        .def("dR_drotvec",
            [](const simcoon::Rotation& self) {
                return simpy::arma_to_numpy::cube_to_arr(self.dR_drotvec());
            },
            simcoon_docs::dR_drotvec)

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
            cube result;
            {
                // parallel loop in libsimcoon: release the GIL (see parallel_nogil.hpp)
                py::gil_scoped_release release;
                result = simcoon::batch_voigt_stress_rotation(q_cpp, active);
            }
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
            cube result;
            {
                py::gil_scoped_release release;   // parallel loop in libsimcoon
                result = simcoon::batch_voigt_strain_rotation(q_cpp, active);
            }
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

    // Batch rotation of Voigt tensors by rotation matrices, on the (6,N) / (6,6,N) layout of umat
    // and objective_rate: no scipy Rotation, no transposes, parallel loop without the GIL.
    m.def("rotate_voigt",
        [](py::array_t<double> values, py::array_t<double> DR, const string& kind, bool active,
           unsigned int n_threads) -> py::array_t<double> {
            static const map<string, simcoon::Tensor2Type> tensor2_kinds = {
                {"stress", simcoon::Tensor2Type::stress}, {"strain", simcoon::Tensor2Type::strain}};
            static const map<string, simcoon::Tensor4Type> tensor4_kinds = {
                {"stiffness", simcoon::Tensor4Type::stiffness},
                {"compliance", simcoon::Tensor4Type::compliance},
                {"strain_concentration", simcoon::Tensor4Type::strain_concentration},
                {"stress_concentration", simcoon::Tensor4Type::stress_concentration}};
            // DR: one (3,3) rotation for every tensor, or a (3,3,N_r) cube
            if ((DR.ndim() != 2 && DR.ndim() != 3) || DR.shape(0) != 3 || DR.shape(1) != 3) {
                throw invalid_argument("rotate_voigt: DR must have shape (3, 3) or (3, 3, N)");
            }
            const cube R = as_cube(DR);

            const auto it2 = tensor2_kinds.find(kind);
            if (it2 != tensor2_kinds.end()) {
                // (6,) one vector, (6,N) one vector per column
                if ((values.ndim() != 1 && values.ndim() != 2) || values.shape(0) != 6) {
                    throw invalid_argument("rotate_voigt: " + kind + " values must have shape (6,) or (6, N)");
                }
                const bool single = values.ndim() == 1;
                const mat v = single ? mat(simpy::numpy_to_arma::arr_to_col(values))
                                     : simpy::numpy_to_arma::arr_to_mat_view(values);
                check_rotation_count(R, v.n_cols);
                mat result;
                {
                    py::gil_scoped_release release;
                    result = simcoon::batch_rotate(v, it2->second, R, active, n_threads);
                }
                if (single) {
                    py::array_t<double> out(6);
                    std::memcpy(out.mutable_data(), result.memptr(), 6 * sizeof(double));
                    return out;
                }
                return simpy::arma_to_numpy::mat_to_arr(std::move(result));
            }

            const auto it4 = tensor4_kinds.find(kind);
            if (it4 == tensor4_kinds.end()) {
                throw invalid_argument("rotate_voigt: unknown kind '" + kind + "'. Valid: stress, strain, "
                                       "stiffness, compliance, strain_concentration, stress_concentration");
            }
            // (6,6) one tensor, (6,6,N) one tensor per slice
            if ((values.ndim() != 2 && values.ndim() != 3) || values.shape(0) != 6 || values.shape(1) != 6) {
                throw invalid_argument("rotate_voigt: " + kind + " values must have shape (6, 6) or (6, 6, N)");
            }
            const bool single = values.ndim() == 2;
            const cube t4 = as_cube(values);
            check_rotation_count(R, t4.n_slices);
            cube result;
            {
                py::gil_scoped_release release;
                result = simcoon::batch_rotate(t4, it4->second, R, active, n_threads);
            }
            if (single) {
                return simpy::arma_to_numpy::mat_to_arr(mat(result.slice(0)));
            }
            return simpy::arma_to_numpy::cube_to_arr(std::move(result));
        },
        py::arg("values"), py::arg("DR"), py::arg("kind"), py::arg("active") = true,
        py::arg("n_threads") = 4,
        simcoon_docs::rotate_voigt);

    m.def("dR_drotvec",
        [](py::array_t<double> rotvec) {
            validate_vector_size(rotvec, 3, "rotvec");
            auto r = rotvec.unchecked<1>();
            vec::fixed<3> omega = {r(0), r(1), r(2)};
            return simpy::arma_to_numpy::cube_to_arr(simcoon::dR_drotvec(omega));
        },
        py::arg("rotvec"),
        simcoon_docs::dR_drotvec_free);
}

} // namespace simpy
