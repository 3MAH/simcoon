#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <string>
#include <carma>
#include <armadillo>

#include <simcoon/Continuum_mechanics/Functions/tensor.hpp>
#include <simcoon/Simulation/Maths/rotation.hpp>
#include <simcoon/python_wrappers/Libraries/Continuum_mechanics/tensor.hpp>

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
    // Helper: numpy (N,6) → arma mat (6,N)  — via carma + transpose
    mat np2d_to_mat6N(const py::array_t<double>& arr) {
        mat m = carma::arr_to_mat(arr);   // numpy (N,6) → arma (N,6)
        return m.t();                      // → arma (6,N)
    }

    // Helper: arma mat (6,N) → numpy (N,6)  — via transpose + carma
    py::array_t<double> mat6N_to_np2d(const mat& m) {
        mat mt = m.t();                    // arma (6,N) → arma (N,6)
        return carma::mat_to_arr(mt);      // → numpy (N,6)
    }
} // anonymous namespace

void register_tensor(py::module_& m) {

    // VoigtType enum
    py::enum_<simcoon::VoigtType>(m, "VoigtType",
        "Type tag for 2nd-order tensors, determines Voigt conversion factors and rotation rules.")
        .value("stress", simcoon::VoigtType::stress, "Voigt = [s11, s22, s33, s12, s13, s23]")
        .value("strain", simcoon::VoigtType::strain, "Voigt = [e11, e22, e33, 2*e12, 2*e13, 2*e23]")
        .value("generic", simcoon::VoigtType::generic, "Same storage as stress (symmetric tensor)")
        .value("none", simcoon::VoigtType::none, "Non-symmetric tensor, voigt() throws")
        .export_values();

    // Tensor4Type enum
    py::enum_<simcoon::Tensor4Type>(m, "Tensor4Type",
        "Type tag for 4th-order tensors, determines rotation rules.")
        .value("stiffness", simcoon::Tensor4Type::stiffness, "QS * L * QS^T")
        .value("compliance", simcoon::Tensor4Type::compliance, "QE * M * QE^T")
        .value("strain_concentration", simcoon::Tensor4Type::strain_concentration, "QE * A * QS^T")
        .value("stress_concentration", simcoon::Tensor4Type::stress_concentration, "QS * B * QE^T")
        .value("generic", simcoon::Tensor4Type::generic, "QS * C * QS^T (default)")
        .export_values();

    // _CppTensor2 class
    py::class_<simcoon::tensor2>(m, "_CppTensor2",
        R"doc(
        Internal C++ tensor2 backend with type-tagged Voigt convention.

        End users should use ``simcoon.Tensor2`` instead.
        )doc")

        // Constructors
        .def(py::init<>())
        .def(py::init<simcoon::VoigtType>(), py::arg("vtype"))

        .def_static("from_mat",
            [](py::array_t<double> m, simcoon::VoigtType vtype) {
                validate_matrix_size(m, 3, 3, "m");
                mat m_cpp = carma::arr_to_mat(m);
                return simcoon::tensor2(mat::fixed<3,3>(m_cpp), vtype);
            },
            py::arg("m"), py::arg("vtype"),
            "Create tensor2 from a 3x3 numpy array")

        .def_static("from_voigt",
            [](py::array_t<double> v, simcoon::VoigtType vtype) {
                validate_vector_size(v, 6, "v");
                vec v_cpp = carma::arr_to_col(v);
                return simcoon::tensor2::from_voigt(vec::fixed<6>(v_cpp.memptr()), vtype);
            },
            py::arg("v"), py::arg("vtype"),
            "Create tensor2 from a 6-element Voigt vector")

        .def_static("zeros",
            static_cast<simcoon::tensor2 (*)(simcoon::VoigtType)>(&simcoon::tensor2::zeros),
            py::arg("vtype") = simcoon::VoigtType::stress)

        .def_static("identity",
            static_cast<simcoon::tensor2 (*)(simcoon::VoigtType)>(&simcoon::tensor2::identity),
            py::arg("vtype") = simcoon::VoigtType::stress)

        // Properties
        .def_property_readonly("mat",
            [](const simcoon::tensor2& self) {
                return carma::mat_to_arr(mat(self.mat()));
            },
            "3x3 matrix representation as numpy array")

        .def_property_readonly("voigt",
            [](const simcoon::tensor2& self) {
                return carma::col_to_arr(vec(self.voigt()));
            },
            "6-element Voigt vector as numpy array")

        .def_property_readonly("vtype", &simcoon::tensor2::vtype)

        // Methods
        .def("is_symmetric", &simcoon::tensor2::is_symmetric,
            py::arg("tol") = 1e-12)

        .def("rotate",
            [](const simcoon::tensor2& self, const simcoon::Rotation& R, bool active) {
                return self.rotate(R, active);
            },
            py::arg("R"), py::arg("active") = true,
            "Rotate the tensor using the Rotation object")

        .def("push_forward",
            [](const simcoon::tensor2& self, py::array_t<double> F, bool metric) {
                validate_matrix_size(F, 3, 3, "F");
                mat F_cpp = carma::arr_to_mat(F);
                return self.push_forward(F_cpp, metric);
            },
            py::arg("F"), py::arg("metric") = true,
            "Push-forward the tensor via deformation gradient F. metric=True includes J factor.")

        .def("pull_back",
            [](const simcoon::tensor2& self, py::array_t<double> F, bool metric) {
                validate_matrix_size(F, 3, 3, "F");
                mat F_cpp = carma::arr_to_mat(F);
                return self.pull_back(F_cpp, metric);
            },
            py::arg("F"), py::arg("metric") = true,
            "Pull-back the tensor via deformation gradient F. metric=True includes J factor.")

        // Arithmetic
        .def("__add__", &simcoon::tensor2::operator+)
        .def("__sub__", static_cast<simcoon::tensor2 (simcoon::tensor2::*)(const simcoon::tensor2&) const>(&simcoon::tensor2::operator-))
        .def("__neg__", static_cast<simcoon::tensor2 (simcoon::tensor2::*)() const>(&simcoon::tensor2::operator-))
        .def("__mul__", static_cast<simcoon::tensor2 (simcoon::tensor2::*)(double) const>(&simcoon::tensor2::operator*))
        .def("__rmul__", [](const simcoon::tensor2& self, double scalar) {
            return scalar * self;
        })
        .def("__truediv__", static_cast<simcoon::tensor2 (simcoon::tensor2::*)(double) const>(&simcoon::tensor2::operator/))
        .def("__mod__", static_cast<simcoon::tensor2 (simcoon::tensor2::*)(const simcoon::tensor2&) const>(&simcoon::tensor2::operator%))
        .def("__eq__", &simcoon::tensor2::operator==)
        .def("__ne__", &simcoon::tensor2::operator!=)

        // String representation
        .def("__repr__", [](const simcoon::tensor2& self) {
            string vtype_str;
            switch (self.vtype()) {
                case simcoon::VoigtType::stress: vtype_str = "stress"; break;
                case simcoon::VoigtType::strain: vtype_str = "strain"; break;
                case simcoon::VoigtType::generic: vtype_str = "generic"; break;
                case simcoon::VoigtType::none: vtype_str = "none"; break;
            }
            return "Tensor2(vtype=" + vtype_str + ")";
        });

    // _CppTensor4 class
    py::class_<simcoon::tensor4>(m, "_CppTensor4",
        R"doc(
        Internal C++ tensor4 backend with type-tagged rotation dispatch.

        End users should use ``simcoon.Tensor4`` instead.
        )doc")

        // Constructors
        .def(py::init<>())
        .def(py::init<simcoon::Tensor4Type>(), py::arg("type"))

        .def_static("from_mat",
            [](py::array_t<double> m, simcoon::Tensor4Type type) {
                validate_matrix_size(m, 6, 6, "m");
                mat m_cpp = carma::arr_to_mat(m);
                return simcoon::tensor4(mat::fixed<6,6>(m_cpp), type);
            },
            py::arg("m"), py::arg("type"),
            "Create tensor4 from a 6x6 Voigt matrix")

        // Static factories
        .def_static("identity", &simcoon::tensor4::identity,
            py::arg("type") = simcoon::Tensor4Type::stiffness)
        .def_static("volumetric", &simcoon::tensor4::volumetric,
            py::arg("type") = simcoon::Tensor4Type::stiffness)
        .def_static("deviatoric", &simcoon::tensor4::deviatoric,
            py::arg("type") = simcoon::Tensor4Type::stiffness)
        .def_static("identity2", &simcoon::tensor4::identity2,
            py::arg("type") = simcoon::Tensor4Type::stiffness)
        .def_static("deviatoric2", &simcoon::tensor4::deviatoric2,
            py::arg("type") = simcoon::Tensor4Type::stiffness)
        .def_static("zeros", &simcoon::tensor4::zeros,
            py::arg("type") = simcoon::Tensor4Type::stiffness)

        // Properties
        .def_property_readonly("mat",
            [](const simcoon::tensor4& self) {
                return carma::mat_to_arr(mat(self.mat()));
            },
            "6x6 Voigt matrix as numpy array")

        .def_property_readonly("type", &simcoon::tensor4::type)

        // Methods
        .def("contract",
            &simcoon::tensor4::contract,
            py::arg("t"),
            "Contract with a tensor2: result = mat * t.voigt()")

        .def("rotate",
            [](const simcoon::tensor4& self, const simcoon::Rotation& R, bool active) {
                return self.rotate(R, active);
            },
            py::arg("R"), py::arg("active") = true,
            "Rotate the tensor using the Rotation object")

        .def("push_forward",
            [](const simcoon::tensor4& self, py::array_t<double> F, bool metric) {
                validate_matrix_size(F, 3, 3, "F");
                mat F_cpp = carma::arr_to_mat(F);
                return self.push_forward(F_cpp, metric);
            },
            py::arg("F"), py::arg("metric") = true,
            "Push-forward via deformation gradient F. metric=True includes J factor.")

        .def("pull_back",
            [](const simcoon::tensor4& self, py::array_t<double> F, bool metric) {
                validate_matrix_size(F, 3, 3, "F");
                mat F_cpp = carma::arr_to_mat(F);
                return self.pull_back(F_cpp, metric);
            },
            py::arg("F"), py::arg("metric") = true,
            "Pull-back via deformation gradient F. metric=True includes J factor.")

        .def("inverse", &simcoon::tensor4::inverse,
            "Invert the 6x6 Voigt matrix. Type: stiffness <-> compliance.")

        // Arithmetic
        .def("__add__", &simcoon::tensor4::operator+)
        .def("__sub__", static_cast<simcoon::tensor4 (simcoon::tensor4::*)(const simcoon::tensor4&) const>(&simcoon::tensor4::operator-))
        .def("__neg__", static_cast<simcoon::tensor4 (simcoon::tensor4::*)() const>(&simcoon::tensor4::operator-))
        .def("__mul__", static_cast<simcoon::tensor4 (simcoon::tensor4::*)(double) const>(&simcoon::tensor4::operator*))
        .def("__rmul__", [](const simcoon::tensor4& self, double scalar) {
            return scalar * self;
        })
        .def("__truediv__", static_cast<simcoon::tensor4 (simcoon::tensor4::*)(double) const>(&simcoon::tensor4::operator/))
        .def("__mod__", &simcoon::tensor4::operator%)
        .def("__matmul__", [](const simcoon::tensor4& self, const simcoon::tensor2& t) {
            return self.contract(t);
        }, "Contract with tensor2 using @ operator")
        .def("__eq__", &simcoon::tensor4::operator==)
        .def("__ne__", &simcoon::tensor4::operator!=)

        // String representation
        .def("__repr__", [](const simcoon::tensor4& self) {
            string type_str;
            switch (self.type()) {
                case simcoon::Tensor4Type::stiffness: type_str = "stiffness"; break;
                case simcoon::Tensor4Type::compliance: type_str = "compliance"; break;
                case simcoon::Tensor4Type::strain_concentration: type_str = "strain_concentration"; break;
                case simcoon::Tensor4Type::stress_concentration: type_str = "stress_concentration"; break;
                case simcoon::Tensor4Type::generic: type_str = "generic"; break;
            }
            return "Tensor4(type=" + type_str + ")";
        });

    // ================================================================
    // Batch functions — thin wrappers calling core C++ batch functions
    // ================================================================
    // Batch operations
    // Python side prepares:
    //   voigt as (N,6) C-order → converted via np2d_to_mat6N (carma + .t())
    //   cubes as (R,C,N) F-order → zero-copy via carma::arr_to_cube

    m.def("_batch_t2_rotate",
        [](
           py::array_t<double> voigt, simcoon::VoigtType vtype,
           py::array_t<double> rot_matrices, bool active) {
            mat v_cpp = np2d_to_mat6N(voigt);
            cube r_cpp = carma::arr_to_cube<double>(rot_matrices);
            mat result = simcoon::batch_rotate(v_cpp, vtype, r_cpp, active);
            return mat6N_to_np2d(result);
        },
        py::arg("voigt"), py::arg("vtype"), py::arg("rot_matrices"), py::arg("active") = true);

    m.def("_batch_t2_push_forward",
        [](
           py::array_t<double> voigt, simcoon::VoigtType vtype,
           py::array_t<double> F_arr, bool metric) {
            mat v_cpp = np2d_to_mat6N(voigt);
            cube f_cpp = carma::arr_to_cube<double>(F_arr);
            mat result = simcoon::batch_push_forward(v_cpp, vtype, f_cpp, metric);
            return mat6N_to_np2d(result);
        },
        py::arg("voigt"), py::arg("vtype"), py::arg("F"), py::arg("metric") = true);

    m.def("_batch_t2_pull_back",
        [](
           py::array_t<double> voigt, simcoon::VoigtType vtype,
           py::array_t<double> F_arr, bool metric) {
            mat v_cpp = np2d_to_mat6N(voigt);
            cube f_cpp = carma::arr_to_cube<double>(F_arr);
            mat result = simcoon::batch_pull_back(v_cpp, vtype, f_cpp, metric);
            return mat6N_to_np2d(result);
        },
        py::arg("voigt"), py::arg("vtype"), py::arg("F"), py::arg("metric") = true);

    m.def("_batch_t2_mises",
        [](py::array_t<double> voigt, simcoon::VoigtType vtype) {
            mat v_cpp = np2d_to_mat6N(voigt);
            vec result = simcoon::batch_mises(v_cpp, vtype);
            return carma::col_to_arr(result);
        },
        py::arg("voigt"), py::arg("vtype"));

    m.def("_batch_t2_trace",
        [](py::array_t<double> voigt, simcoon::VoigtType vtype) {
            mat v_cpp = np2d_to_mat6N(voigt);
            vec result = simcoon::batch_trace(v_cpp, vtype);
            return carma::col_to_arr(result);
        },
        py::arg("voigt"), py::arg("vtype"));

    m.def("_batch_t4_contract",
        [](
           py::array_t<double> t4_arr, simcoon::Tensor4Type t4type,
           py::array_t<double> t2_arr, simcoon::VoigtType t2_vtype) {
            cube t4_cpp = carma::arr_to_cube<double>(t4_arr);
            mat t2_cpp = np2d_to_mat6N(t2_arr);
            mat result = simcoon::batch_contract(t4_cpp, t4type, t2_cpp, t2_vtype);
            simcoon::VoigtType out_vtype = simcoon::infer_contraction_vtype(t4type);
            return py::make_tuple(mat6N_to_np2d(result), out_vtype);
        },
        py::arg("t4"), py::arg("t4type"), py::arg("t2"), py::arg("t2_vtype"));

    m.def("_batch_t4_rotate",
        [](
           py::array_t<double> t4_arr, simcoon::Tensor4Type t4type,
           py::array_t<double> rot_matrices, bool active) {
            cube t4_cpp = carma::arr_to_cube<double>(t4_arr);
            cube r_cpp = carma::arr_to_cube<double>(rot_matrices);
            cube result = simcoon::batch_rotate_t4(t4_cpp, t4type, r_cpp, active);
            return carma::cube_to_arr(result, false);
        },
        py::arg("t4"), py::arg("t4type"), py::arg("rot_matrices"), py::arg("active") = true);

    m.def("_batch_t4_push_forward",
        [](
           py::array_t<double> t4_arr, simcoon::Tensor4Type t4type,
           py::array_t<double> F_arr, bool metric) {
            cube t4_cpp = carma::arr_to_cube<double>(t4_arr);
            cube f_cpp = carma::arr_to_cube<double>(F_arr);
            cube result = simcoon::batch_push_forward_t4(t4_cpp, t4type, f_cpp, metric);
            return carma::cube_to_arr(result, false);
        },
        py::arg("t4"), py::arg("t4type"), py::arg("F"), py::arg("metric") = true);

    m.def("_batch_t4_pull_back",
        [](
           py::array_t<double> t4_arr, simcoon::Tensor4Type t4type,
           py::array_t<double> F_arr, bool metric) {
            cube t4_cpp = carma::arr_to_cube<double>(t4_arr);
            cube f_cpp = carma::arr_to_cube<double>(F_arr);
            cube result = simcoon::batch_pull_back_t4(t4_cpp, t4type, f_cpp, metric);
            return carma::cube_to_arr(result, false);
        },
        py::arg("t4"), py::arg("t4type"), py::arg("F"), py::arg("metric") = true);

    m.def("_batch_t4_inverse",
        [](
           py::array_t<double> t4_arr, simcoon::Tensor4Type t4type) {
            cube t4_cpp = carma::arr_to_cube<double>(t4_arr);
            cube result = simcoon::batch_inverse_t4(t4_cpp, t4type);
            simcoon::Tensor4Type inv_type = simcoon::infer_inverse_type(t4type);
            return py::make_tuple(carma::cube_to_arr(result, false), inv_type);
        },
        py::arg("t4"), py::arg("t4type"));

    // Module-level free functions
    m.def("_dyadic",
        [](const simcoon::tensor2& a, const simcoon::tensor2& b) {
            return simcoon::dyadic(a, b);
        },
        py::arg("a"), py::arg("b"),
        "Dyadic (outer) product of two tensor2 objects, returns a tensor4 (stiffness type).");

    m.def("_auto_dyadic",
        [](const simcoon::tensor2& a) {
            return simcoon::auto_dyadic(a);
        },
        py::arg("a"),
        "Dyadic (outer) product of a tensor2 with itself, returns a tensor4 (stiffness type).");

    m.def("_sym_dyadic",
        [](const simcoon::tensor2& a, const simcoon::tensor2& b) {
            return simcoon::sym_dyadic(a, b);
        },
        py::arg("a"), py::arg("b"),
        "Symmetric Voigt dyadic product of two tensor2 objects, returns a tensor4 (stiffness type).");

    m.def("_auto_sym_dyadic",
        [](const simcoon::tensor2& a) {
            return simcoon::auto_sym_dyadic(a);
        },
        py::arg("a"),
        "Symmetric Voigt dyadic product of a tensor2 with itself, returns a tensor4 (stiffness type).");
}

} // namespace simpy
