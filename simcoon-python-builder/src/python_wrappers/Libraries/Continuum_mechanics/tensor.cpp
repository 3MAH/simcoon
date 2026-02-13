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
            &simcoon::tensor2::zeros,
            py::arg("vtype") = simcoon::VoigtType::stress)

        .def_static("identity",
            &simcoon::tensor2::identity,
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
            [](const simcoon::tensor2& self, py::array_t<double> F) {
                validate_matrix_size(F, 3, 3, "F");
                mat F_cpp = carma::arr_to_mat(F);
                return self.push_forward(mat::fixed<3,3>(F_cpp));
            },
            py::arg("F"),
            "Push-forward the tensor via deformation gradient F")

        .def("pull_back",
            [](const simcoon::tensor2& self, py::array_t<double> F) {
                validate_matrix_size(F, 3, 3, "F");
                mat F_cpp = carma::arr_to_mat(F);
                return self.pull_back(mat::fixed<3,3>(F_cpp));
            },
            py::arg("F"),
            "Pull-back the tensor via deformation gradient F")

        // Arithmetic
        .def("__add__", &simcoon::tensor2::operator+)
        .def("__sub__", &simcoon::tensor2::operator-)
        .def("__mul__", &simcoon::tensor2::operator*)
        .def("__rmul__", [](const simcoon::tensor2& self, double scalar) {
            return scalar * self;
        })
        .def("__eq__", &simcoon::tensor2::operator==)

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
            [](const simcoon::tensor4& self, py::array_t<double> F) {
                validate_matrix_size(F, 3, 3, "F");
                mat F_cpp = carma::arr_to_mat(F);
                return self.push_forward(mat::fixed<3,3>(F_cpp));
            },
            py::arg("F"),
            "Push-forward via deformation gradient F")

        .def("pull_back",
            [](const simcoon::tensor4& self, py::array_t<double> F) {
                validate_matrix_size(F, 3, 3, "F");
                mat F_cpp = carma::arr_to_mat(F);
                return self.pull_back(mat::fixed<3,3>(F_cpp));
            },
            py::arg("F"),
            "Pull-back via deformation gradient F")

        // Arithmetic
        .def("__add__", &simcoon::tensor4::operator+)
        .def("__sub__", &simcoon::tensor4::operator-)
        .def("__mul__", &simcoon::tensor4::operator*)
        .def("__rmul__", [](const simcoon::tensor4& self, double scalar) {
            return scalar * self;
        })
        .def("__matmul__", [](const simcoon::tensor4& self, const simcoon::tensor2& t) {
            return self.contract(t);
        }, "Contract with tensor2 using @ operator")
        .def("__eq__", &simcoon::tensor4::operator==)

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
}

} // namespace simpy
