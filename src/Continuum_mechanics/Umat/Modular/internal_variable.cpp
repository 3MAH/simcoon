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

/**
 * @file internal_variable.cpp
 * @brief Implementation of InternalVariable class
 */

#include <simcoon/Continuum_mechanics/Umat/Modular/internal_variable.hpp>
#include <simcoon/Simulation/Maths/rotation.hpp>

namespace simcoon {

// ========== Constructors ==========

InternalVariable::InternalVariable(const std::string& name, double init)
    : name_(name)
    , type_(IVarType::SCALAR)
    , scalar_value_(init)
    , scalar_start_(init)
    , vec_value_()
    , vec_start_()
    , mat_value_()
    , mat_start_()
    , is_objective_(true)   // scalars are always objective (frame-invariant)
    , vtype_(VoigtType::strain)
    , t4type_(Tensor4Type::stiffness)
    , statev_offset_(0)
{
}

InternalVariable::InternalVariable(const std::string& name, const arma::vec& init,
                                    bool objective, VoigtType vtype)
    : name_(name)
    , type_(IVarType::VECTOR_6)
    , scalar_value_(0.0)
    , scalar_start_(0.0)
    , vec_value_(init.n_elem == 6 ? init : arma::zeros(6))
    , vec_start_(init.n_elem == 6 ? init : arma::zeros(6))
    , mat_value_()
    , mat_start_()
    , is_objective_(objective)
    , vtype_(vtype)
    , t4type_(Tensor4Type::stiffness)
    , statev_offset_(0)
{
    if (init.n_elem != 6) {
        throw std::runtime_error("InternalVariable: vector must have 6 components, got " +
                                std::to_string(init.n_elem));
    }
}

InternalVariable::InternalVariable(const std::string& name, const arma::mat& init, bool objective,
                                    Tensor4Type t4type)
    : name_(name)
    , type_(IVarType::MATRIX_6x6)
    , scalar_value_(0.0)
    , scalar_start_(0.0)
    , vec_value_()
    , vec_start_()
    , mat_value_(init.n_rows == 6 && init.n_cols == 6 ? init : arma::zeros(6, 6))
    , mat_start_(init.n_rows == 6 && init.n_cols == 6 ? init : arma::zeros(6, 6))
    , is_objective_(objective)
    , vtype_(VoigtType::strain)
    , t4type_(t4type)
    , statev_offset_(0)
{
    if (init.n_rows != 6 || init.n_cols != 6) {
        throw std::runtime_error("InternalVariable: matrix must be 6x6, got " +
                                std::to_string(init.n_rows) + "x" + std::to_string(init.n_cols));
    }
}

// ========== Type Information ==========

unsigned int InternalVariable::size() const {
    switch (type_) {
        case IVarType::SCALAR:
            return 1;
        case IVarType::VECTOR_6:
            return 6;
        case IVarType::MATRIX_6x6:
            return 36;
    }
    return 0;  // Should never reach here
}

// ========== Value Accessors ==========

double& InternalVariable::scalar() {
    if (type_ != IVarType::SCALAR) {
        throw std::runtime_error("InternalVariable '" + name_ + "': attempted to access scalar on non-scalar type");
    }
    return scalar_value_;
}

const double& InternalVariable::scalar() const {
    if (type_ != IVarType::SCALAR) {
        throw std::runtime_error("InternalVariable '" + name_ + "': attempted to access scalar on non-scalar type");
    }
    return scalar_value_;
}

arma::vec& InternalVariable::raw_voigt() {
    if (type_ != IVarType::VECTOR_6) {
        throw std::runtime_error("InternalVariable '" + name_ + "': attempted to access vec on non-vector type");
    }
    return vec_value_;
}

const arma::vec& InternalVariable::raw_voigt() const {
    if (type_ != IVarType::VECTOR_6) {
        throw std::runtime_error("InternalVariable '" + name_ + "': attempted to access vec on non-vector type");
    }
    return vec_value_;
}

arma::mat& InternalVariable::raw_mat() {
    if (type_ != IVarType::MATRIX_6x6) {
        throw std::runtime_error("InternalVariable '" + name_ + "': attempted to access mat on non-matrix type");
    }
    return mat_value_;
}

const arma::mat& InternalVariable::raw_mat() const {
    if (type_ != IVarType::MATRIX_6x6) {
        throw std::runtime_error("InternalVariable '" + name_ + "': attempted to access mat on non-matrix type");
    }
    return mat_value_;
}

// ========== Start Value Accessors ==========

double& InternalVariable::scalar_start() {
    if (type_ != IVarType::SCALAR) {
        throw std::runtime_error("InternalVariable '" + name_ + "': attempted to access scalar_start on non-scalar type");
    }
    return scalar_start_;
}

const double& InternalVariable::scalar_start() const {
    if (type_ != IVarType::SCALAR) {
        throw std::runtime_error("InternalVariable '" + name_ + "': attempted to access scalar_start on non-scalar type");
    }
    return scalar_start_;
}

arma::vec& InternalVariable::raw_voigt_start() {
    if (type_ != IVarType::VECTOR_6) {
        throw std::runtime_error("InternalVariable '" + name_ + "': attempted to access vec_start on non-vector type");
    }
    return vec_start_;
}

const arma::vec& InternalVariable::raw_voigt_start() const {
    if (type_ != IVarType::VECTOR_6) {
        throw std::runtime_error("InternalVariable '" + name_ + "': attempted to access vec_start on non-vector type");
    }
    return vec_start_;
}

arma::mat& InternalVariable::raw_mat_start() {
    if (type_ != IVarType::MATRIX_6x6) {
        throw std::runtime_error("InternalVariable '" + name_ + "': attempted to access mat_start on non-matrix type");
    }
    return mat_start_;
}

const arma::mat& InternalVariable::raw_mat_start() const {
    if (type_ != IVarType::MATRIX_6x6) {
        throw std::runtime_error("InternalVariable '" + name_ + "': attempted to access mat_start on non-matrix type");
    }
    return mat_start_;
}

// ========== Increment Computation ==========

double InternalVariable::delta_scalar() const {
    if (type_ != IVarType::SCALAR) {
        throw std::runtime_error("InternalVariable '" + name_ + "': attempted delta_scalar on non-scalar type");
    }
    return scalar_value_ - scalar_start_;
}

arma::vec InternalVariable::delta_vec() const {
    if (type_ != IVarType::VECTOR_6) {
        throw std::runtime_error("InternalVariable '" + name_ + "': attempted delta_vec on non-vector type");
    }
    return vec_value_ - vec_start_;
}

// ========== State Management ==========

void InternalVariable::to_start() {
    switch (type_) {
        case IVarType::SCALAR:
            scalar_start_ = scalar_value_;
            break;
        case IVarType::VECTOR_6:
            vec_start_ = vec_value_;
            break;
        case IVarType::MATRIX_6x6:
            mat_start_ = mat_value_;
            break;
    }
}

void InternalVariable::rotate(const arma::mat& DR) {
    rotate(Rotation::from_matrix(DR));
}

void InternalVariable::rotate(const Rotation& R) {
    if (!is_objective_) {
        return;
    }
    switch (type_) {
        case IVarType::SCALAR:
            break;
        case IVarType::VECTOR_6: {
            // Route through tensor2 so the rotation kernel (apply_strain vs
            // apply_stress, and the future generic path) is chosen once, inside
            // tensor2::rotate, from the variable's authoritative VoigtType.
            const tensor2 t = tensor2::from_voigt(vec_value_, vtype_).rotate(R);
            vec_value_ = t.to_arma_voigt();
            break;
        }
        case IVarType::MATRIX_6x6:
            // Typed congruence: correct for stiffness- AND compliance-like
            // storage (the Mandel rotation is one orthogonal congruence; the
            // Tensor4Type carries the engineering factors).
            mat_value_ = tensor4(mat_value_, t4type_).rotate(R).mat();
            break;
    }
}

// ========== Tensor-typed views ==========

tensor2 InternalVariable::as_tensor2(VoigtType vtype) const {
    if (type_ != IVarType::VECTOR_6) {
        throw std::runtime_error("InternalVariable::as_tensor2: variable '" + name_
                                 + "' is not VECTOR_6");
    }
    return tensor2::from_voigt(vec_value_, vtype);
}

tensor2 InternalVariable::as_tensor2_start(VoigtType vtype) const {
    if (type_ != IVarType::VECTOR_6) {
        throw std::runtime_error("InternalVariable::as_tensor2_start: variable '" + name_
                                 + "' is not VECTOR_6");
    }
    return tensor2::from_voigt(vec_start_, vtype);
}

tensor4 InternalVariable::as_tensor4(Tensor4Type t4type) const {
    if (type_ != IVarType::MATRIX_6x6) {
        throw std::runtime_error("InternalVariable::as_tensor4: variable '" + name_
                                 + "' is not MATRIX_6x6");
    }
    return tensor4(mat_value_, t4type);
}

tensor4 InternalVariable::as_tensor4_start(Tensor4Type t4type) const {
    if (type_ != IVarType::MATRIX_6x6) {
        throw std::runtime_error("InternalVariable::as_tensor4_start: variable '" + name_
                                 + "' is not MATRIX_6x6");
    }
    return tensor4(mat_start_, t4type);
}

void InternalVariable::set_tensor2(const tensor2& t) {
    if (type_ != IVarType::VECTOR_6) {
        throw std::runtime_error("InternalVariable::set_tensor2: variable '" + name_
                                 + "' is not VECTOR_6");
    }
    // Re-express in this variable's own convention via the 3x3 (convention-
    // free): protects against a vtype mismatch between t and the storage.
    vec_value_ = tensor2(t.mat(), vtype_).voigt();
}

void InternalVariable::set_tensor4(const tensor4& t) {
    if (type_ != IVarType::MATRIX_6x6) {
        throw std::runtime_error("InternalVariable::set_tensor4: variable '" + name_
                                 + "' is not MATRIX_6x6");
    }
    // Re-express in this variable's own convention via the Mandel form
    // (convention-free): protects against a Tensor4Type mismatch.
    mat_value_ = tensor4::from_mandel(t.mandel(), t4type_).mat();
}

// ========== Serialization ==========

void InternalVariable::pack(arma::vec& statev) const {
    const unsigned int offset = statev_offset_;
    if (offset + size() > statev.n_elem) {
        throw std::runtime_error("InternalVariable '" + name_ + "': pack out of bounds (offset "
                                 + std::to_string(offset) + " + size " + std::to_string(size())
                                 + " > nstatev " + std::to_string(statev.n_elem) + ")");
    }
    switch (type_) {
        case IVarType::SCALAR:
            statev(offset) = scalar_value_;
            break;
        case IVarType::VECTOR_6:
            statev.subvec(offset, offset + 5) = vec_value_;
            break;
        case IVarType::MATRIX_6x6:
            for (unsigned int i = 0; i < 6; ++i) {
                for (unsigned int j = 0; j < 6; ++j) {
                    statev(offset + i * 6 + j) = mat_value_(i, j);
                }
            }
            break;
    }
}

void InternalVariable::unpack(const arma::vec& statev) {
    const unsigned int offset = statev_offset_;
    if (offset + size() > statev.n_elem) {
        throw std::runtime_error("InternalVariable '" + name_ + "': unpack out of bounds (offset "
                                 + std::to_string(offset) + " + size " + std::to_string(size())
                                 + " > nstatev " + std::to_string(statev.n_elem) + ")");
    }
    switch (type_) {
        case IVarType::SCALAR:
            scalar_value_ = statev(offset);
            scalar_start_ = scalar_value_;
            break;
        case IVarType::VECTOR_6:
            vec_value_ = statev.subvec(offset, offset + 5);
            vec_start_ = vec_value_;
            break;
        case IVarType::MATRIX_6x6:
            for (unsigned int i = 0; i < 6; ++i) {
                for (unsigned int j = 0; j < 6; ++j) {
                    mat_value_(i, j) = statev(offset + i * 6 + j);
                    mat_start_(i, j) = statev(offset + i * 6 + j);
                }
            }
            break;
    }
}

} // namespace simcoon
