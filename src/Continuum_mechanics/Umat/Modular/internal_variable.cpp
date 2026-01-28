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

InternalVariable::InternalVariable(const std::string& name, double init, bool rotate)
    : name_(name)
    , type_(IVarType::SCALAR)
    , scalar_value_(init)
    , scalar_start_(init)
    , vec_value_()
    , vec_start_()
    , mat_value_()
    , mat_start_()
    , requires_rotation_(rotate)
    , statev_offset_(0)
{
}

InternalVariable::InternalVariable(const std::string& name, const arma::vec& init, bool rotate)
    : name_(name)
    , type_(IVarType::VECTOR_6)
    , scalar_value_(0.0)
    , scalar_start_(0.0)
    , vec_value_(init.n_elem == 6 ? init : arma::zeros(6))
    , vec_start_(init.n_elem == 6 ? init : arma::zeros(6))
    , mat_value_()
    , mat_start_()
    , requires_rotation_(rotate)
    , statev_offset_(0)
{
    if (init.n_elem != 6) {
        throw std::runtime_error("InternalVariable: vector must have 6 components, got " +
                                std::to_string(init.n_elem));
    }
}

InternalVariable::InternalVariable(const std::string& name, const arma::mat& init, bool rotate)
    : name_(name)
    , type_(IVarType::MATRIX_6x6)
    , scalar_value_(0.0)
    , scalar_start_(0.0)
    , vec_value_()
    , vec_start_()
    , mat_value_(init.n_rows == 6 && init.n_cols == 6 ? init : arma::zeros(6, 6))
    , mat_start_(init.n_rows == 6 && init.n_cols == 6 ? init : arma::zeros(6, 6))
    , requires_rotation_(rotate)
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

arma::vec& InternalVariable::vec() {
    if (type_ != IVarType::VECTOR_6) {
        throw std::runtime_error("InternalVariable '" + name_ + "': attempted to access vec on non-vector type");
    }
    return vec_value_;
}

const arma::vec& InternalVariable::vec() const {
    if (type_ != IVarType::VECTOR_6) {
        throw std::runtime_error("InternalVariable '" + name_ + "': attempted to access vec on non-vector type");
    }
    return vec_value_;
}

arma::mat& InternalVariable::mat() {
    if (type_ != IVarType::MATRIX_6x6) {
        throw std::runtime_error("InternalVariable '" + name_ + "': attempted to access mat on non-matrix type");
    }
    return mat_value_;
}

const arma::mat& InternalVariable::mat() const {
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

arma::vec& InternalVariable::vec_start() {
    if (type_ != IVarType::VECTOR_6) {
        throw std::runtime_error("InternalVariable '" + name_ + "': attempted to access vec_start on non-vector type");
    }
    return vec_start_;
}

const arma::vec& InternalVariable::vec_start() const {
    if (type_ != IVarType::VECTOR_6) {
        throw std::runtime_error("InternalVariable '" + name_ + "': attempted to access vec_start on non-vector type");
    }
    return vec_start_;
}

arma::mat& InternalVariable::mat_start() {
    if (type_ != IVarType::MATRIX_6x6) {
        throw std::runtime_error("InternalVariable '" + name_ + "': attempted to access mat_start on non-matrix type");
    }
    return mat_start_;
}

const arma::mat& InternalVariable::mat_start() const {
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

arma::mat InternalVariable::delta_mat() const {
    if (type_ != IVarType::MATRIX_6x6) {
        throw std::runtime_error("InternalVariable '" + name_ + "': attempted delta_mat on non-matrix type");
    }
    return mat_value_ - mat_start_;
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

void InternalVariable::set_start() {
    switch (type_) {
        case IVarType::SCALAR:
            scalar_value_ = scalar_start_;
            break;
        case IVarType::VECTOR_6:
            vec_value_ = vec_start_;
            break;
        case IVarType::MATRIX_6x6:
            mat_value_ = mat_start_;
            break;
    }
}

void InternalVariable::rotate(const arma::mat& DR) {
    if (!requires_rotation_) {
        return;
    }

    switch (type_) {
        case IVarType::SCALAR:
            // Scalars don't rotate
            break;
        case IVarType::VECTOR_6:
            // Use rotate_strain for strain-like tensors
            // This applies to plastic strain, backstress, etc.
            vec_value_ = rotate_strain(vec_value_, DR);
            break;
        case IVarType::MATRIX_6x6:
            // Use rotateL for 6x6 matrices (stiffness-like)
            mat_value_ = rotateL(mat_value_, DR);
            break;
    }
}

// ========== Serialization ==========

void InternalVariable::pack(arma::vec& statev) const {
    unsigned int offset = statev_offset_;

    switch (type_) {
        case IVarType::SCALAR:
            if (offset < statev.n_elem) {
                statev(offset) = scalar_value_;
            }
            break;
        case IVarType::VECTOR_6:
            for (unsigned int i = 0; i < 6 && (offset + i) < statev.n_elem; ++i) {
                statev(offset + i) = vec_value_(i);
            }
            break;
        case IVarType::MATRIX_6x6:
            for (unsigned int i = 0; i < 6; ++i) {
                for (unsigned int j = 0; j < 6; ++j) {
                    unsigned int idx = offset + i * 6 + j;
                    if (idx < statev.n_elem) {
                        statev(idx) = mat_value_(i, j);
                    }
                }
            }
            break;
    }
}

void InternalVariable::unpack(const arma::vec& statev) {
    unsigned int offset = statev_offset_;

    switch (type_) {
        case IVarType::SCALAR:
            if (offset < statev.n_elem) {
                scalar_value_ = statev(offset);
                scalar_start_ = statev(offset);
            }
            break;
        case IVarType::VECTOR_6:
            for (unsigned int i = 0; i < 6 && (offset + i) < statev.n_elem; ++i) {
                vec_value_(i) = statev(offset + i);
                vec_start_(i) = statev(offset + i);
            }
            break;
        case IVarType::MATRIX_6x6:
            for (unsigned int i = 0; i < 6; ++i) {
                for (unsigned int j = 0; j < 6; ++j) {
                    unsigned int idx = offset + i * 6 + j;
                    if (idx < statev.n_elem) {
                        mat_value_(i, j) = statev(idx);
                        mat_start_(i, j) = statev(idx);
                    }
                }
            }
            break;
    }
}

} // namespace simcoon
