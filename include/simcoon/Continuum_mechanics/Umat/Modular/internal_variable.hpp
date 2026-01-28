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
 * @file internal_variable.hpp
 * @brief Type-safe internal variable for constitutive models.
 *
 * This class provides storage for internal variables used in constitutive laws
 * (plasticity, viscoelasticity, damage, etc.) with support for:
 * - Scalar variables (accumulated plastic strain, damage, etc.)
 * - 6-component vectors in Voigt notation (plastic strain, backstress)
 * - 6x6 matrices (for complex tensorial internal variables)
 *
 * @version 1.0
 */

#pragma once

#include <string>
#include <stdexcept>
#include <armadillo>

namespace simcoon {

/**
 * @brief Type enumeration for internal variables
 */
enum class IVarType {
    SCALAR,      ///< Single scalar value (double)
    VECTOR_6,    ///< 6-component vector (Voigt notation for symmetric tensors)
    MATRIX_6x6   ///< 6x6 matrix (for fourth-order tensors in Voigt notation)
};

/**
 * @brief Type-safe internal variable with automatic serialization
 *
 * This class encapsulates internal variables used in constitutive models,
 * providing type safety, automatic pack/unpack to statev arrays, and
 * support for objectivity through rotation operations.
 */
class InternalVariable {
private:
    std::string name_;           ///< Identifier for debugging and access
    IVarType type_;              ///< Storage type

    // Storage (only one active based on type_)
    double scalar_value_;        ///< Current scalar value
    double scalar_start_;        ///< Scalar value at start of increment
    arma::vec vec_value_;        ///< Current vector value (6 components)
    arma::vec vec_start_;        ///< Vector value at start of increment
    arma::mat mat_value_;        ///< Current matrix value (6x6)
    arma::mat mat_start_;        ///< Matrix value at start of increment

    bool requires_rotation_;     ///< Whether to apply rotation for objectivity
    unsigned int statev_offset_; ///< Position in flat statev array

public:
    /**
     * @brief Construct a scalar internal variable
     * @param name Identifier for the variable
     * @param init Initial value (default: 0.0)
     * @param rotate Whether rotation should be applied (default: false for scalars)
     */
    InternalVariable(const std::string& name, double init = 0.0, bool rotate = false);

    /**
     * @brief Construct a vector internal variable (6 Voigt components)
     * @param name Identifier for the variable
     * @param init Initial value (default: zeros)
     * @param rotate Whether rotation should be applied (default: true for tensors)
     */
    InternalVariable(const std::string& name, const arma::vec& init, bool rotate = true);

    /**
     * @brief Construct a matrix internal variable (6x6)
     * @param name Identifier for the variable
     * @param init Initial value (default: zeros)
     * @param rotate Whether rotation should be applied (default: true for tensors)
     */
    InternalVariable(const std::string& name, const arma::mat& init, bool rotate = true);

    // Default copy/move/destructor
    InternalVariable(const InternalVariable&) = default;
    InternalVariable(InternalVariable&&) = default;
    InternalVariable& operator=(const InternalVariable&) = default;
    InternalVariable& operator=(InternalVariable&&) = default;
    ~InternalVariable() = default;

    // ========== Type Information ==========

    /**
     * @brief Get the variable type
     * @return The IVarType of this variable
     */
    [[nodiscard]] IVarType type() const noexcept { return type_; }

    /**
     * @brief Get the variable name
     * @return Reference to the name string
     */
    [[nodiscard]] const std::string& name() const noexcept { return name_; }

    /**
     * @brief Get the number of scalar values stored
     * @return 1 for scalar, 6 for vector, 36 for matrix
     */
    unsigned int size() const;

    /**
     * @brief Check if this variable requires rotation for objectivity
     * @return True if rotation should be applied
     */
    [[nodiscard]] bool requires_rotation() const noexcept { return requires_rotation_; }

    // ========== Value Accessors (throw if wrong type) ==========

    /**
     * @brief Access scalar value (mutable)
     * @return Reference to scalar value
     * @throws std::runtime_error if type is not SCALAR
     */
    double& scalar();

    /**
     * @brief Access scalar value (const)
     * @return Const reference to scalar value
     * @throws std::runtime_error if type is not SCALAR
     */
    const double& scalar() const;

    /**
     * @brief Access vector value (mutable)
     * @return Reference to vector value
     * @throws std::runtime_error if type is not VECTOR_6
     */
    arma::vec& vec();

    /**
     * @brief Access vector value (const)
     * @return Const reference to vector value
     * @throws std::runtime_error if type is not VECTOR_6
     */
    const arma::vec& vec() const;

    /**
     * @brief Access matrix value (mutable)
     * @return Reference to matrix value
     * @throws std::runtime_error if type is not MATRIX_6x6
     */
    arma::mat& mat();

    /**
     * @brief Access matrix value (const)
     * @return Const reference to matrix value
     * @throws std::runtime_error if type is not MATRIX_6x6
     */
    const arma::mat& mat() const;

    // ========== Start Value Accessors ==========

    /**
     * @brief Access scalar start value
     * @return Reference to scalar start value
     */
    double& scalar_start();
    const double& scalar_start() const;

    /**
     * @brief Access vector start value
     * @return Reference to vector start value
     */
    arma::vec& vec_start();
    const arma::vec& vec_start() const;

    /**
     * @brief Access matrix start value
     * @return Reference to matrix start value
     */
    arma::mat& mat_start();
    const arma::mat& mat_start() const;

    // ========== Increment Computation ==========

    /**
     * @brief Compute scalar increment (current - start)
     * @return The increment value
     */
    double delta_scalar() const;

    /**
     * @brief Compute vector increment (current - start)
     * @return The increment vector
     */
    arma::vec delta_vec() const;

    /**
     * @brief Compute matrix increment (current - start)
     * @return The increment matrix
     */
    arma::mat delta_mat() const;

    // ========== State Management ==========

    /**
     * @brief Copy current value to start value
     *
     * Call at the beginning of an increment to save the state.
     */
    void to_start();

    /**
     * @brief Copy start value to current value
     *
     * Call to restore state to beginning of increment.
     */
    void set_start();

    /**
     * @brief Apply rotation for objectivity
     * @param DR Rotation increment matrix (3x3)
     *
     * Applies rotation to maintain objectivity during large deformations.
     * - Scalar variables: not rotated (scalars are frame-invariant)
     * - Vector variables (6 Voigt): rotated via rotate_strain(value, DR)
     * - Matrix variables (6x6): not yet supported (placeholder)
     *
     * Only applied if requires_rotation_ is true. By default, scalars
     * have requires_rotation_ = false and tensorial quantities have
     * requires_rotation_ = true.
     *
     * @see rotate_strain() in contimech.hpp
     */
    void rotate(const arma::mat& DR);

    // ========== Serialization ==========

    /**
     * @brief Set the offset in the statev array
     * @param offset Starting index in statev
     */
    void set_offset(unsigned int offset) { statev_offset_ = offset; }

    /**
     * @brief Get the offset in the statev array
     * @return The offset value
     */
    [[nodiscard]] unsigned int offset() const noexcept { return statev_offset_; }

    /**
     * @brief Pack current value into statev array
     * @param statev The state variable vector to pack into
     *
     * Writes size() values starting at statev_offset_.
     */
    void pack(arma::vec& statev) const;

    /**
     * @brief Unpack value from statev array
     * @param statev The state variable vector to unpack from
     *
     * Reads size() values starting at statev_offset_.
     */
    void unpack(const arma::vec& statev);
};

} // namespace simcoon
