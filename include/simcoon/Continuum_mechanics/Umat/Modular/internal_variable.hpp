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
#include <simcoon/Continuum_mechanics/Functions/tensor.hpp>

namespace simcoon {

class Rotation;

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

    bool is_objective_;          ///< Whether the quantity is objective (frame-indifferent)
    VoigtType vtype_;            ///< Authoritative Voigt convention for this variable
                                 ///< (VECTOR_6): strain, stress, or generic. Drives
                                 ///< rotation dispatch (tensor2.rotate) and default
                                 ///< VoigtType for as_tensor2() views. Ignored for
                                 ///< scalars and 6x6 matrices.
    Tensor4Type t4type_;         ///< Authoritative convention for MATRIX_6x6 variables
                                 ///< (stiffness, compliance, ...). Drives the rotation
                                 ///< congruence (tensor4.rotate) and the default type of
                                 ///< as_tensor4() views. Ignored for scalars and vectors.
    unsigned int statev_offset_; ///< Position in flat statev array

public:
    /**
     * @brief Construct a scalar internal variable
     * @param name Identifier for the variable
     * @param init Initial value (default: 0.0)
     *
     * Scalars are ALWAYS objective (frame-invariant) — there is no flag.
     */
    InternalVariable(const std::string& name, double init = 0.0);

    /**
     * @brief Construct a vector internal variable (6 Voigt components)
     * @param name  Identifier for the variable
     * @param init  Initial value (default: zeros)
     * @param objective  Whether the quantity is objective (default: true; see
     *        is_objective() — pass false for velocity/rate-like variables)
     * @param vtype  Voigt convention for this variable:
     *        VoigtType::strain (default) for strain-like variables (EP, EV, α, ...),
     *        VoigtType::stress for stress-like variables (stored generalised forces),
     *        VoigtType::generic for plain 6-component vectors without Voigt semantics.
     *        Controls the rotation kernel (strain-rotation has factor-2 on shear)
     *        and is returned by default from as_tensor2().
     */
    InternalVariable(const std::string& name, const arma::vec& init,
                      bool objective = true,
                      VoigtType vtype = VoigtType::strain);

    /**
     * @brief Construct a matrix internal variable (6x6)
     * @param name Identifier for the variable
     * @param init Initial value (default: zeros)
     * @param objective Whether the quantity is objective (default: true; see
     *        is_objective())
     * @param t4type Tensor4 convention of this variable (default: stiffness).
     *        Controls the rotation congruence and the default type returned by
     *        as_tensor4().
     */
    InternalVariable(const std::string& name, const arma::mat& init, bool objective = true,
                      Tensor4Type t4type = Tensor4Type::stiffness);

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
     * @brief Is the stored quantity objective (frame-indifferent)?
     *
     * If TRUE, a rigid rotation of the components (co-rotation with the
     * material frame, rotate()) is the correct and sufficient objectivity
     * treatment. Scalars are ALWAYS objective (trivially invariant; rotate()
     * is a no-op for them). If FALSE (e.g. a velocity / rate-like variable),
     * a rigid rotation is NOT sufficient — rotate() leaves the variable
     * untouched and the owning mechanism must provide the appropriate
     * transport itself.
     */
    [[nodiscard]] bool is_objective() const noexcept { return is_objective_; }

    /**
     * @brief The Voigt convention of this variable.
     *
     * For VECTOR_6 variables, this tag drives the objectivity rotation kernel
     * and the default VoigtType returned by as_tensor2().
     */
    [[nodiscard]] VoigtType vtype() const noexcept { return vtype_; }

    /**
     * @brief The Tensor4 convention of this variable (MATRIX_6x6 only).
     */
    [[nodiscard]] Tensor4Type t4type() const noexcept { return t4type_; }

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
     * @brief Access the raw Voigt components (mutable).
     *
     * "raw" is deliberate: this is the escape hatch past set_tensor2()'s
     * convention re-expression — the components are engineering Voigt in this
     * variable's own VoigtType, and the CALLER owns keeping them consistent
     * (safe for same-convention linear combinations; use set_tensor2()
     * whenever the source is a typed tensor). To reconstruct the tensorial
     * quantity, do NOT re-tag raw copies — use as_tensor2()/as_tensor2_start(),
     * which pair the components with the variable's own vtype().
     * @throws std::runtime_error if type is not VECTOR_6
     */
    arma::vec& raw_voigt();

    /**
     * @brief Access the raw Voigt components (const)
     * @return Const reference to the 6-component vector (this variable's
     *         VoigtType convention)
     * @throws std::runtime_error if type is not VECTOR_6
     */
    const arma::vec& raw_voigt() const;

    /**
     * @brief Access the raw 6x6 components (mutable). Same escape-hatch
     * semantics as raw_voigt(): components are in this variable's Tensor4Type
     * convention; prefer set_tensor4() for typed sources.
     * @throws std::runtime_error if type is not MATRIX_6x6
     */
    arma::mat& raw_mat();

    /**
     * @brief Access the raw 6x6 components (const)
     * @throws std::runtime_error if type is not MATRIX_6x6
     */
    const arma::mat& raw_mat() const;

    // ========== Start Value Accessors ==========

    /**
     * @brief Access scalar start value
     * @return Reference to scalar start value
     */
    double& scalar_start();
    const double& scalar_start() const;

    /**
     * @brief Access the raw Voigt start value (same convention/escape-hatch
     * semantics as raw_voigt())
     */
    arma::vec& raw_voigt_start();
    const arma::vec& raw_voigt_start() const;

    /**
     * @brief Access the raw 6x6 start value (see raw_mat())
     */
    arma::mat& raw_mat_start();
    const arma::mat& raw_mat_start() const;

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

    // ========== State Management ==========

    /**
     * @brief Copy current value to start value
     *
     * Call at the beginning of an increment to save the state.
     */
    void to_start();

    /**
     * @brief Apply rotation for objectivity
     * @param DR Rotation increment matrix (3x3)
     *
     * Scalars are left unchanged. VECTOR_6 internals rotate through tensor2
     * with the variable's VoigtType (strain-rotation carries the factor-2 on
     * shear); MATRIX_6x6 internals rotate through tensor4 with the variable's
     * Tensor4Type, so the correct congruence is applied for stiffness- AND
     * compliance-like storage. No-op when is_objective() is false: a
     * non-objective quantity must not be blindly co-rotated — its transport
     * is mechanism-specific.
     */
    void rotate(const arma::mat& DR);

    /**
     * @brief Rotation-object overload of rotate(); prefer this when the caller
     * already has a Rotation to avoid rebuilding it from a raw matrix.
     */
    void rotate(const Rotation& R);

    // ========== Tensor-typed views (Tensor2/Tensor4 API) ==========

    /**
     * @brief View the stored vector as a typed tensor2.
     *
     * The caller selects the VoigtType matching the physical quantity:
     * - VoigtType::strain  for plastic strain, viscous strain, strain-form
     *   backstress — assumes the factor-2 shear convention used by rotate()
     *   and pack().
     * - VoigtType::stress  for stress-conjugate variables stored without the
     *   factor-2.
     * - VoigtType::generic when no specific physical convention applies.
     *
     * @throws std::runtime_error if type() != VECTOR_6
     */
    /// When @p vtype_override is left default (none), the variable's own
    /// stored VoigtType is used — callers never need to pass it explicitly
    /// unless they want to reinterpret the raw Voigt components under a
    /// different convention (rare).
    [[nodiscard]] tensor2 as_tensor2() const { return as_tensor2(vtype_); }
    [[nodiscard]] tensor2 as_tensor2(VoigtType vtype) const;

    /// Typed view of the START value (state at the beginning of the
    /// increment) — the tensorial counterpart of raw_voigt_start(), so
    /// backward-Euler closed forms can be written fully typed.
    /// @throws std::runtime_error if type() != VECTOR_6
    [[nodiscard]] tensor2 as_tensor2_start() const { return as_tensor2_start(vtype_); }
    [[nodiscard]] tensor2 as_tensor2_start(VoigtType vtype) const;

    [[nodiscard]] tensor2 as_strain() const { return as_tensor2(VoigtType::strain); }
    [[nodiscard]] tensor2 as_stress() const { return as_tensor2(VoigtType::stress); }

    /// Default view: the variable's own stored Tensor4Type is used (see
    /// as_tensor2() for the same pattern on vectors).
    /// @throws std::runtime_error if type() != MATRIX_6x6
    [[nodiscard]] tensor4 as_tensor4() const { return as_tensor4(t4type_); }
    [[nodiscard]] tensor4 as_tensor4(Tensor4Type t4type) const;

    /// Typed view of the START value (see as_tensor2_start()).
    /// @throws std::runtime_error if type() != MATRIX_6x6
    [[nodiscard]] tensor4 as_tensor4_start() const { return as_tensor4_start(t4type_); }
    [[nodiscard]] tensor4 as_tensor4_start(Tensor4Type t4type) const;

    [[nodiscard]] tensor4 as_stiffness() const { return as_tensor4(Tensor4Type::stiffness); }
    [[nodiscard]] tensor4 as_compliance() const { return as_tensor4(Tensor4Type::compliance); }

    /// Stores the TENSOR, not the raw components: the value is re-expressed
    /// in this variable's own Voigt convention (via the 3x3), so assigning a
    /// stress-typed tensor2 to a strain-typed variable cannot silently drop
    /// the factor-2 shear convention. No-op re-expression when conventions
    /// already match.
    /// @throws std::runtime_error if type() != VECTOR_6
    void set_tensor2(const tensor2& t);

    /// Same contract as set_tensor2: re-expressed in the variable's own
    /// Tensor4Type (via the convention-free Mandel form).
    /// @throws std::runtime_error if type() != MATRIX_6x6
    void set_tensor4(const tensor4& t);

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
     * @throws std::runtime_error if statev_offset_ + size() exceeds statev —
     * a partial write would silently corrupt the state.
     */
    void pack(arma::vec& statev) const;

    /**
     * @brief Unpack value from statev array
     * @param statev The state variable vector to unpack from
     *
     * Reads size() values starting at statev_offset_ into BOTH the current
     * and the start value (state restoration).
     * @throws std::runtime_error if statev_offset_ + size() exceeds statev.
     */
    void unpack(const arma::vec& statev);
};

} // namespace simcoon
