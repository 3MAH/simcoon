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
 * @file internal_variable_collection.hpp
 * @brief Registry of internal variables for modular constitutive models.
 *
 * A name-indexed, registration-ordered registry of InternalVariable —
 * registration order IS the statev serialization layout, which is why a
 * plain map cannot replace it (hash/lexicographic iteration would make the
 * statev ABI nondeterministic). std::deque gives reference stability on
 * registration (add_* returns an InternalVariable& that must survive later
 * add_* calls) without per-element heap allocation.
 *
 * @version 1.0
 */

#pragma once

#include <string>
#include <deque>
#include <unordered_map>
#include <armadillo>
#include <simcoon/Continuum_mechanics/Umat/Modular/internal_variable.hpp>

namespace simcoon {

/**
 * @brief Registry of internal variables for constitutive models
 *
 * Provides named registration and O(1) lookup, automatic offset computation
 * for statev serialization, and the bulk operations the modular UMAT needs
 * per call (pack/unpack, co-rotation, state save).
 */
class InternalVariableCollection {
private:
    std::deque<InternalVariable> variables_;   ///< registration order = statev layout
    std::unordered_map<std::string, size_t> name_to_index_;
    unsigned int total_size_ = 0;
    bool offsets_computed_ = false;

public:
    // ========== Registration ==========

    /**
     * @brief Add a scalar internal variable
     * @param name Unique identifier for the variable
     * @param init Initial value (default: 0.0)
     * @return Reference to the created variable (scalars are always
     *         objective — no flag; see InternalVariable::is_objective())
     * @throws std::runtime_error if name already exists
     */
    InternalVariable& add_scalar(const std::string& name, double init = 0.0);

    /**
     * @brief Add a 6-component vector internal variable
     * @param name Unique identifier for the variable
     * @param init Initial value (default: zeros)
     * @param objective Whether the quantity is objective (default: true —
     *        pass false for velocity/rate-like variables, which rotate_all
     *        must not co-rotate)
     * @return Reference to the created variable
     * @throws std::runtime_error if name already exists
     */
    InternalVariable& add_vec(const std::string& name,
                               const arma::vec& init = arma::zeros(6),
                               bool objective = true,
                               Tensor2Type vtype = Tensor2Type::strain);

    /**
     * @brief Add a 6x6 matrix internal variable
     * @param name Unique identifier for the variable
     * @param init Initial value (default: zeros)
     * @param objective Whether the quantity is objective (default: true)
     * @return Reference to the created variable
     * @throws std::runtime_error if name already exists
     */
    InternalVariable& add_mat(const std::string& name, const arma::mat& init = arma::zeros(6, 6), bool objective = true,
                              Tensor4Type t4type = Tensor4Type::stiffness);

    // ========== Access ==========

    /**
     * @brief Get a variable by name
     * @param name The variable name
     * @return Reference to the variable
     * @throws std::runtime_error if name not found
     */
    InternalVariable& get(const std::string& name);

    /**
     * @brief Get a variable by name (const)
     * @param name The variable name
     * @return Const reference to the variable
     * @throws std::runtime_error if name not found
     */
    const InternalVariable& get(const std::string& name) const;

    /**
     * @brief Check if a variable exists
     * @param name The variable name
     * @return True if the variable exists
     */
    bool has(const std::string& name) const;

    /**
     * @brief Get the number of variables
     * @return Number of registered variables
     */
    [[nodiscard]] size_t size() const noexcept { return variables_.size(); }

    // ========== Offset Management ==========

    /**
     * @brief Get the total size needed in statev
     * @return Total number of scalar values across all variables
     */
    [[nodiscard]] unsigned int total_statev_size() const noexcept { return total_size_; }

    /**
     * @brief Compute offsets for all variables
     * @param base_offset Starting offset in statev (default: 0)
     *
     * This assigns sequential offsets to each variable based on their sizes.
     * Must be called before pack_all/unpack_all.
     */
    void compute_offsets(unsigned int base_offset = 0);

    // ========== Bulk Operations ==========

    /**
     * @brief Pack all variables into statev array
     * @param statev The state variable vector to pack into
     *
     * Requires compute_offsets() to have been called first.
     */
    void pack_all(arma::vec& statev) const;

    /**
     * @brief Unpack all variables from statev array
     * @param statev The state variable vector to unpack from
     *
     * Requires compute_offsets() to have been called first.
     * Also sets start values to current values.
     */
    void unpack_all(const arma::vec& statev);

    /**
     * @brief Co-rotate every objective variable with the rotation increment
     * @param DR Rotation increment matrix (3x3)
     *
     * Non-objective variables (is_objective() == false) are left untouched —
     * their transport is the owning mechanism's responsibility.
     */
    void rotate_all(const arma::mat& DR);

    /**
     * @brief Copy current values to start values for all variables
     */
    void to_start_all();
};

} // namespace simcoon
