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
 * @brief Collection of internal variables for modular constitutive models.
 *
 * This class manages multiple InternalVariable instances, providing:
 * - Named registration and access
 * - Automatic offset computation for statev serialization
 * - Bulk operations (pack/unpack, rotate, state management)
 *
 * @version 1.0
 */

#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include <memory>
#include <armadillo>
#include <simcoon/Continuum_mechanics/Umat/Modular/internal_variable.hpp>

namespace simcoon {

/**
 * @brief Manages a collection of internal variables for constitutive models
 *
 * This class provides a container for InternalVariable instances with
 * named access, automatic offset computation for statev serialization,
 * and bulk operations for common tasks.
 */
class InternalVariableCollection {
private:
    std::vector<std::unique_ptr<InternalVariable>> variables_;
    std::unordered_map<std::string, size_t> name_to_index_;
    unsigned int total_size_;
    bool offsets_computed_;

public:
    /**
     * @brief Default constructor
     */
    InternalVariableCollection();

    // Disable copy (due to unique_ptr)
    InternalVariableCollection(const InternalVariableCollection&) = delete;
    InternalVariableCollection& operator=(const InternalVariableCollection&) = delete;

    // Allow move
    InternalVariableCollection(InternalVariableCollection&&) = default;
    InternalVariableCollection& operator=(InternalVariableCollection&&) = default;

    ~InternalVariableCollection() = default;

    // ========== Registration ==========

    /**
     * @brief Add a scalar internal variable
     * @param name Unique identifier for the variable
     * @param init Initial value (default: 0.0)
     * @param rotate Whether rotation should be applied (default: false)
     * @return Reference to the created variable
     * @throws std::runtime_error if name already exists
     */
    InternalVariable& add_scalar(const std::string& name, double init = 0.0, bool rotate = false);

    /**
     * @brief Add a 6-component vector internal variable
     * @param name Unique identifier for the variable
     * @param init Initial value (default: zeros)
     * @param rotate Whether rotation should be applied (default: true)
     * @return Reference to the created variable
     * @throws std::runtime_error if name already exists
     */
    InternalVariable& add_vec(const std::string& name, const arma::vec& init = arma::zeros(6), bool rotate = true);

    /**
     * @brief Add a 6x6 matrix internal variable
     * @param name Unique identifier for the variable
     * @param init Initial value (default: zeros)
     * @param rotate Whether rotation should be applied (default: true)
     * @return Reference to the created variable
     * @throws std::runtime_error if name already exists
     */
    InternalVariable& add_mat(const std::string& name, const arma::mat& init = arma::zeros(6, 6), bool rotate = true);

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
     * @brief Get variable by index
     * @param i Index in the collection
     * @return Reference to the variable
     */
    InternalVariable& operator[](size_t i);

    /**
     * @brief Get variable by index (const)
     * @param i Index in the collection
     * @return Const reference to the variable
     */
    const InternalVariable& operator[](size_t i) const;

    /**
     * @brief Get the number of variables
     * @return Number of registered variables
     */
    [[nodiscard]] size_t size() const noexcept { return variables_.size(); }

    /**
     * @brief Check if collection is empty
     * @return True if no variables are registered
     */
    [[nodiscard]] bool empty() const noexcept { return variables_.empty(); }

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

    /**
     * @brief Check if offsets have been computed
     * @return True if compute_offsets has been called
     */
    [[nodiscard]] bool offsets_computed() const noexcept { return offsets_computed_; }

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
     * @brief Apply rotation to all variables that require it
     * @param DR Rotation increment matrix (3x3)
     */
    void rotate_all(const arma::mat& DR);

    /**
     * @brief Copy current values to start values for all variables
     */
    void to_start_all();

    /**
     * @brief Copy start values to current values for all variables
     */
    void set_start_all();

    /**
     * @brief Reset all variables to zero
     */
    void reset_all();

    // ========== Iteration ==========

    /**
     * @brief Iterator type for range-based for loops
     */
    using iterator = std::vector<std::unique_ptr<InternalVariable>>::iterator;
    using const_iterator = std::vector<std::unique_ptr<InternalVariable>>::const_iterator;

    iterator begin() { return variables_.begin(); }
    iterator end() { return variables_.end(); }
    const_iterator begin() const { return variables_.begin(); }
    const_iterator end() const { return variables_.end(); }
    const_iterator cbegin() const { return variables_.cbegin(); }
    const_iterator cend() const { return variables_.cend(); }

    // ========== Utility ==========

    /**
     * @brief Get a list of all variable names
     * @return Vector of variable names in registration order
     */
    std::vector<std::string> names() const;

    /**
     * @brief Clear all variables
     */
    void clear();
};

} // namespace simcoon
