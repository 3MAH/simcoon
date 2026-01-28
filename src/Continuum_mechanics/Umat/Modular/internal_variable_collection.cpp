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
 * @file internal_variable_collection.cpp
 * @brief Implementation of InternalVariableCollection class
 */

#include <simcoon/Continuum_mechanics/Umat/Modular/internal_variable_collection.hpp>
#include <stdexcept>

namespace simcoon {

// ========== Constructor ==========

InternalVariableCollection::InternalVariableCollection()
    : variables_()
    , name_to_index_()
    , total_size_(0)
    , offsets_computed_(false)
{
}

// ========== Registration ==========

InternalVariable& InternalVariableCollection::add_scalar(const std::string& name, double init, bool rotate) {
    if (has(name)) {
        throw std::runtime_error("InternalVariableCollection: variable '" + name + "' already exists");
    }

    variables_.push_back(std::make_unique<InternalVariable>(name, init, rotate));
    name_to_index_[name] = variables_.size() - 1;
    total_size_ += 1;
    offsets_computed_ = false;

    return *variables_.back();
}

InternalVariable& InternalVariableCollection::add_vec(const std::string& name, const arma::vec& init, bool rotate) {
    if (has(name)) {
        throw std::runtime_error("InternalVariableCollection: variable '" + name + "' already exists");
    }

    variables_.push_back(std::make_unique<InternalVariable>(name, init, rotate));
    name_to_index_[name] = variables_.size() - 1;
    total_size_ += 6;
    offsets_computed_ = false;

    return *variables_.back();
}

InternalVariable& InternalVariableCollection::add_mat(const std::string& name, const arma::mat& init, bool rotate) {
    if (has(name)) {
        throw std::runtime_error("InternalVariableCollection: variable '" + name + "' already exists");
    }

    variables_.push_back(std::make_unique<InternalVariable>(name, init, rotate));
    name_to_index_[name] = variables_.size() - 1;
    total_size_ += 36;
    offsets_computed_ = false;

    return *variables_.back();
}

// ========== Access ==========

InternalVariable& InternalVariableCollection::get(const std::string& name) {
    auto it = name_to_index_.find(name);
    if (it == name_to_index_.end()) {
        throw std::runtime_error("InternalVariableCollection: variable '" + name + "' not found");
    }
    return *variables_[it->second];
}

const InternalVariable& InternalVariableCollection::get(const std::string& name) const {
    auto it = name_to_index_.find(name);
    if (it == name_to_index_.end()) {
        throw std::runtime_error("InternalVariableCollection: variable '" + name + "' not found");
    }
    return *variables_[it->second];
}

bool InternalVariableCollection::has(const std::string& name) const {
    return name_to_index_.find(name) != name_to_index_.end();
}

InternalVariable& InternalVariableCollection::operator[](size_t i) {
    if (i >= variables_.size()) {
        throw std::out_of_range("InternalVariableCollection: index " + std::to_string(i) +
                               " out of range (size: " + std::to_string(variables_.size()) + ")");
    }
    return *variables_[i];
}

const InternalVariable& InternalVariableCollection::operator[](size_t i) const {
    if (i >= variables_.size()) {
        throw std::out_of_range("InternalVariableCollection: index " + std::to_string(i) +
                               " out of range (size: " + std::to_string(variables_.size()) + ")");
    }
    return *variables_[i];
}

// ========== Offset Management ==========

void InternalVariableCollection::compute_offsets(unsigned int base_offset) {
    unsigned int current_offset = base_offset;

    for (auto& var : variables_) {
        var->set_offset(current_offset);
        current_offset += var->size();
    }

    offsets_computed_ = true;
}

// ========== Bulk Operations ==========

void InternalVariableCollection::pack_all(arma::vec& statev) const {
    if (!offsets_computed_) {
        throw std::runtime_error("InternalVariableCollection: compute_offsets must be called before pack_all");
    }

    for (const auto& var : variables_) {
        var->pack(statev);
    }
}

void InternalVariableCollection::unpack_all(const arma::vec& statev) {
    if (!offsets_computed_) {
        throw std::runtime_error("InternalVariableCollection: compute_offsets must be called before unpack_all");
    }

    for (auto& var : variables_) {
        var->unpack(statev);
    }
}

void InternalVariableCollection::rotate_all(const arma::mat& DR) {
    for (auto& var : variables_) {
        var->rotate(DR);
    }
}

void InternalVariableCollection::to_start_all() {
    for (auto& var : variables_) {
        var->to_start();
    }
}

void InternalVariableCollection::set_start_all() {
    for (auto& var : variables_) {
        var->set_start();
    }
}

void InternalVariableCollection::reset_all() {
    for (auto& var : variables_) {
        switch (var->type()) {
            case IVarType::SCALAR:
                var->scalar() = 0.0;
                var->scalar_start() = 0.0;
                break;
            case IVarType::VECTOR_6:
                var->vec().zeros();
                var->vec_start().zeros();
                break;
            case IVarType::MATRIX_6x6:
                var->mat().zeros();
                var->mat_start().zeros();
                break;
        }
    }
}

// ========== Utility ==========

std::vector<std::string> InternalVariableCollection::names() const {
    std::vector<std::string> result;
    result.reserve(variables_.size());

    for (const auto& var : variables_) {
        result.push_back(var->name());
    }

    return result;
}

void InternalVariableCollection::clear() {
    variables_.clear();
    name_to_index_.clear();
    total_size_ = 0;
    offsets_computed_ = false;
}

} // namespace simcoon
