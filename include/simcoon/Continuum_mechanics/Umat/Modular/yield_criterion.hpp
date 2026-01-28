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
 * @file yield_criterion.hpp
 * @brief Yield criterion module for modular UMAT.
 *
 * This module wraps existing yield criterion functions from criteria.hpp
 * (Eq_stress, dEq_stress) to provide configurable yield surfaces.
 *
 * @version 1.0
 */

#pragma once

#include <string>
#include <armadillo>

namespace simcoon {

/**
 * @brief Types of yield criteria
 */
enum class YieldType {
    VON_MISES = 0,    ///< Isotropic von Mises (J2)
    TRESCA = 1,       ///< Tresca (max shear stress)
    DRUCKER = 2,      ///< Drucker (pressure-dependent)
    HILL = 3,         ///< Hill 1948 anisotropic
    DFA = 4,          ///< Deshpande-Fleck-Ashby
    ANISOTROPIC = 5   ///< Generic anisotropic (9 params)
};

/**
 * @brief Yield criterion module for modular UMAT
 *
 * This class encapsulates yield surface computations, wrapping the
 * existing Eq_stress and dEq_stress functions from criteria.hpp.
 */
class YieldCriterion {
private:
    YieldType type_;
    arma::vec params_;      ///< Criterion-specific parameters
    bool configured_;

    /**
     * @brief Convert YieldType to string for Eq_stress
     * @return String identifier for the criterion
     */
    std::string type_string() const;

public:
    /**
     * @brief Default constructor (von Mises)
     */
    YieldCriterion();

    // Default copy/move/destructor
    YieldCriterion(const YieldCriterion&) = default;
    YieldCriterion(YieldCriterion&&) = default;
    YieldCriterion& operator=(const YieldCriterion&) = default;
    YieldCriterion& operator=(YieldCriterion&&) = default;
    ~YieldCriterion() = default;

    // ========== Configuration ==========

    /**
     * @brief Configure as von Mises criterion
     */
    void configure_von_mises();

    /**
     * @brief Configure as Tresca criterion
     */
    void configure_tresca();

    /**
     * @brief Configure as Drucker criterion
     * @param b J3 influence parameter
     * @param n Exponent parameter
     */
    void configure_drucker(double b, double n);

    /**
     * @brief Configure as Hill 1948 criterion
     * @param F Hill parameter F
     * @param G Hill parameter G
     * @param H Hill parameter H
     * @param L Hill parameter L
     * @param M Hill parameter M
     * @param N Hill parameter N
     */
    void configure_hill(double F, double G, double H, double L, double M, double N);

    /**
     * @brief Configure as DFA criterion
     * @param F DFA parameter F
     * @param G DFA parameter G
     * @param H DFA parameter H
     * @param L DFA parameter L
     * @param M DFA parameter M
     * @param N DFA parameter N
     * @param K Hydrostatic sensitivity parameter
     */
    void configure_dfa(double F, double G, double H, double L, double M, double N, double K);

    /**
     * @brief Configure as generic anisotropic criterion
     * @param P11 Anisotropic tensor component (1,1)
     * @param P22 Anisotropic tensor component (2,2)
     * @param P33 Anisotropic tensor component (3,3)
     * @param P12 Anisotropic tensor component (1,2)
     * @param P13 Anisotropic tensor component (1,3)
     * @param P23 Anisotropic tensor component (2,3)
     * @param P44 Anisotropic tensor component (4,4)
     * @param P55 Anisotropic tensor component (5,5)
     * @param P66 Anisotropic tensor component (6,6)
     */
    void configure_anisotropic(double P11, double P22, double P33,
                               double P12, double P13, double P23,
                               double P44, double P55, double P66);

    /**
     * @brief Configure from props array
     * @param type Yield criterion type
     * @param props Material properties vector
     * @param offset Current offset in props (will be updated)
     */
    void configure(YieldType type, const arma::vec& props, int& offset);

    // ========== Accessors ==========

    /**
     * @brief Get the yield criterion type
     * @return The configured type
     */
    [[nodiscard]] YieldType type() const noexcept { return type_; }

    /**
     * @brief Check if the criterion is configured
     * @return True if configure has been called
     */
    [[nodiscard]] bool is_configured() const noexcept { return configured_; }

    /**
     * @brief Get the criterion parameters
     * @return Const reference to parameters
     */
    [[nodiscard]] const arma::vec& params() const noexcept { return params_; }

    // ========== Yield Function Computations ==========

    /**
     * @brief Compute equivalent stress
     * @param sigma Stress tensor (6 Voigt components)
     * @return Equivalent stress value
     */
    double equivalent_stress(const arma::vec& sigma) const;

    /**
     * @brief Compute equivalent stress with backstress shift
     * @param sigma Stress tensor (6 Voigt components)
     * @param X Backstress tensor (6 Voigt components)
     * @return Equivalent stress of (sigma - X)
     */
    double equivalent_stress(const arma::vec& sigma, const arma::vec& X) const;

    /**
     * @brief Compute flow direction (derivative of equivalent stress)
     * @param sigma Stress tensor (6 Voigt components)
     * @return Flow direction (6 Voigt components)
     */
    arma::vec flow_direction(const arma::vec& sigma) const;

    /**
     * @brief Compute flow direction with backstress shift
     * @param sigma Stress tensor (6 Voigt components)
     * @param X Backstress tensor (6 Voigt components)
     * @return Flow direction at (sigma - X)
     */
    arma::vec flow_direction(const arma::vec& sigma, const arma::vec& X) const;

    /**
     * @brief Compute plastic flow direction (for non-associated flow)
     * @param sigma Stress tensor (6 Voigt components)
     * @return Plastic flow direction
     *
     * For associated flow rules, this is the same as flow_direction.
     * Override in subclasses for non-associated flow.
     */
    arma::vec plastic_flow(const arma::vec& sigma) const;

    /**
     * @brief Compute plastic flow direction with backstress shift
     * @param sigma Stress tensor (6 Voigt components)
     * @param X Backstress tensor (6 Voigt components)
     * @return Plastic flow direction at (sigma - X)
     */
    arma::vec plastic_flow(const arma::vec& sigma, const arma::vec& X) const;

    // ========== Utility ==========

    /**
     * @brief Get number of props consumed by this yield criterion type
     * @param type The yield criterion type
     * @return Number of properties required
     */
    [[nodiscard]] static constexpr int props_count(YieldType type) {
        switch (type) {
            case YieldType::VON_MISES:   return 0;
            case YieldType::TRESCA:      return 0;
            case YieldType::DRUCKER:     return 2;
            case YieldType::HILL:        return 6;
            case YieldType::DFA:         return 7;
            case YieldType::ANISOTROPIC: return 9;
        }
        return 0;
    }
};

} // namespace simcoon
