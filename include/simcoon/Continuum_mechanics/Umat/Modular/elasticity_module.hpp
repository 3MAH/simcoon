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
 * @file elasticity_module.hpp
 * @brief Elasticity module for modular UMAT.
 *
 * This module wraps existing elasticity functions (L_iso, L_ortho, L_isotrans)
 * to provide a configurable elasticity component.
 *
 * @version 1.0
 */

#pragma once

#include <armadillo>

namespace simcoon {

/**
 * @brief Types of linear elasticity
 */
enum class ElasticityType {
    ISOTROPIC = 0,              ///< Isotropic: E, nu, alpha
    CUBIC = 1,                  ///< Cubic: E, nu, G, alpha (3 independent elastic constants)
    TRANSVERSE_ISOTROPIC = 2,   ///< Transverse isotropic: EL, ET, nuTL, nuTT, GLT, alpha_L, alpha_T
    ORTHOTROPIC = 3             ///< Orthotropic: E1, E2, E3, nu12, nu13, nu23, G12, G13, G23, alpha1, alpha2, alpha3
};

/**
 * @brief Elasticity module for modular UMAT
 *
 * This class encapsulates the computation of the elastic stiffness tensor
 * and coefficient of thermal expansion for different types of linear elasticity.
 */
class ElasticityModule {
private:
    ElasticityType type_;
    arma::mat L_;           ///< 6x6 stiffness tensor
    arma::mat M_;           ///< 6x6 compliance tensor
    arma::vec alpha_;       ///< 6-component CTE (Voigt notation)
    bool configured_;

public:
    /**
     * @brief Default constructor
     */
    ElasticityModule();

    // Default copy/move/destructor
    ElasticityModule(const ElasticityModule&) = default;
    ElasticityModule(ElasticityModule&&) = default;
    ElasticityModule& operator=(const ElasticityModule&) = default;
    ElasticityModule& operator=(ElasticityModule&&) = default;
    ~ElasticityModule() = default;

    // ========== Configuration ==========

    /**
     * @brief Configure as isotropic elasticity
     * @param E Young's modulus
     * @param nu Poisson's ratio
     * @param alpha Coefficient of thermal expansion (scalar, isotropic)
     */
    void configure_isotropic(double E, double nu, double alpha);

    /**
     * @brief Configure as cubic elasticity
     * @param E Young's modulus
     * @param nu Poisson's ratio
     * @param G Shear modulus (independent from E and nu for cubic symmetry)
     * @param alpha Coefficient of thermal expansion (scalar, isotropic CTE)
     *
     * Cubic symmetry has 3 independent elastic constants. Unlike isotropic
     * materials where G = E/(2(1+nu)), G is an independent parameter here.
     * The Zener anisotropy ratio A = 2*G*(1+nu)/E quantifies the deviation
     * from isotropy (A=1 for isotropic).
     *
     * Alternative: use configure_cubic_Cii() with C11, C12, C44 directly.
     */
    void configure_cubic(double E, double nu, double G, double alpha);

    /**
     * @brief Configure as cubic elasticity using Cii convention
     * @param C11 Stiffness component C11
     * @param C12 Stiffness component C12
     * @param C44 Stiffness component C44
     * @param alpha Coefficient of thermal expansion (scalar, isotropic CTE)
     */
    void configure_cubic_Cii(double C11, double C12, double C44, double alpha);

    /**
     * @brief Configure as transversely isotropic elasticity
     * @param EL Longitudinal Young's modulus
     * @param ET Transverse Young's modulus
     * @param nuTL Poisson's ratio (transverse-longitudinal)
     * @param nuTT Poisson's ratio (transverse-transverse)
     * @param GLT Shear modulus
     * @param alpha_L Longitudinal CTE
     * @param alpha_T Transverse CTE
     * @param axis Axis of symmetry (1=x, 2=y, 3=z)
     */
    void configure_transverse_isotropic(double EL, double ET, double nuTL, double nuTT,
                                        double GLT, double alpha_L, double alpha_T, int axis = 3);

    /**
     * @brief Configure as orthotropic elasticity
     * @param E1 Young's modulus in direction 1
     * @param E2 Young's modulus in direction 2
     * @param E3 Young's modulus in direction 3
     * @param nu12 Poisson's ratio 1-2
     * @param nu13 Poisson's ratio 1-3
     * @param nu23 Poisson's ratio 2-3
     * @param G12 Shear modulus 1-2
     * @param G13 Shear modulus 1-3
     * @param G23 Shear modulus 2-3
     * @param alpha1 CTE in direction 1
     * @param alpha2 CTE in direction 2
     * @param alpha3 CTE in direction 3
     */
    void configure_orthotropic(double E1, double E2, double E3,
                               double nu12, double nu13, double nu23,
                               double G12, double G13, double G23,
                               double alpha1, double alpha2, double alpha3);

    /**
     * @brief Configure from props array
     * @param type Elasticity type
     * @param props Material properties vector
     * @param offset Current offset in props (will be updated)
     */
    void configure(ElasticityType type, const arma::vec& props, int& offset);

    // ========== Accessors ==========

    /**
     * @brief Get the elasticity type
     * @return The configured type
     */
    [[nodiscard]] ElasticityType type() const noexcept { return type_; }

    /**
     * @brief Check if the module is configured
     * @return True if configure has been called
     */
    [[nodiscard]] bool is_configured() const noexcept { return configured_; }

    /**
     * @brief Get the 6x6 stiffness tensor
     * @return Const reference to L
     */
    [[nodiscard]] const arma::mat& L() const noexcept { return L_; }

    /**
     * @brief Get the 6x6 compliance tensor
     * @return Const reference to M
     */
    [[nodiscard]] const arma::mat& M() const noexcept { return M_; }

    /**
     * @brief Get the CTE tensor
     * @return Const reference to alpha (6 components)
     */
    [[nodiscard]] const arma::vec& alpha() const noexcept { return alpha_; }

    // ========== Derived Quantities ==========

    /**
     * @brief Get damaged stiffness tensor
     * @param d Damage variable (0 = undamaged, 1 = fully damaged)
     * @return Degraded stiffness (1-d)*L
     */
    arma::mat damaged_L(double d) const;

    /**
     * @brief Get thermal strain
     * @param DT Temperature increment
     * @return Thermal strain vector (alpha * DT)
     */
    arma::vec thermal_strain(double DT) const;

    /**
     * @brief Get number of props consumed by this elasticity type
     * @param type The elasticity type
     * @return Number of properties required
     */
    [[nodiscard]] static constexpr int props_count(ElasticityType type) {
        switch (type) {
            case ElasticityType::ISOTROPIC:            return 3;
            case ElasticityType::CUBIC:                return 4;
            case ElasticityType::TRANSVERSE_ISOTROPIC: return 8;
            case ElasticityType::ORTHOTROPIC:          return 12;
        }
        return 0;
    }
};

} // namespace simcoon
