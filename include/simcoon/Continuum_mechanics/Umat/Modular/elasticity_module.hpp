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
 * This module wraps existing elasticity functions (L_iso, L_cubic, L_ortho,
 * L_isotrans) to provide a configurable elasticity component.
 *
 * Elastic constants are passed as ordinal slots (C1, C2, ...) whose meaning
 * is selected by a per-symmetry convention enum (IsoConv, CubicConv, ...).
 * The convention codes are stable integers so they can travel in the flat
 * props stream: every elasticity block starts with one convention slot,
 * followed by the constants and the CTE values. A convention changes the
 * INTERPRETATION of the slots, never their count — the props layout is
 * invariant per ElasticityType. All parameterization conversions live in
 * the classical builders (constitutive.cpp); this module only maps the
 * enum to the builders' convention strings.
 *
 * @version 1.0
 */

#pragma once

#include <armadillo>
#include <stdexcept>
#include <simcoon/Continuum_mechanics/Functions/tensor.hpp>

namespace simcoon {

/**
 * @brief Types of linear elasticity
 */
enum class ElasticityType {
    ISOTROPIC = 0,              ///< Isotropic: conv, C1, C2, alpha
    CUBIC = 1,                  ///< Cubic: conv, C1, C2, C3, alpha (3 independent elastic constants)
    TRANSVERSE_ISOTROPIC = 2,   ///< Transverse isotropic: conv, EL, ET, nuTL, nuTT, GLT, alpha_L, alpha_T, axis
    ORTHOTROPIC = 3             ///< Orthotropic: conv, C1..C9, alpha1, alpha2, alpha3
};

/**
 * @brief Parameterization of the two isotropic elastic constants (C1, C2).
 *
 * Mirrors the convention strings of L_iso/M_iso:
 * | code | convention | C1      | C2      |
 * |------|------------|---------|---------|
 * | 0    | Enu        | E       | nu      |
 * | 1    | nuE        | nu      | E       |
 * | 2    | Kmu        | K       | mu (=G) |
 * | 3    | muK        | mu (=G) | K       |
 * | 4    | lambdamu   | lambda  | mu (=G) |
 * | 5    | mulambda   | mu (=G) | lambda  |
 */
enum class IsoConv : int {
    Enu = 0,
    nuE = 1,
    Kmu = 2,
    muK = 3,
    lambdamu = 4,
    mulambda = 5
};

/**
 * @brief Parameterization of the three cubic elastic constants (C1, C2, C3).
 *
 * Mirrors the convention strings of L_cubic/M_cubic:
 * | code | convention | C1  | C2  | C3  |
 * |------|------------|-----|-----|-----|
 * | 0    | EnuG       | E   | nu  | G   |
 * | 1    | Cii        | C11 | C12 | C44 |
 */
enum class CubicConv : int {
    EnuG = 0,
    Cii = 1
};

/**
 * @brief Parameterization of the transversely isotropic constants.
 *
 * L_isotrans/M_isotrans implement a single parameterization
 * (EL, ET, nuTL, nuTT, GLT); the enum exists so the props layout keeps a
 * convention slot uniformly across all elasticity types, and so further
 * parameterizations remain additive.
 */
enum class IsotransConv : int {
    EnuG = 0
};

/**
 * @brief Parameterization of the nine orthotropic elastic constants (C1..C9).
 *
 * Mirrors the convention strings of L_ortho/M_ortho:
 * | code | convention | C1  | C2  | C3  | C4   | C5   | C6   | C7  | C8  | C9  |
 * |------|------------|-----|-----|-----|------|------|------|-----|-----|-----|
 * | 0    | EnuG       | E1  | E2  | E3  | nu12 | nu13 | nu23 | G12 | G13 | G23 |
 * | 1    | Cii        | C11 | C12 | C13 | C22  | C23  | C33  | C44 | C55 | C66 |
 */
enum class OrthoConv : int {
    EnuG = 0,
    Cii = 1
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
    tensor4 L_t_;           ///< Typed stiffness, rebuilt by configure_* (eng→Mandel once)
    tensor4 M_t_;           ///< Typed compliance, rebuilt by configure_*
    bool configured_;

    /// Rebuild the typed L/M mirrors after L_/M_ change; every configure_*
    /// path must call this so L_tensor()/M_tensor() stay cheap const-refs.
    /// Also validates the material: L must be strictly positive definite —
    /// a wrong parameter order or inadmissible constants (e.g. Poisson
    /// combinations) otherwise reach the return mapping as an indefinite
    /// stiffness and fail far from the cause.
    void refresh_tensors() {
        // Cholesky is the cheap positive-definiteness test (succeeds iff SPD,
        // ~2-3x faster than a full eigendecomposition) — matters because
        // umat_modular reconstructs the module on every integration-point call.
        arma::mat R;
        if (!arma::chol(R, L_)) {
            throw std::runtime_error(
                "ElasticityModule: stiffness is not positive definite "
                "(check the elastic constants and their order)");
        }
        L_t_ = tensor4(L_, Tensor4Type::stiffness);
        M_t_ = tensor4(M_, Tensor4Type::compliance);
    }

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
     * @param C1 First elastic constant, interpreted per @p conv (see IsoConv)
     * @param C2 Second elastic constant, interpreted per @p conv
     * @param alpha Coefficient of thermal expansion (scalar, isotropic)
     * @param conv Parameterization of (C1, C2) — default Enu: C1 = E, C2 = nu
     */
    void configure_isotropic(double C1, double C2, double alpha,
                             IsoConv conv = IsoConv::Enu);

    /**
     * @brief Configure as cubic elasticity
     * @param C1 First elastic constant, interpreted per @p conv (see CubicConv)
     * @param C2 Second elastic constant, interpreted per @p conv
     * @param C3 Third elastic constant, interpreted per @p conv
     * @param alpha Coefficient of thermal expansion (scalar, isotropic CTE)
     * @param conv Parameterization of (C1, C2, C3) — default EnuG:
     *        C1 = E, C2 = nu, C3 = G; Cii: C1 = C11, C2 = C12, C3 = C44
     *
     * Cubic symmetry has 3 independent elastic constants. Unlike isotropic
     * materials where G = E/(2(1+nu)), G is an independent parameter here.
     * The Zener anisotropy ratio A = 2*G*(1+nu)/E quantifies the deviation
     * from isotropy (A=1 for isotropic).
     */
    void configure_cubic(double C1, double C2, double C3, double alpha,
                         CubicConv conv = CubicConv::EnuG);

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
     * @param conv Parameterization — a single one exists (see IsotransConv)
     */
    void configure_transverse_isotropic(double EL, double ET, double nuTL, double nuTT,
                                        double GLT, double alpha_L, double alpha_T, int axis = 3,
                                        IsotransConv conv = IsotransConv::EnuG);

    /**
     * @brief Configure as orthotropic elasticity
     * @param C1,C2,C3,C4,C5,C6,C7,C8,C9 The nine elastic constants,
     *        interpreted per @p conv (see OrthoConv). Default EnuG:
     *        (E1, E2, E3, nu12, nu13, nu23, G12, G13, G23);
     *        Cii: (C11, C12, C13, C22, C23, C33, C44, C55, C66).
     * @param alpha1 CTE in direction 1
     * @param alpha2 CTE in direction 2
     * @param alpha3 CTE in direction 3
     * @param conv Parameterization of the nine constants
     */
    void configure_orthotropic(double C1, double C2, double C3,
                               double C4, double C5, double C6,
                               double C7, double C8, double C9,
                               double alpha1, double alpha2, double alpha3,
                               OrthoConv conv = OrthoConv::EnuG);

    /**
     * @brief Configure from props array
     * @param type Elasticity type
     * @param props Material properties vector
     * @param offset Current offset in props (will be updated)
     *
     * Block layout (uniform across types): [conv, constants..., alphas...]
     * (+ axis for TRANSVERSE_ISOTROPIC). The conv slot is validated against
     * the type's convention enum and selects the interpretation of the
     * constant slots; see props_count() for per-type totals.
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

    // ========== Tensor-typed accessors (Tensor2/Tensor4 API) ==========

    /**
     * @brief Stiffness as a typed Tensor4 (cached — built once per configure).
     *
     * Use `.contract(strain_tensor)` to obtain the elastic stress tensor2 with
     * the correct VoigtType automatically inferred. Returned by const-ref so
     * hot loops can contract without paying the eng→Mandel congruence per call.
     */
    [[nodiscard]] const tensor4& L_tensor() const noexcept { return L_t_; }

    [[nodiscard]] const tensor4& M_tensor() const noexcept { return M_t_; }

    [[nodiscard]] tensor2 alpha_tensor() const {
        return strain(alpha_);
    }

    // ========== Derived Quantities ==========

    /**
     * @brief Get thermal strain
     * @param DT Temperature increment
     * @return Thermal strain vector (alpha * DT)
     */
    [[nodiscard]] arma::vec thermal_strain(double DT) const;

    [[nodiscard]] tensor2 thermal_strain_tensor(double DT) const {
        return strain(thermal_strain(DT));
    }

    /**
     * @brief Get number of props consumed by this elasticity type
     * @param type The elasticity type
     * @return Number of properties required (including the convention slot)
     */
    [[nodiscard]] static constexpr int props_count(ElasticityType type) {
        switch (type) {
            case ElasticityType::ISOTROPIC:            return 4;   // conv, C1, C2, alpha
            case ElasticityType::CUBIC:                return 5;   // conv, C1..C3, alpha
            case ElasticityType::TRANSVERSE_ISOTROPIC: return 9;   // conv, EL..GLT, alpha_L, alpha_T, axis
            case ElasticityType::ORTHOTROPIC:          return 13;  // conv, C1..C9, alpha1..3
        }
        return 0;
    }
};

} // namespace simcoon
