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
 * @file hardening.hpp
 * @brief Hardening modules for modular UMAT.
 *
 * Provides isotropic and kinematic hardening models:
 * - Isotropic: Linear, Power-law, Voce, Combined Voce
 * - Kinematic: Prager (linear), Armstrong-Frederick, Chaboche (multi-AF)
 *
 * @see PlasticityMechanism, ModularUMAT
 * @version 1.0
 */

#pragma once

#include <memory>
#include <vector>
#include <string>
#include <armadillo>
#include <simcoon/Continuum_mechanics/Umat/Modular/internal_variable_collection.hpp>

namespace simcoon {

// ============================================================================
// ISOTROPIC HARDENING
// ============================================================================

/**
 * @brief Types of isotropic hardening
 */
enum class IsoHardType {
    NONE = 0,           ///< No isotropic hardening
    LINEAR = 1,         ///< Linear: R = H * p
    POWER_LAW = 2,      ///< Power law: R = k * p^m
    VOCE = 3,           ///< Voce saturation: R = Q * (1 - exp(-b*p))
    COMBINED_VOCE = 4   ///< Sum of N Voce terms
};

/**
 * @brief Base class for isotropic hardening
 */
class IsotropicHardening {
public:
    virtual ~IsotropicHardening() = default;

    /**
     * @brief Configure from props array
     * @param props Material properties vector
     * @param offset Current offset in props (will be updated)
     */
    virtual void configure(const arma::vec& props, int& offset) = 0;

    /**
     * @brief Compute hardening stress R(p)
     * @param p Accumulated plastic strain
     * @return Hardening stress
     */
    [[nodiscard]] virtual double R(double p) const = 0;

    /**
     * @brief Compute hardening modulus dR/dp
     * @param p Accumulated plastic strain
     * @return Derivative of hardening stress
     */
    [[nodiscard]] virtual double dR_dp(double p) const = 0;

    /**
     * @brief Get the hardening type
     * @return Hardening type enum
     */
    [[nodiscard]] virtual IsoHardType type() const = 0;

    /**
     * @brief Create an isotropic hardening instance
     * @param type Hardening type
     * @param N Number of terms (for COMBINED_VOCE)
     * @return Unique pointer to hardening instance
     */
    static std::unique_ptr<IsotropicHardening> create(IsoHardType type, int N = 1);

    /**
     * @brief Get number of props consumed by this hardening type
     * @param type The hardening type
     * @param N Number of terms (for COMBINED_VOCE)
     * @return Number of properties required
     */
    [[nodiscard]] static constexpr int props_count(IsoHardType type, int N = 1) {
        switch (type) {
            case IsoHardType::NONE:          return 0;
            case IsoHardType::LINEAR:        return 1;
            case IsoHardType::POWER_LAW:     return 2;
            case IsoHardType::VOCE:          return 2;
            case IsoHardType::COMBINED_VOCE: return 2 * N;
        }
        return 0;
    }
};

/**
 * @brief No isotropic hardening
 */
class NoIsotropicHardening final : public IsotropicHardening {
public:
    void configure(const arma::vec& props, int& offset) override {}
    double R(double p) const override { return 0.0; }
    double dR_dp(double p) const override { return 0.0; }
    IsoHardType type() const override { return IsoHardType::NONE; }
};

/**
 * @brief Linear isotropic hardening: R = H * p
 */
class LinearHardening final : public IsotropicHardening {
private:
    double H_;  ///< Hardening modulus
public:
    LinearHardening() : H_(0.0) {}
    void configure(const arma::vec& props, int& offset) override;
    double R(double p) const override { return H_ * p; }
    double dR_dp(double p) const override { return H_; }
    IsoHardType type() const override { return IsoHardType::LINEAR; }
};

/**
 * @brief Power-law isotropic hardening: R = k * p^m
 */
class PowerLawHardening final : public IsotropicHardening {
private:
    double k_;  ///< Hardening coefficient
    double m_;  ///< Hardening exponent
public:
    PowerLawHardening() : k_(0.0), m_(1.0) {}
    void configure(const arma::vec& props, int& offset) override;
    double R(double p) const override;
    double dR_dp(double p) const override;
    IsoHardType type() const override { return IsoHardType::POWER_LAW; }
};

/**
 * @brief Voce saturation hardening: R = Q * (1 - exp(-b*p))
 */
class VoceHardening final : public IsotropicHardening {
private:
    double Q_;  ///< Saturation stress
    double b_;  ///< Hardening rate
public:
    VoceHardening() : Q_(0.0), b_(0.0) {}
    void configure(const arma::vec& props, int& offset) override;
    double R(double p) const override;
    double dR_dp(double p) const override;
    IsoHardType type() const override { return IsoHardType::VOCE; }
};

/**
 * @brief Combined Voce hardening: R = sum_i Q_i * (1 - exp(-b_i*p))
 */
class CombinedVoceHardening final : public IsotropicHardening {
private:
    int N_;           ///< Number of Voce terms
    arma::vec Q_;     ///< Saturation stresses
    arma::vec b_;     ///< Hardening rates
public:
    explicit CombinedVoceHardening(int N = 1) : N_(N), Q_(arma::zeros(N)), b_(arma::zeros(N)) {}
    void configure(const arma::vec& props, int& offset) override;
    double R(double p) const override;
    double dR_dp(double p) const override;
    IsoHardType type() const override { return IsoHardType::COMBINED_VOCE; }
    int num_terms() const { return N_; }
};

// ============================================================================
// KINEMATIC HARDENING
// ============================================================================

/**
 * @brief Types of kinematic hardening
 */
enum class KinHardType {
    NONE = 0,                   ///< No kinematic hardening
    PRAGER = 1,                 ///< Linear Prager: X = (2/3)*C*alpha
    ARMSTRONG_FREDERICK = 2,    ///< Single AF: dX = (2/3)*C*dep - D*X*dp
    CHABOCHE = 3                ///< Multiple AF terms
};

/**
 * @brief Base class for kinematic hardening
 */
class KinematicHardening {
public:
    virtual ~KinematicHardening() = default;

    /**
     * @brief Configure from props array
     * @param props Material properties vector
     * @param offset Current offset in props (will be updated)
     */
    virtual void configure(const arma::vec& props, int& offset) = 0;

    /**
     * @brief Register internal variables for this hardening model
     * @param ivc Internal variable collection
     */
    virtual void register_variables(InternalVariableCollection& ivc) = 0;

    /**
     * @brief Get total backstress X = sum(X_i)
     * @param ivc Internal variable collection
     * @return Total backstress tensor (6 components)
     */
    virtual arma::vec total_backstress(const InternalVariableCollection& ivc) const = 0;

    /**
     * @brief Compute flow direction for kinematic variable alpha_i
     * @param i Index of the backstress (0 for single term)
     * @param n Plastic flow direction
     * @param ivc Internal variable collection
     * @return Flow direction for alpha_i
     */
    virtual arma::vec alpha_flow(int i, const arma::vec& n, const InternalVariableCollection& ivc) const = 0;

    /**
     * @brief Compute total kinematic hardening modulus contribution
     * @param n Plastic flow direction
     * @param ivc Internal variable collection
     * @return Hardening modulus contribution to K
     */
    virtual double hardening_modulus(const arma::vec& n, const InternalVariableCollection& ivc) const = 0;

    /**
     * @brief Update internal variables given plastic multiplier increment
     * @param dp Plastic multiplier increment
     * @param n Plastic flow direction
     * @param ivc Internal variable collection
     */
    virtual void update(double dp, const arma::vec& n, InternalVariableCollection& ivc) = 0;

    /**
     * @brief Get the hardening type
     * @return Hardening type enum
     */
    [[nodiscard]] virtual KinHardType type() const = 0;

    /**
     * @brief Get number of backstress terms
     * @return Number of backstress tensors
     */
    [[nodiscard]] virtual int num_backstresses() const = 0;

    /**
     * @brief Create a kinematic hardening instance
     * @param type Hardening type
     * @param N Number of backstress terms (for CHABOCHE)
     * @return Unique pointer to hardening instance
     */
    static std::unique_ptr<KinematicHardening> create(KinHardType type, int N = 1);

    /**
     * @brief Get number of props consumed by this hardening type
     * @param type The hardening type
     * @param N Number of backstress terms (for CHABOCHE)
     * @return Number of properties required
     */
    [[nodiscard]] static constexpr int props_count(KinHardType type, int N = 1) {
        switch (type) {
            case KinHardType::NONE:                 return 0;
            case KinHardType::PRAGER:               return 1;
            case KinHardType::ARMSTRONG_FREDERICK:   return 2;
            case KinHardType::CHABOCHE:             return 2 * N;
        }
        return 0;
    }
};

/**
 * @brief No kinematic hardening
 */
class NoKinematicHardening final : public KinematicHardening {
public:
    void configure(const arma::vec& props, int& offset) override {}
    void register_variables(InternalVariableCollection& ivc) override {}
    arma::vec total_backstress(const InternalVariableCollection& ivc) const override {
        return arma::zeros(6);
    }
    arma::vec alpha_flow(int i, const arma::vec& n, const InternalVariableCollection& ivc) const override {
        return arma::zeros(6);
    }
    double hardening_modulus(const arma::vec& n, const InternalVariableCollection& ivc) const override {
        return 0.0;
    }
    void update(double dp, const arma::vec& n, InternalVariableCollection& ivc) override {}
    KinHardType type() const override { return KinHardType::NONE; }
    int num_backstresses() const override { return 0; }
};

/**
 * @brief Linear Prager kinematic hardening: X = (2/3)*C*alpha
 */
class PragerHardening final : public KinematicHardening {
private:
    double C_;  ///< Kinematic hardening modulus
public:
    PragerHardening() : C_(0.0) {}
    void configure(const arma::vec& props, int& offset) override;
    void register_variables(InternalVariableCollection& ivc) override;
    arma::vec total_backstress(const InternalVariableCollection& ivc) const override;
    arma::vec alpha_flow(int i, const arma::vec& n, const InternalVariableCollection& ivc) const override;
    double hardening_modulus(const arma::vec& n, const InternalVariableCollection& ivc) const override;
    void update(double dp, const arma::vec& n, InternalVariableCollection& ivc) override;
    KinHardType type() const override { return KinHardType::PRAGER; }
    int num_backstresses() const override { return 1; }
};

/**
 * @brief Armstrong-Frederick kinematic hardening: dX = (2/3)*C*dep - D*X*dp
 */
class ArmstrongFrederickHardening final : public KinematicHardening {
private:
    double C_;  ///< Hardening parameter
    double D_;  ///< Dynamic recovery parameter
public:
    ArmstrongFrederickHardening() : C_(0.0), D_(0.0) {}
    void configure(const arma::vec& props, int& offset) override;
    void register_variables(InternalVariableCollection& ivc) override;
    arma::vec total_backstress(const InternalVariableCollection& ivc) const override;
    arma::vec alpha_flow(int i, const arma::vec& n, const InternalVariableCollection& ivc) const override;
    double hardening_modulus(const arma::vec& n, const InternalVariableCollection& ivc) const override;
    void update(double dp, const arma::vec& n, InternalVariableCollection& ivc) override;
    KinHardType type() const override { return KinHardType::ARMSTRONG_FREDERICK; }
    int num_backstresses() const override { return 1; }
};

/**
 * @brief Chaboche kinematic hardening: dX_i = (2/3)*C_i*dep - D_i*X_i*dp
 *
 * Multiple Armstrong-Frederick terms for capturing different hardening time scales.
 */
class ChabocheHardening final : public KinematicHardening {
private:
    int N_;           ///< Number of backstress terms
    arma::vec C_;     ///< Hardening parameters
    arma::vec D_;     ///< Dynamic recovery parameters
public:
    explicit ChabocheHardening(int N = 1) : N_(N), C_(arma::zeros(N)), D_(arma::zeros(N)) {}
    void configure(const arma::vec& props, int& offset) override;
    void register_variables(InternalVariableCollection& ivc) override;
    arma::vec total_backstress(const InternalVariableCollection& ivc) const override;
    arma::vec alpha_flow(int i, const arma::vec& n, const InternalVariableCollection& ivc) const override;
    double hardening_modulus(const arma::vec& n, const InternalVariableCollection& ivc) const override;
    void update(double dp, const arma::vec& n, InternalVariableCollection& ivc) override;
    KinHardType type() const override { return KinHardType::CHABOCHE; }
    int num_backstresses() const override { return N_; }
};

} // namespace simcoon
