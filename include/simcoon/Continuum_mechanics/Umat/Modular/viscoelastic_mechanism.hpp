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
 * @file viscoelastic_mechanism.hpp
 * @brief Viscoelastic strain mechanism using Prony series (generalized Maxwell model)
 *
 * This mechanism implements a generalized Maxwell model with N Prony terms.
 * Each branch i has:
 * - E_i: Branch modulus (or g_i relative to E_0)
 * - tau_i: Relaxation time
 * - EV_i: Viscous strain tensor (internal variable)
 *
 * The viscoelastic strain evolves according to:
 *   dEV_i/dt = (1/tau_i) * (C_i^{-1} : sigma - EV_i)
 *
 * Or equivalently using the hereditary integral formulation.
 *
 * @see StrainMechanism, ModularUMAT
 * @version 1.0
 */

#pragma once

#include <vector>
#include <string>
#include <armadillo>
#include <simcoon/Continuum_mechanics/Umat/Modular/strain_mechanism.hpp>

namespace simcoon {

/**
 * @brief Viscoelastic mechanism using Prony series
 *
 * Implements a generalized Maxwell model with N branches.
 * Each branch contributes viscous strain that relaxes over time.
 */
class ViscoelasticMechanism final : public StrainMechanism {
private:
    int N_prony_;                       ///< Number of Prony terms
    std::vector<double> g_i_;           ///< Relative moduli (g_i = E_i/E_0)
    std::vector<double> tau_i_;         ///< Relaxation times
    arma::mat L_0_;                     ///< Reference stiffness (for computing branch stiffnesses)
    arma::mat M_0_;                     ///< Reference compliance (inverse of L_0_, stored to avoid recomputation)
    bool use_deviatoric_only_;          ///< If true, only apply viscoelasticity to deviatoric part

    // Cached values for tangent computation
    mutable std::vector<arma::vec> EV_n_;       ///< Viscous strains at start of increment
    mutable std::vector<arma::vec> flow_i_;     ///< Flow directions for each branch
    mutable std::vector<double> factor_i_;      ///< Update factors exp(-DTime/tau_i)

public:
    /**
     * @brief Constructor
     * @param N_prony Number of Prony terms
     */
    explicit ViscoelasticMechanism(int N_prony);

    // StrainMechanism interface
    void configure(const arma::vec& props, int& offset) override;
    void register_variables(InternalVariableCollection& ivc) override;

    [[nodiscard]] int num_constraints() const override { return N_prony_; }
    [[nodiscard]] MechanismType type() const override { return MechanismType::VISCOELASTICITY; }

    void compute_constraints(
        const arma::vec& sigma,
        const arma::mat& L,
        double DTime,
        const InternalVariableCollection& ivc,
        arma::vec& Phi,
        arma::vec& Y_crit
    ) const override;

    void compute_flow_directions(
        const arma::vec& sigma,
        const InternalVariableCollection& ivc,
        std::map<std::string, arma::vec>& Lambda_map
    ) const override;

    void compute_jacobian_contribution(
        const arma::vec& sigma,
        const arma::mat& L,
        const InternalVariableCollection& ivc,
        arma::mat& B,
        int row_offset
    ) const override;

    arma::vec inelastic_strain(const InternalVariableCollection& ivc) const override;

    void update(
        const arma::vec& ds,
        int offset,
        InternalVariableCollection& ivc
    ) override;

    void tangent_contribution(
        const arma::vec& sigma,
        const arma::mat& L,
        const arma::vec& Ds,
        int offset,
        const InternalVariableCollection& ivc,
        arma::mat& Lt
    ) const override;

    void compute_work(
        const arma::vec& sigma_start,
        const arma::vec& sigma,
        const InternalVariableCollection& ivc,
        double& Wm_r,
        double& Wm_ir,
        double& Wm_d
    ) const override;

    // Accessors
    [[nodiscard]] int num_prony_terms() const noexcept { return N_prony_; }
    [[nodiscard]] double g(int i) const { return g_i_[i]; }
    [[nodiscard]] double tau(int i) const { return tau_i_[i]; }

    /**
     * @brief Set reference stiffness for computing branch contributions
     * @param L_0 Reference elastic stiffness
     */
    void set_reference_stiffness(const arma::mat& L_0) { L_0_ = L_0; M_0_ = arma::inv(L_0); }

    /**
     * @brief Get the effective instantaneous stiffness reduction
     *
     * For a generalized Maxwell model, the instantaneous modulus is E_0,
     * and the long-term modulus is E_inf = E_0 * (1 - sum(g_i)).
     *
     * @return Long-term stiffness factor (1 - sum(g_i))
     */
    double long_term_factor() const;
};

} // namespace simcoon
