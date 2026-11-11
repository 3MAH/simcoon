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
 * @brief Generalized Maxwell viscoelasticity, porting the Prony_Nfast kernel.
 *
 * Each Prony branch i is a Maxwell element with its own stiffness L_i and
 * viscosity tensor H_i (bulk/shear split): driving stress L_i (eps - EV_i)
 * produces a strain rate invH_i . L_i (eps - EV_i). EV_i is the branch
 * viscous strain tensor, v_i accumulates || Delta EV_i ||.
 *
 * The branches are driven by the total strain alone and are linear, so their
 * backward-Euler step is taken in closed form (maxwell_parallel_step,
 * linear_viscoelastic.hpp), exactly as the PRONK kernel: the mechanism keeps
 * its rows in the Fischer-Burmeister system, but they are always satisfied
 * (Phi = -Y_crit, Delta s = 0), and predict() sets EV_i to the closed-form
 * state before the elastic prediction, so the other mechanisms' return
 * mapping starts from it. The stress sees it through inelastic_strain(); the tangent through
 * total_strain_map(), I - sum_i M_0 L_i C_i L_i, applied by the orchestrator
 * after every other contribution (exact, whatever the other mechanisms).
 *
 * Props per branch (4 values): E_i, nu_i, etaB_i, etaS_i.
 * Statev per branch (7 values): v_i (scalar), EV_i (6 Voigt).
 *
 * @see StrainMechanism, ModularUMAT, Prony_Nfast.cpp (original kernel)
 * @version 2.0
 */

#pragma once

#include <vector>
#include <string>
#include <armadillo>
#include <simcoon/Continuum_mechanics/Umat/Modular/strain_mechanism.hpp>

namespace simcoon {

class ViscoelasticMechanism final : public StrainMechanism {
private:
    int N_prony_;
    std::vector<double> E_i_;       ///< Branch Young's moduli
    std::vector<double> nu_i_;      ///< Branch Poisson ratios
    std::vector<double> etaB_i_;    ///< Branch bulk viscosities
    std::vector<double> etaS_i_;    ///< Branch shear viscosities

    // Cached per-branch tensors (set in configure; constant over the step)
    std::vector<arma::mat> L_i_;        ///< Branch stiffnesses
    std::vector<arma::mat> H_i_;        ///< Branch viscosity tensors
    std::vector<arma::mat> M0_L_i_;     ///< Cached M_0 · L_i (for inelastic_strain)

    arma::mat L_0_;                     ///< Reference stiffness, set by orchestrator (zero until then)

    // IVC keys of register_variables; the hot loop reads the variables by index
    std::vector<std::string> ev_key_;   ///< "EV_<i>" per branch
    std::vector<std::string> v_key_;    ///< "v_<i>" per branch
    std::vector<size_t> ev_index_;      ///< ivc_ index of "EV_<i>", resolved by register_variables
    std::vector<size_t> v_index_;       ///< ivc_ index of "v_<i>", resolved by register_variables

    // Per-iteration caches (CCP: frozen flow direction from previous iteration)
    arma::mat kappa_t_;                           ///< Zero fluxes (6 x N): the rows carry no multiplier
    mutable arma::vec::fixed<6> EV_tilde_;        ///< inelastic_strain() buffer, valid for the branch strains of EV_at_
    mutable std::vector<double> EV_at_;           ///< the EV_i (6 each) EV_tilde_ was summed from; empty: not yet
    arma::mat dEVtilde_dE_;                       ///< d(sum_i M_0 L_i EV_i)/d(eps), set by predict()
    arma::mat map_;                               ///< I - dEVtilde_dE_, set by predict()
    arma::vec alpha_;                             ///< Thermal expansion of the block (zero until set_thermal_expansion)
    std::vector<arma::mat> dEVdE_i_;              ///< d(EV_i)/d(eps) per branch, set by predict()

public:
    explicit ViscoelasticMechanism(int N_prony);

    void configure(const arma::vec& props, int& offset) override;
    void register_variables() override;

    [[nodiscard]] int num_constraints() const override { return N_prony_; }
    [[nodiscard]] MechanismType type() const override { return MechanismType::VISCOELASTICITY; }

    /// Strain-form rate rows, stress-decoupled — excluded from the
    /// drift-guard ARMING count (see StrainMechanism::guarded_constraints);
    /// the commit checks themselves are no-ops here through the default
    /// multiplier_cap()/state_drift().
    [[nodiscard]] bool guarded_constraints() const override { return false; }
    [[nodiscard]] bool carries_multipliers() const override { return false; }   // rows integrated in predict()

    void compute_constraints(
        const arma::vec& sigma,
        const arma::vec& E_total,
        const arma::mat& L,
        double DTime,
        arma::vec& Phi,
        arma::vec& Y_crit
    ) const override;
    void compute_jacobian_contribution(
        const arma::vec& sigma,
        const arma::mat& L,
        arma::mat& B,
        int row_offset
    ) const override;

    // Viscoelastic Phi is written in strain form (Prony_Nfast) and does not
    // depend on sigma, so dPhi_dsigma returns empty → no cross-mechanism
    // coupling inbound from stress. The mechanism still contributes kappa to
    // outbound stress perturbations seen by other mechanisms.
    [[nodiscard]] const arma::mat& kappa(
        const arma::vec& sigma,
        double DT,
        const arma::mat& L_ref) const override;

    [[nodiscard]] const arma::vec& inelastic_strain() const override;

    void update(
        const arma::vec& ds,
        int offset
    ) override;

    [[nodiscard]] const arma::mat& total_strain_map() const override;
    void predict(const arma::vec& E_total_end, double DTime) override;
    /**
     * @brief Derivatives of the viscous dissipation of the increment, the trapezoid
     * \f$ \Delta\gamma = \sum_i \tfrac12 (\mathbf{A}_i^n + \mathbf{A}_i^{n+1}) : \Delta\mathbf{E}^v_i \f$ on the
     * branch forces \f$ \mathbf{A}_i = \mathbf{L}_i (\boldsymbol{\varepsilon} - \mathbf{E}^v_i) \f$ (compute_work),
     * exact for the backward-Euler branches (PRONK_T):
     * \f$ \partial\Delta\gamma/\partial\boldsymbol{\varepsilon} = \sum_i \tfrac12 [(\mathbf{L}_i(\mathbf{I} -
     * \partial\mathbf{E}^v_i/\partial\boldsymbol{\varepsilon}))^T \Delta\mathbf{E}^v_i + (\partial\mathbf{E}^v_i/\partial\boldsymbol{\varepsilon})^T
     * (\mathbf{A}_i^n + \mathbf{A}_i^{n+1})] \f$; every branch quantity depends on
     * \f$ \boldsymbol{\varepsilon} - \boldsymbol{\alpha}\,\Delta\theta \f$ alone, so the temperature one is
     * \f$ -\boldsymbol{\alpha} \cdot \partial\Delta\gamma/\partial\boldsymbol{\varepsilon} \f$.
     * The rows carry no multiplier: the multiplier sensitivities are not read.
     */
    void dissipation_sensitivity(const DissipationInputs& in, arma::vec& dDgamma_dE, double& dDgamma_dT) const override;

    void tangent_contribution(
        const arma::vec& sigma,
        const arma::mat& L,
        const arma::vec& Ds,
        int offset,
        arma::mat& Lt
    ) const override;

    void compute_work(
        const arma::vec& sigma_start,
        const arma::vec& sigma,
        const arma::vec& E_start,
        const arma::vec& E_end,
        double& Wm_r,
        double& Wm_ir,
        double& Wm_d
    ) const override;

    // Accessors
    [[nodiscard]] int num_prony_terms() const noexcept { return N_prony_; }
    [[nodiscard]] double E(int i)    const { return E_i_[i]; }
    [[nodiscard]] double nu(int i)   const { return nu_i_[i]; }
    [[nodiscard]] double etaB(int i) const { return etaB_i_[i]; }
    [[nodiscard]] double etaS(int i) const { return etaS_i_[i]; }
    [[nodiscard]] const arma::mat& L_branch(int i) const { return L_i_[i]; }

    /**
     * @brief Set reference compliance used to convert branch-stored viscous
     * strain into the mechanism's contribution to the total inelastic strain.
     * Called by ModularUMAT after the elasticity module is configured.
     */
    void set_reference_stiffness(const arma::mat& L_0);
    /**
     * @brief Thermal expansion of the block, for \f$ \partial\mathbf{E}^v_i/\partial\theta \f$ of the
     * thermomechanical face. Called by ModularUMAT after the elasticity module is configured.
     */
    void set_thermal_expansion(const arma::vec& alpha);
};

} // namespace simcoon
