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
 * @file viscoelastic_mechanism.cpp
 * @brief Implementation of ViscoelasticMechanism class
 */

#include <simcoon/Continuum_mechanics/Umat/Modular/viscoelastic_mechanism.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/parameter.hpp>
#include <stdexcept>
#include <cmath>

namespace simcoon {

// ========== Constructor ==========

ViscoelasticMechanism::ViscoelasticMechanism(int N_prony)
    : StrainMechanism("viscoelastic")
    , N_prony_(N_prony)
    , g_i_(N_prony, 0.0)
    , tau_i_(N_prony, 1.0)
    , L_0_(arma::eye(6, 6))
    , M_0_(arma::eye(6, 6))
    , use_deviatoric_only_(false)
    , EV_n_(N_prony)
    , flow_i_(N_prony)
    , factor_i_(N_prony, 0.0)
{
}

// ========== Configuration ==========

void ViscoelasticMechanism::configure(const arma::vec& props, int& offset) {
    // Props layout for viscoelasticity:
    // For each Prony term i:
    //   props[offset + 2*i]: g_i (relative modulus)
    //   props[offset + 2*i + 1]: tau_i (relaxation time)

    for (int i = 0; i < N_prony_; ++i) {
        g_i_[i] = props(offset + 2 * i);
        tau_i_[i] = props(offset + 2 * i + 1);

        if (tau_i_[i] <= 0.0) {
            throw std::runtime_error("ViscoelasticMechanism: relaxation time must be positive");
        }
    }
    offset += 2 * N_prony_;

    // Initialize cached vectors
    for (int i = 0; i < N_prony_; ++i) {
        EV_n_[i] = arma::zeros(6);
        flow_i_[i] = arma::zeros(6);
    }
}

void ViscoelasticMechanism::register_variables(InternalVariableCollection& ivc) {
    // Register viscous strain tensor for each Prony term
    for (int i = 0; i < N_prony_; ++i) {
        ivc.add_vec("EV_" + std::to_string(i), arma::zeros(6), true);
    }
}

// ========== Constitutive Computations ==========

void ViscoelasticMechanism::compute_constraints(
    const arma::vec& sigma,
    const arma::mat& L,
    double DTime,
    const InternalVariableCollection& ivc,
    arma::vec& Phi,
    arma::vec& Y_crit
) const {
    Phi.set_size(N_prony_);
    Y_crit.set_size(N_prony_);

    for (int i = 0; i < N_prony_; ++i) {
        // Get current viscous strain
        arma::vec EV_i = ivc.get("EV_" + std::to_string(i)).vec();
        EV_n_[i] = ivc.get("EV_" + std::to_string(i)).vec_start();

        // Compute exponential factor
        double exp_factor = std::exp(-DTime / tau_i_[i]);
        factor_i_[i] = exp_factor;

        // For generalized Maxwell model with implicit integration:
        // EV_i^{n+1} = exp(-dt/tau) * EV_i^n + (1 - exp(-dt/tau)) * C_i^{-1} : sigma
        //
        // The constraint equation is:
        // Phi_i = EV_i - [exp(-dt/tau) * EV_i^n + (1 - exp(-dt/tau)) * g_i * S^{-1} : sigma]
        //
        // where S^{-1} is the compliance. For isotropic elasticity with only deviatoric
        // viscoelasticity, C_i^{-1} : sigma = (1/2G) * dev(sigma)

        // Compliance (inverse of L_0)
        arma::mat S = M_0_;

        // Target viscous strain (equilibrium value scaled by g_i)
        arma::vec EV_target = g_i_[i] * (S * sigma);

        // Predicted viscous strain using implicit integration
        arma::vec EV_pred = exp_factor * EV_n_[i] + (1.0 - exp_factor) * EV_target;

        // Constraint: actual - predicted should be zero
        // We use a relaxation formulation where ds represents the strain rate
        Phi(i) = arma::norm(EV_i - EV_pred);

        // Flow direction (normalized difference)
        if (Phi(i) > sim_iota) {
            flow_i_[i] = (EV_target - EV_n_[i]) / (Phi(i) + sim_iota);
        } else {
            flow_i_[i] = arma::zeros(6);
        }

        // Critical value for convergence
        Y_crit(i) = std::max(arma::norm(EV_target), 1e-6);
    }
}

void ViscoelasticMechanism::compute_flow_directions(
    const arma::vec& sigma,
    const InternalVariableCollection& ivc,
    std::map<std::string, arma::vec>& Lambda_map
) const {
    // Compliance
    arma::mat S = M_0_;

    for (int i = 0; i < N_prony_; ++i) {
        // Flow direction for viscous strain is towards the target
        arma::vec EV_target = g_i_[i] * (S * sigma);
        arma::vec EV_i = ivc.get("EV_" + std::to_string(i)).vec();

        arma::vec direction = EV_target - EV_i;
        double norm = arma::norm(direction);
        if (norm > sim_iota) {
            Lambda_map["EV_" + std::to_string(i)] = direction / norm;
        } else {
            Lambda_map["EV_" + std::to_string(i)] = arma::zeros(6);
        }
    }
}

void ViscoelasticMechanism::compute_jacobian_contribution(
    const arma::vec& sigma,
    const arma::mat& L,
    const InternalVariableCollection& ivc,
    arma::mat& B,
    int row_offset
) const {
    // For viscoelasticity, the Jacobian contribution comes from the
    // coupling between stress and viscous strain evolution

    for (int i = 0; i < N_prony_; ++i) {
        // The diagonal term represents the "stiffness" of the constraint
        // For implicit integration: dPhi/dDs = 1 (approximately)
        B(row_offset + i, row_offset + i) = 1.0;
    }
}

arma::vec ViscoelasticMechanism::inelastic_strain(const InternalVariableCollection& ivc) const {
    arma::vec EV_total = arma::zeros(6);

    for (int i = 0; i < N_prony_; ++i) {
        EV_total += ivc.get("EV_" + std::to_string(i)).vec();
    }

    return EV_total;
}

void ViscoelasticMechanism::update(
    const arma::vec& ds,
    int offset,
    InternalVariableCollection& ivc
) {
    // For viscoelasticity, we update each branch using implicit integration
    // The viscous strain is computed directly from the implicit formula

    // Note: In the standard return mapping, ds represents multiplier increments.
    // For viscoelasticity with implicit integration, the update is typically
    // done directly based on the stress and time increment.

    // This update function is called during the Newton iteration.
    // Since viscoelasticity is often treated differently (exponential integrator),
    // we may need to adapt this approach.

    // For now, use a simple update based on the flow direction
    for (int i = 0; i < N_prony_; ++i) {
        double ds_i = ds(offset + i);

        if (std::abs(ds_i) > sim_iota) {
            arma::vec& EV_i = ivc.get("EV_" + std::to_string(i)).vec();
            EV_i += ds_i * flow_i_[i];
        }
    }
}

void ViscoelasticMechanism::tangent_contribution(
    const arma::vec& sigma,
    const arma::mat& L,
    const arma::vec& Ds,
    int offset,
    const InternalVariableCollection& ivc,
    arma::mat& Lt
) const {
    // The viscoelastic tangent modification
    // For a generalized Maxwell model, the tangent at a given time increment is:
    //
    // L_t = L_inf + sum_i [ g_i * (1 - exp(-dt/tau_i)) / (dt/tau_i) * L_0 ]
    //
    // However, for the algorithmic tangent, we need to account for the
    // implicit integration scheme.

    // For a generalized Maxwell model, the algorithmic tangent accounts
    // for the viscous strain evolution over the time increment.
    // The effective tangent reduction is:
    //   Lt = L * (1 - sum_i [ g_i * (1 - exp(-dt/tau_i)) ])
    //
    // This reduces to L_inf = (1 - sum(g_i)) * L for dt -> infinity (fully relaxed)
    // and to L for dt -> 0 (instantaneous response).

    double reduction = 0.0;
    for (int i = 0; i < N_prony_; ++i) {
        reduction += g_i_[i] * (1.0 - factor_i_[i]);
    }

    // Apply reduction to the tangent
    Lt *= (1.0 - reduction);
}

void ViscoelasticMechanism::compute_work(
    const arma::vec& sigma_start,
    const arma::vec& sigma,
    const InternalVariableCollection& ivc,
    double& Wm_r,
    double& Wm_ir,
    double& Wm_d
) const {
    Wm_r = 0.0;
    Wm_ir = 0.0;
    Wm_d = 0.0;

    // Average stress
    arma::vec sigma_avg = 0.5 * (sigma_start + sigma);

    for (int i = 0; i < N_prony_; ++i) {
        // Get viscous strain increment
        arma::vec DEV_i = ivc.get("EV_" + std::to_string(i)).delta_vec();

        // Dissipated work: sigma : dEV
        Wm_d += arma::dot(sigma_avg, DEV_i);
    }

    // For linear viscoelasticity, there's also stored energy in the springs
    // This is a simplified computation
    Wm_ir = 0.0;  // Stored in Maxwell elements (approximation)
}

double ViscoelasticMechanism::long_term_factor() const {
    double sum_g = 0.0;
    for (int i = 0; i < N_prony_; ++i) {
        sum_g += g_i_[i];
    }
    return 1.0 - sum_g;
}

} // namespace simcoon
