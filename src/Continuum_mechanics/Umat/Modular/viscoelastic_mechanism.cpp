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
 * @brief Port of the Prony_Nfast generalized-Maxwell kernel to the modular
 *        framework. See viscoelastic_mechanism.hpp for the physics summary.
 */

#include <simcoon/Continuum_mechanics/Umat/Modular/viscoelastic_mechanism.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/parameter.hpp>
#include <stdexcept>
#include <cmath>

namespace simcoon {

ViscoelasticMechanism::ViscoelasticMechanism(int N_prony)
    : StrainMechanism("viscoelastic")
    , N_prony_(N_prony)
    , E_i_(N_prony, 0.0)
    , nu_i_(N_prony, 0.0)
    , etaB_i_(N_prony, 0.0)
    , etaS_i_(N_prony, 0.0)
    , L_i_(N_prony)
    , H_i_(N_prony)
    , invH_i_(N_prony)
    , M0_L_i_(N_prony)
    , M_0_(arma::eye(6, 6))
    , ev_key_(N_prony)
    , v_key_(N_prony)
    , flow_i_(N_prony)
    , Lambda_i_(N_prony)
    , kappa_i_(N_prony)
    , dPhi_i_dv_(N_prony)
    , K_diag_(arma::zeros(N_prony))
{
    for (int i = 0; i < N_prony_; ++i) {
        flow_i_[i]     = arma::zeros(6);
        Lambda_i_[i]   = arma::zeros(6);
        kappa_i_[i]    = arma::zeros(6);
        dPhi_i_dv_[i]  = arma::zeros(6);
    }
}

void ViscoelasticMechanism::configure(const arma::vec& props, int& offset) {
    // Four scalars per Prony branch: E_i, nu_i, etaB_i, etaS_i.
    for (int i = 0; i < N_prony_; ++i) {
        E_i_[i]    = props(offset + 4 * i);
        nu_i_[i]   = props(offset + 4 * i + 1);
        etaB_i_[i] = props(offset + 4 * i + 2);
        etaS_i_[i] = props(offset + 4 * i + 3);

        if (etaS_i_[i] <= 0.0 || etaB_i_[i] <= 0.0) {
            throw std::runtime_error(
                "ViscoelasticMechanism: bulk and shear viscosities must be > 0");
        }

        L_i_[i]    = L_iso(E_i_[i], nu_i_[i], "Enu");
        H_i_[i]    = H_iso(etaB_i_[i], etaS_i_[i]);
        invH_i_[i] = arma::inv(H_i_[i]);
    }
    offset += 4 * N_prony_;
}

void ViscoelasticMechanism::register_variables(InternalVariableCollection& ivc) {
    // Resolve and cache the full prefixed keys once — they're used every FB
    // iteration in compute_constraints / inelastic_strain / update.
    for (int i = 0; i < N_prony_; ++i) {
        v_key_[i]  = key("v_"  + std::to_string(i));
        ev_key_[i] = key("EV_" + std::to_string(i));
        ivc.add_scalar(v_key_[i],  0.0, false);
        ivc.add_vec   (ev_key_[i], arma::zeros(6), true);
    }
}

void ViscoelasticMechanism::set_reference_stiffness(const arma::mat& L_0) {
    M_0_ = arma::inv(L_0);
    // Pre-multiply (M_0 · L_i) once — both factors are frozen for the step.
    for (int i = 0; i < N_prony_; ++i) {
        M0_L_i_[i] = M_0_ * L_i_[i];
    }
}

void ViscoelasticMechanism::compute_constraints(
    const arma::vec& /*sigma*/,
    const arma::vec& E_total,
    const arma::mat& /*L*/,
    double DTime,
    const InternalVariableCollection& ivc,
    arma::vec& Phi,
    arma::vec& Y_crit
) const {
    Phi.set_size(N_prony_);
    Y_crit.set_size(N_prony_);

    for (int i = 0; i < N_prony_; ++i) {
        const arma::vec& EV_i = ivc.get(ev_key_[i]).vec();

        // Branch flow rate: driving stress through branch viscosity.
        flow_i_[i] = invH_i_[i] * (L_i_[i] * (E_total - EV_i));
        Lambda_i_[i] = eta_norm_strain(flow_i_[i]);
        kappa_i_[i] = L_i_[i] * Lambda_i_[i];
        // dPhi_i/dv_i stress-type term (Prony_Nfast line 188)
        dPhi_i_dv_[i] = invH_i_[i] * (Lambda_i_[i] % Ir05());

        const double flow_mag = norm_strain(flow_i_[i]);
        const double Delta_v_i = ivc.get(v_key_[i]).delta_scalar();

        if (DTime > simcoon::iota) {
            Phi(i) = flow_mag - Delta_v_i / DTime;
            K_diag_(i) = -arma::dot(dPhi_i_dv_[i], kappa_i_[i]) - 1.0 / DTime;
        } else {
            Phi(i) = flow_mag;
            K_diag_(i) = -arma::dot(dPhi_i_dv_[i], kappa_i_[i]);
        }

        Y_crit(i) = std::max(flow_mag, simcoon::precision_umat);
    }
}

void ViscoelasticMechanism::compute_flow_directions(
    const arma::vec& /*sigma*/,
    const InternalVariableCollection& /*ivc*/,
    std::map<std::string, arma::vec>& Lambda_map
) const {
    // Lambda_i_ is populated by compute_constraints (CCP: frozen from previous iterate).
    for (int i = 0; i < N_prony_; ++i) {
        Lambda_map["EV_" + std::to_string(i)] = Lambda_i_[i];
    }
}

void ViscoelasticMechanism::compute_jacobian_contribution(
    const arma::vec& /*sigma*/,
    const arma::mat& /*L*/,
    const InternalVariableCollection& /*ivc*/,
    arma::mat& B,
    int row_offset
) const {
    // Diagonal entry K(i,i) cached by compute_constraints (it has DTime).
    // Cross-branch off-diagonals dPhi_i/dv_j are filled in the cross-coupling
    // pass (separate commit).
    for (int i = 0; i < N_prony_; ++i) {
        B(row_offset + i, row_offset + i) = K_diag_(i);
    }
}

const std::vector<arma::vec>& ViscoelasticMechanism::kappa(
    const arma::vec& /*sigma*/, double /*DT*/, const arma::mat& /*L_ref*/,
    const InternalVariableCollection& /*ivc*/) const {
    // kappa_i_ is the cached [L_i · Lambda_i] per branch; return it directly.
    return kappa_i_;
}

arma::vec ViscoelasticMechanism::inelastic_strain(const InternalVariableCollection& ivc) const {
    // EV_tilde = sum_i (M_0 · L_i) · EV_i  (Prony_Nfast form).
    // M0_L_i_ is pre-multiplied in set_reference_stiffness to avoid a 6x6·6x6
    // matmul on every FB iteration.
    arma::vec EV_tilde = arma::zeros(6);
    for (int i = 0; i < N_prony_; ++i) {
        EV_tilde += M0_L_i_[i] * ivc.get(ev_key_[i]).vec();
    }
    return EV_tilde;
}

void ViscoelasticMechanism::update(
    const arma::vec& ds,
    int offset,
    InternalVariableCollection& ivc
) {
    // Apply the multiplier increment to each branch:
    //   v_i   += ds_i      (lead scalar / accumulated flow length)
    //   EV_i  += ds_i * Lambda_i
    for (int i = 0; i < N_prony_; ++i) {
        const double ds_i = ds(offset + i);
        double& v_i = ivc.get(v_key_[i]).scalar();
        v_i += ds_i;
        arma::vec& EV_i = ivc.get(ev_key_[i]).vec();
        EV_i += ds_i * Lambda_i_[i];
    }
}

void ViscoelasticMechanism::tangent_contribution(
    const arma::vec& /*sigma*/,
    const arma::mat& /*L*/,
    const arma::vec& Ds,
    int offset,
    const InternalVariableCollection& /*ivc*/,
    arma::mat& Lt
) const {
    // Consistent tangent contribution, Prony_Nfast lines 231-269. Active
    // branches (Ds_i > 0) contribute -kappa_i ⊗ P_eps_i where P_eps_i is built
    // from the inverse of the active block of (-K_{ij}). For a single branch
    // or weakly-coupled multi-branch case, the diagonal approximation used
    // here is the original Prony_Nfast behavior with op(i) = 0/1 activation
    // masking.
    std::vector<double> op(N_prony_, 0.0);
    arma::mat Bhat = arma::zeros(N_prony_, N_prony_);
    arma::mat Bbar = arma::eye(N_prony_, N_prony_);

    for (int i = 0; i < N_prony_; ++i) {
        if (Ds(offset + i) > simcoon::iota) {
            op[i] = 1.0;
        }
        const double K_ii = -arma::dot(dPhi_i_dv_[i], kappa_i_[i]);
        Bhat(i, i) = -K_ii;
    }

    for (int i = 0; i < N_prony_; ++i) {
        for (int j = 0; j < N_prony_; ++j) {
            Bbar(i, j) = op[i] * op[j] * Bhat(i, j)
                       + (i == j ? 1.0 - op[i] * op[j] : 0.0);
        }
    }

    const arma::mat invBbar = arma::inv(Bbar);
    arma::mat invBhat = arma::zeros(N_prony_, N_prony_);
    for (int i = 0; i < N_prony_; ++i) {
        for (int j = 0; j < N_prony_; ++j) {
            invBhat(i, j) = op[i] * op[j] * invBbar(i, j);
        }
    }

    for (int i = 0; i < N_prony_; ++i) {
        arma::vec P_eps_i = arma::zeros(6);
        for (int j = 0; j < N_prony_; ++j) {
            P_eps_i += invBhat(j, i) * (L_i_[j] * dPhi_i_dv_[j]);
        }
        Lt -= kappa_i_[i] * P_eps_i.t();
    }
}

void ViscoelasticMechanism::compute_work(
    const arma::vec& sigma_start,
    const arma::vec& sigma,
    const InternalVariableCollection& ivc,
    double& Wm_r,
    double& Wm_ir,
    double& Wm_d
) const {
    // Dissipated work per Prony_Nfast lines 274-285:
    //   W_d = sum_i 0.5 (A_v_i_start + A_v_i) . DEV_i
    // with A_v_i = L_i . (eps - EV_i). The reversible contribution stored in
    // the branch springs is picked up by the modular orchestrator's elastic
    // work accounting (sigma . D(E - E_inel)).
    Wm_r  = 0.0;
    Wm_ir = 0.0;
    Wm_d  = 0.0;

    const arma::vec sigma_avg = 0.5 * (sigma_start + sigma);
    for (int i = 0; i < N_prony_; ++i) {
        const arma::vec DEV_i = ivc.get(ev_key_[i]).delta_vec();
        Wm_d += arma::dot(sigma_avg, DEV_i);
    }
}

} // namespace simcoon
