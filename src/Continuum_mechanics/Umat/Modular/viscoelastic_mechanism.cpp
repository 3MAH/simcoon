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
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Viscoelasticity/linear_viscoelastic.hpp>
#include <stdexcept>

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
    , M0_L_i_(N_prony)
    , ev_key_(N_prony)
    , v_key_(N_prony)
    , ev_index_(N_prony, 0)
    , v_index_(N_prony, 0)
    , kappa_t_(N_prony, tensor2::zeros(Tensor2Type::stress))
    , dEVtilde_dE_(arma::zeros(6, 6))
{
    L_0_.zeros(6, 6);   // until set_reference_stiffness: an unset reference spring is zero
}

void ViscoelasticMechanism::configure(const arma::vec& props, int& offset) {
    // Four scalars per Prony branch: E_i, nu_i, etaB_i, etaS_i.
    for (int i = 0; i < N_prony_; ++i) {
        E_i_[i]    = props(offset + 4 * i);
        nu_i_[i]   = props(offset + 4 * i + 1);
        etaB_i_[i] = props(offset + 4 * i + 2);
        etaS_i_[i] = props(offset + 4 * i + 3);

        if (etaS_i_[i] <= 0.0 || etaB_i_[i] <= 0.0) {
            throw std::invalid_argument(
                "ViscoelasticMechanism: bulk and shear viscosities must be > 0");
        }

        L_i_[i]    = L_iso(E_i_[i], nu_i_[i], "Enu");
        H_i_[i]    = H_iso(etaB_i_[i], etaS_i_[i]);
    }
    offset += 4 * N_prony_;
}

void ViscoelasticMechanism::register_variables() {
    for (int i = 0; i < N_prony_; ++i) {
        v_key_[i]  = "v_"  + std::to_string(i);
        ev_key_[i] = "EV_" + std::to_string(i);
        ivc_.add_scalar(v_key_[i],  0.0);
        ivc_.add_vec   (ev_key_[i], arma::zeros(6), true);
        v_index_[i]  = ivc_.index_of(v_key_[i]);
        ev_index_[i] = ivc_.index_of(ev_key_[i]);
    }
}

void ViscoelasticMechanism::set_reference_stiffness(const arma::mat& L_0) {
    L_0_ = L_0;
    const arma::mat M_0 = arma::inv(L_0);
    // Pre-multiply (M_0 · L_i) once — both factors are frozen for the step.
    for (int i = 0; i < N_prony_; ++i) {
        M0_L_i_[i] = M_0 * L_i_[i];
    }
}

void ViscoelasticMechanism::compute_constraints(
    const arma::vec& /*sigma*/,
    const arma::vec& /*E_total*/,
    const arma::mat& /*L*/,
    double /*DTime*/,
    arma::vec& Phi,
    arma::vec& Y_crit
) const {
    // branches solved in closed form by predict(): rows kept but always satisfied (Phi = -Y_crit, Delta s = 0)
    Phi.set_size(N_prony_);
    Y_crit.set_size(N_prony_);
    Y_crit.fill(1.0);
    Phi.fill(-1.0);
}

void ViscoelasticMechanism::compute_jacobian_contribution(
    const arma::vec& /*sigma*/,
    const arma::mat& /*L*/,
    arma::mat& B,
    int row_offset
) const {
    // Unit diagonal for the always-satisfied rows (no coupling: the fluxes are zero).
    for (int i = 0; i < N_prony_; ++i) {
        B(row_offset + i, row_offset + i) = -1.0;
    }
}

const std::vector<tensor2>& ViscoelasticMechanism::kappa(
    const arma::vec& /*sigma*/, double /*DT*/, const arma::mat& /*L_ref*/) const {
    return kappa_t_;   // zero: the rows carry no multiplier
}

arma::vec ViscoelasticMechanism::inelastic_strain() const {
    // EV_tilde = sum_i (M_0 L_i) EV_i (Prony_Nfast form); M0_L_i_ precomputed in set_reference_stiffness
    arma::vec EV_tilde = arma::zeros(6);
    for (int i = 0; i < N_prony_; ++i) {
        EV_tilde += M0_L_i_[i] * ivc_.at(ev_index_[i]).raw_voigt();
    }
    return EV_tilde;
}

void ViscoelasticMechanism::update(
    const arma::vec& /*ds*/,
    int /*offset*/
) {
    // Nothing: the branches took their step in predict(), before the elastic prediction.
}

void ViscoelasticMechanism::predict(const arma::vec& E_total_end, double DTime) {
    // closed-form backward Euler from the rotated start, final before the other mechanisms; DTime = 0 -> inactive
    std::vector<arma::vec> EV_start(N_prony_);
    for (int i = 0; i < N_prony_; ++i) {
        EV_start[i] = ivc_.at(ev_index_[i]).raw_voigt_start();
    }
    const LinearViscoStep st = maxwell_parallel_step(L_0_, L_i_, H_i_, EV_start, E_total_end, arma::zeros(6), DTime);
    dEVtilde_dE_.zeros();
    for (int i = 0; i < N_prony_; ++i) {
        InternalVariable& ev = ivc_.at(ev_index_[i]);
        ev.raw_voigt() = st.EV_i[i];
        InternalVariable& v = ivc_.at(v_index_[i]);
        v.scalar() = v.scalar_start() + norm_strain(st.EV_i[i] - EV_start[i]);
        dEVtilde_dE_ += M0_L_i_[i] * st.dEVdE_i[i];
    }
}

arma::mat ViscoelasticMechanism::total_strain_map() const {
    return arma::eye(6, 6) - dEVtilde_dE_;
}

void ViscoelasticMechanism::tangent_contribution(
    const arma::vec& /*sigma*/,
    const arma::mat& /*L*/,
    const arma::vec& /*Ds*/,
    int /*offset*/,
    arma::mat& /*Lt*/
) const {
    // Nothing here: the branches enter the tangent through total_strain_map(), applied by
    // the orchestrator after every other contribution.
}

void ViscoelasticMechanism::compute_work(
    const arma::vec& /*sigma_start*/,
    const arma::vec& /*sigma*/,
    const arma::vec& E_start,
    const arma::vec& E_end,
    double& Wm_r,
    double& Wm_ir,
    double& Wm_d
) const {
    // W_d = sum_i 1/2 (A_i_start + A_i_end) . DEV_i on the BRANCH stress A_i = L_i (eps - EV_i)
    // (the total stress overcounts by L_0/L_i); Wm_r is closed by the orchestrator
    Wm_r  = 0.0;
    Wm_ir = 0.0;
    Wm_d  = 0.0;

    for (int i = 0; i < N_prony_; ++i) {
        const InternalVariable& ev = ivc_.at(ev_index_[i]);
        const arma::vec A_start = L_i_[i] * (E_start - ev.raw_voigt_start());
        const arma::vec A_end = L_i_[i] * (E_end - ev.raw_voigt());
        Wm_d += 0.5 * arma::dot(A_start + A_end, ev.delta_vec());
    }
}

} // namespace simcoon
