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
    , kappa_t_(N_prony, tensor2::zeros(Tensor2Type::stress))
    , dEVtilde_dE_(arma::zeros(6, 6))
{
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

void ViscoelasticMechanism::register_variables() {
    // Resolve and cache the keys once — they're used every FB iteration in
    // compute_constraints / inelastic_strain / update.
    for (int i = 0; i < N_prony_; ++i) {
        v_key_[i]  = "v_"  + std::to_string(i);
        ev_key_[i] = "EV_" + std::to_string(i);
        ivc_.add_scalar(v_key_[i],  0.0);
        ivc_.add_vec   (ev_key_[i], arma::zeros(6), true);
    }
}

void ViscoelasticMechanism::set_reference_stiffness(const arma::mat& L_0) {
    L_0_ = L_0;
    M_0_ = arma::inv(L_0);
    // Pre-multiply (M_0 · L_i) once — both factors are frozen for the step.
    for (int i = 0; i < N_prony_; ++i) {
        M0_L_i_[i] = M_0_ * L_i_[i];
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
    // The branches are not solved by the FB system: predict() took their backward-Euler step
    // in closed form. Their rows stay in the system (uniform bookkeeping) but are always
    // satisfied: Phi = -Y_crit gives Delta s = 0.
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
    // EV_tilde = sum_i (M_0 · L_i) · EV_i  (Prony_Nfast form).
    // M0_L_i_ is pre-multiplied in set_reference_stiffness to avoid a 6x6·6x6
    // matmul on every FB iteration.
    arma::vec EV_tilde = arma::zeros(6);
    for (int i = 0; i < N_prony_; ++i) {
        EV_tilde += M0_L_i_[i] * ivc_.get(ev_key_[i]).raw_voigt();
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
    // Closed-form backward-Euler step from the (rotated) start state: the branches see the
    // total strain only, so it is final before the other mechanisms' return mapping starts.
    // A zero time increment leaves them inactive (C_i = 0).
    std::vector<arma::vec> EV_start(N_prony_);
    for (int i = 0; i < N_prony_; ++i) {
        EV_start[i] = ivc_.get(ev_key_[i]).raw_voigt_start();
    }
    const LinearViscoStep st = maxwell_parallel_step(L_0_.is_empty() ? arma::mat(arma::zeros(6, 6)) : L_0_,
                                                     L_i_, H_i_, EV_start, E_total_end, arma::zeros(6), DTime);
    dEVtilde_dE_.zeros();
    for (int i = 0; i < N_prony_; ++i) {
        InternalVariable& ev = ivc_.get(ev_key_[i]);
        ev.raw_voigt() = st.EV_i[i];
        InternalVariable& v = ivc_.get(v_key_[i]);
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
    // Dissipation of the dashpots (Prony_Nfast lines 274-288): trapezoidal
    //   W_d = sum_i 0.5 (A_i_start + A_i_end) . DEV_i
    // on the BRANCH stress A_i = L_i (eps - EV_i), the one the dashpot carries,
    // not the total stress (that overcounts by L_0/L_i and made Wm_d exceed Wm).
    // The recoverable part is closed by the orchestrator (Wm - Wm_ir - Wm_d).
    Wm_r  = 0.0;
    Wm_ir = 0.0;
    Wm_d  = 0.0;

    for (int i = 0; i < N_prony_; ++i) {
        const InternalVariable& ev = ivc_.get(ev_key_[i]);
        const arma::vec A_start = L_i_[i] * (E_start - ev.raw_voigt_start());
        const arma::vec A_end = L_i_[i] * (E_end - ev.raw_voigt());
        Wm_d += 0.5 * arma::dot(A_start + A_end, ev.delta_vec());
    }
}

} // namespace simcoon
