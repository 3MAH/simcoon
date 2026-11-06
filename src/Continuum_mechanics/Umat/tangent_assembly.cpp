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

///@file tangent_assembly.cpp
///@brief Generic tangent-modulus assembly for dissipative UMATs.
///@version 1.0

#include <armadillo>
#include <vector>
#include <cassert>

#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Umat/tangent_assembly.hpp>
#include <stdexcept>

namespace simcoon {

namespace {
// Lt_dPhi = L^T maps dPhi to the strain sensitivity: Delta s = invBhat ((L^T dPhi) . Delta eps);
// the algorithmic form passes (M^-1 L)^T, not symmetric when dLambda/dsigma is not
ContinuumTangent assemble_tangent(
    const arma::mat& Bhat,
    const std::vector<arma::vec>& kappa_j,
    const std::vector<arma::vec>& dPhidsigma_l,
    const arma::vec& Ds_j,
    const arma::mat& L,
    const arma::mat& Lt_dPhi)
{
    const arma::uword Nmech = Bhat.n_rows;
    assert(Bhat.n_cols == Nmech);
    assert(kappa_j.size() == Nmech);
    assert(dPhidsigma_l.size() == Nmech);
    assert(Ds_j.n_elem == Nmech);

    arma::vec op = arma::zeros<arma::vec>(Nmech);
    for (arma::uword i = 0; i < Nmech; ++i) {
        if (Ds_j(i) > simcoon::iota) op(i) = 1.;
    }

    // Bbar = mask % Bhat + diag(1 - op): identity on the inactive rows keeps it invertible
    const arma::mat mask    = op * op.t();
    const arma::mat Bbar    = mask % Bhat + arma::diagmat(1.0 - op);
    arma::mat       invBbar;
    if (!arma::inv(invBbar, Bbar) || !invBbar.is_finite()) {
        // degenerate Bhat (e.g. SMA plateau): elastic operator rather than inf in Lt/P, which would abort the solve
        ContinuumTangent out;
        out.Lt = L;
        out.P_epsilon.assign(Nmech, arma::zeros(6));
        out.invBhat.zeros(Nmech, Nmech);
        return out;
    }
    arma::mat invBhat = mask % invBbar;   // zero out inactive rows/cols

    // elastic increment: every product below sums signed zeros, so Lt = L bit for bit; skipped when no term can be inf/NaN
    if (!any(op) && Lt_dPhi.is_finite()) {
        double max_dPhi = 0.;
        bool finite = true;
        for (arma::uword m = 0; m < Nmech; ++m) {
            finite = finite && dPhidsigma_l[m].is_finite() && kappa_j[m].is_finite();
            if (finite) max_dPhi = std::max(max_dPhi, arma::abs(dPhidsigma_l[m]).max());
        }
        if (finite && arma::abs(Lt_dPhi).max()*max_dPhi < 1.E300) {
            ContinuumTangent out;
            out.Lt = L;
            out.P_epsilon.assign(Nmech, arma::zeros(6));
            out.invBhat = std::move(invBhat);
            return out;
        }
    }

    // P = (L dPhi) invBhat^T: Bhat rows = criteria, columns = mechanisms
    arma::mat DP(6, Nmech);
    for (arma::uword m = 0; m < Nmech; ++m) DP.col(m) = dPhidsigma_l[m];
    const arma::mat P = (Lt_dPhi * DP) * invBhat.t();   // 6 \times Nmech

    // Stack \kappa^l into a 6\times Nmech matrix for the rank-Nmech update Lt = L - \kappa P^T.
    arma::mat K6(6, Nmech);
    for (arma::uword l = 0; l < Nmech; ++l) K6.col(l) = kappa_j[l];

    const arma::mat::fixed<6,6> KP = K6 * P.t();
    ContinuumTangent out;
    out.Lt = L - KP;
    out.P_epsilon.resize(Nmech);
    for (arma::uword l = 0; l < Nmech; ++l) out.P_epsilon[l] = P.col(l);
    out.invBhat = std::move(invBhat);
    return out;
}
} // namespace

ContinuumTangent assemble_continuum_tangent(
    const arma::mat& Bhat,
    const std::vector<arma::vec>& kappa_j,
    const std::vector<arma::vec>& dPhidsigma_l,
    const arma::vec& Ds_j,
    const arma::mat& L)
{
    return assemble_tangent(Bhat, kappa_j, dPhidsigma_l, Ds_j, L, L);
}

ContinuumTangent assemble_continuum_tangent(
    double Bhat_scalar,
    const arma::vec& kappa,
    const arma::vec& dPhidsigma,
    double Ds,
    const arma::mat& L)
{
    // 1-mech specialisation: skip the 1\times 1 arma::inv and the std::vector
    // brace-init that the N-mech form would otherwise pay per UMAT call.
    ContinuumTangent out;
    out.P_epsilon.assign(1, arma::zeros(6));
    out.invBhat.zeros(1, 1);

    if (Ds <= simcoon::iota || std::abs(Bhat_scalar) < simcoon::iota ||
        !std::isfinite(Bhat_scalar)) {
        // inactive or degenerate (Bhat -> 0): Lt = L, 1/Bhat would put inf in Lt
        out.Lt = L;
        return out;
    }
    const double invBhat = 1.0 / Bhat_scalar;
    out.invBhat(0, 0) = invBhat;
    out.P_epsilon[0] = invBhat * (L * dPhidsigma);
    out.Lt = L - kappa * out.P_epsilon[0].t();
    return out;
}

ContinuumTangent assemble_algorithmic_tangent(
    const arma::mat& Bhat_continuum,
    const std::vector<arma::vec>& kappa_j,
    const std::vector<arma::vec>& dPhidsigma_l,
    const arma::vec& Ds_j,
    const arma::mat& L,
    const std::vector<arma::mat>& dLambda_dsigma_l)
{
    const arma::uword Nmech = Bhat_continuum.n_rows;
    assert(Bhat_continuum.n_cols == Nmech);
    assert(kappa_j.size() == Nmech);
    assert(dPhidsigma_l.size() == Nmech);
    assert(Ds_j.n_elem == Nmech);
    assert(dLambda_dsigma_l.size() == Nmech);

    // M = I + L \cdot \sum_j (\Delta s^j \cdot \partial \Lambda^j/\partial \sigma ), summed over the active set only.
    const arma::mat::fixed<6,6> I6 = arma::eye(6, 6);
    arma::mat::fixed<6,6> sumDsLambdaSigma = arma::zeros(6, 6);
    for (arma::uword j = 0; j < Nmech; ++j) {
        if (Ds_j(j) > simcoon::iota) {
            sumDsLambdaSigma += Ds_j(j) * dLambda_dsigma_l[j];
        }
    }

    // no dLambda/dsigma correction: M = I, the algorithmic tangent is the continuum one
    if (arma::approx_equal(sumDsLambdaSigma, arma::zeros(6, 6), "absdiff", simcoon::iota)) {
        return assemble_continuum_tangent(Bhat_continuum, kappa_j, dPhidsigma_l, Ds_j, L);
    }

    // singular M (FD dLambda/dsigma near activation, LAPACK-dependent): continuum fallback
    arma::mat::fixed<6,6> M = L * sumDsLambdaSigma;
    M += I6;
    arma::mat::fixed<6,6> Minv;
    if (!arma::inv(Minv, M)) {
        return assemble_continuum_tangent(Bhat_continuum, kappa_j, dPhidsigma_l, Ds_j, L);
    }

    // Algorithmic operators:  \tilde{L} = M^{-1} L,  \tilde{\kappa}^j = M^{-1} \kappa^j,
    // \tilde{B}^{lj} = Bhat_continuum^{lj} + \partial \Phi^l \cdot (M^{-1} - I) \cdot \kappa^j  (Bhat_continuum = \partial \Phi \cdot \kappa - K).
    const arma::mat::fixed<6,6> L_tilde = Minv * L;
    std::vector<arma::vec> kappa_tilde(Nmech);
    for (arma::uword j = 0; j < Nmech; ++j) {
        kappa_tilde[j] = Minv * kappa_j[j];
    }
    arma::mat Bhat_tilde = Bhat_continuum;
    const arma::mat::fixed<6,6> Minv_minus_I = Minv - I6;
    for (arma::uword j = 0; j < Nmech; ++j) {
        const arma::vec::fixed<6> dkappa_j = Minv_minus_I * kappa_j[j];   // shared by every row l
        for (arma::uword l = 0; l < Nmech; ++l) {
            Bhat_tilde(l, j) += arma::dot(dPhidsigma_l[l], dkappa_j);
        }
    }

    // L~ = M^-1 L is not symmetric: its transpose acts on dPhi
    const arma::mat::fixed<6,6> L_tilde_t = L_tilde.t();
    return assemble_tangent(Bhat_tilde, kappa_tilde, dPhidsigma_l, Ds_j, L_tilde, L_tilde_t);
}

ContinuumTangent assemble_algorithmic_tangent(
    double Bhat_scalar,
    const arma::vec& kappa,
    const arma::vec& dPhidsigma,
    double Ds,
    const arma::mat& L,
    const arma::mat& dLambda_dsigma)
{
    // Inactive mechanism OR no \sigma -correction -> continuum result.
    if (Ds <= simcoon::iota ||
        arma::approx_equal(dLambda_dsigma, arma::zeros(6, 6), "absdiff", simcoon::iota))
    {
        return assemble_continuum_tangent(Bhat_scalar, kappa, dPhidsigma, Ds, L);
    }
    arma::mat Minv;
    if (!arma::inv(Minv, arma::eye(6, 6) + L * (Ds * dLambda_dsigma))) {
        // Ill-conditioned FD Hessian (backend-dependent detection): fall back
        // to the continuum operator instead of aborting the solve.
        return assemble_continuum_tangent(Bhat_scalar, kappa, dPhidsigma, Ds, L);
    }
    const arma::mat L_tilde  = Minv * L;
    const arma::vec kappa_t  = Minv * kappa;
    const double    Bhat_t   = Bhat_scalar + arma::dot(dPhidsigma, (Minv - arma::eye(6, 6)) * kappa);
    // 1-mech form of the assembly with the transpose of L~ on dPhi (L~ is not symmetric when
    // dLambda/dsigma is not); same masking as the continuum overload.
    ContinuumTangent out;
    out.P_epsilon.assign(1, arma::zeros(6));
    out.invBhat.zeros(1, 1);
    if (std::abs(Bhat_t) < simcoon::iota || !std::isfinite(Bhat_t)) {
        out.Lt = L_tilde;
        return out;
    }
    out.invBhat(0, 0) = 1.0 / Bhat_t;
    out.P_epsilon[0] = out.invBhat(0, 0) * (L_tilde.t() * dPhidsigma);
    out.Lt = L_tilde - kappa_t * out.P_epsilon[0].t();
    return out;
}



ContinuumTangent compute_tangent_operator(
    int tangent_mode,
    const arma::mat& Bhat,
    const std::vector<arma::vec>& kappa_j,
    const std::vector<arma::vec>& dPhidsigma_l,
    const arma::vec& Ds_j,
    const arma::mat& L,
    const HessianProvider& dLambda_dsigma) {

    if (tangent_mode == simcoon::tangent_none) {
        // Explicit integration: elastic operator, zeroed sensitivities so the
        // thermomechanical P_theta algebra degrades consistently.
        ContinuumTangent ct;
        ct.Lt = L;
        ct.P_epsilon.assign(kappa_j.size(), arma::zeros(6));
        ct.invBhat = arma::zeros(Bhat.n_rows, Bhat.n_cols);
        return ct;
    }
    if (tangent_mode == simcoon::tangent_continuum) {
        return assemble_continuum_tangent(Bhat, kappa_j, dPhidsigma_l, Ds_j, L);
    }
    if (tangent_mode == simcoon::tangent_algorithmic || tangent_mode == simcoon::tangent_closest_point) {
        // closest-point kernels use cpp_consistent_tangent() and never come here; for the others mode 3 = algorithmic
        bool active = false;
        for (arma::uword j = 0; j < Ds_j.n_elem; ++j) {
            active = active || (Ds_j(j) > simcoon::iota);
        }
        if (!dLambda_dsigma || !active) {
            // stress-independent flow or no active mechanism: the algorithmic operator IS the continuum one
            return assemble_continuum_tangent(Bhat, kappa_j, dPhidsigma_l, Ds_j, L);
        }
        return assemble_algorithmic_tangent(Bhat, kappa_j, dPhidsigma_l, Ds_j, L,
                                            dLambda_dsigma());
    }
    throw std::invalid_argument("compute_tangent_operator: unknown tangent_mode "
                                + std::to_string(tangent_mode));
}

ContinuumTangent compute_tangent_operator(
    int tangent_mode,
    double Bhat_scalar,
    const arma::vec& kappa,
    const arma::vec& dPhidsigma,
    double Ds,
    const arma::mat& L,
    const std::function<arma::mat()>& dLambda_dsigma) {

    HessianProvider provider = nullptr;
    if (dLambda_dsigma) {
        provider = [&dLambda_dsigma]() {
            return std::vector<arma::mat>{dLambda_dsigma()};
        };
    }
    arma::mat Bhat(1, 1);
    Bhat(0, 0) = Bhat_scalar;
    return compute_tangent_operator(tangent_mode, Bhat, {kappa}, {dPhidsigma},
                                    arma::vec{Ds}, L, provider);
}

void heat_source_increment(const double &Dgamma, const arma::vec &dDgamma_dE, const double &dDgamma_dT,
                           const double &Deta, const arma::vec &deta_dE, const double &deta_dT,
                           const double &T1, double &r, arma::mat &drdE, arma::mat &drdT)
{
    r = Dgamma - T1*Deta;
    drdE = dDgamma_dE - T1*deta_dE;
    drdT = dDgamma_dT - Deta - T1*deta_dT;   // 1x1
}

} // namespace simcoon
