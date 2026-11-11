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
#include <simcoon/Simulation/Maths/small_lu.hpp>
#include <simcoon/Simulation/Maths/small_dense.hpp>
#include <stdexcept>

namespace simcoon {

namespace {
// Lt_dPhi = L^T maps dPhi to the strain sensitivity: Delta s = invBhat ((L^T dPhi) . Delta eps);
// the algorithmic form passes (M^-1 L)^T, not symmetric when dLambda/dsigma is not
ContinuumTangent assemble_tangent(
    const arma::mat& Bhat,
    const arma::mat& kappa_j,
    const arma::mat& dPhidsigma_l,
    const arma::vec& Ds_j,
    const arma::mat& L,
    const arma::mat& Lt_dPhi)
{
    const arma::uword Nmech = Bhat.n_rows;
    assert(Bhat.n_cols == Nmech);
    assert(kappa_j.n_rows == 6 && kappa_j.n_cols == Nmech);
    assert(dPhidsigma_l.n_rows == 6 && dPhidsigma_l.n_cols == Nmech);
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
        out.L_tilde = L;
        out.P_epsilon.zeros(6, Nmech);
        out.invBhat.zeros(Nmech, Nmech);
        return out;
    }
    arma::mat invBhat = mask % invBbar;   // zero out inactive rows/cols

    // elastic increment: every product below sums signed zeros, so Lt = L bit for bit; skipped when no term can be inf/NaN
    if (!any(op) && Lt_dPhi.is_finite()) {
        const bool finite = dPhidsigma_l.is_finite() && kappa_j.is_finite();
        const double max_dPhi = finite ? arma::abs(dPhidsigma_l).max() : 0.;
        if (finite && arma::abs(Lt_dPhi).max()*max_dPhi < 1.E300) {
            ContinuumTangent out;
            out.Lt = L;
            out.L_tilde = L;
            out.P_epsilon.zeros(6, Nmech);
            out.invBhat = std::move(invBhat);
            return out;
        }
    }

    // P^l = sum_m invBhat(l,m) Lt_dPhi dPhi^m (Bhat rows = criteria, columns = mechanisms), then Lt = L - sum_l kappa^l P^l^T
    assert(L.n_rows == 6 && L.n_cols == 6 && Lt_dPhi.n_rows == 6 && Lt_dPhi.n_cols == 6);
    arma::mat LD(6, Nmech);
    for (arma::uword m = 0; m < Nmech; ++m) mul6(Lt_dPhi.memptr(), dPhidsigma_l.colptr(m), LD.colptr(m));
    ContinuumTangent out;
    out.L_tilde = L;
    out.P_epsilon.zeros(6, Nmech);
    for (arma::uword l = 0; l < Nmech; ++l) {
        double* p = out.P_epsilon.colptr(l);
        for (arma::uword m = 0; m < Nmech; ++m) {
            const double b = invBhat(l, m);
            const double* d = LD.colptr(m);
            for (int i = 0; i < 6; ++i) p[i] += d[i]*b;
        }
    }
    out.Lt = L;
    for (arma::uword j = 0; j < 6; ++j)
        for (arma::uword i = 0; i < 6; ++i) {
            double s = 0.;
            for (arma::uword l = 0; l < Nmech; ++l) s += kappa_j(i, l)*out.P_epsilon(j, l);
            out.Lt(i, j) -= s;
        }
    out.invBhat = std::move(invBhat);
    return out;
}
} // namespace

ContinuumTangent assemble_continuum_tangent(
    const arma::mat& Bhat,
    const arma::mat& kappa_j,
    const arma::mat& dPhidsigma_l,
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
    // 1-mech specialisation: no 1x1 inverse
    ContinuumTangent out;
    out.L_tilde = L;
    out.P_epsilon.zeros(6, 1);
    out.invBhat.zeros(1, 1);

    if (Ds <= simcoon::iota || std::abs(Bhat_scalar) < simcoon::iota ||
        !std::isfinite(Bhat_scalar)) {
        // inactive or degenerate (Bhat -> 0): Lt = L, 1/Bhat would put inf in Lt
        out.Lt = L;
        return out;
    }
    const double invBhat = 1.0 / Bhat_scalar;
    out.invBhat(0, 0) = invBhat;
    out.P_epsilon = invBhat * (L * dPhidsigma);
    out.Lt = L - kappa * out.P_epsilon.t();
    return out;
}

ContinuumTangent assemble_algorithmic_tangent(
    const arma::mat& Bhat_continuum,
    const arma::mat& kappa_j,
    const arma::mat& dPhidsigma_l,
    const arma::vec& Ds_j,
    const arma::mat& L,
    const arma::mat& dLambda_dsigma_l)
{
    const arma::uword Nmech = Bhat_continuum.n_rows;
    assert(Bhat_continuum.n_cols == Nmech);
    assert(kappa_j.n_cols == Nmech);
    assert(dPhidsigma_l.n_cols == Nmech);
    assert(Ds_j.n_elem == Nmech);
    assert(dLambda_dsigma_l.n_rows == 6 && dLambda_dsigma_l.n_cols == 6 * Nmech);

    // M = I + L \cdot \sum_j (\Delta s^j \cdot \partial \Lambda^j/\partial \sigma ), summed over the active set only.
    const arma::mat::fixed<6,6> I6 = arma::eye(6, 6);
    arma::mat::fixed<6,6> sumDsLambdaSigma = arma::zeros(6, 6);
    for (arma::uword j = 0; j < Nmech; ++j) {
        if (Ds_j(j) > simcoon::iota) {
            sumDsLambdaSigma += Ds_j(j) * dLambda_dsigma_l.cols(6 * j, 6 * j + 5);
        }
    }

    // no dLambda/dsigma correction: M = I, the algorithmic tangent is the continuum one
    if (sumDsLambdaSigma.is_zero(simcoon::iota)) {
        return assemble_continuum_tangent(Bhat_continuum, kappa_j, dPhidsigma_l, Ds_j, L);
    }

    // singular M (FD dLambda/dsigma near activation, LAPACK-dependent): continuum fallback
    arma::mat::fixed<6,6> M = mul66(L, sumDsLambdaSigma);
    M += I6;
    small_lu<6> lu;
    if (!lu.factor(M)) {
        return assemble_continuum_tangent(Bhat_continuum, kappa_j, dPhidsigma_l, Ds_j, L);
    }

    // Algorithmic operators:  \tilde{L} = M^{-1} L,  \tilde{\kappa}^j = M^{-1} \kappa^j,
    // \tilde{B}^{lj} = Bhat_continuum^{lj} + \partial \Phi^l \cdot (M^{-1} - I) \cdot \kappa^j  (Bhat_continuum = \partial \Phi \cdot \kappa - K).
    arma::mat::fixed<6,6> L_tilde = L;
    lu.solve_columns_in_place(L_tilde);
    arma::mat kappa_tilde = kappa_j;
    lu.solve_columns_in_place(kappa_tilde);
    arma::mat Bhat_tilde = Bhat_continuum;
    for (arma::uword j = 0; j < Nmech; ++j) {
        const arma::vec::fixed<6> dkappa_j = kappa_tilde.col(j) - kappa_j.col(j);   // (M^-1 - I) kappa_j, shared by every row l
        for (arma::uword l = 0; l < Nmech; ++l) {
            const arma::vec::fixed<6> dPhi_l = dPhidsigma_l.col(l);
            Bhat_tilde(l, j) += arma::dot(dPhi_l, dkappa_j);
        }
    }

    // L~ = M^-1 L is not symmetric: its transpose acts on dPhi
    const arma::mat::fixed<6,6> L_tilde_t = L_tilde.t();
    ContinuumTangent out = assemble_tangent(Bhat_tilde, kappa_tilde, dPhidsigma_l, Ds_j, L_tilde, L_tilde_t);
    out.kappa_tilde = std::move(kappa_tilde);
    return out;
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
        dLambda_dsigma.is_zero(simcoon::iota))
    {
        return assemble_continuum_tangent(Bhat_scalar, kappa, dPhidsigma, Ds, L);
    }
    // singular M (ill-conditioned FD Hessian): continuum fallback
    const arma::mat::fixed<6,6> DsdLambda = Ds * dLambda_dsigma;
    arma::mat::fixed<6,6> M = mul66(L, DsdLambda);
    M += arma::eye(6, 6);
    small_lu<6> lu;
    if (!lu.factor(M)) {
        return assemble_continuum_tangent(Bhat_scalar, kappa, dPhidsigma, Ds, L);
    }
    arma::mat::fixed<6,6> L_tilde = L;
    lu.solve_columns_in_place(L_tilde);
    arma::vec::fixed<6> kappa_t = kappa;
    lu.solve_in_place(kappa_t.memptr());
    const double Bhat_t = Bhat_scalar + arma::dot(dPhidsigma, kappa_t - kappa);   // dPhi . (M^-1 - I) kappa
    // 1-mech form of the assembly with the transpose of L~ on dPhi (L~ is not symmetric when
    // dLambda/dsigma is not); same masking as the continuum overload.
    ContinuumTangent out;
    out.L_tilde = L_tilde;
    out.kappa_tilde = kappa_t;
    out.P_epsilon.zeros(6, 1);
    out.invBhat.zeros(1, 1);
    if (std::abs(Bhat_t) < simcoon::iota || !std::isfinite(Bhat_t)) {
        out.Lt = L_tilde;
        return out;
    }
    out.invBhat(0, 0) = 1.0 / Bhat_t;
    out.P_epsilon = out.invBhat(0, 0) * (L_tilde.t() * dPhidsigma);
    out.Lt = L_tilde - kappa_t * out.P_epsilon.t();
    return out;
}

ContinuumTangent compute_tangent_operator(
    int tangent_mode,
    const arma::mat& Bhat,
    const arma::mat& kappa_j,
    const arma::mat& dPhidsigma_l,
    const arma::vec& Ds_j,
    const arma::mat& L,
    const HessianProvider& dLambda_dsigma) {

    if (tangent_mode == simcoon::tangent_none) {
        // Explicit integration: elastic operator, zeroed sensitivities so the
        // thermomechanical P_theta algebra degrades consistently.
        ContinuumTangent ct;
        ct.Lt = L;
        ct.L_tilde = L;
        ct.P_epsilon.zeros(6, kappa_j.n_cols);
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

    arma::mat Bhat(1, 1);
    Bhat(0, 0) = Bhat_scalar;
    return compute_tangent_operator(tangent_mode, Bhat, kappa, dPhidsigma,
                                    arma::vec{Ds}, L, dLambda_dsigma);
}

void heat_source_increment(const double &Dgamma, const arma::vec &dDgamma_dE, const double &dDgamma_dT,
                           const double &Deta, const arma::vec &deta_dE, const double &deta_dT,
                           const double &T1, double &r, arma::vec &drdE, double &drdT)
{
    r = Dgamma - T1*Deta;
    drdE = dDgamma_dE - T1*deta_dE;
    drdT = dDgamma_dT - Deta - T1*deta_dT;
}

} // namespace simcoon
