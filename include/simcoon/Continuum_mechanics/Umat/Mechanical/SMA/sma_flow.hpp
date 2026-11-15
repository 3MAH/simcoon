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
 * @file sma_flow.hpp
 * @brief Shared pieces of the unified SMA models: forward flow direction, complementarity option,
 *        the increment both faces build on.
 */

#pragma once

#include <cmath>
#include <stdexcept>
#include <string>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Continuum_mechanics/Functions/criteria.hpp>
#include <simcoon/Simulation/Maths/num_solve.hpp>
#include <simcoon/Simulation/Maths/lagrange.hpp>
#include <simcoon/Continuum_mechanics/Umat/return_mapping.hpp>
#include <simcoon/Continuum_mechanics/Umat/tangent_assembly.hpp>

namespace simcoon {

/**
 * @brief Forward-transformation strain-flow direction of the unified SMA
 * models: \f$ \boldsymbol{\Lambda}^{TF}(\boldsymbol{\sigma}) =
 * H_{cur}(\boldsymbol{\sigma})\, \partial_\sigma \sigma^{Drucker} \f$ with the
 * exponential \f$ H_{cur} \f$ envelope.
 *
 * Single definition shared by unified_T (mechanical and thermomechanical) and
 * unified_TR — in particular their finite-difference flow Hessians
 * differentiate exactly this function, so the three kernels stay consistent.
 */
inline arma::vec sma_transformation_flow(const arma::vec& s, double sigmacrit,
                                         double Hmin, double Hmax, double k1,
                                         bool aniso_criteria,
                                         const arma::vec& DFA_params,
                                         double prager_b, double prager_n) {
    double sstar = Mises_stress(s) - sigmacrit;
    if (sstar < 0.) {
        sstar = 0.;
    }
    const double Hc = Hmin + (Hmax - Hmin) * (1. - std::exp(-1. * k1 * sstar));
    return aniso_criteria ? Hc * dDrucker_ani_stress(s, DFA_params, prager_b, prager_n)
                          : Hc * dDrucker_stress(s, prager_b, prager_n);
}

/**
 * @brief The complementarity function of a unified SMA kernel, checked against its integrator.
 *
 * The semi-smooth functions (ncp_is_semismooth: @c min, @c active_set) need the closest point:
 * their step sends a multiplier back to zero, which the cutting plane cannot undo since it adds
 * the transformation strain of every iterate to \f$ \boldsymbol{\varepsilon}^T \f$ and the
 * dissipated work to the increment (the reason MODUL refuses them with an Armstrong-Frederick
 * back stress). The closest point rebuilds the state from the start values at every iterate.
 *
 * @param ncp NCPMethod code (ncp_method_of)
 * @param tangent_mode kernel tangent code (parameter.hpp)
 * @param where the kernel, for the error message
 * @return the method
 * @throw std::invalid_argument for an unknown code, or a semi-smooth function under the cutting plane
 */
inline NCPMethod sma_ncp_method(const int &ncp, const int &tangent_mode, const std::string &where) {
    const NCPMethod method = ncp_method_of(ncp, where);
    if (ncp_is_semismooth(method) && !integrates_closest_point(tangent_mode)) {
        throw std::invalid_argument(where + ": ncp='min' / 'active_set' need integrator='closest_point' "
                                    "(the cutting plane keeps the transformation strain of its iterates)");
    }
    return method;
}

/// The transformation hardening functions of a unified SMA law and their slopes (sma_hardening()).
struct SMAHardening {
    double F = 0.;    ///< \f$ H_f^F \f$
    double R = 0.;    ///< \f$ H_f^R \f$
    double dF = 0.;   ///< \f$ \mathrm{d}H_f^F/\mathrm{d}\xi \f$
    double dR = 0.;   ///< \f$ \mathrm{d}H_f^R/\mathrm{d}\xi \f$
};

/**
 * @brief Smooth hardening functions of the unified SMA laws at \f$ \xi \f$ clamped to \f$ [0,1] \f$,
 *        slopes zero outside \f$ (0,1) \f$: one definition for every loop of both faces.
 *
 * \f[
 * H_f^F = \tfrac12 a_1\left(1 + \bar\xi^{n_1} - (1-\bar\xi)^{n_2}\right) + a_3, \qquad
 * H_f^R = \tfrac12 a_2\left(1 + \bar\xi^{n_3} - (1-\bar\xi)^{n_4}\right) - a_3, \qquad
 * \bar\xi = \min(\max(\xi, 0), 1)
 * \f]
 * (\f$ a_1\bar\xi + a_3 \f$, \f$ a_2\bar\xi - a_3 \f$ when both exponents are 1). The Lagrange
 * penalties alone hold \f$ \xi \f$ in \f$ [0,1] \f$; outside, the hardening is the C0 continuation of its
 * value at the bound (the hardening energy of a full cycle closes on \f$ [0,1] \f$ only, through
 * \f$ a_3 \f$), and a fractional power is NaN there. The kinks at the bounds make the local residual
 * piecewise C1.
 */
inline SMAHardening sma_hardening(double xi, double a1, double a2, double a3, double n1, double n2, double n3, double n4) {
    const double x = std::min(std::max(xi, 0.), 1.);
    const bool inside = (xi > 0.) && (xi < 1.);
    SMAHardening h;
    if ((n1 == 1.) && (n2 == 1.)) {
        h.F = a1*x + a3;
        h.dF = inside ? a1 : 0.;
    }
    else {
        h.F = 0.5*a1*(1. + std::pow(x, n1) - std::pow(1. - x, n2)) + a3;
        h.dF = inside ? 0.5*a1*(n1*std::pow(xi, n1 - 1.) + n2*std::pow(1. - xi, n2 - 1.)) : 0.;
    }
    if ((n3 == 1.) && (n4 == 1.)) {
        h.R = a2*x - a3;
        h.dR = inside ? a2 : 0.;
    }
    else {
        h.R = 0.5*a2*(1. + std::pow(x, n3) - std::pow((1. - x), n4)) - a3;
        h.dR = inside ? 0.5*a2*(n3*std::pow(xi, n3 - 1.) + n4*std::pow(1. - xi, n4 - 1.)) : 0.;
    }
    return h;
}

/**
 * @brief One increment of a unified SMA law (umat_sma_unified_T, umat_sma_unified_TR), as the
 *        mechanical and the thermomechanical faces read it after the local solve.
 *
 * The mechanisms are, in this order, the forward transformation \f$ \xi_F \f$, the reverse
 * transformation \f$ \xi_R \f$ and, for the SMR* laws, the reorientation \f$ p^{TR} \f$ (N = 2 or
 * 3). The columns of the 6xN members follow that order. Under the cutting plane the tangent
 * ingredients are those of the last iterate, under the closest point those of the converged
 * state (ReturnMappingResult). The thermodynamic forces \f$ A_j \f$ are those of the transformation
 * rows, the bound penalties of \f$ \xi \f$ included (add_bound_penalties); the reorientation row has
 * none (\f$ A_{Reo} = 0 \f$). Energy split of an increment:
 * \f[
 * \Delta\gamma = \tfrac12(\boldsymbol{\sigma}_n + \boldsymbol{\sigma}_{n+1}) : \Delta\boldsymbol{\varepsilon}^T
 *   + \sum_j \tfrac12\left(A_j^n + A_j^{n+1}\right)\Delta s_j - \Delta\Psi^{Reo}, \qquad
 * \Delta\Psi^{Reo} = \tfrac12\left(\mathbf{S}_n + \mathbf{S}_{n+1}\right) : \Delta\boldsymbol{\varepsilon}^{Reo},
 * \quad \mathbf{S} = (1 + \lambda^{Reo})\,\mathbf{X}
 * \f]
 * \f$ \Delta\gamma \f$ goes to Wm_d (and to the heat source), the rest of the stress power to Wm_r:
 * the transformation hardening, the bound penalties and the back stress of the reorientation are
 * recovered at unloading.
 */
struct SMAIncrement {
    bool closest_point = false;      ///< the closest point ran (rm, L_tan frozen at \f$ \xi_n \f$)
    ReturnMappingResult rm;          ///< closest point: the converged solve
    arma::mat L_tan;                 ///< stiffness the tangent is built on (6x6)
    arma::mat Bhat;                  ///< cutting plane: \f$ \hat{B}^{lj} \f$ of the last iterate (NxN)
    arma::mat kappa;                 ///< cutting plane: \f$ \boldsymbol{\kappa}^j \f$ (6xN)
    arma::mat dPhi;                  ///< cutting plane: \f$ \partial\Phi^l/\partial\boldsymbol{\sigma} \f$ (6xN)
    arma::vec Ds;                    ///< multipliers of the increment (N)
    arma::mat dflow_dsigma;          ///< closest point: total \f$ \mathrm{d}\boldsymbol{\Lambda}^j/\mathrm{d}\boldsymbol{\sigma} \f$ of the flows below, side by side (6x6N); empty otherwise
    arma::mat flow_chain;            ///< closest point: \f$ \sum_k \Delta s_k\,\partial\boldsymbol{\Lambda}^k/\partial\Delta s_j \f$ (6xN, ReturnStateHooks::flow_state_coupling); empty otherwise
    arma::mat flows;                 ///< strain flows \f$ \boldsymbol{\Lambda}^F, \boldsymbol{\Lambda}^R = -\bar{\boldsymbol{\varepsilon}}^T, \boldsymbol{\Lambda}^{Reo} \f$ (6xN)
    arma::vec A_start;               ///< \f$ A_j \f$ at the start of the increment (N)
    arma::vec A_end;                 ///< \f$ A_j \f$ at its end (N)
    arma::mat dA_dsigma;             ///< \f$ \partial A_j/\partial\boldsymbol{\sigma} \f$ (6xN)
    arma::mat dA_ds;                 ///< \f$ \partial A_j/\partial s_k \f$ (NxN)
    arma::vec dA_dtheta;             ///< \f$ \partial A_j/\partial\theta \f$ (N)
    arma::vec dPhi_dtheta;           ///< \f$ \partial\Phi^l/\partial\theta \f$ (N)
    arma::vec sigma_start;           ///< stress at the start of the increment
    arma::vec DET;                   ///< transformation strain increment (reorientation included)
    arma::vec DETF, DETR, DEReo;     ///< its parts: forward, reverse (stored positive), reorientation
    arma::vec alpha;                 ///< thermal expansion of the elastic prediction (\f$ \xi_n \f$)
    arma::vec Dalpha;                ///< \f$ \boldsymbol{\alpha}_M - \boldsymbol{\alpha}_A \f$
    double xi_start = 0.;            ///< \f$ \xi_n \f$
    double xi = 0.;                  ///< \f$ \xi_{n+1} \f$
    double rhoDs0 = 0.;              ///< \f$ \rho\Delta s_0 \f$
    double T_init = 0.;              ///< reference temperature (statev 0)
    double stored = 0.;              ///< \f$ \Delta\Psi^{Reo} \f$, the work of the reorientation back stress (0 for SMA*)
    arma::vec dstored_dsigma;        ///< \f$ \partial\Delta\Psi^{Reo}/\partial\boldsymbol{\sigma}_{n+1} \f$ at fixed multipliers (6); empty: zero
    arma::vec dstored_ds;            ///< \f$ \partial\Delta\Psi^{Reo}/\partial\Delta s_j \f$ at fixed stress (N); empty: zero

    /**
     * @brief Moves the bound penalties of \f$ \xi \f$ into the forces of the transformation rows:
     *        \f$ A_F \leftarrow A_F - \lambda_1(\xi) \f$, \f$ A_R \leftarrow A_R + \lambda_0(\xi) \f$
     *        (\f$ \lambda_0 = -\f$ lagrange_pow_0), at \f$ \xi_n \f$ and \f$ \xi_{n+1} \f$, with their
     *        slopes in dA_ds. The penalty stores energy like a hardening: its work is recovered.
     *        Call once A_start, A_end, dA_ds, xi_start and xi are set.
     */
    void add_bound_penalties(const double &c, const double &p0, const double &n, const double &alpha) {
        A_start(0) -= lagrange_pow_1(xi_start, c, p0, n, alpha);
        A_start(1) -= lagrange_pow_0(xi_start, c, p0, n, alpha);
        A_end(0) -= lagrange_pow_1(xi, c, p0, n, alpha);
        A_end(1) -= lagrange_pow_0(xi, c, p0, n, alpha);
        const double dl1 = dlagrange_pow_1(xi, c, p0, n, alpha);
        const double dl0 = -1.*dlagrange_pow_0(xi, c, p0, n, alpha);
        dA_ds(0, 0) -= dl1;
        dA_ds(0, 1) += dl1;
        dA_ds(1, 0) += dl0;
        dA_ds(1, 1) -= dl0;
    }

    /**
     * @brief The operator of the solve for a kernel tangent code: cpp_tangent() under the closest
     *        point; under the cutting plane the continuum operator (algorithmic clamped to it: the
     *        finite-difference flow Hessian is near-singular on the plateau), or the elastic one.
     */
    ContinuumTangent tangent(const int &tangent_mode) const {
        if (closest_point) return cpp_tangent(rm, L_tan, tangent_mode);
        const int mode = (tangent_mode == tangent_none) ? tangent_none : tangent_continuum;
        return compute_tangent_operator(mode, Bhat, kappa, dPhi, Ds, L_tan);
    }
    /// \f$ \boldsymbol{\kappa}^j \f$ of the operator (6xN)
    const arma::mat &flux() const { return closest_point ? rm.kappa_j : kappa; }
    /// \f$ \partial\Phi^l/\partial\boldsymbol{\sigma} \f$ of the operator (6xN)
    const arma::mat &gradients() const { return closest_point ? rm.dPhidsigma_l : dPhi; }
};

} // namespace simcoon
