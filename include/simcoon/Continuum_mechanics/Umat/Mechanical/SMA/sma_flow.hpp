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

/**
 * @brief One increment of a unified SMA law (umat_sma_unified_T, umat_sma_unified_TR), as the
 *        mechanical and the thermomechanical faces read it after the local solve.
 *
 * The mechanisms are, in this order, the forward transformation \f$ \xi_F \f$, the reverse
 * transformation \f$ \xi_R \f$ and, for the SMR* laws, the reorientation \f$ p^{TR} \f$ (N = 2 or
 * 3). The columns of the 6xN members follow that order. Under the cutting plane the tangent
 * ingredients are those of the last iterate, under the closest point those of the converged
 * state (ReturnMappingResult). The thermodynamic forces \f$ A_j \f$ are those of the transformation
 * rows; the reorientation row has none (\f$ A_{Reo} = 0 \f$, its work is the stress power on
 * \f$ \Delta\boldsymbol{\varepsilon}^{Reo} \f$).
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
