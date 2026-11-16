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
 * @file sma.hpp
 * @brief Constitutive functions shared by the unified SMA laws (SMA*, SMR*, mechanical and
 *        thermomechanical): forward flow direction, transformation hardening, complementarity option.
 */

#pragma once

#include <string>
#include <armadillo>
#include <simcoon/Simulation/Maths/num_solve.hpp>

namespace simcoon {

/**
 * @brief Forward-transformation strain-flow direction of the unified SMA
 * models: \f$ \boldsymbol{\Lambda}^{TF}(\boldsymbol{\sigma}) =
 * H_{cur}(\boldsymbol{\sigma})\, \partial_\sigma \sigma^{Drucker} \f$ with the
 * exponential \f$ H_{cur} \f$ envelope.
 *
 * Shared by the unified SMA kernels (SMA*, SMR*, mechanical and thermomechanical).
 */
arma::vec sma_transformation_flow(const arma::vec &s, const double &sigmacrit, const double &Hmin,
                                  const double &Hmax, const double &k1, const bool &aniso_criteria,
                                  const arma::vec &DFA_params, const double &prager_b, const double &prager_n);

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
NCPMethod sma_ncp_method(const int &ncp, const int &tangent_mode, const std::string &where);

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
SMAHardening sma_hardening(const double &xi, const double &a1, const double &a2, const double &a3,
                           const double &n1, const double &n2, const double &n3, const double &n4);

} // namespace simcoon
