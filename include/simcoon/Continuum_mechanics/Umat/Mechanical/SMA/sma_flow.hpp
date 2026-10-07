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
 * @brief Shared transformation-flow direction of the unified SMA models.
 */

#pragma once

#include <cmath>
#include <armadillo>
#include <simcoon/Continuum_mechanics/Functions/criteria.hpp>

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

} // namespace simcoon
