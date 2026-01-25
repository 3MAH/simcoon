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


#pragma once
#include <iostream>
#include <armadillo>

namespace simcoon{

/**
 * @file damage.hpp
 * @author Yves Chemisky 
 * @brief The damage library that computes damage evolution laws.
 */

/** @addtogroup damage
 *  @{
 */

/**
 * @brief Provides the damage evolution \f$ \delta D \f$ considering a Weibull damage law.
 * @param stress The stress vector \f$ \sigma \f$
 * @param damage The old damage \f$ D_{old} \f$
 * @param alpha The shape parameter \f$ \alpha \f$
 * @param beta The scale parameter \f$ \beta \f$
 * @param DTime The time increment \f$ \Delta T \f$
 * @param criterion The criterion (default is "vonmises")
 * @return The damage evolution \f$ \delta D \f$
 * @details It is given by :
    \f$ \Delta D = (1-D_{old})*\Big(1-exp\big(-1(\frac{crit}{\beta})^{\alpha}\big)\Big) \f$
 *   Parameters of this function are: the stress vector \f$ \sigma \f$, the old damage \f$ D_{old} \f$, the shape parameter \f$ \alpha \f$, the scale parameter \f$ \beta \f$, the time increment \f$ \Delta T \f$ and the criterion (which is a string).
 *   The criterion possibilities are :
 *   “vonmises” : \f$ crit = \sigma_{Mises} \f$
 *   “hydro” : \f$ crit = tr(\sigma) \f$
 *   “J3” : \f$ crit = J3(\sigma) \f$
 *   Default value of the criterion is “vonmises”.
 * @code 
        double varD = damage_weibull(stress, damage, alpha, beta, DTime, criterion);
 * @endcode
*/
double damage_weibull(const arma::vec &stress, const double &damage, const double &alpha, const double &beta, const double &DTime, const std::string&criterion = "vonmises");

/**
 * @brief Provides the damage evolution \f$ \delta D \f$ considering a Weibull damage law.
 * @param stress The stress vector \f$ \sigma \f$
 * @param strain The strain vector \f$ \epsilon \f$
 * @param damage The old damage \f$ D_{old} \f$
 * @param A0 The material properties characteristic of creep damage \f$ A_0 \f$
 * @param r The parameter \f$ r \f$
 * @param criterion The criterion (which is a string)
 * @return The damage evolution \f$ \delta D \f$
 * @details It is given by :
 *    \f$ \delta D = \Big(\frac{crit}{A_0(1-D_{old})}\Big)^r \f$
 *   the stress vector \f$ \sigma \f$, the strain vector \f$ \epsilon \f$, the old damage \f$ D_{old} \f$, the material properties characteristic of creep damage \f$ A_0 \f$, \f$ r \f$ and the criterion (which is a string).
 *   The criterion possibilities are :
 *   “vonmises” : \f$ crit = \sigma_{Mises} \f$
 *   “hydro” : \f$ crit = tr(\sigma) \f$
 *   “J3” : \f$ crit = J3(\sigma) \f$
 *   Here, the criterion has no default value.
 * @code 
        varD = damage_kachanov(stress, strain, damage, A0, r, criterion);
 * @endcode
*/
double damage_kachanov(const arma::vec &stress, const arma::vec &strain, const double &damage, const double &A0, const double &r, const std::string &criterion);

/**
 * @brief  Provides the constant damage evolution \f$ \Delta D \f$ considering a Woehler- Miner’s damage law.
 * @param S_max The max stress value \f$ \sigma_{Max} \f$
 * @param S_mean The mean stress value \f$\ sigma_{Mean} \f$
 * @param S_ult The “ultimate” stress value \f$ \sigma_{ult} \f$
 * @param b The parameter \f$ b \f$
 * @param B0 The parameter \f$ B_0 \f$
 * @param beta The parameter \f$ \beta \f$
 * @param Sl_0 The parameter \f$ Sl_0 \f$ (default is 0.0)
 * @return The damage evolution \f$ \delta D \f$
 * @details It is given by :
    \f$ \Delta D = \big(\frac{S_{Max}-S_{Mean}+Sl_0*(1-b*S_{Mean})}{S_{ult}-S_{Max}}\big)*\big(\frac{S_{Max}-S_{Mean}}{B_0*(1-b*S_{Mean})}\big)^\beta \f$
 *   Parameters of this function are: the max stress value \f$ \sigma_{Max} \f$, the mean stress value \f$\ sigma_{Mean} \f$, the “ultimate” stress value \f$ \sigma_{ult} \f$, the parameter \f$ b \f$, the parameter \f$ B_0 \f$, the parameter \f$ \beta \f$ and the parameter \f$ Sl_0 \f$.
 *   Default value of :math:`Sl_0` is 0.0.
 * @code 
        double varD = damage_miner(S_max, S_mean, S_ult, b, B0, beta, Sl_0);
 * @endcode
*/
double damage_miner(const double &S_max, const double &S_mean, const double &S_ult, const double &b, const double &B0, const double &beta, const double &Sl_0 = 0.);

/**
 * @brief Provides the constant damage evolution \f$ \Delta D \f$ considering a Coffin-Manson’s damage law.
 * @param S_amp The stress amplitude \f$ \sigma_{Amp} \f$
 * @param C2 The parameter \f$ C_2 \f$
 * @param gamma2 The parameter \f$ \gamma_2 \f$
 * @return The damage evolution \f$ \delta D \f$
 * @details It is given by :
    \f$ \Delta D = \big(\frac{\sigma_{Amp}}{C_{2}}\big)^{\gamma_2} \f$
 *   Parameters of this function are: the stress amplitude \f$ \sigma_{Amp} \f$, the parameter \f$ C_2 \f$ and the parameter \f$ \gamma_2 \f$.
 * @code 
        double varD = damage_manson(S_amp, C2, gamma2);
 * @endcode
*/
double damage_manson(const double &S_amp, const double &C2, const double &gamma2);

/** @} */ // end of damage group

} //namespace simcoon
