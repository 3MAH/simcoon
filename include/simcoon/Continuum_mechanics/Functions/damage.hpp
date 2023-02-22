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
* @section The damage library that computes damage evolution laws
*/

//This function returns damage evolution (/dt) considering a Weibull damage law

/**
 * @brief Provides the damage evolution :math:`\delta D` considering a Weibull damage law.
 *    It is given by : :math:`\delta D = (1-D_{old})*\Big(1-exp\big(-1(\frac{crit}{\beta})^{\alpha}\big)\Big)`
 *   Parameters of this function are: the stress vector :math:`\sigma`, the old damage :math:`D_{old}`, the shape parameter :math:`\alpha`, the scale parameter :math:`\beta`, the time increment :math:`\Delta T` and the criterion (which is a string).
 *   The criterion possibilities are :
 *   “vonmises” : :math:`crit = \sigma_{Mises}`
 *   “hydro” : :math:`crit = tr(\sigma)`
 *   “J3” : :math:`crit = J3(\sigma)`
 *   Default value of the criterion is “vonmises”.
 * 
 * 
 * @param m
 * @return The 3x3 deviatoric part of the matrix input (arma::mat)
 * @details Example: 
 * @code 
        mat m = randu(3,3;)
        mat m_dev = dev(m);
 * @endcode
*/
double damage_weibull(const arma::vec &stress, const double &damage, const double &alpha, const double &beta, const double &DTime, const std::string&criterion = "vonmises");

//This function returns damage evolution (/dt) considering Kachanov's creep damage law
double damage_kachanov(const arma::vec &, const arma::vec &, const double &, const double &, const double &, const std::string &);

//This function returns the constant damage evolution (/dN) considering Woehler- Miner's damage law
double damage_miner(const double &, const double &, const double &, const double &, const double &, const double &, const double & =0.);

//This function returns the constant damage evolution (/dN) considering Coffin-Manson's damage law
double damage_manson(const double &, const double &, const double &);

} //namespace simcoon
