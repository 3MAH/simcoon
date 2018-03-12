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

///@file damage.hpp
///@brief Functions that computes damage evolution laws
///@version 1.0
#pragma once

#include <iostream>
#include <armadillo>

namespace simcoon{

//This function returns damage evolution (/dt) considering a Weibull damage law
double damage_weibull(const arma::vec &, const double &, const double &, const double &, const double &, const std::string& = "vonmises");

//This function returns damage evolution (/dt) considering Kachanov's creep damage law
double damage_kachanov(const arma::vec &, const arma::vec &, const double &, const double &, const double &, const std::string &);

//This function returns the constant damage evolution (/dN) considering Woehler- Miner's damage law
double damage_miner(const double &, const double &, const double &, const double &, const double &, const double &, const double & =0.);

//This function returns the constant damage evolution (/dN) considering Coffin-Manson's damage law
double damage_manson(const double &, const double &, const double &);

} //namespace simcoon
