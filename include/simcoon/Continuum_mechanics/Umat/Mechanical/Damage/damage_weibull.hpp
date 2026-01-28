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

///@file damage_LLD_0.hpp
///@brief Ladeveze-Le Dantec model that accounts for anisotropic damage and plasticity
///@brief Damage evolution is considered as a function of time

#pragma once

#include <string>
#include <iostream>
#include <armadillo>

namespace simcoon{

/**
 * @file damage_weibull.hpp
 * @brief Mechanical damage model.
 */

/** @addtogroup umat_mechanical
 *  @{
 */

    
void umat_damage_weibull(const std::string &umat_name, const arma::vec &Etot, const arma::vec &DEtot, arma::vec &sigma, arma::mat &Lt, arma::mat &L, const arma::mat &DR, const int &nprops, const arma::vec &props, const int &nstatev, arma::vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, double &tnew_dt);
    

/** @} */ // end of umat_mechanical group

} //namespace simcoon
