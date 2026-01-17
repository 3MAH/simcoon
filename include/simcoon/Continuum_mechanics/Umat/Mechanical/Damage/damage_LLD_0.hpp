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

#include <iostream>
#include <armadillo>

namespace simcoon{

/**
 * @file damage_LLD_0.hpp
 * @brief Mechanical damage model.
 */

/** @addtogroup umat_mechanical
 *  @{
 */

    
void umat_damage_LLD_0(const arma::vec &, const arma::vec &, arma::vec &, arma::mat &, arma::mat &, arma::vec &, const arma::mat &, const int &, const arma::vec &, const int &, arma::vec &, const double &, const double &,const double &,const double &, double &, double &, double &, double &, const int &, const int &, const bool &, const int &, double &);
    

/** @} */ // end of umat_mechanical group

} //namespace simcoon
