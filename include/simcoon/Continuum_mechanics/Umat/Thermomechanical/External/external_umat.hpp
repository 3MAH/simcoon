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

///@file external_umat.hpp
///@brief User thermomechanical subroutine for Isotropic elastic materials in 3D case
///@version 1.0

#pragma once

#include <armadillo>

/**
 * @brief User thermomechanical UMAT, same argument list and contract as the built-in kernels.
 * @details r is the heat of the increment per unit reference volume, storage included,
 * \f$ r = \Delta\gamma - \theta_{n+1} (\eta_{n+1} - \eta_n) \f$ (heat_source_increment()), and
 * drdE / drdT its derivatives: not a rate, and defined at a zero time increment.
 */
extern void umat_external_T(const arma::vec &, const arma::vec &, arma::vec &, double &, arma::mat &, arma::mat &, arma::mat &, arma::mat &, const arma::mat &, const int &, const arma::vec &, const int &, arma::vec &, const double &, const double &,const double &,const double &, double &, double &, double &, double &, double &, double &, double &, const int &, const int &, const bool &, double &);