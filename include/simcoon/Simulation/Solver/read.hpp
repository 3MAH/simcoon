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

///@file read.hpp
///@brief Solver utility functions for mixed boundary conditions
///@version 2.0
///
///@note The primary configuration method for v2.0 is JSON via the Python API.
///      Legacy functions (read_matprops, read_path) are retained for internal
///      test compatibility but deprecated for new code.

#pragma once
#include <armadillo>
#include <string>
#include "block.hpp"
#include "output.hpp"

namespace simcoon{

/**
 * @file read.hpp
 * @brief Solver utility functions.
 */

/** @addtogroup solver
 *  @{
 */

/**
 * @brief Read material properties from a legacy .dat file.
 * @deprecated Use Python JSON API instead for new code.
 * @param umat_name Output: UMAT name
 * @param nprops Output: Number of properties
 * @param props Output: Properties vector
 * @param nstatev Output: Number of state variables
 * @param psi_rve Output: First Euler angle
 * @param theta_rve Output: Second Euler angle
 * @param phi_rve Output: Third Euler angle
 * @param path_data Directory containing the file
 * @param materialfile Filename
 */
void read_matprops(std::string &umat_name, unsigned int &nprops, arma::vec &props, unsigned int &nstatev,
                   double &psi_rve, double &theta_rve, double &phi_rve,
                   const std::string &path_data, const std::string &materialfile);

arma::Col<int> subdiag2vec();

/// Function that fills the matrix K for mixed strain/stress boundary conditions
void Lt_2_K(const arma::mat &, arma::mat &, const arma::Col<int> &, const double &);

/// Function that fills the matrix K for mixed strain/stress/thermal boundary conditions
void Lth_2_K(const arma::mat &, arma::mat &, arma::mat &, arma::mat &, arma::mat &, const arma::Col<int> &, const int &, const double &);

/// Function that checks the coherency between the path and the step increments provided
void check_path_output(const std::vector<block> &, const solver_output &);


/** @} */ // end of solver group

} //namespace simcoon
