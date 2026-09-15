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

///@file solver_assembly.hpp
///@brief Jacobian assembly for mixed boundary conditions, and the path/output coherency check
///@version 1.0

#pragma once
#include <armadillo>
#include <vector>
#include "block.hpp"
#include "output.hpp"

namespace simcoon{

/** @addtogroup solver
 *  @{
 */

/**
 * @brief Assemble the Newton Jacobian \f$ \mathbf{K} \f$ under mixed strain/stress control.
 *
 * A stress-controlled component keeps its tangent row; a strain-controlled one is replaced by
 * \f$ \lambda \f$ on the diagonal, so the increment it asks for is the prescribed one.
 *
 * @param[in] Lt mechanical tangent operator
 * @param[out] K the 6x6 Jacobian
 * @param[in] cBC_meca per-component control flags, Voigt order [11,22,33,12,13,23]
 * @param[in] lambda diagonal stiffness given to the prescribed components
 */
void Lt_2_K(const arma::mat &, arma::mat &, const arma::Col<int> &, const double &);

/**
 * @brief Thermomechanical counterpart of Lt_2_K: the 7x7 Jacobian with the heat equation.
 *
 * @param[in] dSdE, dSdT, dQdE, dQdT the four coupled tangent blocks
 * @param[out] K the 7x7 Jacobian
 * @param[in] cBC_meca per-component mechanical control flags
 * @param[in] cBC_T temperature control flag
 * @param[in] lambda diagonal stiffness given to the prescribed components
 */
void Lth_2_K(const arma::mat &, arma::mat &, arma::mat &, arma::mat &, arma::mat &, const arma::Col<int> &, const int &, const double &);

/**
 * @brief Check a loading programme against its output settings.
 *
 * @param[in] blocks the loading blocks
 * @param[in] so the output settings
 */
void check_path_output(const std::vector<block> &, const solver_output &);

/** @} */ // end of solver group

} //namespace simcoon
