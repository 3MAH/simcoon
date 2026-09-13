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

/**
 * @file solver_assembly.hpp
 * @brief The file-free part of what used to live in Solver/read.hpp.
 *
 * The text readers (path.txt, material.dat, output.dat, solver_*.inp) left the library with
 * the JSON-only migration and now live in test/support/file_readers.hpp, which only the C++
 * test suite compiles. These three functions never touched a file: two assemble the Newton
 * Jacobian under mixed strain/stress control, the third validates a loading programme against
 * its output settings.
 */

/** @addtogroup solver
 *  @{
 */

/// Function that fills the matrix Tdsde for mix strain/stress conditions
void Lt_2_K(const arma::mat &, arma::mat &, const arma::Col<int> &, const double &);

/// Function that fills the matrix Tdsde for mix strain/stress conditions
void Lth_2_K(const arma::mat &, arma::mat &, arma::mat &, arma::mat &, arma::mat &, const arma::Col<int> &, const int &, const double &);

/// Function that checks the coherency between the path and the step increments provided
void check_path_output(const std::vector<block> &, const solver_output &);

/** @} */ // end of solver group

} //namespace simcoon
