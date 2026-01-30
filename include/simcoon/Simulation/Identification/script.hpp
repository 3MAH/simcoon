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

///@file script.hpp
///@brief Identification script functions - DEPRECATED in v2.0
///@version 2.0
///
///@note The C++ identification workflow that used legacy text files (path.txt,
///      material.dat, solver_essentials.inp, solver_control.inp) has been deprecated
///      in v2.0. Use Python-based identification workflows instead.
///
/// For parameter identification, use scipy.optimize with the Python solver:
/// @code
/// from simcoon.solver import Solver, Block, StepMeca
/// from scipy.optimize import minimize
///
/// def objective(params):
///     props = np.array([params[0], params[1]])  # E, nu
///     block = Block(steps=[...], umat_name="ELISO", props=props, nstatev=1)
///     solver = Solver(blocks=[block])
///     history = solver.solve()
///     # Compare with experimental data
///     return error
///
/// result = minimize(objective, x0=[210000, 0.3], method='Nelder-Mead')
/// @endcode

#pragma once

#include <iostream>
#include <armadillo>
#include "constants.hpp"
#include "parameters.hpp"
#include "opti_data.hpp"
#include "individual.hpp"
#include "generation.hpp"

namespace simcoon{

/**
 * @file script.hpp
 * @brief Parameter identification functions - DEPRECATED in v2.0.
 *
 * @note Use Python-based identification workflows instead.
 *       The legacy C++ identification used text file formats that are no longer supported.
 */

/** @addtogroup identification
 *  @{
 */


// File copy/apply utilities - still available for general use
void copy_parameters(const std::vector<parameters> &, const std::string &, const std::string &);
void copy_constants(const std::vector<constants> &, const std::string &, const std::string &);
void apply_parameters(const std::vector<parameters> &, const std::string &);
void apply_constants(const std::vector<constants> &, const std::string &);

// DEPRECATED: The following functions are deprecated in v2.0.
// They will throw std::runtime_error when called.
// Use Python-based workflows instead.

[[deprecated("Use Python solver and scipy.optimize instead")]]
void launch_solver(const generation &, const int &, std::vector<parameters> &, std::vector<constants> &, const std::string &, const std::string &, const std::string &, const std::string &, const std::string&);

[[deprecated("Use Python solver and scipy.optimize instead")]]
void launch_odf(const generation &, std::vector<parameters> &, const std::string &, const std::string &, const std::string &, const std::string &, const std::string &, const std::string&);

[[deprecated("Use Python solver and scipy.optimize instead")]]
void launch_func_N(const generation &, const int &, std::vector<parameters> &, std::vector<constants> &, const std::string &, const std::string &, const std::string &, const std::string &, const std::string &, const std::string&);

[[deprecated("Use Python solver and scipy.optimize instead")]]
void run_simulation(const std::string &, const individual &, const int &, std::vector<parameters> &, std::vector<constants> &, std::vector<opti_data> &, const std::string &, const std::string &, const std::string &, const std::string &, const std::string&);

double calc_cost(const arma::vec &, arma::vec &, const arma::vec &, const std::vector<opti_data> &, const std::vector<opti_data> &, const int &, const int &);

[[deprecated("Use Python solver and scipy.optimize instead")]]
arma::mat calc_sensi(const individual &, generation &, const std::string &, const int &, const int &, std::vector<parameters> &, std::vector<constants> &, arma::vec &, std::vector<opti_data> &, std::vector<opti_data> &, const std::string &, const std::string &, const std::string &, const std::string &, const int &, const arma::vec &, const std::string&);


/** @} */ // end of identification group

} //namespace simcoon
