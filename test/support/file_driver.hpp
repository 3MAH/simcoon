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

///@file file_driver.hpp
///@brief The historical file-driven solver entry point, kept for the C++ test suite
///@version 1.0

#pragma once
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <string>

namespace simcoon{

/**
 * @file file_driver.hpp
 * @brief solver(): read a case from its files, run it, write the results files.
 *
 * This used to be the shipped entry point of the library. With the 2.0 JSON-only migration the
 * library stopped reading files altogether — a loading programme is built in Python and handed
 * to solver_run() in memory — so the wrapper moved here, where the reference cases of
 * test/Umats still drive the engine from their historical fixtures.
 *
 * It keeps its name and its namespace on purpose: the test files call solver(...) unchanged.
 */

/** @addtogroup solver
 *  @{
 */

void solver(const std::string &umat_name, const arma::vec &props, const unsigned int &nstatev, const double &psi_rve, const double &theta_rve, const double &phi_rve, const int &solver_type, const int &corate_type, const double &div = 0.5, const double &mul = 2., const int &miniter = 10, const int &maxiter = 100, const int &inforce_solver = 1, const double &precision = 1.E-6, const double &lambda_eff = 10000., const std::string &path_data = "data", const std::string &path_results = "results", const std::string &pathfile = "path.txt", const std::string &outputfile = "result_job.txt", const int &tangent_mode = tangent_default);

/** @} */ // end of solver group

} //namespace simcoon
