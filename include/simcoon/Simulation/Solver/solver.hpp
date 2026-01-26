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

///@file solver.hpp
///@brief To solver an homogeneous thermomechanical problem
///@version 1.0

#pragma once
#include <armadillo>
#include <string>

namespace simcoon{

/**
 * @file solver.hpp
 * @brief Solver functions and classes.
 */

/** @addtogroup solver
 *  @{
 */


//function that solves a
/**
 * @brief Main solver function for homogeneous thermomechanical problems.
 * 
 * @param umat_name Name of the constitutive model (UMAT)
 * @param props Vector of material properties
 * @param nstatev Number of internal state variables
 * @param psi_rve First Euler angle of RVE orientation (rad)
 * @param theta_rve Second Euler angle of RVE orientation (rad)
 * @param phi_rve Third Euler angle of RVE orientation (rad)
 * @param solver_type Type of solver (0: small strain, 1: finite strain)
 * @param corate_type Type of corotational formulation
 * @param div Divisor for time stepping (default: 0.5)
 * @param mul Multiplier for time stepping (default: 2.0)
 * @param miniter Minimum iterations per increment (default: 10)
 * @param maxiter Maximum iterations per increment (default: 100)
 * @param inforce_solver Enforce solver convergence (default: 1)
 * @param precision Convergence tolerance (default: 1e-6)
 * @param lambda_eff Effective stiffness estimate for mixed control (default: 10000)
 * @param path_data Path to data directory (default: "data")
 * @param path_results Path to results directory (default: "results")
 * @param pathfile Name of loading path file (default: "path.txt")
 * @param outputfile Name of output file (default: "result_job.txt")
 * 
 * @details This function drives the simulation by:
 * - Reading the loading path from input files
 * - Managing time stepping with adaptive incrementation
 * - Calling the UMAT for constitutive updates
 * - Writing results to output files
 * 
 */
void solver(const std::string &umat_name, const arma::vec &props, const unsigned int &nstatev, const double &psi_rve, const double &theta_rve, const double &phi_rve, const int &solver_type, const int &corate_type, const double &div = 0.5, const double &mul = 2., const int &miniter = 10, const int &maxiter = 100, const int &inforce_solver = 1, const double &precision = 1.E-6, const double &lambda_eff = 10000., const std::string &path_data = "data", const std::string &path_results = "results", const std::string &pathfile = "path.txt", const std::string &outputfile = "result_job.txt");


/** @} */ // end of solver group

} //namespace simcoon
