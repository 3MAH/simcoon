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


#pragma once
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <string>

namespace simcoon{

/**
 * @file solver.hpp
 * @brief Solver functions and classes.
 */

/** @addtogroup solver
 *  @{
 */


/**
 * The file-driven entry point solver() left the library with the 2.0 JSON-only migration: the
 * loading programme is built in Python and handed to solver_run() (solver_sink.hpp) in memory,
 * and nothing here reads a file: JSON inputs are read in Python (simcoon.solver), the
 * reference cases of testBin run through pytest, and the legacy text formats are converted
 * once by scripts/legacy_to_json.py.
 */

/** @} */ // end of solver group

} //namespace simcoon
