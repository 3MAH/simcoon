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
///@brief Solver header - see Python API for solver functionality
///@version 2.0
///
///@note The legacy C++ solver function that read path.txt/material.dat files
///      has been removed in v2.0. Use the Python simcoon.solver.Solver class instead.
///      For direct UMAT calls, use the umat functions in umat_smart.hpp.

#pragma once
#include <armadillo>
#include <string>

namespace simcoon{

/**
 * @file solver.hpp
 * @brief Solver header - functionality moved to Python API.
 *
 * @note The legacy solver() function has been removed in v2.0.
 *       Use the Python simcoon.solver.Solver class for material point simulations.
 *
 * Example Python usage:
 * @code
 * from simcoon.solver import Solver, Block, StepMeca
 * import numpy as np
 *
 * props = np.array([210000.0, 0.3])  # E, nu for ELISO
 * step = StepMeca(
 *     DEtot_end=np.array([0.01, 0, 0, 0, 0, 0]),
 *     control=['strain', 'stress', 'stress', 'stress', 'stress', 'stress']
 * )
 * block = Block(steps=[step], umat_name="ELISO", props=props, nstatev=1)
 * solver = Solver(blocks=[block])
 * history = solver.solve()
 * @endcode
 */

/** @addtogroup solver
 *  @{
 */

// NOTE: The legacy solver() function signature has been removed.
// Use the Python API (simcoon.solver.Solver) for solver workflows.
// For C++ UMAT integration (Abaqus/Ansys), use the umat functions directly.

/** @} */ // end of solver group

} //namespace simcoon
