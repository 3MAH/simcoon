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

///@file multiphase.hpp
///@brief User subroutine for non-linear N-phases heterogeneous materials using the method
///@version 1.0

#pragma once

#include <armadillo>
#include <simcoon/Simulation/Phase/phase_characteristics.hpp>

namespace simcoon{

/**
 * @file multiphase.hpp
 * @brief Multiphase material modeling functions.
 */

/** @addtogroup micromechanics
 *  @{
 */


// The multiphase function works with the following material properties
///@brief props[0] : Number of phases
///@brief props[1] : Number of the file NPhase[i].dat utilized
///@brief props[2] : Number of integration points in the 1 direction
///@brief props[3] : Number of integration points in the 2 direction

/**
 * @brief User material subroutine for non-linear N-phases heterogeneous materials.
 * 
 * @param rve Phase characteristics of the RVE containing sub-phases
 * @param DR Rotation increment tensor
 * @param Time Current simulation time
 * @param DTime Time increment
 * @param ndi Number of direct stress components
 * @param nshr Number of shear stress components
 * @param start Flag indicating if this is the first call (updated on output)
 * @param solver_type Type of solver algorithm
 * @param tnew_dt Suggested new time increment ratio
 * @param control Control flag for time stepping
 * 
 * @details This function performs the constitutive update for a multi-phase material
 * using a specified homogenization scheme. It handles:
 * - Localization of macroscopic strain to each phase
 * - Call to constituent UMATs for each phase
 * - Homogenization of phase responses to macroscopic level
 * 
 * @code
 *     phase_characteristics rve;
 *     // ... initialize rve with sub-phases ...
 *     mat DR = eye(3,3);
 *     double Time = 0.0, DTime = 0.01;
 *     bool start = true;
 *     double tnew_dt = 1.0;
 *     umat_multi(rve, DR, Time, DTime, 3, 3, start, 0, tnew_dt, 1);
 * @endcode
 */
void umat_multi(phase_characteristics &, const arma::mat &, const double &,const double &, const int &, const int &, bool &, const unsigned int &, double &, const int &);


/** @} */ // end of micromechanics group

} //namespace simcoon
