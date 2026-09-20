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
void umat_multi(phase_characteristics &rve, const arma::mat &DR, const double &Time, const double &DTime, const int &ndi, const int &nshr, bool &start, const unsigned int &solver_type, double &tnew_dt, const int &control);

/**
 * @brief Geometry of the sub-phases a mean-field model is built from.
 *
 * @param umat_name name of the model
 * @return the code phase_characteristics::construct takes: 2 (ellipsoids) for MIHEN, MIMTN
 * and MISCN, 1 (layers) for MIPLN, 0 for a homogeneous model, which has no sub-phases
 */
int sub_phase_shape(const std::string &umat_name);

/**
 * @brief Check the sub-phases and the props of a mean-field phase.
 *
 * The props hold the scheme's settings only: MIHEN [mp, np], MIMTN [mp, np, n_matrix],
 * MISCN [mp, np, n_matrix, (start)], MIPLN []. Checked: at least one sub-phase, the
 * exact props length, mp and np >= 1, n_matrix among the phases given.
 * @param phase the mean-field phase
 * @throws std::invalid_argument on any of the above
 */
void check_sub_phases(const phase_characteristics &phase);

/**
 * @brief First guess of the self-consistent scheme: props(3) when given, else 1.
 * @return 1 Mori-Tanaka (default), 0 homogeneous strain (n_matrix < 0)
 */
int self_consistent_start(const phase_characteristics &phase);


/** @} */ // end of micromechanics group

} //namespace simcoon
