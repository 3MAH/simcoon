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

///@file external_umat.hpp
///@brief External user material subroutine interface for custom mechanical models
///@version 1.0

#pragma once

#include <armadillo>

/** @addtogroup umat_mechanical
 *  @{
 */

/**
 * @brief External mechanical UMAT subroutine wrapper
 *
 * @details Wrapper function for external user-defined material subroutines (UMATs)
 * loaded as plugins at runtime. The actual constitutive model is implemented by the
 * user through the plugin API (see umat_plugin_api.hpp).
 *
 * Material parameters and state variables are user-defined and depend on
 * the specific constitutive model implemented in the plugin.
 *
 * @param Etot Total strain tensor at beginning of increment (Voigt notation: 6x1 vector)
 * @param DEtot Strain increment tensor (Voigt notation: 6x1 vector)
 * @param sigma Stress tensor (Voigt notation: 6x1 vector) [output]
 * @param Lt Consistent tangent modulus (6x6 matrix) [output]
 * @param L Elastic stiffness tensor (6x6 matrix) [output]
 * @param sigma_in Internal stress contribution for explicit solvers (6x1 vector) [output]
 * @param DR Rotation increment matrix (3x3) for objective integration
 * @param nprops Number of material properties
 * @param props Material properties vector
 * @param nstatev Number of state variables
 * @param statev State variables vector [input/output]
 * @param T Temperature at beginning of increment
 * @param DT Temperature increment
 * @param Time Time at beginning of increment
 * @param DTime Time increment
 * @param Wm Total mechanical work [output]
 * @param Wm_r Recoverable (elastic) work [output]
 * @param Wm_ir Irrecoverable work [output]
 * @param Wm_d Dissipated work [output]
 * @param ndi Number of direct stress components (typically 3)
 * @param nshr Number of shear stress components (typically 3)
 * @param start Flag indicating first increment (true) or continuation (false)
 * @param solver_type Solver type: 0=implicit (Newton), 1=explicit (RNL), 2=dynamic implicit
 * @param tnew_dt Suggested new time step ratio for adaptive time stepping [output]
 *
 * @see umat_plugin_ext_api for the native simcoon plugin base class
 * @see umat_plugin_aba_api for Abaqus-compatible plugins
 * @see umat_plugin_ans_api for ANSYS-compatible plugins
 */
extern void umat_external(const arma::vec &Etot, const arma::vec &DEtot, arma::vec &sigma, arma::mat &Lt, arma::mat &L, arma::vec &sigma_in, const arma::mat &DR, const int &nprops, const arma::vec &props, const int &nstatev, arma::vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, const int &solver_type, double &tnew_dt);

/** @} */ // end of umat_mechanical group
