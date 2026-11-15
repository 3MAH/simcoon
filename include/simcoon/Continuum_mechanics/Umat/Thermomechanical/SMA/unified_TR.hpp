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
/**
 * @file unified_TR.hpp
 * @brief Thermomechanical face of the unified SMA model with martensite reorientation (SMR*)
 */

#pragma once
#include <string>
#include <armadillo>
#include <simcoon/parameter.hpp>

namespace simcoon{

/** @addtogroup umat_thermomechanical
 *  @{
 */

/**
 * @brief Unified SMA model with transformation and reorientation, coupled to the heat equation
 *        (SMRDI, SMRDC, SMRAI, SMRAC).
 *
 * The mechanical model is umat_sma_unified_TR (sma_unified_TR_increment), on either integrator; the
 * thermal face is the one of umat_sma_unified_T_T (sma_thermal_face) with the reorientation row as a
 * third mechanism. The reorientation criterion does not depend on the temperature
 * (\f$ \partial\Phi^{Reo}/\partial\theta = 0 \f$), its work
 * \f$ \tfrac12(\boldsymbol{\sigma}_n + \boldsymbol{\sigma}):\Delta\boldsymbol{\varepsilon}^{Reo} \f$ is
 * dissipated (the back-strain stores no energy, as on the mechanical face) and it does not enter
 * the entropy.
 *
 * **Props:** \f$ \rho \f$, \f$ c_p^A \f$, \f$ c_p^M \f$, then the props of umat_sma_unified_TR (the
 * SMR* layout, from flagT on): 38 (SMRDI), 40 (SMRDC), 45 (SMRAI), 47 (SMRAC).
 *
 * **State variables:** the 30 of umat_sma_unified_TR.
 *
 * The tangent modes are those of umat_sma_unified_T_T.
 *
 * @param umat_name SMRDI, SMRDC, SMRAI or SMRAC
 * @param Etot, DEtot total strain at the start of the increment and its increment (6)
 * @param sigma [in,out] stress (6)
 * @param r [out] heat of the increment per reference volume (heat_source_increment())
 * @param dSdE, dSdT, drdE, drdT [out] tangents (sma_thermal_face())
 * @param DR rotation increment (3x3)
 * @param nprops, nstatev numbers of props and state variables (unused: read from the name)
 * @param props material properties (above)
 * @param statev [in,out] state variables (30)
 * @param T, DT temperature at the start of the increment and its increment
 * @param Time, DTime time at the start of the increment and its increment (unused: rate independent)
 * @param Wm, Wm_r, Wm_ir, Wm_d, Wt, Wt_r, Wt_ir [in,out] work terms
 * @param ndi, nshr direct and shear components (3D only)
 * @param start first increment (initialisation)
 * @param tnew_dt [out] step-cut request
 * @param tangent_mode kernel tangent code (parameter.hpp)
 * @param ncp NCPMethod code of the local step (0 Fischer-Burmeister, the default; see sma_ncp_method):
 *        read by the closest point, the cutting plane takes Fischer-Burmeister
 */
void umat_sma_unified_TR_T(const std::string &umat_name, const arma::vec &Etot, const arma::vec &DEtot, arma::vec &sigma, double &r, arma::mat &dSdE, arma::vec &dSdT, arma::vec &drdE, double &drdT, const arma::mat &DR, const int &nprops, const arma::vec &props, const int &nstatev, arma::vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, double &Wt, double &Wt_r, double &Wt_ir, const int &ndi, const int &nshr, const bool &start, double &tnew_dt, const int &tangent_mode = tangent_default, const int &ncp = 0);

/** @} */ // end of umat_thermomechanical group

} //namespace simcoon
