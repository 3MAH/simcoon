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
 * @file Zener_Nfast.hpp
 * @brief Thermomechanical generalized Kelvin chain (Zener) viscoelastic model with N units
 * @author Yves Chemisky, George Chatzigeorgiou
 * @version 2.0
 */

#pragma once

#include <iostream>
#include <armadillo>

namespace simcoon {

/** @addtogroup umat_thermomechanical
 *  @{
 */

/**
 * @brief Thermomechanical generalized Kelvin chain (Zener) viscoelastic model with N units (``ZENNK``, thermomechanical)
 *
 * @details Fully coupled thermomechanical version of umat_zener_Nfast(): the same rheology,
 * stress update and state variables (see that kernel for the constitutive equations and
 * the instantaneous-stiffness convention), plus the heat equation terms. The volumetric
 * heat capacity is \f$ c_0 = \rho\, c_p \f$ and the heat source \f$ r \f$ collects the
 * thermoelastic coupling \f$ -T\, \boldsymbol{\alpha} : \dot{\boldsymbol{\sigma}} \f$ and the
 * viscous dissipation of the dashpots; its sensitivities `drdE` and `drdT` and the
 * stress sensitivity `dSdT` are returned for the monolithic thermomechanical solver.
 * At constant temperature this kernel reproduces umat_zener_Nfast() exactly.
 *
 * **Material parameters (props), \f$ n_{props} = 6 + 4N \f$:**
 *
 * | Index | Symbol | Description | Units |
 * |-------|--------|-------------|-------|
 * | props[0] | \f$ \rho \f$ | Density | Mass/Volume |
 * | props[1] | \f$ c_p \f$ | Specific heat capacity | Energy/(Mass x Temperature) |
 * | props[2] | \f$ E_0 \f$ | Instantaneous Young's modulus | Stress |
 * | props[3] | \f$ \nu_0 \f$ | Instantaneous Poisson's ratio | - |
 * | props[4] | \f$ \alpha \f$ | Isotropic thermal expansion coefficient | 1/Temperature |
 * | props[5] | \f$ N \f$ | Number of Kelvin units | - |
 * | props[6+4(i-1)] | \f$ E_i \f$ | Young's modulus of branch i | Stress |
 * | props[7+4(i-1)] | \f$ \nu_i \f$ | Poisson's ratio of branch i | - |
 * | props[8+4(i-1)] | \f$ \eta_{B,i} \f$ | Bulk viscosity of branch i | Stress x Time |
 * | props[9+4(i-1)] | \f$ \eta_{S,i} \f$ | Shear viscosity of branch i | Stress x Time |
 *
 * **State variables (statev), \f$ n_{statev} = 7 + 7N \f$:**
 *
 * | Index | Symbol | Description | Units |
 * |-------|--------|-------------|-------|
 * | statev[0] | \f$ T_{init} \f$ | Initial temperature | Temperature |
 * | statev[1:6] | \f$ \boldsymbol{\varepsilon}^{v} \f$ | Total viscous strain \f$ \sum_i \boldsymbol{\varepsilon}^{v}_i \f$ (Voigt) | - |
 * | statev[7+7(i-1)] | \f$ v_i \f$ | Accumulated flow length of branch i | - |
 * | statev[8+7(i-1):13+7(i-1)] | \f$ \boldsymbol{\varepsilon}^{v}_i \f$ | Viscous strain of branch i (Voigt) | - |
 *
 * @param Etot Total strain at the beginning of the increment (Voigt notation: \f$6 \times 1\f$ vector)
 * @param DEtot Strain increment (Voigt notation: \f$6 \times 1\f$ vector)
 * @param sigma Stress (Voigt notation: \f$6 \times 1\f$ vector) [input/output]
 * @param r Heat source per unit volume (thermoelastic coupling + viscous dissipation) [output]
 * @param dSdE Consistent mechanical tangent \f$ \partial \boldsymbol{\sigma} / \partial \boldsymbol{\varepsilon} \f$ (\f$6 \times 6\f$) [output]
 * @param dSdT Stress sensitivity to temperature \f$ \partial \boldsymbol{\sigma} / \partial T \f$ (\f$6 \times 1\f$) [output]
 * @param drdE Heat source sensitivity to strain \f$ \partial r / \partial \boldsymbol{\varepsilon} \f$ (\f$1 \times 6\f$) [output]
 * @param drdT Heat source sensitivity to temperature \f$ \partial r / \partial T \f$ (\f$1 \times 1\f$) [output]
 * @param DR Rotation increment matrix (\f$3 \times 3\f$) for objective integration
 * @param nprops Number of material properties
 * @param props Material properties vector (see table above)
 * @param nstatev Number of state variables
 * @param statev State variables vector (see table above) [input/output]
 * @param T Temperature at the beginning of the increment
 * @param DT Temperature increment
 * @param Time Time at the beginning of the increment
 * @param DTime Time increment
 * @param Wm Total mechanical work [input/output, cumulative]
 * @param Wm_r Recoverable (stored) work [input/output, cumulative]
 * @param Wm_ir Irrecoverable stored work, always 0 for this model [input/output]
 * @param Wm_d Dissipated (viscous) work [input/output, cumulative]
 * @param Wt Total thermal work \f$ \int T\, d\eta \f$ [input/output, cumulative]
 * @param Wt_r Reversible thermal work [input/output, cumulative]
 * @param Wt_ir Irreversible thermal work [output]
 * @param ndi Number of direct stress components (typically 3)
 * @param nshr Number of shear stress components (typically 3)
 * @param start Flag indicating the first increment (true) or a continuation (false)
 * @param tnew_dt Suggested new time step size for adaptive time stepping [output]
 * @param tangent_mode Tangent operator selection (see parameter.hpp)
 *
 * @note The mechanical props of umat_zener_Nfast() are shifted by two slots (\f$ \rho, c_p \f$ first)
 *
 * @see umat_zener_Nfast() the mechanical twin
 */
void umat_zener_Nfast_T(const arma::vec &Etot, const arma::vec &DEtot, arma::vec &sigma, double &r, arma::mat &dSdE, arma::mat &dSdT, arma::mat &drdE, arma::mat &drdT, const arma::mat &DR, const int &nprops, const arma::vec &props, const int &nstatev, arma::vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, double &Wt, double &Wt_r, double &Wt_ir, const int &ndi, const int &nshr, const bool &start, double &tnew_dt, const int &tangent_mode = 0);

/** @} */ // end of umat_thermomechanical group

} //namespace simcoon
