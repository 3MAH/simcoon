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

///@file Zener_fast.hpp
///@brief User subroutine for Zener viscoelastic model in 3D case, with thermoelastic effect(Poynting-Thomson model)
///@author Chemisky, Chatzigeorgiou
///@version 1.0

#pragma once

#include <iostream>
#include <armadillo>

namespace simcoon {

/** @addtogroup umat_thermomechanical
 *  @{
 */

/**
 * @brief Thermomechanical Zener (Poynting--Thomson) viscoelastic model with one Kelvin unit (``ZENER``, thermomechanical)
 *
 * @details Single-unit case of umat_zener_Nfast_T(): an instantaneous isotropic spring
 * \f$ \mathbf{L}_0 = \mathbf{L}_{iso}(E_0, \nu_0) \f$ in series with one Kelvin--Voigt unit,
 * spring \f$ \mathbf{L}_1 = \mathbf{L}_{iso}(E_1, \nu_1) \f$ in parallel with the dashpot
 * \f$ \mathbf{H}_1 = 3\,\eta_{B,1}\,\mathbf{I}_{vol} + 2\,\eta_{S,1}\,\mathbf{I}_{dev} \f$:
 * \f[
 * \boldsymbol{\sigma} = \mathbf{L}_0 : \left( \boldsymbol{\varepsilon} - \boldsymbol{\alpha}\,(T - T_{init}) - \boldsymbol{\varepsilon}^{v}_1 \right),
 * \qquad
 * \dot{\boldsymbol{\varepsilon}}^{v}_1 = \mathbf{H}_1^{-1} : \left( \boldsymbol{\sigma} - \mathbf{L}_1 : \boldsymbol{\varepsilon}^{v}_1 \right),
 * \f]
 * with \f$ \boldsymbol{\alpha} = \alpha\,\mathbf{I} \f$. Backward Euler in closed form
 * (kelvin_series_step, linear_viscoelastic.hpp), exact consistent tangents `dSdE` and `dSdT`
 * (`tangent_none` returns \f$ \mathbf{L}_0 \f$ and \f$ -\mathbf{L}_0\,\boldsymbol{\alpha} \f$).
 * The volumetric heat capacity is \f$ c_0 = \rho\, c_p \f$, the entropy
 * \f$ \eta = c_0 \ln(T/T_{init}) + \boldsymbol{\alpha}:\boldsymbol{\sigma} \f$, and the dissipation
 * is the trapezoidal work of the viscous stress \f$ \boldsymbol{\sigma} - \mathbf{L}_1 : \boldsymbol{\varepsilon}^{v}_1 \f$
 * on \f$ \Delta\boldsymbol{\varepsilon}^{v}_1 \f$; the heat of the increment \f$ r \f$ and its
 * derivatives come from heat_source_increment(). The arguments are those of umat_zener_Nfast_T().
 *
 * **Material parameters (props), \f$ n_{props} = 9 \f$:**
 *
 * | Index | Symbol | Description | Units |
 * |-------|--------|-------------|-------|
 * | props[0] | \f$ \rho \f$ | Density | Mass/Volume |
 * | props[1] | \f$ c_p \f$ | Specific heat capacity | Energy/(Mass x Temperature) |
 * | props[2] | \f$ E_0 \f$ | Instantaneous Young's modulus | Stress |
 * | props[3] | \f$ \nu_0 \f$ | Instantaneous Poisson's ratio | - |
 * | props[4] | \f$ \alpha \f$ | Isotropic thermal expansion coefficient | 1/Temperature |
 * | props[5] | \f$ E_1 \f$ | Young's modulus of the Kelvin unit | Stress |
 * | props[6] | \f$ \nu_1 \f$ | Poisson's ratio of the Kelvin unit | - |
 * | props[7] | \f$ \eta_{B,1} \f$ | Bulk viscosity of the Kelvin unit | Stress x Time |
 * | props[8] | \f$ \eta_{S,1} \f$ | Shear viscosity of the Kelvin unit | Stress x Time |
 *
 * **State variables (statev), \f$ n_{statev} = 8 \f$:**
 *
 * | Index | Symbol | Description | Units |
 * |-------|--------|-------------|-------|
 * | statev[0] | \f$ T_{init} \f$ | Initial temperature | Temperature |
 * | statev[1] | \f$ v \f$ | Accumulated flow length \f$ \sum \| \Delta\boldsymbol{\varepsilon}^{v}_1 \| \f$ | - |
 * | statev[2:7] | \f$ \boldsymbol{\varepsilon}^{v}_1 \f$ | Viscous strain of the Kelvin unit (Voigt [11,22,33,12,13,23], engineering shears) | - |
 *
 * @see umat_zener_Nfast_T() the N-unit version (props shifted by the N slot, statev 7 + 7N)
 */
void umat_zener_fast_T(const arma::vec &, const arma::vec &, arma::vec &, double &, arma::mat &, arma::mat &, arma::mat &, arma::mat &, const arma::mat &, const int &, const arma::vec &, const int &, arma::vec &, const double &, const double &,const double &,const double &, double &, double &, double &, double &, double &, double &, double &, const int &, const int &, const bool &, double &, const int & = 0);

/** @} */ // end of umat_thermomechanical group

} //namespace smart
