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

///@file plastic_isotropic.hpp
///@brief User subroutine for elastic-plastic materials in 1D-2D-3D case
///@brief Convex cutting plane return mapping (tangent_mode 0-2) or closest-point projection (tangent_mode 3)
///@brief Isotropic hardening with a power-law hardening is considered (PowerLawHardening: exact for
///       p >= 1e-6, C1 onset blend below, see the mechanical header)
///@version 1.0

#pragma once
#include <armadillo>

namespace simcoon{

/**
 * @file plastic_isotropic.hpp
 * @brief Thermomechanical plasticity model.
 *
 * @note **Integrator.** Under tangent_mode 0-2 the return mapping is the convex cutting plane
 * loop; under tangent_closest_point (3, full 3D states) the kernel runs the modular
 * PlasticityMechanism (von Mises + PowerLawHardening) on its own props and state through
 * closest_point_return_mapping() and returns cpp_consistent_tangent() for dSdE; the thermal
 * cross-tangents keep the continuum form in every mode. Non-convergence requests tnew_dt = 0.5.
 */

/** @addtogroup umat_thermomechanical
 *  @{
 */


/**
 * @brief Thermomechanical elastic-plastic model with isotropic hardening (J2, power-law hardening)
 *
 * @details Yield function \f$ \Phi = \sigma_{eq} - k\,p^m - \sigma_Y \leq 0 \f$ with associated
 * J2 flow. The Gibbs free energy of the lead mechanism (internal variables \f$ p \f$ and
 * \f$ \boldsymbol{\varepsilon}^p \f$) reads
 * \f[
 * G = -\frac{1}{2}\,\boldsymbol{\sigma} : \mathbf{S} : \boldsymbol{\sigma}
 * - \boldsymbol{\sigma} : \boldsymbol{\varepsilon}^p
 * - \boldsymbol{\sigma} : \boldsymbol{\alpha}\,(\theta - \theta_0)
 * + c_0\left[(\theta - \theta_0) - \theta \ln\frac{\theta}{\theta_0}\right]
 * - \eta_0\,\theta + E_0 + F(\theta, p), \qquad c_0 = \rho\, c_p .
 * \f]
 * Under tangent_none the kernel still assembles the continuum operator: the sensitivities
 * \f$ P_{\varepsilon} \f$ and \f$ P_\theta \f$ feed the physical heat source \f$ r \f$ and its
 * linearization, not only the Newton operator. Under tangent_algorithmic only dSdE carries the
 * consistent correction; dSdT, drdE and drdT keep the continuum form. drdE and drdT
 * differentiate the dissipation through \f$ \mathbf{P}_\varepsilon \f$, \f$ P_\theta \f$ and the
 * turning of the flow direction over the increment
 * (\f$ \tfrac12\Delta p\,\partial\boldsymbol{\Lambda}/\partial\boldsymbol{\sigma}:\boldsymbol{\sigma}^n \f$):
 * exact to round-off against finite differences.
 *
 * **Material Parameters (props, 8):**
 * | Index | Symbol | Description | Units |
 * |-------|--------|-------------|-------|
 * | props[0] | \f$ \rho \f$ | Density | Mass/Volume |
 * | props[1] | \f$ c_p \f$ | Specific heat capacity | Energy/(Mass.Temperature) |
 * | props[2] | \f$ E \f$ | Young's modulus | Stress (MPa) |
 * | props[3] | \f$ \nu \f$ | Poisson's ratio | - |
 * | props[4] | \f$ \alpha \f$ | Isotropic thermal expansion coefficient | 1/Temperature |
 * | props[5] | \f$ \sigma_Y \f$ | J2 equivalent yield stress | Stress (MPa) |
 * | props[6] | \f$ k \f$ | Hardening parameter | Stress (MPa) |
 * | props[7] | \f$ m \f$ | Hardening exponent | - |
 *
 * **State Variables (statev, 8):**
 * | Index | Symbol | Description |
 * |-------|--------|-------------|
 * | statev[0] | \f$ T_{init} \f$ | Initial temperature |
 * | statev[1] | \f$ p \f$ | Accumulated plastic strain |
 * | statev[2..7] | \f$ \boldsymbol{\varepsilon}^p \f$ | Plastic strain, Voigt [11,22,33,12,13,23], engineering shear (\f$\times 2\f$) |
 */
void umat_plasticity_iso_T(const arma::vec &, const arma::vec &, arma::vec &, double &, arma::mat &, arma::vec &, arma::vec &, double &, const arma::mat &, const int &, const arma::vec &, const int &, arma::vec &, const double &, const double &,const double &,const double &, double &, double &, double &, double &, double &, double &, double &, const int &, const int &, const bool &, double &, const int & = 0);
    

/** @} */ // end of umat_thermomechanical group

} //namespace simcoon
