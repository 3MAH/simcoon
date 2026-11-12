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

///@file plastic_kin_iso.hpp
///@brief User subroutine for elastic-plastic materials in 1D-2D-3D case
///@brief Convex cutting plane return mapping, or closest-point projection (integrator='closest_point')
///@brief Isotropic hardening with a power-law hardenig is considered
///@version 1.0

#pragma once
#include <armadillo>

namespace simcoon{

/**
 * @file plastic_kin_iso.hpp
 * @brief Thermomechanical plasticity model.
 *
 * @note **Integrator.** By default the return mapping is the convex cutting plane loop; under
 * integrator='closest_point' (kernel code tangent_closest_point, full 3D states) the kernel runs the modular
 * PlasticityMechanism (von Mises + PowerLawHardening + PragerHardening with \f$ C = 3 k_X / 2 \f$
 * for the kernel's \f$ \mathbf{X} = k_X \mathbf{T} \mathbf{a} \f$) on its own props and state through
 * closest_point_return_mapping() and returns cpp_consistent_tangent() for dSdE; the thermal
 * cross-tangents keep the continuum form in every mode. Non-convergence requests tnew_dt = 0.5.
 */

/** @addtogroup umat_thermomechanical
 *  @{
 */


/**
 * @brief Thermomechanical elastic-plastic UMAT with isotropic (power-law) and linear kinematic
 *        (Prager) hardening, von Mises criterion (EPKCP).
 *
 * The backstress is not stored: it follows from the back-strain \f$ \mathbf{a} \f$ (engineering
 * shear components, rotated as a strain) as
 * \f$ \mathbf{X} = k_X\,(\mathbf{a} \circ \mathbf{I}_{r05}) \f$, i.e. Prager hardening with
 * \f$ C = 3 k_X / 2 \f$. \f$ W_m^{ir} \f$ books the hardening energy
 * \f$ \tfrac12(R_n+R_{n+1})\Delta p + \tfrac12(\mathbf{X}_n+\mathbf{X}_{n+1}):\Delta\mathbf{a} \f$,
 * \f$ W_m^d \f$ the rest of the plastic work; drdE / drdT carry the turning of the flow direction
 * on \f$ \boldsymbol{\sigma}-\mathbf{X} \f$ and hold the back-stress direction fixed (1e-6 on a
 * non-proportional probe).
 *
 * **Material properties** (9; MPa, K):
 * | Index | Symbol | Description |
 * |-------|--------|-------------|
 * | props[0] | \f$ \rho \f$ | density |
 * | props[1] | \f$ c_p \f$ | specific heat capacity |
 * | props[2] | \f$ E \f$ | Young modulus |
 * | props[3] | \f$ \nu \f$ | Poisson ratio |
 * | props[4] | \f$ \alpha \f$ | coefficient of thermal expansion |
 * | props[5] | \f$ \sigma_Y \f$ | J2 equivalent yield stress |
 * | props[6] | \f$ k \f$ | isotropic hardening modulus, \f$ H_p = k\,p^m \f$ |
 * | props[7] | \f$ m \f$ | isotropic hardening exponent |
 * | props[8] | \f$ k_X \f$ | linear kinematic hardening modulus |
 *
 * **State variables** (14; Voigt order [11, 22, 33, 12, 13, 23], engineering shear strains):
 * | Index | Symbol | Description |
 * |-------|--------|-------------|
 * | statev[0] | \f$ T_{init} \f$ | initial temperature |
 * | statev[1] | \f$ p \f$ | accumulated plastic strain |
 * | statev[2..7] | \f$ \boldsymbol{\varepsilon}^p \f$ | plastic strain |
 * | statev[8..13] | \f$ \mathbf{a} \f$ | kinematic back-strain (the backstress is \f$ \mathbf{X} = k_X\,(\mathbf{a} \circ \mathbf{I}_{r05}) \f$) |
 */
void umat_plasticity_kin_iso_T(const arma::vec &, const arma::vec &, arma::vec &, double &, arma::mat &, arma::vec &, arma::vec &, double &, const arma::mat &, const int &, const arma::vec &, const int &, arma::vec &, const double &, const double &,const double &,const double &, double &, double &, double &, double &, double &, double &, double &, const int &, const int &, const bool &, double &, const int & = 0);
    

/** @} */ // end of umat_thermomechanical group

} //namespace simcoon
