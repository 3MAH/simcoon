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

///@file generic_hyper_invariants.hpp
///@brief User subroutine for generic hyperelastic materials using invariants
///@version 1.0

#pragma once

#include <armadillo>
#include <simcoon/parameter.hpp>

namespace simcoon{

/** @addtogroup umat_finite
 *  @{
 */

/**
 * @brief Generic hyperelastic UMAT for potentials expressed in the isochoric invariants.
 *
 * Serves every HyperPotential of the form
 * \f[
    W\left(\bar{I}_1, \bar{I}_2\right) + U(J)
 * \f]
 * and, for an anisotropic potential, the fibre pseudo-invariants
 * \f$ \bar{I}^{*}_{4,i} \f$ as well. @p umat_name selects the potential: NEOHC, MOORI,
 * YEOHH, ISHAH, GETHH, SWANH and HOLZA. The kernel evaluates the potential's
 * derivatives (hyper_potential_derivatives) and hands them to
 * hyper_invariants_response, which returns the Cauchy stress and the canonical box
 * tangent \f$ \partial \hat{\boldsymbol{\tau}} / \partial \mathbf{D}_e \f$.
 *
 * **props** are the selected potential's own parameters, in the order HyperPotential
 * documents for it, optionally followed by one more entry selecting the volumetric term
 * \f$ U(J) \f$ (see VolumetricPotential: absent or 0 for
 * \f$ \kappa (J \ln J - J + 1) \f$, 1 for \f$ \frac{\kappa}{2} (J-1)^2 \f$).
 *
 * **statev**: one is required, @c statev(0), which stores the initial temperature.
 * The law is hyperelastic, so \f$ W_m = W_{m,r} \f$ and
 * \f$ W_{m,ir} = W_{m,d} = 0 \f$ throughout.
 *
 * @param[in] umat_name the 5-letter name selecting the potential
 * @param[in] etot,Detot total strain and its increment
 * @param[in] F0,F1 deformation gradient at the start and the end of the increment
 * @param[in,out] sigma Cauchy stress, 6-Voigt
 * @param[out] Lt canonical box tangent, 6x6
 * @param[out] L the ground-state stiffness, set on the first call
 * @param[in] DR the increment of rigid-body rotation
 * @param[in] nprops,props the potential's parameters
 * @param[in] nstatev,statev the state variables (one, the initial temperature)
 * @param[in] T,DT temperature and its increment
 * @param[in] Time,DTime time and its increment
 * @param[in,out] Wm_0,Wm_1,Wm_2,Wm_3 the cumulative work terms \f$ W_m, W_{m,r}, W_{m,ir}, W_{m,d} \f$
 * @param[in] ndi,nshr the number of direct and shear components
 * @param[in] start true on the first call of a block
 * @param[out] tnew_dt the suggested time-step scaling
 * @param[in] tangent_mode unused: a hyperelastic law returns its exact tangent
 */
void umat_generic_hyper_invariants(const std::string &umat_name, const arma::vec &etot, const arma::vec &Detot, const arma::mat &F0, const arma::mat &F1, arma::vec &sigma, arma::mat &Lt, arma::mat &L, const arma::mat &DR, const int &nprops, const arma::vec &props, const int &nstatev, arma::vec &statev, const double &T, const double &DT,const double &Time,const double &DTime, double &Wm_0, double &Wm_1, double &Wm_2, double &Wm_3, const int &ndi, const int &nshr, const bool &start, double &tnew_dt, const int &corate_type, const int &tangent_mode = tangent_default);
                        

/** @} */ // end of umat_finite group

} //namespace simcoon
