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
 * @file Prony_Nfast.hpp
 * @brief Linear viscoelastic material model using Prony series representation
 * @author Yves Chemisky
 * @version 1.0
 */

#pragma once

#include <string>
#include <iostream>
#include <armadillo>

namespace simcoon {

/** @addtogroup umat_mechanical
 *  @{
 */

/**
 * @brief Linear viscoelastic constitutive model with Prony series representation (Generalized Maxwell model)
 *
 * @details This function implements a linear viscoelastic material model using the Prony series
 * decomposition, also known as the Generalized Maxwell model. The model consists of multiple
 * Maxwell elements in parallel, each characterized by a relaxation time and modulus.
 *
 * **Constitutive Equations:**
 *
 * The relaxation modulus is expressed as a Prony series:
 * \f[
 * E(t) = E_\infty + \sum_{i=1}^N E_i e^{-t/\tau_i}
 * \f]
 * where:
 * - \f$ E_\infty \f$ is the long-term (equilibrium) modulus
 * - \f$ E_i \f$ is the modulus of the i-th Maxwell element
 * - \f$ \tau_i \f$ is the relaxation time of the i-th element
 * - \f$ N \f$ is the number of Maxwell elements
 *
 * **Incremental Form:**
 *
 * The stress at time \f$ t + \Delta t \f$ is computed as:
 * \f[
 * \boldsymbol{\sigma}(t + \Delta t) = \mathbf{L}_\infty : \boldsymbol{\varepsilon}(t + \Delta t) + \sum_{i=1}^N \mathbf{q}_i(t + \Delta t)
 * \f]
 *
 * where \f$ \mathbf{q}_i \f$ are the internal state variables (stress-like quantities) that evolve according to:
 * \f[
 * \mathbf{q}_i(t + \Delta t) = e^{-\Delta t/\tau_i} \mathbf{q}_i(t) + \frac{E_i}{E_\infty} \mathbf{L}_\infty : \left( e^{-\Delta t/\tau_i} - 1 \right) \Delta \boldsymbol{\varepsilon}
 * \f]
 *
 * **Bulk and Shear Decomposition:**
 *
 * The model handles volumetric and deviatoric responses independently:
 * - Volumetric: \f$ K(t) = K_\infty + \sum_{i=1}^N K_i e^{-t/\tau_i^K} \f$
 * - Deviatoric: \f$ G(t) = G_\infty + \sum_{i=1}^N G_i e^{-t/\tau_i^G} \f$
 *
 * **Consistent Tangent Modulus:**
 *
 * For implicit finite element analysis:
 * \f[
 * \mathbf{L}_t = \mathbf{L}_\infty + \sum_{i=1}^N \left( 1 - e^{-\Delta t/\tau_i} \right) \mathbf{L}_i
 * \f]
 *
 * **Material Parameters (props):**
 *
 * For N Maxwell elements, the props vector contains:
 * | Index | Symbol | Description | Units |
 * |-------|--------|-------------|-------|
 * | props[0] | \f$ N \f$ | Number of Maxwell elements | - |
 * | props[1] | \f$ E_\infty \f$ | Equilibrium Young's modulus | Stress |
 * | props[2] | \f$ \nu_\infty \f$ | Equilibrium Poisson's ratio | - |
 * | props[3] | \f$ \alpha \f$ | Thermal expansion coefficient | 1/Temperature |
 * | props[4] | \f$ E_1 \f$ | Modulus of 1st Maxwell element | Stress |
 * | props[5] | \f$ \tau_1 \f$ | Relaxation time of 1st element | Time |
 * | ... | ... | ... | ... |
 * | props[4+2(i-1)] | \f$ E_i \f$ | Modulus of i-th element | Stress |
 * | props[5+2(i-1)] | \f$ \tau_i \f$ | Relaxation time of i-th element | Time |
 *
 * **State Variables (statev):**
 *
 * For each Maxwell element, 6 internal state variables store the stress-like quantities \f$ \mathbf{q}_i \f$:
 * | Index | Symbol | Description | Units |
 * |-------|--------|-------------|-------|
 * | statev[0] | \f$ T_{init} \f$ | Initial temperature | Temperature |
 * | statev[1:6] | \f$ \mathbf{q}_1 \f$ | Internal stresses of 1st element (Voigt) | Stress |
 * | statev[7:12] | \f$ \mathbf{q}_2 \f$ | Internal stresses of 2nd element (Voigt) | Stress |
 * | ... | ... | ... | ... |
 * | statev[1+6(i-1):6i] | \f$ \mathbf{q}_i \f$ | Internal stresses of i-th element (Voigt) | Stress |
 *
 * Total state variables required: \f$ n_{statev} = 1 + 6N \f$
 *
 * @param Etot Total strain tensor at beginning of increment (Voigt notation: 6×1 vector)
 * @param DEtot Strain increment tensor (Voigt notation: 6×1 vector)
 * @param sigma Stress tensor (Voigt notation: 6×1 vector) [output]
 * @param Lt Consistent tangent modulus (6×6 matrix) [output]
 * @param DR Rotation increment matrix (3×3) for objective integration
 * @param nprops Number of material properties
 * @param props Material properties vector (see table above)
 * @param nstatev Number of state variables
 * @param statev State variables vector (see table above) [input/output]
 * @param T Temperature at beginning of increment
 * @param DT Temperature increment
 * @param Time Time at beginning of increment
 * @param DTime Time increment
 * @param Wm Total mechanical work [output]
 * @param Wm_r Recoverable (elastic) work [output]
 * @param Wm_ir Irrecoverable work stored in viscous elements [output]
 * @param Wm_d Dissipated (viscous) work [output]
 * @param ndi Number of direct stress components (typically 3)
 * @param nshr Number of shear stress components (typically 3)
 * @param start Flag indicating first increment (true) or continuation (false)
 * @param tnew_dt Suggested new time step size for adaptive time stepping [output]
 *
 * @note The Prony series provides an efficient representation for viscoelastic relaxation
 * @note Relaxation times should span the expected loading time scales
 * @note For frequency-domain data, use fitting algorithms to extract Prony parameters
 * @note The model assumes small strains and linear viscoelasticity
 * @note Time step should be small relative to shortest relaxation time for accuracy
 *
 * @see L_iso() for isotropic elastic stiffness construction
 * @see Ivol() for volumetric projection tensor
 * @see Idev() for deviatoric projection tensor
 *
 * @code
 * // Example usage: 2 Maxwell elements
 * int N = 2;
 * vec props(8);
 * props(0) = N;          // 2 Maxwell elements
 * props(1) = 1000;       // E_inf = 1000 MPa
 * props(2) = 0.3;        // nu_inf = 0.3
 * props(3) = 1e-5;       // alpha = 1e-5 /K
 * props(4) = 500;        // E_1 = 500 MPa
 * props(5) = 0.1;        // tau_1 = 0.1 s
 * props(6) = 200;        // E_2 = 200 MPa
 * props(7) = 1.0;        // tau_2 = 1.0 s
 *
 * vec statev = zeros(1 + 6*N);  // 1 + 12 state variables
 * vec Etot = {0.001, 0, 0, 0, 0, 0};
 * vec DEtot = {0.0001, 0, 0, 0, 0, 0};
 * vec sigma = zeros(6);
 * mat Lt = zeros(6,6);
 * mat DR = eye(3,3);
 *
 * umat_prony_Nfast(Etot, DEtot, sigma, Lt, DR, 8, props, 13, statev,
 *                  20.0, 0.0, 0.0, 0.01, Wm, Wm_r, Wm_ir, Wm_d,
 *                  3, 3, false, tnew_dt);
 * @endcode
 *
 * **References:**
 * - Ferry, J. D. (1980). *Viscoelastic Properties of Polymers* (3rd ed.). Wiley.
 * - Simo, J. C., & Hughes, T. J. R. (1998). *Computational Inelasticity*. Springer.
 * - Park, S. W., & Schapery, R. A. (1999). "Methods of interconversion between linear viscoelastic material functions. Part I—A numerical method based on Prony series." *International Journal of Solids and Structures*, 36(11), 1653-1675.
 */
void umat_prony_Nfast(const std::string &umat_name, const arma::vec &Etot, const arma::vec &DEtot, arma::vec &sigma, arma::mat &Lt, arma::mat &L, const arma::mat &DR, const int &nprops, const arma::vec &props, const int &nstatev, arma::vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, double &tnew_dt);

/** @} */ // end of umat_mechanical group

} //namespace simcoon