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
 * @brief Generalized Zener viscoelastic model with N Maxwell elements in parallel
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
 * @brief Generalized Zener (Standard Linear Solid) viscoelastic model with N parallel branches
 *
 * @details This function implements a generalized linear viscoelastic material model using
 * N Zener (Maxwell) elements in parallel with a long-term equilibrium spring. This is
 * equivalent to a Generalized Maxwell model and provides a discrete approximation of the
 * continuous relaxation spectrum.
 *
 * **Rheological Representation:**
 *
 * ```
 *                    E_∞ (equilibrium spring)
 *   ────────────────────/\/\/\────────────────
 *            │                       │
 *            ├──── E_1 ─[η_1]────────┤   Branch 1
 *            │                       │
 *            ├──── E_2 ─[η_2]────────┤   Branch 2
 *            │                       │
 *            │         ...           │
 *            │                       │
 *            └──── E_N ─[η_N]────────┘   Branch N
 * ```
 *
 * **Constitutive Equations:**
 *
 * The relaxation modulus is expressed as a sum of exponentials (Prony series):
 * \f[
 * E(t) = E_\infty + \sum_{i=1}^N E_i e^{-t/\tau_i}
 * \f]
 * where:
 * - \f$ E_\infty \f$ is the long-term (equilibrium) modulus
 * - \f$ E_i \f$ is the modulus of the i-th Maxwell element
 * - \f$ \tau_i = \eta_i / E_i \f$ is the relaxation time of the i-th element
 * - \f$ N \f$ is the number of Maxwell elements
 *
 * **Instantaneous Modulus:**
 *
 * At \f$ t = 0 \f$:
 * \f[
 * E_0 = E(0) = E_\infty + \sum_{i=1}^N E_i
 * \f]
 *
 * **Stress Decomposition:**
 *
 * The total stress is the sum of the equilibrium stress and internal stresses:
 * \f[
 * \boldsymbol{\sigma} = \mathbf{L}_\infty : \boldsymbol{\varepsilon} + \sum_{i=1}^N \mathbf{q}_i
 * \f]
 * where \f$ \mathbf{q}_i \f$ are the internal stress-like variables for each Maxwell element.
 *
 * **Internal Variable Evolution:**
 *
 * Each internal variable evolves according to:
 * \f[
 * \dot{\mathbf{q}}_i + \frac{1}{\tau_i} \mathbf{q}_i = \frac{E_i}{E_\infty} \mathbf{L}_\infty : \dot{\boldsymbol{\varepsilon}}
 * \f]
 *
 * **Incremental Form:**
 *
 * For numerical integration over a time step \f$ \Delta t \f$:
 * \f[
 * \mathbf{q}_{i,n+1} = e^{-\Delta t/\tau_i} \mathbf{q}_{i,n} + \frac{E_i}{E_\infty} \mathbf{L}_\infty : \frac{\tau_i}{\Delta t} \left( 1 - e^{-\Delta t/\tau_i} \right) \Delta \boldsymbol{\varepsilon}
 * \f]
 *
 * **Consistent Tangent Modulus:**
 *
 * For implicit finite element analysis:
 * \f[
 * \mathbf{L}_t = \mathbf{L}_\infty \left( 1 + \sum_{i=1}^N \frac{E_i}{E_\infty} \frac{\tau_i}{\Delta t} \left( 1 - e^{-\Delta t/\tau_i} \right) \right)
 * \f]
 *
 * **Material Parameters (props):**
 *
 * For N Maxwell elements, the props vector contains:
 *
 * | Index | Symbol | Description | Units |
 * |-------|--------|-------------|-------|
 * | props[0] | \f$ N \f$ | Number of Maxwell elements | - |
 * | props[1] | \f$ E_\infty \f$ | Equilibrium Young's modulus | Stress |
 * | props[2] | \f$ \nu \f$ | Poisson's ratio (assumed constant) | - |
 * | props[3] | \f$ \alpha \f$ | Thermal expansion coefficient | 1/Temperature |
 * | props[4] | \f$ E_1 \f$ | Modulus of 1st Maxwell element | Stress |
 * | props[5] | \f$ \tau_1 \f$ | Relaxation time of 1st element | Time |
 * | ... | ... | ... | ... |
 * | props[4+2(i-1)] | \f$ E_i \f$ | Modulus of i-th element | Stress |
 * | props[5+2(i-1)] | \f$ \tau_i \f$ | Relaxation time of i-th element | Time |
 *
 * **State Variables (statev):**
 *
 * For each Maxwell element, 6 internal state variables store the stress-like quantities:
 *
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
 * @note Relaxation times should span the expected loading time scales (decades in log-time)
 * @note Typically 5-10 Maxwell elements are sufficient for most polymers
 * @note For frequency-domain data, use fitting algorithms to extract Prony parameters
 * @note The model assumes small strains and linear viscoelasticity
 * @note Time step should be small relative to shortest relaxation time for accuracy
 *
 * @see umat_zener_fast() for single Maxwell element version
 * @see umat_prony_Nfast() for equivalent Prony series formulation
 *
 * **References:**
 * - Zener, C. (1948). *Elasticity and Anelasticity of Metals*. University of Chicago Press.
 * - Ferry, J. D. (1980). *Viscoelastic Properties of Polymers* (3rd ed.). Wiley.
 * - Simo, J. C., & Hughes, T. J. R. (1998). *Computational Inelasticity*. Springer.
 * - Park, S. W., & Schapery, R. A. (1999). "Methods of interconversion between linear viscoelastic material functions." *Int. J. Solids Struct.*, 36(11), 1653-1675.
 */
void umat_zener_Nfast(const std::string &umat_name, const arma::vec &Etot, const arma::vec &DEtot, arma::vec &sigma, arma::mat &Lt, arma::mat &L, const arma::mat &DR, const int &nprops, const arma::vec &props, const int &nstatev, arma::vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, double &tnew_dt);

/** @} */ // end of umat_mechanical group

} //namespace simcoon
