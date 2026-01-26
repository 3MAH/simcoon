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
 * @file Zener_fast.hpp
 * @brief Linear viscoelastic material model using the Zener (Standard Linear Solid) representation
 * @author Yves Chemisky
 * @version 1.0
 */

#pragma once

#include <iostream>
#include <armadillo>

namespace simcoon {

/** @addtogroup umat_mechanical
 *  @{
 */

/**
 * @brief Linear viscoelastic constitutive model using the Zener (Standard Linear Solid) model
 *
 * @details This function implements a linear viscoelastic material model using the Zener model,
 * also known as the Standard Linear Solid (SLS). The Zener model consists of a spring in series
 * with a Kelvin-Voigt element (spring and dashpot in parallel), providing both instantaneous
 * elastic response and time-dependent relaxation.
 *
 * **Rheological Representation:**
 *
 * ```
 *     E_0 (instantaneous)
 *   ─────────┬─────────
 *            │
 *     ┌──────┴──────┐
 *     │   E_1      η_1   │  (Kelvin-Voigt element)
 *     │ ──/\/\──┬──[=]──│
 *     │         │        │
 *     └─────────┴────────┘
 * ```
 *
 * **Constitutive Equations:**
 *
 * The Zener model stress-strain relationship is given by:
 * \f[
 * \boldsymbol{\sigma} + \tau \dot{\boldsymbol{\sigma}} = E_0 \boldsymbol{\varepsilon} + (E_0 + E_1) \tau \dot{\boldsymbol{\varepsilon}}
 * \f]
 * where:
 * - \f$ E_0 \f$ is the instantaneous (glassy) modulus
 * - \f$ E_1 \f$ is the modulus of the Kelvin-Voigt element
 * - \f$ \tau = \eta_1 / E_1 \f$ is the relaxation time
 * - \f$ \eta_1 \f$ is the viscosity of the dashpot
 *
 * **Relaxation Modulus:**
 *
 * The relaxation modulus for the Zener model is:
 * \f[
 * E(t) = E_\infty + (E_0 - E_\infty) e^{-t/\tau}
 * \f]
 * where:
 * - \f$ E_\infty = \frac{E_0 E_1}{E_0 + E_1} \f$ is the long-term (rubbery) modulus
 * - \f$ E_0 \f$ is the instantaneous modulus
 * - \f$ \tau \f$ is the relaxation time
 *
 * **Creep Compliance:**
 *
 * The creep compliance is:
 * \f[
 * J(t) = \frac{1}{E_0} + \frac{1}{E_1} \left( 1 - e^{-t/\tau_c} \right)
 * \f]
 * where \f$ \tau_c = \eta_1 / E_1 \f$ is the retardation time.
 *
 * **Incremental Form:**
 *
 * For numerical integration over a time step \f$ \Delta t \f$:
 * \f[
 * \boldsymbol{\sigma}_{n+1} = \mathbf{L}_\infty : \boldsymbol{\varepsilon}_{n+1} + \mathbf{q}_{n+1}
 * \f]
 * where the internal variable \f$ \mathbf{q} \f$ evolves as:
 * \f[
 * \mathbf{q}_{n+1} = e^{-\Delta t/\tau} \mathbf{q}_n + \frac{E_1 - E_\infty}{E_\infty} \mathbf{L}_\infty : \left( e^{-\Delta t/\tau} - 1 \right) \Delta \boldsymbol{\varepsilon}
 * \f]
 *
 * **Material Parameters (props):**
 *
 * | Index | Symbol | Description | Units |
 * |-------|--------|-------------|-------|
 * | props[0] | \f$ E_\infty \f$ | Equilibrium (long-term) Young's modulus | Stress |
 * | props[1] | \f$ \nu \f$ | Poisson's ratio (assumed constant) | - |
 * | props[2] | \f$ \alpha \f$ | Thermal expansion coefficient | 1/Temperature |
 * | props[3] | \f$ E_1 \f$ | Modulus of Kelvin-Voigt element | Stress |
 * | props[4] | \f$ \tau \f$ | Relaxation time | Time |
 *
 * **State Variables (statev):**
 *
 * Total state variables required: \f$ n_{statev} = 7 \f$
 *
 * | Index | Symbol | Description | Units |
 * |-------|--------|-------------|-------|
 * | statev[0] | \f$ T_{init} \f$ | Initial temperature | Temperature |
 * | statev[1:6] | \f$ \mathbf{q} \f$ | Internal stress-like variable (Voigt 6×1) | Stress |
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
 * @note The Zener model is the simplest model capturing both creep and relaxation
 * @note For multiple relaxation times, use the Prony series model (umat_prony_Nfast)
 * @note Time step should be small relative to relaxation time for accuracy
 * @note The model assumes small strains and linear viscoelasticity
 *
 * **References:**
 * - Zener, C. (1948). *Elasticity and Anelasticity of Metals*. University of Chicago Press.
 * - Ferry, J. D. (1980). *Viscoelastic Properties of Polymers* (3rd ed.). Wiley.
 * - Simo, J. C., & Hughes, T. J. R. (1998). *Computational Inelasticity*. Springer.
 */
void umat_zener_fast(const arma::vec &Etot, const arma::vec &DEtot, arma::vec &sigma, arma::mat &Lt, const arma::mat &DR, const int &nprops, const arma::vec &props, const int &nstatev, arma::vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, double &tnew_dt);

/** @} */ // end of umat_mechanical group

} //namespace simcoon
