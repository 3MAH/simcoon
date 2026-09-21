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
 * @brief Generalized Maxwell (Prony series) linear viscoelastic model with N branches
 * @author Yves Chemisky, George Chatzigeorgiou
 * @version 2.0
 */

#pragma once

#include <string>
#include <iostream>
#include <armadillo>
#include <simcoon/parameter.hpp>

namespace simcoon {

/** @addtogroup umat_mechanical
 *  @{
 */

/**
 * @brief Generalized Maxwell (Prony series) linear viscoelastic model with N branches (``PRONK``)
 *
 * @details Isotropic, rate-dependent, linear viscoelastic model with thermal expansion:
 * an equilibrium spring in parallel with \f$ N \f$ Maxwell branches, each a spring
 * \f$ \mathbf{L}_i \f$ in series with a dashpot \f$ \mathbf{H}_i \f$.
 *
 * **Rheological representation:**
 *
 * ```
 *              L_inf = L_0 - sum_i L_i   (equilibrium spring)
 *   ─────────────────────/\/\/\─────────────────────
 *            │                                 │
 *            ├──── L_1 /\/\/\──[ H_1 ]──────────┤   branch 1
 *            │                                 │
 *            │                ...              │
 *            │                                 │
 *            └──── L_N /\/\/\──[ H_N ]──────────┘   branch N
 * ```
 *
 * **Constitutive equations** (Voigt notation, engineering shear strains):
 *
 * The instantaneous stiffness acts on the elastic part of the strain,
 * \f[
 * \boldsymbol{\sigma} = \mathbf{L}_0 : \left( \boldsymbol{\varepsilon} - \boldsymbol{\alpha}\,(T - T_{init}) - \tilde{\boldsymbol{\varepsilon}}^{v} \right),
 * \qquad
 * \tilde{\boldsymbol{\varepsilon}}^{v} = \sum_{i=1}^N \mathbf{M}_0 : \mathbf{L}_i : \boldsymbol{\varepsilon}^{v}_i ,
 * \f]
 * with \f$ \mathbf{L}_0 = \mathbf{L}_{iso}(E_0, \nu_0) \f$, \f$ \mathbf{M}_0 = \mathbf{L}_0^{-1} \f$
 * and \f$ \mathbf{L}_i = \mathbf{L}_{iso}(E_i, \nu_i) \f$; this is the parallel assembly
 * \f$ \boldsymbol{\sigma} = \mathbf{L}_0 : (\boldsymbol{\varepsilon} - \boldsymbol{\varepsilon}^{th}) - \sum_i \mathbf{L}_i : \boldsymbol{\varepsilon}^{v}_i \f$.
 * The dashpot of branch \f$ i \f$ carries the branch stress
 * \f$ \mathbf{A}_i = \mathbf{L}_i : (\boldsymbol{\varepsilon} - \boldsymbol{\varepsilon}^{v}_i) \f$
 * and flows as
 * \f[
 * \dot{\boldsymbol{\varepsilon}}^{v}_i = \mathbf{H}_i^{-1} : \mathbf{L}_i : \left( \boldsymbol{\varepsilon} - \boldsymbol{\varepsilon}^{v}_i \right),
 * \qquad
 * \mathbf{H}_i = 3\,\eta_{B,i}\,\mathbf{I}_{vol} + 2\,\eta_{S,i}\,\mathbf{I}_{dev} ,
 * \f]
 * so each branch has a bulk and a shear relaxation time
 * \f$ \tau_{B,i} = \eta_{B,i} / K_i \f$ and \f$ \tau_{S,i} = \eta_{S,i} / \mu_i \f$, with
 * \f$ K_i = E_i / (3 (1 - 2\nu_i)) \f$ and \f$ \mu_i = E_i / (2 (1 + \nu_i)) \f$.
 *
 * **Instantaneous versus long-term stiffness (convention):**
 *
 * \f$ (E_0, \nu_0) \f$ define the **instantaneous** (glassy) stiffness: at \f$ t = 0 \f$
 * every \f$ \boldsymbol{\varepsilon}^{v}_i \f$ vanishes and
 * \f$ \boldsymbol{\sigma} = \mathbf{L}_0 : \boldsymbol{\varepsilon} \f$. As \f$ t \to \infty \f$
 * each \f$ \boldsymbol{\varepsilon}^{v}_i \f$ tends to \f$ \boldsymbol{\varepsilon} \f$ and the
 * response relaxes to
 * \f[
 * \mathbf{L}_\infty = \mathbf{L}_0 - \sum_{i=1}^N \mathbf{L}_i ,
 * \f]
 * i.e. the uniaxial relaxation modulus is the Prony series
 * \f$ E(t) = E_\infty + \sum_i E_i\, e^{-t/\tau_i} \f$ written with
 * \f$ E(0) = E_0 = E_\infty + \sum_i E_i \f$. The branch moduli must satisfy
 * \f$ \sum_i \mathbf{L}_i < \mathbf{L}_0 \f$ (\f$ \sum_i E_i < E_0 \f$ for equal Poisson
 * ratios); nothing validates it, and a violation gives a negative long-term stiffness (the
 * stress crosses zero during a hold). A calibration that provides the long-term modulus must
 * pass \f$ E_0 = E_\infty + \sum_i E_i \f$.
 *
 * **Integration:**
 *
 * Backward Euler on the branch multipliers \f$ v_i \f$ (accumulated flow length), the flow
 * direction \f$ \boldsymbol{\Lambda}_i = \dot{\boldsymbol{\varepsilon}}^{v}_i / \| \dot{\boldsymbol{\varepsilon}}^{v}_i \| \f$
 * refreshed at every iteration: the residuals
 * \f$ \Phi_i = \| \dot{\boldsymbol{\varepsilon}}^{v}_i \| - \Delta v_i / \Delta t \f$ are solved by
 * Newton--Raphson and the consistent tangent
 * \f[
 * \mathbf{L}_t = \mathbf{L}_0 - \sum_{i=1}^N \frac{(\mathbf{L}_i : \boldsymbol{\Lambda}_i) \otimes \partial_{\boldsymbol{\varepsilon}} \Phi_i}{K_{ii}}
 * \f]
 * follows.
 * When \f$ \Delta t \le \f$ `iota` the branches are **inactive** (\f$ \Phi_i = 0 \f$): the
 * solver probes the tangent with a zero time increment at every block start and commits
 * that answer, so a kernel returning its stationary condition there would relax one
 * branch per block boundary.
 *
 * @note **Stress measure.** The stress returned by this model (the `stress` argument,
 * written \f$ \boldsymbol{\sigma} \f$ in the relations above) is the Cauchy stress under
 * infinitesimal strain; under finite strain the update runs in a corotational frame, so it
 * is the rotated Kirchhoff stress
 * \f$ \hat{\boldsymbol{\tau}} = \boldsymbol{Q}^{T}\boldsymbol{\tau}\,\boldsymbol{Q} \f$ on the
 * frame fixed by the chosen objective rate (\f$ \boldsymbol{Q} = \boldsymbol{R} \f$ for
 * Green--Naghdi and \f$ \log_R \f$, the logarithmic frame for the XBM/log rate,
 * \f$ \boldsymbol{F} \f$ for \f$ \log_F \f$).
 *
 * **Material parameters (props), \f$ n_{props} = 4 + 4N \f$:**
 *
 * | Index | Symbol | Description | Units |
 * |-------|--------|-------------|-------|
 * | props[0] | \f$ E_0 \f$ | Instantaneous Young's modulus | Stress |
 * | props[1] | \f$ \nu_0 \f$ | Instantaneous Poisson's ratio | - |
 * | props[2] | \f$ \alpha \f$ | Isotropic thermal expansion coefficient | 1/Temperature |
 * | props[3] | \f$ N \f$ | Number of Maxwell branches | - |
 * | props[4+4(i-1)] | \f$ E_i \f$ | Young's modulus of branch i | Stress |
 * | props[5+4(i-1)] | \f$ \nu_i \f$ | Poisson's ratio of branch i | - |
 * | props[6+4(i-1)] | \f$ \eta_{B,i} \f$ | Bulk viscosity of branch i | Stress x Time |
 * | props[7+4(i-1)] | \f$ \eta_{S,i} \f$ | Shear viscosity of branch i | Stress x Time |
 *
 * **State variables (statev), \f$ n_{statev} = 7 + 7N \f$:**
 *
 * | Index | Symbol | Description | Units |
 * |-------|--------|-------------|-------|
 * | statev[0] | \f$ T_{init} \f$ | Initial temperature | Temperature |
 * | statev[1:6] | \f$ \tilde{\boldsymbol{\varepsilon}}^{v} \f$ | Total viscous strain \f$ \sum_i \mathbf{M}_0 : \mathbf{L}_i : \boldsymbol{\varepsilon}^{v}_i \f$ (Voigt) | - |
 * | statev[7+7(i-1)] | \f$ v_i \f$ | Accumulated flow length of branch i | - |
 * | statev[8+7(i-1):13+7(i-1)] | \f$ \boldsymbol{\varepsilon}^{v}_i \f$ | Viscous strain of branch i (Voigt) | - |
 *
 * @param umat_name Name of the constitutive law (unused)
 * @param Etot Total strain at the beginning of the increment (Voigt notation: \f$6 \times 1\f$ vector)
 * @param DEtot Strain increment (Voigt notation: \f$6 \times 1\f$ vector)
 * @param stress Stress (Voigt notation: \f$6 \times 1\f$ vector) [input/output]
 * @param Lt Consistent tangent modulus (\f$6 \times 6\f$ matrix) [output]
 * @param L Instantaneous elastic stiffness \f$ \mathbf{L}_0 \f$ (\f$6 \times 6\f$ matrix) [output]
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
 * @param ndi Number of direct stress components (typically 3)
 * @param nshr Number of shear stress components (typically 3)
 * @param start Flag indicating the first increment (true) or a continuation (false)
 * @param tnew_dt Suggested new time step size for adaptive time stepping [output]
 * @param tangent_mode Tangent operator selection (see parameter.hpp)
 *
 * @note Not equivalent to umat_zener_Nfast(): identical props layout, but that kernel is a
 * generalized Kelvin chain (branches in series with the spring), not a generalized Maxwell model
 * @note The thermal strain enters the stress, not the branch flow: branch i is driven by
 * \f$ \boldsymbol{\varepsilon} - \boldsymbol{\varepsilon}^{v}_i \f$
 * @note Relaxation times should span the loading time scales; the time step should stay small
 * against the shortest one for accuracy
 *
 * @see umat_prony_Nfast_T() thermomechanical twin (props prefixed by \f$ \rho, c_p \f$)
 * @see umat_zener_Nfast() same props layout, different rheology
 * @see ViscoelasticMechanism the modular (``MODUL``) twin of this kernel
 * @see H_iso() for the bulk/shear viscosity tensor
 *
 * @code
 * // One Maxwell branch: E_0 = 3000 MPa, nu_0 = 0.35, no thermal expansion,
 * // E_1 = 1500 MPa, nu_1 = 0.35, etaB_1 = 3000 MPa.s, etaS_1 = 1200 MPa.s
 * // (long-term modulus E_0 - E_1 = 1500 MPa).
 * vec props = {3000., 0.35, 0., 1., 1500., 0.35, 3000., 1200.};
 * vec statev = zeros(7 + 7*1);
 * @endcode
 */
void umat_prony_Nfast(const std::string &umat_name, const arma::vec &Etot, const arma::vec &DEtot, arma::vec &stress, arma::mat &Lt, arma::mat &L, const arma::mat &DR, const int &nprops, const arma::vec &props, const int &nstatev, arma::vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, double &tnew_dt, const int &tangent_mode = tangent_default);

/** @} */ // end of umat_mechanical group

} //namespace simcoon
