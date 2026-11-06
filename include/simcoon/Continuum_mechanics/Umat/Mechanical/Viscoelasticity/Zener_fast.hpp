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
 * @brief Standard linear solid, Poynting--Thomson (Kelvin) form (ZENER)
 * @author Yves Chemisky
 * @version 1.0
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
 * @brief Standard linear solid in its Poynting--Thomson (Kelvin) form (``ZENER``): a spring
 *        \f$ \mathbf{L}_0 \f$ in series with one Kelvin--Voigt unit \f$ (\mathbf{L}_1, \mathbf{H}_1) \f$.
 *
 * @details The instantaneous spring and the Kelvin unit carry the same stress:
 * \f[
 * \boldsymbol{\sigma} = \mathbf{L}_0 : \left( \boldsymbol{\varepsilon} - \boldsymbol{\alpha}\,(T - T_{init})
 *   - \boldsymbol{\varepsilon}^{v} \right)
 *   = \mathbf{L}_1 : \boldsymbol{\varepsilon}^{v} + \mathbf{H}_1 : \dot{\boldsymbol{\varepsilon}}^{v},
 * \f]
 * with \f$ \mathbf{L}_0 = \mathbf{L}_{iso}(E_0, \nu_0) \f$, \f$ \mathbf{L}_1 = \mathbf{L}_{iso}(E_1, \nu_1) \f$
 * and the viscosity tensor \f$ \mathbf{H}_1 = 3\eta_{B,1}\,\mathbf{I}_{vol} + 2\eta_{S,1}\,\mathbf{I}_{dev} \f$
 * (H_iso). ``ZENNK`` (umat_zener_Nfast()) chains \f$ N \f$ such units in series.
 *
 * **One-dimensional response.** With the dashpot viscosity \f$ \eta_1 \f$:
 * \f[
 * \sigma + \tau_r\,\dot{\sigma} = E_\infty\,\varepsilon + E_0\,\tau_r\,\dot{\varepsilon},\qquad
 * J(t) = \frac{1}{E_0} + \frac{1}{E_1}\left(1 - e^{-t/\tau_c}\right),\qquad
 * E(t) = E_\infty + (E_0 - E_\infty)\,e^{-t/\tau_r},
 * \f]
 * with the retardation time \f$ \tau_c = \eta_1 / E_1 \f$, the relaxation time
 * \f$ \tau_r = \eta_1 / (E_0 + E_1) \f$ and the relaxed modulus \f$ E_\infty = E_0 E_1 / (E_0 + E_1) \f$.
 * The same solid in its Zener (Maxwell) form, a spring \f$ E_\infty \f$ in parallel with a Maxwell
 * branch \f$ (E_M, \eta_M) \f$, has \f$ E_M = E_0^2 / (E_0 + E_1) \f$ and
 * \f$ \eta_M = \eta_1 E_0^2 / (E_0 + E_1)^2 \f$. In 3D the props of the two forms do not map one to
 * one (each branch has its own Poisson ratio and bulk and shear viscosities): this law is not the
 * one-branch umat_prony_Nfast() (``PRONK``, generalized Maxwell).
 *
 * **Discrete update.** Backward Euler on \f$ \boldsymbol{\varepsilon}^{v} \f$, solved in closed form
 * (kelvin_series_step, linear_viscoelastic.hpp); the tangent is the exact derivative of that update
 * (\f$ \mathbf{L}_0 \f$ under tangent_none).
 *
 * **Work split.** \f$ \Delta W_d = \tfrac{1}{2}\left(\mathbf{A}^{v}_n + \mathbf{A}^{v}_{n+1}\right) : \Delta\boldsymbol{\varepsilon}^{v} \f$
 * with \f$ \mathbf{A}^{v} = \boldsymbol{\sigma} - \mathbf{L}_1 : \boldsymbol{\varepsilon}^{v} \f$,
 * \f$ W_r = W - W_d \f$, \f$ W_{ir} = 0 \f$.
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
 * **Material properties** (7; MPa, MPa.s):
 *
 * | Index | Symbol | Description |
 * |-------|--------|-------------|
 * | props[0] | \f$ E_0 \f$ | Young modulus of the instantaneous spring |
 * | props[1] | \f$ \nu_0 \f$ | Poisson ratio of the instantaneous spring |
 * | props[2] | \f$ \alpha \f$ | isotropic thermal expansion coefficient |
 * | props[3] | \f$ E_1 \f$ | Young modulus of the Kelvin unit |
 * | props[4] | \f$ \nu_1 \f$ | Poisson ratio of the Kelvin unit |
 * | props[5] | \f$ \eta_{B,1} \f$ | bulk viscosity of the Kelvin unit |
 * | props[6] | \f$ \eta_{S,1} \f$ | shear viscosity of the Kelvin unit |
 *
 * **State variables** (8; Voigt order [11, 22, 33, 12, 13, 23], engineering shear strains):
 *
 * | Index | Symbol | Description |
 * |-------|--------|-------------|
 * | statev[0] | \f$ T_{init} \f$ | initial temperature |
 * | statev[1] | \f$ v \f$ | accumulated viscous strain \f$ \sum \|\Delta\boldsymbol{\varepsilon}^{v}\| \f$ |
 * | statev[2..7] | \f$ \boldsymbol{\varepsilon}^{v} \f$ | viscous strain of the Kelvin unit |
 *
 * @param Etot Total strain tensor at beginning of increment (Voigt notation: \f$6 \times 1\f$ vector)
 * @param DEtot Strain increment tensor (Voigt notation: \f$6 \times 1\f$ vector)
 * @param stress Stress tensor (Voigt notation: \f$6 \times 1\f$ vector) [output]
 * @param Lt Consistent tangent modulus (\f$6 \times 6\f$ matrix) [output]
 * @param DR Rotation increment matrix (\f$3 \times 3\f$) for objective integration
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
 * **References:**
 * - Zener, C. (1948). *Elasticity and Anelasticity of Metals*. University of Chicago Press.
 * - Tschoegl, N. W. (1989). *The Phenomenological Theory of Linear Viscoelastic Behavior*. Springer.
 */
void umat_zener_fast(const std::string &umat_name, const arma::vec &Etot, const arma::vec &DEtot, arma::vec &stress, arma::mat &Lt, arma::mat &L, const arma::mat &DR, const int &nprops, const arma::vec &props, const int &nstatev, arma::vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, double &tnew_dt, const int &tangent_mode = tangent_default);

/** @} */ // end of umat_mechanical group

} //namespace simcoon
