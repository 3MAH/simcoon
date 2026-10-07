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

///@file plastic_johnson_cook_ccp.hpp
///@brief Johnson-Cook elastic-viscoplastic UMAT (EPJCK): J2 flow with strain-rate and thermal-softening
///       dependent yield stress, integrated by the convex cutting plane (CCP) algorithm
///@version 1.0

#pragma once
#include <string>
#include <armadillo>
#include <simcoon/parameter.hpp>

namespace simcoon{

/** @addtogroup umat_mechanical
 *  @{
 */

/**
 * @brief Johnson-Cook elastic-viscoplastic constitutive law (EPJCK), mechanical version.
 *
 * @details J2 (von Mises) plasticity with associated flow and the Johnson-Cook yield stress
 * (Johnson & Cook, 1983), which multiplies a power-law strain hardening by a logarithmic
 * strain-rate factor and a thermal-softening factor:
 * \f[
 *   \sigma_Y\left(p,\dot p,\theta\right) = \left(A + B\,p^{n}\right)
 *   \left(1 + C\,\ln\frac{\dot p}{\dot\varepsilon_0}\right)
 *   \left(1 - \theta^{*m}\right), \qquad
 *   \theta^{*} = \frac{\theta - \theta_{\mathrm{ref}}}{\theta_{\mathrm{melt}} - \theta_{\mathrm{ref}}} .
 * \f]
 * The yield function, flow rule and loading/unloading conditions are
 * \f[
 *   \Phi = \bar{\sigma} - \sigma_Y\left(p,\dot p,\theta\right) \leq 0, \qquad
 *   \dot{\boldsymbol{\varepsilon}}^{p} = \dot p\,\boldsymbol{\Lambda}, \quad
 *   \boldsymbol{\Lambda} = \frac{3}{2}\frac{\boldsymbol{\sigma}'}{\bar{\sigma}}, \qquad
 *   \dot p \geq 0, \quad \dot p\,\Phi = 0 ,
 * \f]
 * with the stress given by \f$ \boldsymbol{\sigma} = \mathbf{L} : \left(\boldsymbol{\varepsilon} -
 * \boldsymbol{\alpha}\,(\theta - \theta_0) - \boldsymbol{\varepsilon}^{p}\right) \f$.
 *
 * **Regularization of the Johnson-Cook law.** The logarithmic rate factor is clamped at its
 * reference value, \f$ \ln(\dot p/\dot\varepsilon_0) \to \max\left(\ln(\dot p/\dot\varepsilon_0), 0\right) \f$:
 * below \f$ \dot\varepsilon_0 \f$ (and in particular at the onset of yielding, \f$ \dot p = 0 \f$)
 * the material responds with its quasi-static, reference-rate yield stress. The homologous
 * temperature is clamped to \f$ [0, 1) \f$ so that the law is defined below
 * \f$ \theta_{\mathrm{ref}} \f$ and never reaches the singular melting point.
 *
 * **Discrete update.** The plastic strain rate is treated fully implicitly over the increment,
 * \f$ \dot p = \Delta p / \Delta t \f$. The rate factor therefore depends on the unknown
 * \f$ \Delta p \f$ and its derivative enters the scalar Jacobian of the CCP iteration:
 * \f[
 *   K = \frac{\partial \Phi}{\partial p} + \frac{\partial \Phi}{\partial \Delta p}
 *     = -\,n B p^{n-1}\, f_{\dot\varepsilon}\, f_\theta
 *       - \left(A + B p^{n}\right) \frac{C}{\dot p\,\Delta t}\, f_\theta ,
 * \f]
 * so the iterates converge on the rate-dependent yield surface. With \f$ \Delta t = 0 \f$ the
 * law degenerates to its rate-independent form.
 *
 * **Tangent operator.** Assembled by compute_tangent_operator() from the converged
 * \f$ \hat{B} \f$, \f$ \boldsymbol{\kappa} = \mathbf{L}:\boldsymbol{\Lambda} \f$ and
 * \f$ \partial\Phi/\partial\boldsymbol{\sigma} \f$, in the three modes of parameter.hpp:
 * - `tangent_none` (0): \f$ \mathbf{L}_t = \mathbf{L} \f$, for explicit integration schemes;
 * - `tangent_continuum` (1): \f$ \mathbf{L}_t = \mathbf{L} - \hat{B}^{-1}\,
 *   \boldsymbol{\kappa}\otimes\left(\mathbf{L}:\partial\Phi/\partial\boldsymbol{\sigma}\right) \f$;
 * - `tangent_algorithmic` (2, default): the Simo-Hughes consistent operator; the flow
 *   Hessian is \f$ \partial\boldsymbol{\Lambda}/\partial\boldsymbol{\sigma} \f$ (deta_stress) and
 *   the rate term is already part of \f$ \hat{B} \f$, so the operator is the exact Jacobian of
 *   the discrete map for a fixed \f$ \Delta t \f$.
 *
 * **Energy split.** The hardening force is \f$ A_p = -B p^{n} \f$ (reference-rate, isothermal
 * hardening stress), as in EPICP; the rate and thermal multipliers act on the dissipation:
 * \f$ \Delta W_m^{ir} = -\tfrac12 (A_p^n + A_p^{n+1})\,\Delta p \f$,
 * \f$ \Delta W_m^{d} = \tfrac12 (\boldsymbol{\sigma}^n + \boldsymbol{\sigma}^{n+1}):\Delta\boldsymbol{\varepsilon}^{p}
 * + \tfrac12 (A_p^n + A_p^{n+1})\,\Delta p \f$.
 *
 * **Material properties (props, 11):**
 * | Index | Symbol                      | Description                                   |
 * |-------|-----------------------------|-----------------------------------------------|
 * | 0     | \f$ E \f$                   | Young's modulus                               |
 * | 1     | \f$ \nu \f$                 | Poisson's ratio                               |
 * | 2     | \f$ \alpha \f$              | Coefficient of thermal expansion (isotropic)  |
 * | 3     | \f$ A \f$                   | Initial yield stress                          |
 * | 4     | \f$ B \f$                   | Hardening coefficient                         |
 * | 5     | \f$ n \f$                   | Hardening exponent                            |
 * | 6     | \f$ C \f$                   | Strain-rate sensitivity                       |
 * | 7     | \f$ \dot\varepsilon_0 \f$   | Reference strain rate (1/s)                   |
 * | 8     | \f$ m \f$                   | Thermal-softening exponent                    |
 * | 9     | \f$ \theta_{\mathrm{ref}} \f$  | Reference temperature (K)                  |
 * | 10    | \f$ \theta_{\mathrm{melt}} \f$ | Melting temperature (K)                    |
 *
 * **State variables (statev, 9):**
 * | Index | Symbol                             | Description                                      |
 * |-------|------------------------------------|--------------------------------------------------|
 * | 0     | \f$ \theta_0 \f$                   | Initial temperature                              |
 * | 1     | \f$ p \f$                          | Accumulated plastic strain                       |
 * | 2-7   | \f$ \boldsymbol{\varepsilon}^{p} \f$ | Plastic strain (Voigt 11 22 33 12 13 23)       |
 * | 8     | \f$ \dot p \f$                     | Plastic strain rate over the last increment (output) |
 *
 * The same statev layout as EPICP (plastic strain at offset 2) is declared in the finite-strain
 * convention table of umat_smart.cpp, so EPJCK is served by the log-strain / Kirchhoff box
 * route under every objective rate.
 *
 * @param umat_name Name of the constitutive law (EPJCK)
 * @param Etot Total strain at \f$ t_n \f$ (6x1 Voigt)
 * @param DEtot Strain increment (6x1 Voigt)
 * @param sigma Stress [input/output] (6x1 Voigt)
 * @param Lt Tangent modulus [output] (6x6), see @p tangent_mode
 * @param L Elastic stiffness tensor [output] (6x6)
 * @param DR Rotation increment (3x3)
 * @param nprops Number of material properties (11)
 * @param props Material properties vector
 * @param nstatev Number of state variables (9)
 * @param statev State variables vector [input/output]
 * @param T Temperature at \f$ t_n \f$
 * @param DT Temperature increment
 * @param Time Time at \f$ t_n \f$
 * @param DTime Time increment (drives the plastic strain rate \f$ \dot p = \Delta p/\Delta t \f$)
 * @param Wm Total mechanical work [input/output]
 * @param Wm_r Recoverable mechanical work [input/output]
 * @param Wm_ir Irrecoverable (stored) mechanical work [input/output]
 * @param Wm_d Dissipated mechanical work [input/output]
 * @param ndi Number of direct stress components
 * @param nshr Number of shear stress components
 * @param start Flag for the first increment
 * @param tnew_dt Suggested time step ratio [output]
 * @param tangent_mode Tangent assembly mode: tangent_none (0, elastic L for explicit
 *        integration), tangent_continuum (1) or tangent_algorithmic (2, default)
 */
void umat_plasticity_johnson_cook_CCP(const std::string &umat_name, const arma::vec &Etot, const arma::vec &DEtot, arma::vec &sigma, arma::mat &Lt, arma::mat &L, const arma::mat &DR, const int &nprops, const arma::vec &props, const int &nstatev, arma::vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, double &tnew_dt, const int &tangent_mode = tangent_default);

/** @} */ // end of umat_mechanical group

} //namespace simcoon
