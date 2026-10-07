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
///@brief Thermomechanical Johnson-Cook elastic-viscoplastic UMAT (EPJCK): strain-rate and
///       thermal-softening dependent yield stress, exact heat source from the Gibbs framework
///@version 1.0

#pragma once
#include <armadillo>

namespace simcoon{

/** @addtogroup umat_thermomechanical
 *  @{
 */

/**
 * @brief Johnson-Cook elastic-viscoplastic constitutive law (EPJCK), thermomechanical version.
 *
 * @details Same constitutive model and CCP integration as the mechanical kernel
 * umat_plasticity_johnson_cook_CCP() (see that header for the yield stress, the regularization
 * and the discrete rate treatment), extended with the coupled thermal outputs of the
 * thermomechanical framework: the heat source \f$ r \f$ and the four linearizations
 * \f$ \partial\boldsymbol{\sigma}/\partial\boldsymbol{\varepsilon} \f$,
 * \f$ \partial\boldsymbol{\sigma}/\partial\theta \f$, \f$ \partial r/\partial\boldsymbol{\varepsilon} \f$,
 * \f$ \partial r/\partial\theta \f$ consumed by the coupled Newton-Raphson solver.
 *
 * **Thermal coupling.** Unlike EPICP, the yield function depends explicitly on the temperature
 * through the thermal-softening factor \f$ f_\theta = 1 - \theta^{*m} \f$:
 * \f[
 *   \frac{\partial\Phi}{\partial\theta} = \left(A + B p^{n}\right) f_{\dot\varepsilon}\,
 *   \frac{m\,\theta^{*\,m-1}}{\theta_{\mathrm{melt}} - \theta_{\mathrm{ref}}} > 0 ,
 * \f]
 * so a temperature rise lowers the yield stress and feeds back on the plastic flow through
 * \f$ \mathbf{P}_\theta = \hat{B}^{-1}\left(\partial\Phi/\partial\theta -
 * \partial\Phi/\partial\boldsymbol{\sigma} : \mathbf{L} : \boldsymbol{\alpha}\right) \f$ and
 * \f$ \partial\boldsymbol{\sigma}/\partial\theta = -\mathbf{L}:\boldsymbol{\alpha} -
 * \boldsymbol{\kappa}\,P_\theta \f$.
 *
 * **Heat source.** No Taylor-Quinney coefficient is introduced: \f$ r \f$ is the sum of the
 * thermoelastic term \f$ N \f$ of the Gibbs framework, linearized with respect to the strain
 * and temperature increments with the entropy
 * \f$ \eta = c_0 \ln(\theta/\theta_0) + \boldsymbol{\alpha}:\boldsymbol{\sigma} \f$, and of
 * the intrinsic dissipation evaluated on the converged increment,
 * \f$ \Gamma\,\Delta t = \tfrac12 (\boldsymbol{\sigma}^n + \boldsymbol{\sigma}^{n+1}):\Delta\boldsymbol{\varepsilon}^{p}
 * + \tfrac12 (A_p^n + A_p^{n+1})\,\Delta p \f$ with \f$ A_p = -B p^{n} \f$ (the increment
 * \f$ W_m^d \f$ accumulates, same energy split as the mechanical kernel). EPICP's
 * thermomechanical kernel reconstructs \f$ \Gamma \f$ from its linearization
 * \f$ \boldsymbol{\Gamma}_\varepsilon : \Delta\boldsymbol{\varepsilon} + \Gamma_\theta \Delta\theta \f$,
 * which is exact to first order for a rate-independent law. It is not for Johnson-Cook: at
 * fixed \f$ \Delta t \f$ the plastic increment is a logarithmic function of the strain
 * increment, the rate term \f$ C (A + B p^n) / \Delta p \f$ in \f$ \hat{B} \f$ drives
 * \f$ \mathbf{P}_\varepsilon \to 0 \f$ under step refinement and the linearized dissipation
 * vanishes while the physical one stays finite. The linearizations
 * \f$ \partial r/\partial\boldsymbol{\varepsilon} \f$ and \f$ \partial r/\partial\theta \f$
 * keep the \f$ \boldsymbol{\Gamma}_\varepsilon \f$, \f$ \Gamma_\theta \f$ terms as the Newton
 * derivatives.
 *
 * **Tangent modes.** The mechanical block \f$ \partial\boldsymbol{\sigma}/\partial\boldsymbol{\varepsilon} \f$
 * follows @p tangent_mode (continuum or algorithmic). As in the other thermomechanical kernels,
 * `tangent_none` is promoted to the continuum operator: the sensitivities
 * \f$ \mathbf{P}_\varepsilon \f$ and \f$ \hat{B}^{-1} \f$ feed the physical heat source and
 * its linearization, not only the Newton operator.
 *
 * **Material properties (props, 13):**
 * | Index | Symbol                         | Description                                  |
 * |-------|--------------------------------|----------------------------------------------|
 * | 0     | \f$ \rho \f$                   | Density                                      |
 * | 1     | \f$ c_p \f$                    | Specific heat capacity                       |
 * | 2     | \f$ E \f$                      | Young's modulus                              |
 * | 3     | \f$ \nu \f$                    | Poisson's ratio                              |
 * | 4     | \f$ \alpha \f$                 | Coefficient of thermal expansion (isotropic) |
 * | 5     | \f$ A \f$                      | Initial yield stress                         |
 * | 6     | \f$ B \f$                      | Hardening coefficient                        |
 * | 7     | \f$ n \f$                      | Hardening exponent                           |
 * | 8     | \f$ C \f$                      | Strain-rate sensitivity                      |
 * | 9     | \f$ \dot\varepsilon_0 \f$      | Reference strain rate (1/s)                  |
 * | 10    | \f$ m \f$                      | Thermal-softening exponent                   |
 * | 11    | \f$ \theta_{\mathrm{ref}} \f$  | Reference temperature (K)                    |
 * | 12    | \f$ \theta_{\mathrm{melt}} \f$ | Melting temperature (K)                      |
 *
 * **State variables (statev, 9):** \f$ \theta_0 \f$, \f$ p \f$, \f$ \boldsymbol{\varepsilon}^{p} \f$
 * (6, Voigt), \f$ \dot p \f$ (output), the layout of the mechanical kernel.
 *
 * @param Etot Total strain at \f$ t_n \f$ (6x1 Voigt)
 * @param DEtot Strain increment (6x1 Voigt)
 * @param sigma Stress [input/output] (6x1 Voigt)
 * @param r Internal heat production [output]
 * @param dSdE Mechanical tangent \f$ \partial\boldsymbol{\sigma}/\partial\boldsymbol{\varepsilon} \f$ [output] (6x6)
 * @param dSdT Thermal stress sensitivity \f$ \partial\boldsymbol{\sigma}/\partial\theta \f$ [output] (6x1)
 * @param drdE Heat source strain sensitivity \f$ \partial r/\partial\boldsymbol{\varepsilon} \f$ [output] (6x1)
 * @param drdT Heat source temperature sensitivity \f$ \partial r/\partial\theta \f$ [output] (1x1)
 * @param DR Rotation increment (3x3)
 * @param nprops Number of material properties (13)
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
 * @param Wt Total thermal work [input/output]
 * @param Wt_r Reversible thermal work [input/output]
 * @param Wt_ir Irreversible thermal work [input/output]
 * @param ndi Number of direct stress components
 * @param nshr Number of shear stress components
 * @param start Flag for the first increment
 * @param tnew_dt Suggested time step ratio [output]
 * @param tangent_mode Tangent assembly mode (parameter.hpp); tangent_none is promoted to
 *        tangent_continuum, see above
 */
void umat_plasticity_johnson_cook_CCP_T(const arma::vec &Etot, const arma::vec &DEtot, arma::vec &sigma, double &r, arma::mat &dSdE, arma::mat &dSdT, arma::mat &drdE, arma::mat &drdT, const arma::mat &DR, const int &nprops, const arma::vec &props, const int &nstatev, arma::vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, double &Wt, double &Wt_r, double &Wt_ir, const int &ndi, const int &nshr, const bool &start, double &tnew_dt, const int &tangent_mode = 0);

/** @} */ // end of umat_thermomechanical group

} //namespace simcoon
