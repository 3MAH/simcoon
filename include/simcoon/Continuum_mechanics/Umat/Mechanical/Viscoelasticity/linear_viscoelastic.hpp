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
 * @file linear_viscoelastic.hpp
 * @brief Closed-form backward-Euler update of the linear viscoelastic rheologies (Zener, Prony),
 *        with its exact consistent tangent.
 *
 * The branches are linear springs \f$ \mathbf{L}_i \f$ and dashpots \f$ \mathbf{H}_i \f$
 * (H_iso), so the implicit Euler step is a linear system in the viscous strains
 * \f$ \boldsymbol{\varepsilon}^v_i \f$: it is solved exactly, without iteration, and its
 * derivatives with respect to the strain and the temperature are the consistent tangent.
 * With \f$ \mathbf{B}_i = \mathbf{I} + \Delta t\,\mathbf{H}_i^{-1}\mathbf{L}_i \f$ and
 * \f$ \mathbf{C}_i = \Delta t\,(\mathbf{H}_i + \Delta t\,\mathbf{L}_i)^{-1} = \Delta t\,\mathbf{B}_i^{-1}\mathbf{H}_i^{-1} \f$
 * (symmetric), and \f$ \boldsymbol{\varepsilon}^e = \boldsymbol{\varepsilon}_{n+1}
 * - \boldsymbol{\alpha}(T_{n+1} - T_{init}) \f$. A zero time increment leaves the branches
 * inactive (\f$ \mathbf{C}_i = \mathbf{0} \f$). Engineering Voigt, MPa.
 */

#pragma once

#include <vector>
#include <armadillo>

namespace simcoon {

/**
 * @brief Backward-Euler state and consistent derivatives of a linear viscoelastic step.
 */
struct LinearViscoStep {
    arma::vec sigma;                    ///< Stress at the end of the increment (6)
    std::vector<arma::vec> EV_i;        ///< Viscous strain of each branch (6 each)
    arma::mat dSdE;                     ///< \f$ \partial\boldsymbol{\sigma}/\partial\boldsymbol{\varepsilon} \f$ (6x6, symmetric)
    arma::vec dSdT;                     ///< \f$ \partial\boldsymbol{\sigma}/\partial T \f$ (6)
    std::vector<arma::mat> dEVdE_i;     ///< \f$ \partial\boldsymbol{\varepsilon}^v_i/\partial\boldsymbol{\varepsilon} \f$ (6x6 each)
    std::vector<arma::vec> dEVdT_i;     ///< \f$ \partial\boldsymbol{\varepsilon}^v_i/\partial T \f$ (6 each)
};

/**
 * @brief Kelvin branches in series with the spring \f$ \mathbf{L}_0 \f$ (Zener, generalised Kelvin).
 *
 * \f$ \boldsymbol{\sigma} = \mathbf{L}_0 (\boldsymbol{\varepsilon}^e - \sum_i \boldsymbol{\varepsilon}^v_i) \f$ and,
 * per branch, \f$ \boldsymbol{\sigma} = \mathbf{L}_i \boldsymbol{\varepsilon}^v_i + \mathbf{H}_i \dot{\boldsymbol{\varepsilon}}^v_i \f$.
 * The implicit step gives \f$ \boldsymbol{\varepsilon}^v_i = \mathbf{B}_i^{-1}\boldsymbol{\varepsilon}^v_{i,n} + \mathbf{C}_i \boldsymbol{\sigma} \f$ and
 * \f[ (\mathbf{I} + \mathbf{L}_0 \textstyle\sum_i \mathbf{C}_i)\,\boldsymbol{\sigma}
 *     = \mathbf{L}_0 (\boldsymbol{\varepsilon}^e - \textstyle\sum_i \mathbf{B}_i^{-1}\boldsymbol{\varepsilon}^v_{i,n}), \qquad
 *     \frac{\partial\boldsymbol{\sigma}}{\partial\boldsymbol{\varepsilon}} = (\mathbf{L}_0^{-1} + \textstyle\sum_i \mathbf{C}_i)^{-1}. \f]
 *
 * @param L0 instantaneous spring (6x6)
 * @param L_i branch springs
 * @param H_i branch viscosities (H_iso)
 * @param EV_i_start branch viscous strains at the start of the increment
 * @param eps_e \f$ \boldsymbol{\varepsilon}_{n+1} - \boldsymbol{\alpha}(T_{n+1} - T_{init}) \f$
 * @param alpha thermal expansion tensor (6)
 * @param DTime time increment
 * @return the end state and its consistent derivatives
 */
LinearViscoStep kelvin_series_step(const arma::mat &L0, const std::vector<arma::mat> &L_i,
                                   const std::vector<arma::mat> &H_i,
                                   const std::vector<arma::vec> &EV_i_start, const arma::vec &eps_e,
                                   const arma::vec &alpha, const double &DTime);

/**
 * @brief Maxwell branches in parallel with the spring \f$ \mathbf{L}_0 \f$ (Prony series).
 *
 * \f$ \boldsymbol{\sigma} = \mathbf{L}_0 \boldsymbol{\varepsilon}^e - \sum_i \mathbf{L}_i \boldsymbol{\varepsilon}^v_i \f$
 * (\f$ \mathbf{L}_0 \f$ is the instantaneous stiffness) and, per branch,
 * \f$ \mathbf{L}_i(\boldsymbol{\varepsilon}^e - \boldsymbol{\varepsilon}^v_i) = \mathbf{H}_i \dot{\boldsymbol{\varepsilon}}^v_i \f$:
 * the branches decouple, \f$ \boldsymbol{\varepsilon}^v_i = \mathbf{B}_i^{-1}\boldsymbol{\varepsilon}^v_{i,n}
 * + \mathbf{C}_i\mathbf{L}_i\boldsymbol{\varepsilon}^e \f$ and
 * \f$ \partial\boldsymbol{\sigma}/\partial\boldsymbol{\varepsilon} = \mathbf{L}_0 - \sum_i \mathbf{L}_i\mathbf{C}_i\mathbf{L}_i \f$.
 *
 * @param L0 instantaneous spring (6x6)
 * @param L_i branch springs
 * @param H_i branch viscosities (H_iso)
 * @param EV_i_start branch viscous strains at the start of the increment
 * @param eps_e \f$ \boldsymbol{\varepsilon}_{n+1} - \boldsymbol{\alpha}(T_{n+1} - T_{init}) \f$
 * @param alpha thermal expansion tensor (6)
 * @param DTime time increment
 * @return the end state and its consistent derivatives
 */
LinearViscoStep maxwell_parallel_step(const arma::mat &L0, const std::vector<arma::mat> &L_i,
                                      const std::vector<arma::mat> &H_i,
                                      const std::vector<arma::vec> &EV_i_start, const arma::vec &eps_e,
                                      const arma::vec &alpha, const double &DTime);

/**
 * @brief Heat source of a viscous thermomechanical step and its exact derivatives.
 *
 * \f[ r = \frac{\Delta\gamma - T_m\,\boldsymbol{\alpha} : (\boldsymbol{\sigma} - \boldsymbol{\sigma}_n)
 *     - \rho c_p\,\Delta T}{\Delta t}, \qquad T_m = T + \tfrac12\Delta T \f]
 * (the midpoint temperature, as \f$ W_t \f$), with \f$ \Delta\gamma \f$ the dissipation of the branches:
 * the heat of the actual increments, which keeps flowing while the branches relax under a strain
 * hold. \f$ \partial r/\partial\boldsymbol{\varepsilon} \f$ and \f$ \partial r/\partial T \f$ follow
 * from the step's consistent derivatives. A zero time increment gives \f$ r = 0 \f$.
 *
 * @param Dgamma dissipation of the increment, \f$ \Delta\gamma \f$
 * @param dDgamma_dE \f$ \partial\Delta\gamma/\partial\boldsymbol{\varepsilon} \f$ (6)
 * @param dDgamma_dT \f$ \partial\Delta\gamma/\partial T \f$
 * @param st the step (its dSdE, dSdT)
 * @param alpha thermal expansion tensor (6)
 * @param sigma stress at the end of the increment (6)
 * @param sigma_start stress at the start of the increment (6)
 * @param T temperature at the start of the increment
 * @param DT temperature increment
 * @param rho_cp \f$ \rho c_p \f$
 * @param DTime time increment
 * @param[out] r heat source
 * @param[out] drdE \f$ \partial r/\partial\boldsymbol{\varepsilon} \f$ (6x1, the thermomechanical UMAT layout)
 * @param[out] drdT \f$ \partial r/\partial T \f$ (1x1)
 */
void viscous_heat_source(const double &Dgamma, const arma::vec &dDgamma_dE, const double &dDgamma_dT,
                         const LinearViscoStep &st, const arma::vec &alpha, const arma::vec &sigma,
                         const arma::vec &sigma_start, const double &T, const double &DT,
                         const double &rho_cp, const double &DTime,
                         double &r, arma::mat &drdE, arma::mat &drdT);

} //namespace simcoon
