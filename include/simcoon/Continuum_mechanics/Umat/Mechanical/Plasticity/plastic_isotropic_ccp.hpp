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

///@file plastic_isotropic_ccp.hpp
///@brief User subroutine for elastic-plastic materials in 1D-2D-3D case
///@brief This subroutines uses a convex cutting plane algorithm
///@brief Isotropic hardening with a power-law hardenig is considered
///@version 1.0

#pragma once
#include <string>
#include <armadillo>

namespace simcoon{

/**
 * @file plastic_isotropic_ccp.hpp
 * @brief Elastic-plastic material model with isotropic hardening using the Convex Cutting Plane algorithm
 * @author Yves Chemisky
 * @version 1.0
 */

/** @addtogroup umat_mechanical
 *  @{
 */

/**
 * @brief Elastic-plastic constitutive model with isotropic hardening solved by the Convex Cutting Plane (CCP) algorithm
 *
 * @details This function implements an elastic-plastic material model for small and finite strain analysis.
 * The model features:
 * - J2 (von Mises) plasticity with associative flow rule
 * - Isotropic hardening following a power law: \f$ H_p = k \cdot p^m \f$
 * - Convex Cutting Plane algorithm for return mapping
 * - Thermal expansion effects
 * - Consistent tangent modulus for implicit FE analysis
 *
 * **Constitutive Equations:**
 *
 * The yield function is defined as:
 * \f[
 * \Phi(\boldsymbol{\sigma}, p) = \sigma_{eq} - H_p(p) - \sigma_Y \leq 0
 * \f]
 * where:
 * - \f$ \sigma_{eq} = \sqrt{\frac{3}{2} \mathbf{s} : \mathbf{s}} \f$ is the von Mises equivalent stress
 * - \f$ \mathbf{s} \f$ is the deviatoric stress tensor
 * - \f$ H_p(p) = k \cdot p^m \f$ is the isotropic hardening function
 * - \f$ p = \int \sqrt{\frac{2}{3} \dot{\boldsymbol{\varepsilon}}^p : \dot{\boldsymbol{\varepsilon}}^p} \, dt \f$ is the accumulated plastic strain
 * - \f$ \sigma_Y \f$ is the initial yield stress
 *
 * The plastic flow rule (associative plasticity):
 * \f[
 * \dot{\boldsymbol{\varepsilon}}^p = \dot{\lambda} \frac{\partial \Phi}{\partial \boldsymbol{\sigma}} = \dot{\lambda} \frac{3}{2} \frac{\mathbf{s}}{\sigma_{eq}}
 * \f]
 *
 * Evolution of accumulated plastic strain:
 * \f[
 * \dot{p} = \sqrt{\frac{2}{3} \dot{\boldsymbol{\varepsilon}}^p : \dot{\boldsymbol{\varepsilon}}^p} = \dot{\lambda}
 * \f]
 *
 * **Convex Cutting Plane Algorithm:**
 *
 * The CCP algorithm solves the return mapping problem by reformulating it as a complementarity problem:
 * - Find \f$ \Delta p \geq 0 \f$ such that \f$ \Phi(\boldsymbol{\sigma}, p) \leq 0 \f$ and \f$ \Delta p \cdot \Phi = 0 \f$
 * - This is solved using the Fischer-Burmeister function for robust convergence
 * - The method provides a consistent tangent modulus for implicit finite element analysis
 *
 * **Material Parameters (props):**
 *
 * The material properties vector must contain 6 constants:
 * | Index | Symbol | Description | Units |
 * |-------|--------|-------------|-------|
 * | props[0] | \f$ E \f$ | Young's modulus | Stress |
 * | props[1] | \f$ \nu \f$ | Poisson's ratio | - |
 * | props[2] | \f$ \alpha \f$ | Isotropic thermal expansion coefficient | 1/Temperature |
 * | props[3] | \f$ \sigma_Y \f$ | Initial yield stress | Stress |
 * | props[4] | \f$ k \f$ | Hardening parameter | Stress |
 * | props[5] | \f$ m \f$ | Hardening exponent | - |
 *
 * **State Variables (statev):**
 *
 * The state variables vector contains 8 internal variables:
 * | Index | Symbol | Description | Units |
 * |-------|--------|-------------|-------|
 * | statev[0] | \f$ T_{init} \f$ | Initial temperature | Temperature |
 * | statev[1] | \f$ p \f$ | Accumulated plastic strain | - |
 * | statev[2] | \f$ \varepsilon^p_{11} \f$ | Plastic strain component 11 | - |
 * | statev[3] | \f$ \varepsilon^p_{22} \f$ | Plastic strain component 22 | - |
 * | statev[4] | \f$ \varepsilon^p_{33} \f$ | Plastic strain component 33 | - |
 * | statev[5] | \f$ \varepsilon^p_{12} \f$ | Plastic strain component 12 (×2 in Voigt) | - |
 * | statev[6] | \f$ \varepsilon^p_{13} \f$ | Plastic strain component 13 (×2 in Voigt) | - |
 * | statev[7] | \f$ \varepsilon^p_{23} \f$ | Plastic strain component 23 (×2 in Voigt) | - |
 *
 * @param umat_name Name of the constitutive model (EPICP)
 * @param Etot Total strain tensor at beginning of increment (Voigt notation: 6×1 vector)
 * @param DEtot Strain increment tensor (Voigt notation: 6×1 vector)
 * @param sigma Stress tensor (Voigt notation: 6×1 vector) [output]
 * @param Lt Consistent tangent modulus \f$ \mathbf{L}_t = \frac{\partial \boldsymbol{\sigma}}{\partial \boldsymbol{\varepsilon}} \f$ (6×6 matrix) [output]
 * @param L Elastic stiffness tensor (6×6 matrix) [output]
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
 * @param Wm_ir Irrecoverable work stored in hardening [output]
 * @param Wm_d Dissipated (plastic) work [output]
 * @param ndi Number of direct stress components (typically 3)
 * @param nshr Number of shear stress components (typically 3)
 * @param start Flag indicating first increment (true) or continuation (false)
 * @param tnew_dt Suggested new time step size for adaptive time stepping [output]
 *
 * @note Voigt notation convention: [11, 22, 33, 12, 13, 23] with engineering shear strains (γ = 2ε)
 * @note The consistent tangent modulus Lt ensures quadratic convergence in implicit Newton-Raphson schemes
 * @note The tangent modulus Lt is always computed
 *
 * @see Fischer_Burmeister_m() for the complementarity solver
 * @see L_iso() for isotropic elastic stiffness construction
 * @see eta_stress() for plastic flow direction computation
 *
 * @code
 * // Example usage:
 * vec Etot = {0.001, 0.0, 0.0, 0.0, 0.0, 0.0};
 * vec DEtot = {0.0001, 0.0, 0.0, 0.0, 0.0, 0.0};
 * vec sigma = zeros(6);
 * mat Lt = zeros(6,6);
 * mat L = zeros(6,6);
 * mat DR = eye(3,3);
 * vec props = {70000, 0.3, 1e-5, 200, 500, 0.2};
 * vec statev = zeros(8);
 *
 * umat_plasticity_iso_CCP("EPICP", Etot, DEtot, sigma, Lt, L, DR, 6, props, 8, statev,
 *                         20.0, 0.0, 0.0, 1.0, Wm, Wm_r, Wm_ir, Wm_d,
 *                         3, 3, true, tnew_dt);
 * @endcode
 */
void umat_plasticity_iso_CCP(const std::string &umat_name, const arma::vec &Etot, const arma::vec &DEtot, arma::vec &sigma, arma::mat &Lt, arma::mat &L, const arma::mat &DR, const int &nprops, const arma::vec &props, const int &nstatev, arma::vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, double &tnew_dt);
    

/** @} */ // end of umat_mechanical group

} //namespace simcoon
