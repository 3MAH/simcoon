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

///@file elastic_isotropic.hpp
///@brief User subroutine for Isotropic elastic materials in 3D case
///@version 1.0

#pragma once

#include <string>
#include <armadillo>

namespace simcoon{

/**
 * @file elastic_isotropic.hpp
 * @brief Linear elastic isotropic material model for small strain analysis
 * @author Yves Chemisky
 * @version 1.0
 */

/** @addtogroup umat_mechanical
 *  @{
 */

/**
 * @brief Linear elastic isotropic constitutive model with thermal expansion
 *
 * @details This function implements the fundamental linear elastic isotropic material model
 * based on Hooke's law, suitable for small strain analysis of materials exhibiting
 * reversible elastic deformation.
 *
 * **Constitutive Equations:**
 *
 * The stress-strain relationship follows the generalized Hooke's law:
 * \f[
 * \boldsymbol{\sigma} = \mathbf{L} : (\boldsymbol{\varepsilon} - \boldsymbol{\varepsilon}^{th})
 * \f]
 *
 * where:
 * - \f$ \boldsymbol{\sigma} \f$ is the stress tensor (Voigt notation)
 * - \f$ \mathbf{L} \f$ is the fourth-order elastic stiffness tensor
 * - \f$ \boldsymbol{\varepsilon} \f$ is the total strain tensor
 * - \f$ \boldsymbol{\varepsilon}^{th} = \alpha (T - T_0) \mathbf{I} \f$ is the thermal strain
 * - \f$ \alpha \f$ is the coefficient of thermal expansion
 * - \f$ T \f$ is the current temperature
 * - \f$ T_0 \f$ is the reference temperature
 *
 * **Elastic Stiffness Tensor:**
 *
 * For an isotropic material, the stiffness tensor is constructed from two independent constants:
 * \f[
 * \mathbf{L} = 3K \mathbf{I}_{vol} + 2\mu \mathbf{I}_{dev}
 * \f]
 *
 * where:
 * - \f$ K = \frac{E}{3(1-2\nu)} \f$ is the bulk modulus
 * - \f$ \mu = \frac{E}{2(1+\nu)} \f$ is the shear modulus
 * - \f$ \mathbf{I}_{vol} \f$ is the volumetric projection tensor
 * - \f$ \mathbf{I}_{dev} \f$ is the deviatoric projection tensor
 *
 * In component form (Voigt notation):
 * \f[
 * \mathbf{L} = \frac{E}{(1+\nu)(1-2\nu)}
 * \begin{bmatrix}
 * 1-\nu & \nu & \nu & 0 & 0 & 0 \\
 * \nu & 1-\nu & \nu & 0 & 0 & 0 \\
 * \nu & \nu & 1-\nu & 0 & 0 & 0 \\
 * 0 & 0 & 0 & \frac{1-2\nu}{2} & 0 & 0 \\
 * 0 & 0 & 0 & 0 & \frac{1-2\nu}{2} & 0 \\
 * 0 & 0 & 0 & 0 & 0 & \frac{1-2\nu}{2}
 * \end{bmatrix}
 * \f]
 *
 * **Stress Update:**
 *
 * For an increment of strain \f$ \Delta \boldsymbol{\varepsilon} \f$:
 * \f[
 * \boldsymbol{\sigma}_{n+1} = \boldsymbol{\sigma}_n + \mathbf{L} : \left( \Delta \boldsymbol{\varepsilon} - \alpha \Delta T \mathbf{I} \right)
 * \f]
 *
 * **Tangent Modulus:**
 *
 * For linear elasticity, the tangent modulus is constant and equal to the elastic stiffness:
 * \f[
 * \mathbf{L}_t = \mathbf{L}
 * \f]
 *
 * This ensures optimal convergence in Newton-Raphson iterations.
 *
 * **Material Parameters (props):**
 *
 * | Index | Symbol | Description | Units | Typical Range |
 * |-------|--------|-------------|-------|---------------|
 * | props[0] | \f$ E \f$ | Young's modulus | Stress | 1-500 GPa |
 * | props[1] | \f$ \nu \f$ | Poisson's ratio | - | 0.0-0.5 |
 * | props[2] | \f$ \alpha \f$ | Coefficient of thermal expansion | 1/Temperature | 1e-6 to 1e-4 /K |
 *
 * **Constraint on Poisson's ratio:**
 * - For 3D: \f$ -1 < \nu < 0.5 \f$
 * - \f$ \nu = 0.5 \f$ implies incompressibility (infinite bulk modulus)
 * - \f$ \nu = 0 \f$ implies no lateral strain under uniaxial loading
 *
 * **State Variables (statev):**
 *
 * | Index | Symbol | Description | Units |
 * |-------|--------|-------------|-------|
 * | statev[0] | \f$ T_{init} \f$ | Initial/reference temperature | Temperature |
 *
 * Total state variables required: \f$ n_{statev} = 1 \f$
 *
 * @param umat_name Name of the constitutive model (ELISO)
 * @param Etot Total strain tensor at beginning of increment (Voigt notation: 6×1 vector)
 * @param DEtot Strain increment tensor (Voigt notation: 6×1 vector)
 * @param sigma Stress tensor (Voigt notation: 6×1 vector) [output]
 * @param Lt Tangent modulus \f$ \mathbf{L}_t = \mathbf{L} \f$ (6×6 matrix) [output]
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
 * @param Wm_ir Irrecoverable work (0 for elastic model) [output]
 * @param Wm_d Dissipated work (0 for elastic model) [output]
 * @param ndi Number of direct stress components (typically 3)
 * @param nshr Number of shear stress components (typically 3)
 * @param start Flag indicating first increment (true) or continuation (false)
 * @param tnew_dt Suggested new time step size for adaptive time stepping [output]
 *
 * @note This model is restricted to small strains (< 1%)
 * @note For large strains, use hyperelastic models (Neo-Hookean, Mooney-Rivlin)
 * @note The model is path-independent (no history dependence)
 * @note No energy dissipation occurs (Wm_ir = Wm_d = 0, Wm = Wm_r)
 * @note For plane stress/strain, set ndi=2 or ndi=1 accordingly
 *
 * @see L_iso() for elastic stiffness tensor construction
 * @see Ivol() for volumetric projection tensor
 * @see Idev() for deviatoric projection tensor
 * @see Ith() for thermal expansion vector
 *
 * @code
 * // Example usage: Steel
 * vec props = {210000, 0.3, 1.2e-5};  // E=210 GPa, nu=0.3, alpha=12e-6 /K
 * vec statev = zeros(1);
 * statev(0) = 20.0;  // Reference temperature 20°C
 *
 * vec Etot = {0.001, 0.0, 0.0, 0.0, 0.0, 0.0};  // 0.1% strain in direction 1
 * vec DEtot = {0.0001, 0.0, 0.0, 0.0, 0.0, 0.0};
 * vec sigma = zeros(6);
 * mat Lt = zeros(6,6);
 * mat L = zeros(6,6);
 * mat DR = eye(3,3);
 *
 * umat_elasticity_iso("ELISO", Etot, DEtot, sigma, Lt, L, DR,
 *                     3, props, 1, statev, 25.0, 5.0, 0.0, 1.0,
 *                     Wm, Wm_r, Wm_ir, Wm_d, 3, 3, false, tnew_dt);
 *
 * // Expected stress: sigma(0) ≈ E*Etot(0)*(1-nu)/((1+nu)(1-2nu)) - E*alpha*DT/(1-2nu)
 * @endcode
 *
 * **References:**
 * - Timoshenko, S. P., & Goodier, J. N. (1970). *Theory of Elasticity* (3rd ed.). McGraw-Hill.
 * - Bower, A. F. (2009). *Applied Mechanics of Solids*. CRC Press.
 * - Gurtin, M. E. (1981). *An Introduction to Continuum Mechanics*. Academic Press.
 */
void umat_elasticity_iso(const std::string &umat_name, const arma::vec &Etot, const arma::vec &DEtot, arma::vec &sigma, arma::mat &Lt, arma::mat &L, const arma::mat &DR, const int &nprops, const arma::vec &props, const int &nstatev, arma::vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, double &tnew_dt);
                            

/** @} */ // end of umat_mechanical group

} //namespace simcoon
