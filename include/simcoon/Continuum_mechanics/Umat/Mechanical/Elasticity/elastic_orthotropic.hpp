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

///@file elastic_orthotropic.hpp
///@brief User subroutine for orthotropic elastic materials in 3D case
///@version 1.0

#pragma once

#include <armadillo>

namespace simcoon{

/**
 * @file elastic_orthotropic.hpp
 * @brief Linear elastic orthotropic material model
 * @author Yves Chemisky
 * @version 1.0
 *
 * This file implements a linear elastic orthotropic material model with
 * thermal expansion. Orthotropic materials have three mutually perpendicular
 * planes of symmetry, resulting in 9 independent elastic constants.
 */

/** @addtogroup umat_mechanical
 *  @{
 */

/**
 * @brief Linear elastic orthotropic constitutive model with thermal expansion
 *
 * @details This function implements a linear elastic orthotropic material model
 * suitable for materials with three perpendicular planes of symmetry, such as:
 * - Wood and timber
 * - Rolled metals
 * - Unidirectional composites (approximation)
 * - Crystals with orthorhombic symmetry
 *
 * **Constitutive Equations:**
 *
 * The stress-strain relationship in the material principal axes:
 * \f[
 * \begin{bmatrix} \varepsilon_{11} \\ \varepsilon_{22} \\ \varepsilon_{33} \\ \gamma_{12} \\ \gamma_{13} \\ \gamma_{23} \end{bmatrix}
 * = \begin{bmatrix}
 * 1/E_1 & -\nu_{12}/E_1 & -\nu_{13}/E_1 & 0 & 0 & 0 \\
 * -\nu_{12}/E_1 & 1/E_2 & -\nu_{23}/E_2 & 0 & 0 & 0 \\
 * -\nu_{13}/E_1 & -\nu_{23}/E_2 & 1/E_3 & 0 & 0 & 0 \\
 * 0 & 0 & 0 & 1/G_{12} & 0 & 0 \\
 * 0 & 0 & 0 & 0 & 1/G_{13} & 0 \\
 * 0 & 0 & 0 & 0 & 0 & 1/G_{23}
 * \end{bmatrix}
 * \begin{bmatrix} \sigma_{11} \\ \sigma_{22} \\ \sigma_{33} \\ \tau_{12} \\ \tau_{13} \\ \tau_{23} \end{bmatrix}
 * + \begin{bmatrix} \alpha_1 \\ \alpha_2 \\ \alpha_3 \\ 0 \\ 0 \\ 0 \end{bmatrix} \Delta T
 * \f]
 *
 * **Symmetry Requirements:**
 *
 * The compliance matrix must be symmetric, which imposes:
 * \f[
 * \frac{\nu_{ij}}{E_i} = \frac{\nu_{ji}}{E_j}
 * \f]
 *
 * **Material Parameters (props):**
 *
 * | Index | Symbol | Description | Units |
 * |-------|--------|-------------|-------|
 * | props[0] | \f$ E_1 \f$ | Young's modulus in direction 1 | Stress |
 * | props[1] | \f$ E_2 \f$ | Young's modulus in direction 2 | Stress |
 * | props[2] | \f$ E_3 \f$ | Young's modulus in direction 3 | Stress |
 * | props[3] | \f$ \nu_{12} \f$ | Poisson's ratio (strain in 2 due to stress in 1) | - |
 * | props[4] | \f$ \nu_{13} \f$ | Poisson's ratio (strain in 3 due to stress in 1) | - |
 * | props[5] | \f$ \nu_{23} \f$ | Poisson's ratio (strain in 3 due to stress in 2) | - |
 * | props[6] | \f$ G_{12} \f$ | Shear modulus in 1-2 plane | Stress |
 * | props[7] | \f$ G_{13} \f$ | Shear modulus in 1-3 plane | Stress |
 * | props[8] | \f$ G_{23} \f$ | Shear modulus in 2-3 plane | Stress |
 * | props[9] | \f$ \alpha_1 \f$ | CTE in direction 1 | 1/Temperature |
 * | props[10] | \f$ \alpha_2 \f$ | CTE in direction 2 | 1/Temperature |
 * | props[11] | \f$ \alpha_3 \f$ | CTE in direction 3 | 1/Temperature |
 *
 * **State Variables (statev):**
 *
 * | Index | Symbol | Description | Units |
 * |-------|--------|-------------|-------|
 * | statev[0] | \f$ T_{init} \f$ | Initial/reference temperature | Temperature |
 *
 * @param Etot Total strain tensor at beginning of increment (Voigt notation: 6×1)
 * @param DEtot Strain increment tensor (Voigt notation: 6×1)
 * @param sigma Stress tensor [output] (Voigt notation: 6×1)
 * @param Lt Tangent modulus (= L for linear elasticity) [output] (6×6)
 * @param L Elastic stiffness tensor [output] (6×6)
 * @param sigma_in Internal stress for explicit solvers [output] (6×1)
 * @param DR Rotation increment matrix (3×3)
 * @param nprops Number of material properties (12)
 * @param props Material properties vector
 * @param nstatev Number of state variables (1)
 * @param statev State variables vector [input/output]
 * @param T Temperature at beginning of increment
 * @param DT Temperature increment
 * @param Time Time at beginning of increment
 * @param DTime Time increment
 * @param Wm Total mechanical work [output]
 * @param Wm_r Recoverable (elastic) work [output]
 * @param Wm_ir Irrecoverable work (0 for elastic) [output]
 * @param Wm_d Dissipated work (0 for elastic) [output]
 * @param ndi Number of direct stress components
 * @param nshr Number of shear stress components
 * @param start Flag indicating first increment
 * @param solver_type Solver type (0=implicit, 1=explicit)
 * @param tnew_dt Suggested new time step [output]
 *
 * @note Material axes must be aligned with global axes. For rotated materials,
 *       use local orientation definitions in the FE software.
 * @note The stiffness matrix must be positive definite for physical stability
 *
 * @see L_ortho() for orthotropic stiffness tensor construction
 */
void umat_elasticity_ortho(const arma::vec &Etot, const arma::vec &DEtot, arma::vec &sigma, arma::mat &Lt, arma::mat &L, arma::vec &sigma_in, const arma::mat &DR, const int &nprops, const arma::vec &props, const int &nstatev, arma::vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, const int &solver_type, double &tnew_dt);


/** @} */ // end of umat_mechanical group

} //namespace simcoon
