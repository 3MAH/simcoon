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

///@file elastic_transverse_isotropic.hpp
///@brief User subroutine for transversely isotropic elastic materials in 3D case
///@version 1.0

#pragma once

#include <armadillo>

namespace simcoon{

/**
 * @file elastic_transverse_isotropic.hpp
 * @brief Linear elastic transversely isotropic material model
 * @author Yves Chemisky
 * @version 1.0
 *
 * This file implements a linear elastic transversely isotropic material model
 * with thermal expansion. Transversely isotropic materials have one axis of
 * symmetry (fiber direction) and isotropic properties in the transverse plane.
 */

/** @addtogroup umat_mechanical
 *  @{
 */

/**
 * @brief Linear elastic transversely isotropic constitutive model with thermal expansion
 *
 * @details This function implements a linear elastic transversely isotropic material
 * model suitable for materials with one preferred direction, such as:
 * - Unidirectional fiber composites
 * - Hexagonal crystals
 * - Biological tissues (muscles, tendons)
 * - Extruded or drawn materials
 *
 * **Material Symmetry:**
 *
 * The material has rotational symmetry about axis 1 (fiber direction).
 * Properties are isotropic in the 2-3 plane (transverse plane).
 *
 * **Constitutive Equations:**
 *
 * The compliance matrix in the material principal axes:
 * \f[
 * \begin{bmatrix} \varepsilon_{11} \\ \varepsilon_{22} \\ \varepsilon_{33} \\ \gamma_{12} \\ \gamma_{13} \\ \gamma_{23} \end{bmatrix}
 * = \begin{bmatrix}
 * 1/E_L & -\nu_{LT}/E_L & -\nu_{LT}/E_L & 0 & 0 & 0 \\
 * -\nu_{LT}/E_L & 1/E_T & -\nu_{TT}/E_T & 0 & 0 & 0 \\
 * -\nu_{LT}/E_L & -\nu_{TT}/E_T & 1/E_T & 0 & 0 & 0 \\
 * 0 & 0 & 0 & 1/G_{LT} & 0 & 0 \\
 * 0 & 0 & 0 & 0 & 1/G_{LT} & 0 \\
 * 0 & 0 & 0 & 0 & 0 & 2(1+\nu_{TT})/E_T
 * \end{bmatrix}
 * \begin{bmatrix} \sigma_{11} \\ \sigma_{22} \\ \sigma_{33} \\ \tau_{12} \\ \tau_{13} \\ \tau_{23} \end{bmatrix}
 * \f]
 *
 * Note: The transverse shear modulus is not independent:
 * \f[
 * G_{TT} = \frac{E_T}{2(1 + \nu_{TT})}
 * \f]
 *
 * **Material Parameters (props):**
 *
 * | Index | Symbol | Description | Units |
 * |-------|--------|-------------|-------|
 * | props[0] | \f$ E_L \f$ | Longitudinal Young's modulus (fiber direction) | Stress |
 * | props[1] | \f$ E_T \f$ | Transverse Young's modulus | Stress |
 * | props[2] | \f$ \nu_{LT} \f$ | Poisson's ratio (transverse strain due to longitudinal stress) | - |
 * | props[3] | \f$ \nu_{TT} \f$ | Transverse Poisson's ratio (in the isotropic plane) | - |
 * | props[4] | \f$ G_{LT} \f$ | Longitudinal shear modulus | Stress |
 * | props[5] | \f$ \alpha_L \f$ | Longitudinal CTE | 1/Temperature |
 * | props[6] | \f$ \alpha_T \f$ | Transverse CTE | 1/Temperature |
 *
 * **Typical Values (UD Carbon/Epoxy):**
 * - \f$ E_L \f$ ≈ 140 GPa
 * - \f$ E_T \f$ ≈ 10 GPa
 * - \f$ \nu_{LT} \f$ ≈ 0.3
 * - \f$ \nu_{TT} \f$ ≈ 0.4
 * - \f$ G_{LT} \f$ ≈ 5 GPa
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
 * @param nprops Number of material properties (7)
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
 * @note The fiber direction (axis 1) must be aligned with the global x-axis.
 *       For rotated fibers, use local orientation definitions.
 * @note The 2-3 plane is the plane of isotropy (transverse plane)
 *
 * @see L_isotrans() for transversely isotropic stiffness tensor construction
 */
void umat_elasticity_trans_iso(const arma::vec &Etot, const arma::vec &DEtot, arma::vec &sigma, arma::mat &Lt, arma::mat &L, arma::vec &sigma_in, const arma::mat &DR, const int &nprops, const arma::vec &props, const int &nstatev, arma::vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, const int &solver_type, double &tnew_dt);


/** @} */ // end of umat_mechanical group

} //namespace simcoon
