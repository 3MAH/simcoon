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
 * @file Hill_isoh_Nfast.hpp
 * @brief Elastic-plastic material with Hill anisotropic yield and multiple isotropic hardening terms
 * @author Yves Chemisky
 * @version 1.0
 */

#pragma once
#include <string>
#include <armadillo>

namespace simcoon {

/** @addtogroup umat_mechanical
 *  @{
 */

/**
 * @brief Elastic-plastic model with Hill anisotropic criterion and N isotropic hardening terms
 *
 * @details This function implements the Hill anisotropic plasticity model with multiple
 * isotropic hardening terms for enhanced flexibility in capturing complex hardening behavior.
 * The model features:
 * - Hill (1948) quadratic anisotropic yield criterion
 * - Multiple Voce-type isotropic hardening terms: \f$ R(p) = \sum_{i=1}^{N} Q_i (1 - e^{-b_i p}) \f$
 * - Associative flow rule
 * - Convex Cutting Plane algorithm for return mapping
 * - Consistent tangent modulus for implicit FE analysis
 *
 * **Constitutive Equations:**
 *
 * The Hill anisotropic yield function is defined as:
 * \f[
 * \Phi(\boldsymbol{\sigma}, p) = \sigma_{eq}^{Hill} - R(p) - \sigma_Y \leq 0
 * \f]
 * where:
 * \f[
 * \sigma_{eq}^{Hill} = \sqrt{F(\sigma_{22} - \sigma_{33})^2 + G(\sigma_{33} - \sigma_{11})^2 + H(\sigma_{11} - \sigma_{22})^2 + 2L\sigma_{23}^2 + 2M\sigma_{13}^2 + 2N\sigma_{12}^2}
 * \f]
 *
 * **Multiple Isotropic Hardening (Voce Law):**
 *
 * The isotropic hardening is the sum of N exponential terms:
 * \f[
 * R(p) = \sum_{i=1}^{N_{iso}} Q_i \left( 1 - e^{-b_i p} \right)
 * \f]
 * where:
 * - \f$ Q_i \f$ is the saturation stress of the i-th hardening term
 * - \f$ b_i \f$ is the hardening rate parameter of the i-th term
 * - \f$ p \f$ is the accumulated plastic strain
 *
 * **Hardening Modulus:**
 *
 * The derivative of the isotropic hardening function:
 * \f[
 * H_{iso} = \frac{dR}{dp} = \sum_{i=1}^{N_{iso}} Q_i b_i e^{-b_i p}
 * \f]
 *
 * **Advantages of Multiple Hardening Terms:**
 *
 * - Better fit to experimental stress-strain curves
 * - Captures both rapid initial hardening and gradual saturation
 * - Different terms can represent different physical mechanisms:
 *   - Fast term (high b): dislocation pile-up, initial work hardening
 *   - Slow term (low b): long-range backstress development, texture evolution
 *
 * **Material Parameters (props):**
 *
 * | Index | Symbol | Description | Units |
 * |-------|--------|-------------|-------|
 * | props[0] | \f$ E \f$ | Young's modulus | Stress |
 * | props[1] | \f$ \nu \f$ | Poisson's ratio | - |
 * | props[2] | \f$ \alpha \f$ | Thermal expansion coefficient | 1/Temperature |
 * | props[3] | \f$ \sigma_Y \f$ | Initial yield stress | Stress |
 * | props[4] | \f$ N_{iso} \f$ | Number of isotropic hardening terms | - |
 * | props[5] | \f$ F \f$ | Hill parameter F | 1/Stress² |
 * | props[6] | \f$ G \f$ | Hill parameter G | 1/Stress² |
 * | props[7] | \f$ H \f$ | Hill parameter H | 1/Stress² |
 * | props[8] | \f$ L \f$ | Hill parameter L | 1/Stress² |
 * | props[9] | \f$ M \f$ | Hill parameter M | 1/Stress² |
 * | props[10] | \f$ N \f$ | Hill parameter N | 1/Stress² |
 * | props[11+2j] | \f$ Q_j \f$ | Saturation stress of j-th term | Stress |
 * | props[12+2j] | \f$ b_j \f$ | Hardening rate of j-th term | 1/Strain |
 *
 * **State Variables (statev):**
 *
 * Total state variables required: \f$ n_{statev} = 8 \f$
 *
 * | Index | Symbol | Description | Units |
 * |-------|--------|-------------|-------|
 * | statev[0] | \f$ T_{init} \f$ | Initial temperature | Temperature |
 * | statev[1] | \f$ p \f$ | Accumulated plastic strain | Strain |
 * | statev[2] | \f$ \varepsilon^p_{11} \f$ | Plastic strain component 11 | Strain |
 * | statev[3] | \f$ \varepsilon^p_{22} \f$ | Plastic strain component 22 | Strain |
 * | statev[4] | \f$ \varepsilon^p_{33} \f$ | Plastic strain component 33 | Strain |
 * | statev[5] | \f$ \varepsilon^p_{12} \f$ | Plastic strain component 12 | Strain |
 * | statev[6] | \f$ \varepsilon^p_{13} \f$ | Plastic strain component 13 | Strain |
 * | statev[7] | \f$ \varepsilon^p_{23} \f$ | Plastic strain component 23 | Strain |
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
 * @param Wm_ir Irrecoverable work stored in plastic deformation [output]
 * @param Wm_d Dissipated (plastic) work [output]
 * @param ndi Number of direct stress components (typically 3)
 * @param nshr Number of shear stress components (typically 3)
 * @param start Flag indicating first increment (true) or continuation (false)
 * @param tnew_dt Suggested new time step size for adaptive time stepping [output]
 *
 * @note Typically 2-3 hardening terms are sufficient for most metals
 * @note For single power-law hardening, use umat_plasticity_hill_isoh_CCP instead
 * @note The Hill parameters must satisfy physical constraints (positive definite)
 *
 * @see umat_plasticity_hill_isoh_CCP() for single power-law hardening version
 * @see umat_hill_chaboche_CCP() for combined kinematic and isotropic hardening
 *
 * **References:**
 * - Hill, R. (1948). "A theory of the yielding and plastic flow of anisotropic metals." *Proc. Roy. Soc. A*, 193(1033), 281-297.
 * - Voce, E. (1948). "The relationship between stress and strain for homogeneous deformation." *J. Inst. Met.*, 74, 537-562.
 * - Barlat, F., et al. (2003). "Plane stress yield function for aluminum alloy sheets." *Int. J. Plasticity*, 19(9), 1297-1319.
 */
void umat_plasticity_hill_isoh_CCP_N(const std::string &umat_name, const arma::vec &Etot, const arma::vec &DEtot, arma::vec &sigma, arma::mat &Lt, arma::mat &L, const arma::mat &DR, const int &nprops, const arma::vec &props, const int &nstatev, arma::vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, double &tnew_dt);

/** @} */ // end of umat_mechanical group

} //namespace simcoon
