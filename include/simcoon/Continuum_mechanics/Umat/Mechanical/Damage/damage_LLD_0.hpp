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
 * @file damage_LLD_0.hpp
 * @brief Ladevèze-Le Dantec anisotropic damage model for composite materials
 * @author Yves Chemisky
 * @version 1.0
 */

#pragma once

#include <string>
#include <iostream>
#include <armadillo>

namespace simcoon{

/** @addtogroup umat_mechanical
 *  @{
 */

/**
 * @brief Ladevèze-Le Dantec (LLD) anisotropic damage model with coupled plasticity for composite materials
 *
 * @details This function implements the Ladevèze-Le Dantec damage model specifically designed for
 * unidirectional fiber-reinforced composite materials. The model features:
 * - Transversely isotropic elastic behavior
 * - Anisotropic damage in matrix-dominated modes (transverse tension, in-plane shear)
 * - Coupled transverse-shear plasticity
 * - Thermodynamically consistent damage evolution based on energy release rate
 * - Separate damage variables for different failure modes
 * - No fiber failure (matrix-dominated damage only)
 *
 * **Material Symmetry:**
 *
 * The model assumes transverse isotropy with the fiber direction as the axis of symmetry:
 * - Direction 1: Fiber direction (longitudinal, no damage)
 * - Directions 2-3: Transverse plane (isotropic, subject to damage)
 *
 * **Damage Variables:**
 *
 * Two independent scalar damage variables characterize the material state:
 * - \f$ d_{22} \in [0,1] \f$: Transverse damage (matrix cracking perpendicular to fibers)
 * - \f$ d_{12} \in [0,1] \f$: Shear damage (matrix/interface damage in fiber-matrix plane)
 *
 * where \f$ d = 0 \f$ is undamaged and \f$ d = 1 \f$ is fully damaged.
 *
 * **Effective Stress Concept:**
 *
 * The effective stress acting on the undamaged material configuration is:
 * \f[
 * \tilde{\boldsymbol{\sigma}} = \mathbf{M}^{-1} : \boldsymbol{\sigma}
 * \f]
 * where \f$ \mathbf{M} \f$ is the damage effect tensor:
 * \f[
 * \mathbf{M} = \text{diag}(1, 1-d_{22}, 1-d_{22}, 1-d_{12}, 1-d_{12}, 1)
 * \f]
 *
 * **Damaged Elastic Stiffness:**
 *
 * The elastic stiffness degrades with damage:
 * \f[
 * \mathbf{L}(d_{22}, d_{12}) = \mathbf{L}_0 : \mathbf{M}
 * \f]
 * where \f$ \mathbf{L}_0 \f$ is the undamaged transversely isotropic stiffness.
 *
 * **Specific moduli degradation:**
 * - \f$ E_2 = E_{2,0} (1 - d_{22}) \f$ (transverse Young's modulus)
 * - \f$ E_3 = E_{3,0} (1 - d_{22}) \f$ (out-of-plane Young's modulus)
 * - \f$ G_{12} = G_{12,0} (1 - d_{12}) \f$ (in-plane shear modulus)
 * - \f$ G_{13} = G_{13,0} (1 - d_{12}) \f$ (out-of-plane shear modulus)
 * - \f$ E_1 \f$ (fiber direction) remains constant (no fiber damage)
 *
 * **Thermodynamic Forces:**
 *
 * The energy release rates driving damage evolution are:
 * \f[
 * Y_{22} = \frac{1}{2} \boldsymbol{\varepsilon}^e : \frac{\partial \mathbf{L}}{\partial d_{22}} : \boldsymbol{\varepsilon}^e
 * \f]
 * \f[
 * Y_{12} = \frac{1}{2} \boldsymbol{\varepsilon}^e : \frac{\partial \mathbf{L}}{\partial d_{12}} : \boldsymbol{\varepsilon}^e
 * \f]
 *
 * **Damage Evolution Laws:**
 *
 * **Transverse Damage (d₂₂):**
 * \f[
 * \dot{d}_{22} = \begin{cases}
 * 0 & \text{if } Y_{22} < Y_{22,0} \\
 * \left( \frac{Y_{22} - Y_{22,c}}{Y_{22,u} - Y_{22,c}} \right)^b & \text{if } Y_{22,0} \leq Y_{22} < Y_{22,u} \\
 * \infty & \text{if } Y_{22} \geq Y_{22,u}
 * \end{cases}
 * \f]
 * where:
 * - \f$ Y_{22,0} \f$ is the damage initiation threshold
 * - \f$ Y_{22,c} \f$ is the characteristic energy release rate
 * - \f$ Y_{22,u} \f$ is the ultimate energy release rate (failure)
 * - \f$ b \f$ is the damage evolution exponent
 *
 * **Shear Damage (d₁₂):**
 * \f[
 * \dot{d}_{12} = \begin{cases}
 * 0 & \text{if } Y_{12} < Y_{12,0} \\
 * \frac{Y_{12} - Y_{12,c}}{Y_{12,c}} & \text{if } Y_{12} \geq Y_{12,0}
 * \end{cases}
 * \f]
 * where:
 * - \f$ Y_{12,0} \f$ is the shear damage initiation threshold
 * - \f$ Y_{12,c} \f$ is the characteristic shear energy release rate
 *
 * **Coupled Transverse-Shear Plasticity:**
 *
 * A Hill-type yield criterion couples transverse and shear stresses:
 * \f[
 * \Phi = \sqrt{\left( \frac{\sigma_{22}}{A_{ts}} \right)^2 + \sigma_{12}^2} - \sigma_{ts,0} - \alpha_{ts} p_{ts} - \beta_{ts} p_{ts}^2 \leq 0
 * \f]
 * where:
 * - \f$ A_{ts} \f$ is the transverse-shear coupling parameter
 * - \f$ \sigma_{ts,0} \f$ is the initial yield stress
 * - \f$ \alpha_{ts}, \beta_{ts} \f$ are hardening parameters
 * - \f$ p_{ts} \f$ is the accumulated plastic strain
 *
 * **Plastic Flow Rule:**
 *
 * Associative flow in the transverse-shear plane:
 * \f[
 * \dot{\boldsymbol{\varepsilon}}^p = \dot{p}_{ts} \frac{\partial \Phi}{\partial \boldsymbol{\sigma}}
 * \f]
 *
 * **Material Parameters (props):**
 *
 * | Index | Symbol | Description | Units | Typical Range (CFRP) |
 * |-------|--------|-------------|-------|----------------------|
 * | props[0] | axis | Fiber orientation axis (1, 2, or 3) | - | 1 |
 * | props[1] | \f$ E_L \f$ | Longitudinal Young's modulus | Stress | 100-200 GPa |
 * | props[2] | \f$ E_T \f$ | Transverse Young's modulus | Stress | 5-15 GPa |
 * | props[3] | \f$ \nu_{TL} \f$ | Major Poisson's ratio | - | 0.25-0.35 |
 * | props[4] | \f$ \nu_{TT} \f$ | Transverse Poisson's ratio | - | 0.35-0.45 |
 * | props[5] | \f$ G_{LT} \f$ | In-plane shear modulus | Stress | 3-8 GPa |
 * | props[6] | \f$ \alpha_L \f$ | Longitudinal CTE | 1/Temperature | -0.5 to 0 e-6 /K |
 * | props[7] | \f$ \alpha_T \f$ | Transverse CTE | 1/Temperature | 20-40 e-6 /K |
 * | props[8] | \f$ Y_{12,0} \f$ | Shear damage initiation threshold | Energy/Volume | 0.05-0.2 MPa |
 * | props[9] | \f$ Y_{12,c} \f$ | Characteristic shear energy release | Energy/Volume | 0.1-0.5 MPa |
 * | props[10] | \f$ Y_{22,0} \f$ | Transverse damage initiation | Energy/Volume | 0.1-0.3 MPa |
 * | props[11] | \f$ Y_{22,c} \f$ | Characteristic transverse energy | Energy/Volume | 0.2-0.6 MPa |
 * | props[12] | \f$ Y_{22,u} \f$ | Ultimate transverse energy | Energy/Volume | 0.5-2.0 MPa |
 * | props[13] | \f$ b \f$ | Transverse damage exponent | - | 1-5 |
 * | props[14] | \f$ A_{ts} \f$ | Transverse-shear coupling | - | 1-3 |
 * | props[15] | \f$ \sigma_{ts,0} \f$ | Initial yield stress | Stress | 30-80 MPa |
 * | props[16] | \f$ \alpha_{ts} \f$ | Linear hardening parameter | Stress | 0-500 MPa |
 * | props[17] | \f$ \beta_{ts} \f$ | Quadratic hardening parameter | Stress | 0-5000 MPa |
 *
 * **State Variables (statev):**
 *
 * Total state variables required: \f$ n_{statev} = 10 \f$
 *
 * | Index | Symbol | Description | Units |
 * |-------|--------|-------------|-------|
 * | statev[0] | \f$ T_{init} \f$ | Initial/reference temperature | Temperature |
 * | statev[1] | \f$ d_{22} \f$ | Transverse damage variable | - |
 * | statev[2] | \f$ d_{12} \f$ | Shear damage variable | - |
 * | statev[3] | \f$ p_{ts} \f$ | Accumulated plastic strain | Strain |
 * | statev[4] | \f$ \varepsilon^p_{11} \f$ | Plastic strain component 11 | Strain |
 * | statev[5] | \f$ \varepsilon^p_{22} \f$ | Plastic strain component 22 | Strain |
 * | statev[6] | \f$ \varepsilon^p_{33} \f$ | Plastic strain component 33 | Strain |
 * | statev[7] | \f$ \varepsilon^p_{12} \f$ | Plastic strain component 12 (engineering) | Strain |
 * | statev[8] | \f$ \varepsilon^p_{13} \f$ | Plastic strain component 13 (engineering) | Strain |
 * | statev[9] | \f$ \varepsilon^p_{23} \f$ | Plastic strain component 23 (engineering) | Strain |
 *
 * @param Etot Total strain tensor at beginning of increment (Voigt notation: 6×1 vector)
 * @param DEtot Strain increment tensor (Voigt notation: 6×1 vector)
 * @param sigma Stress tensor (Voigt notation: 6×1 vector) [output]
 * @param Lt Consistent tangent modulus (6×6 matrix) [output]
 * @param L Damaged elastic stiffness tensor (6×6 matrix) [output]
 * @param sigma_in Internal stress contribution for explicit solvers (6×1 vector) [output]
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
 * @param Wm_ir Irrecoverable work stored in damage [output]
 * @param Wm_d Dissipated work (damage + plasticity) [output]
 * @param ndi Number of direct stress components (typically 3)
 * @param nshr Number of shear stress components (typically 3)
 * @param start Flag indicating first increment (true) or continuation (false)
 * @param solver_type Solver type: 0=implicit, 1=explicit, 2=dynamic implicit
 * @param tnew_dt Suggested new time step size for adaptive time stepping [output]
 *
 * @note The LLD model is specifically designed for unidirectional fiber composites
 * @note Damage is irreversible and monotonically increasing
 * @note The model does NOT account for fiber failure (compression/tension in direction 1)
 * @note Suitable for matrix-dominated failure modes: transverse cracking, delamination
 * @note Parameter identification requires multiple test configurations:
 * @note - Transverse tension for d₂₂ parameters
 * @note - In-plane shear for d₁₂ parameters
 * @note - Off-axis tests for plasticity coupling
 * @note Convergence requires small load steps once damage initiates
 * @note Material axes must be properly oriented relative to global coordinates
 *
 * @see L_isotrans() for transversely isotropic stiffness construction
 * @see Lagrange_exp() for penalty function (damage bounds enforcement)
 * @see rotate_strain() for strain tensor rotation
 * @see rotate_stress() for stress tensor rotation
 *
 * @code
 * // Example usage: Carbon/Epoxy unidirectional composite (T300/914)
 * vec props(18);
 * props(0) = 1;           // axis = 1 (fibers along direction 1)
 * props(1) = 138000;      // EL = 138 GPa (fiber-dominated)
 * props(2) = 11000;       // ET = 11 GPa (matrix-dominated)
 * props(3) = 0.28;        // nuTL = 0.28
 * props(4) = 0.40;        // nuTT = 0.40
 * props(5) = 5500;        // GLT = 5.5 GPa
 * props(6) = -0.3e-6;     // alphaL = -0.3e-6 /K (negative for carbon fibers)
 * props(7) = 30e-6;       // alphaT = 30e-6 /K (resin dominated)
 * props(8) = 0.10;        // Y_12_0 = 0.10 MPa (shear damage threshold)
 * props(9) = 0.30;        // Y_12_c = 0.30 MPa
 * props(10) = 0.15;       // Y_22_0 = 0.15 MPa (transverse damage threshold)
 * props(11) = 0.40;       // Y_22_c = 0.40 MPa
 * props(12) = 1.20;       // Y_22_u = 1.20 MPa (ultimate failure)
 * props(13) = 2.5;        // b = 2.5
 * props(14) = 2.0;        // A_ts = 2.0
 * props(15) = 50;         // sigma_ts_0 = 50 MPa
 * props(16) = 300;        // alpha_ts = 300 MPa
 * props(17) = 1000;       // beta_ts = 1000 MPa
 *
 * vec statev = zeros(10);
 * statev(0) = 20.0;  // Reference temperature 20°C
 *
 * vec Etot = {0.0, 0.005, 0.0, 0.0, 0.0, 0.0};  // 0.5% transverse strain (induces damage)
 * vec DEtot = {0.0, 0.0001, 0.0, 0.0, 0.0, 0.0};
 * vec sigma = zeros(6);
 * mat Lt = zeros(6,6);
 * mat L = zeros(6,6);
 * vec sigma_in = zeros(6);
 * mat DR = eye(3,3);
 *
 * umat_damage_LLD_0(Etot, DEtot, sigma, Lt, L, sigma_in, DR,
 *                   18, props, 10, statev, 20.0, 0.0, 0.0, 1.0,
 *                   Wm, Wm_r, Wm_ir, Wm_d, 3, 3, false, 0, tnew_dt);
 *
 * // Check damage state
 * double d22 = statev(1);  // Transverse damage
 * double d12 = statev(2);  // Shear damage
 * cout << "Transverse damage: " << d22 << ", Shear damage: " << d12 << endl;
 * @endcode
 *
 * **References:**
 * - Ladevèze, P., & Le Dantec, E. (1992). "Damage modelling of the elementary ply for laminated composites." *Composites Science and Technology*, 43(3), 257-267.
 * - Allix, O., & Ladevèze, P. (1992). "Interlaminar interface modelling for the prediction of delamination." *Composite Structures*, 22(4), 235-242.
 * - Ladevèze, P. (1992). "A damage computational method for composite structures." *Computers & Structures*, 44(1-2), 79-87.
 * - Matzenmiller, A., Lubliner, J., & Taylor, R. L. (1995). "A constitutive model for anisotropic damage in fiber-composites." *Mechanics of Materials*, 20(2), 125-152.
 * - Pinho, S. T., et al. (2012). "Material and structural response of polymer-matrix fibre-reinforced composites." *Journal of Composite Materials*, 46(19-20), 2313-2341.
 */
void umat_damage_LLD_0(const std::string &umat_name, const arma::vec &Etot, const arma::vec &DEtot, arma::vec &sigma, arma::mat &Lt, arma::mat &L, const arma::mat &DR, const int &nprops, const arma::vec &props, const int &nstatev, arma::vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, double &tnew_dt);


/** @} */ // end of umat_mechanical group

} //namespace simcoon
