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
 * @file Hill_chaboche_ccp.hpp
 * @brief Elastic-plastic material with Hill anisotropic yield and Chaboche kinematic hardening
 * @author Yves Chemisky
 * @version 1.0
 */

#pragma once
#include <string>
#include <armadillo>

namespace simcoon{

/** @addtogroup umat_mechanical
 *  @{
 */

/**
 * @brief Elastic-plastic model combining Hill anisotropic yield criterion with Chaboche kinematic hardening
 *
 * @details This function implements an advanced plasticity model combining the Hill (1948)
 * anisotropic yield criterion with Chaboche nonlinear kinematic hardening. This allows
 * modeling of anisotropic materials under cyclic loading conditions. The model features:
 * - Hill (1948) quadratic anisotropic yield criterion
 * - Multiple Armstrong-Frederick backstresses for nonlinear kinematic hardening
 * - Optional Voce-type isotropic hardening
 * - Convex Cutting Plane algorithm for return mapping
 * - Consistent tangent modulus for implicit FE analysis
 *
 * **Constitutive Equations:**
 *
 * The anisotropic yield function with kinematic hardening is defined as:
 * \f[
 * \Phi(\boldsymbol{\sigma}, \mathbf{X}, p) = \sigma_{eq}^{Hill}(\boldsymbol{\sigma} - \mathbf{X}) - R(p) - \sigma_Y \leq 0
 * \f]
 * where the Hill equivalent stress of the shifted stress tensor is:
 * \f[
 * \sigma_{eq}^{Hill}(\boldsymbol{\eta}) = \sqrt{\boldsymbol{\eta} : \mathbf{P}^{Hill} : \boldsymbol{\eta}}
 * \f]
 * with \f$ \boldsymbol{\eta} = \boldsymbol{\sigma} - \mathbf{X} \f$ being the effective (shifted) stress.
 *
 * **Hill Anisotropic Yield Surface:**
 *
 * The Hill projection tensor \f$ \mathbf{P}^{Hill} \f$ is defined by six parameters \f$ (F, G, H, L, M, N) \f$:
 * \f[
 * (\sigma_{eq}^{Hill})^2 = F(\eta_{22} - \eta_{33})^2 + G(\eta_{33} - \eta_{11})^2 + H(\eta_{11} - \eta_{22})^2 + 2L\eta_{23}^2 + 2M\eta_{13}^2 + 2N\eta_{12}^2
 * \f]
 *
 * **Armstrong-Frederick Backstress Evolution:**
 *
 * Each backstress component evolves with dynamic recovery:
 * \f[
 * \dot{\mathbf{X}}_i = \frac{2}{3} C_i \dot{\boldsymbol{\varepsilon}}^p - D_i \mathbf{X}_i \dot{p}
 * \f]
 * where the total backstress is: \f$ \mathbf{X} = \sum_{i=1}^{N_{kin}} \mathbf{X}_i \f$
 *
 * **Flow Rule:**
 *
 * The plastic strain rate follows the associative flow rule:
 * \f[
 * \dot{\boldsymbol{\varepsilon}}^p = \dot{p} \frac{\partial \Phi}{\partial \boldsymbol{\sigma}} = \dot{p} \frac{\mathbf{P}^{Hill} : \boldsymbol{\eta}}{\sigma_{eq}^{Hill}(\boldsymbol{\eta})}
 * \f]
 *
 * **Applications:**
 *
 * This model is suited for:
 * - Rolled sheet metals under cyclic loading
 * - Textured polycrystalline materials with Bauschinger effect
 * - Drawn tubes and wires with anisotropic yield and kinematic hardening
 * - Any orthotropic material exhibiting cyclic plasticity
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
 * | props[5] | \f$ N_{kin} \f$ | Number of kinematic hardening terms | - |
 * | props[6] | \f$ F \f$ | Hill parameter F | 1/Stress² |
 * | props[7] | \f$ G \f$ | Hill parameter G | 1/Stress² |
 * | props[8] | \f$ H \f$ | Hill parameter H | 1/Stress² |
 * | props[9] | \f$ L \f$ | Hill parameter L | 1/Stress² |
 * | props[10] | \f$ M \f$ | Hill parameter M | 1/Stress² |
 * | props[11] | \f$ N \f$ | Hill parameter N | 1/Stress² |
 * | props[12+2j] | \f$ Q_j \f$ | Isotropic saturation stress (j-th term) | Stress |
 * | props[13+2j] | \f$ b_j \f$ | Isotropic hardening rate (j-th term) | 1/Strain |
 * | props[12+2N_{iso}+2i] | \f$ C_i \f$ | Kinematic modulus (i-th backstress) | Stress |
 * | props[13+2N_{iso}+2i] | \f$ D_i \f$ | Dynamic recovery (i-th backstress) | - |
 *
 * **State Variables (statev):**
 *
 * Total: \f$ n_{statev} = 1 + 1 + 6 + 6 \times N_{kin} \f$
 *
 * | Index | Symbol | Description | Units |
 * |-------|--------|-------------|-------|
 * | statev[0] | \f$ T_{init} \f$ | Initial temperature | Temperature |
 * | statev[1] | \f$ p \f$ | Accumulated plastic strain | Strain |
 * | statev[2:7] | \f$ \boldsymbol{\varepsilon}^p \f$ | Plastic strain tensor (Voigt) | Strain |
 * | statev[8+6i:13+6i] | \f$ \mathbf{X}_i \f$ | Backstress i (Voigt) | Stress |
 *
 * @param Etot Total strain tensor at beginning of increment (Voigt notation: 6×1 vector)
 * @param DEtot Strain increment tensor (Voigt notation: 6×1 vector)
 * @param sigma Stress tensor (Voigt notation: 6×1 vector) [output]
 * @param Lt Consistent tangent modulus (6×6 matrix) [output]
 * @param L Elastic stiffness tensor (6×6 matrix) [output]
 * @param sigma_in Internal stress for explicit solvers (6×1 vector) [output]
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
 * @param solver_type Solver type: 0=implicit, 1=explicit, 2=dynamic implicit
 * @param tnew_dt Suggested new time step size for adaptive time stepping [output]
 *
 * @note Material axes must align with global coordinates or use rotation tensors
 * @note For isotropic yield with Chaboche hardening, use umat_plasticity_chaboche_CCP
 * @note Captures both anisotropic yielding and Bauschinger effect
 *
 * @see umat_plasticity_hill_isoh_CCP() for Hill with isotropic hardening only
 * @see umat_plasticity_chaboche_CCP() for isotropic (von Mises) Chaboche model
 *
 * **References:**
 * - Hill, R. (1948). "A theory of the yielding and plastic flow of anisotropic metals." *Proc. Roy. Soc. A*, 193(1033), 281-297.
 * - Chaboche, J. L. (1986). "Time-independent constitutive theories for cyclic plasticity." *Int. J. Plasticity*, 2(2), 149-188.
 * - Barlat, F., et al. (2005). "Linear transformation-based anisotropic yield functions." *Int. J. Plasticity*, 21(5), 1009-1039.
 */
void umat_hill_chaboche_CCP(const std::string &umat_name, const arma::vec &Etot, const arma::vec &DEtot, arma::vec &sigma, arma::mat &Lt, arma::mat &L, const arma::mat &DR, const int &nprops, const arma::vec &props, const int &nstatev, arma::vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, double &tnew_dt);


/** @} */ // end of umat_mechanical group

} //namespace simcoon
