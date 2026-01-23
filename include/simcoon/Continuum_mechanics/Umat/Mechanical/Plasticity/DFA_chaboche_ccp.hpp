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
 * @file DFA_chaboche_ccp.hpp
 * @brief Elastic-plastic material with Distortional/Directional hardening and Chaboche kinematic hardening
 * @author Yves Chemisky
 * @version 1.0
 */

#pragma once
#include <armadillo>

namespace simcoon{

/** @addtogroup umat_mechanical
 *  @{
 */

/**
 * @brief Elastic-plastic model with Distortional/Directional Forming Anisotropy (DFA) and Chaboche hardening
 *
 * @details This function implements an advanced plasticity model that captures the evolution of
 * plastic anisotropy during deformation (distortional hardening). Unlike classical models where
 * the yield surface only translates (kinematic) or expands (isotropic), the DFA model allows
 * the yield surface shape to evolve with plastic straining.
 *
 * **Key Features:**
 * - Distortional hardening: yield surface shape evolves with deformation
 * - Multiple Armstrong-Frederick backstresses for kinematic hardening
 * - Optional Voce-type isotropic hardening
 * - Captures cross-hardening and directional effects
 * - Convex Cutting Plane algorithm for return mapping
 *
 * **Physical Motivation:**
 *
 * Classical hardening models (isotropic + kinematic) cannot capture:
 * - **Cross-hardening**: Yield stress in one direction affected by prior loading in another
 * - **Permanent softening**: Reduced yield in tension after prior compression (or vice versa)
 * - **Yield surface distortion**: Changes in yield locus shape observed in metals after pre-straining
 *
 * The DFA model addresses these limitations by allowing the anisotropy tensor \f$ \mathbf{H} \f$
 * to evolve with plastic deformation.
 *
 * **Yield Function with Evolving Anisotropy:**
 *
 * \f[
 * \Phi = \sqrt{(\boldsymbol{\sigma} - \mathbf{X}) : \mathbf{H}(p, \boldsymbol{\varepsilon}^p) : (\boldsymbol{\sigma} - \mathbf{X})} - R(p) - \sigma_Y \leq 0
 * \f]
 *
 * where the anisotropy tensor \f$ \mathbf{H} \f$ evolves according to:
 * \f[
 * \dot{\mathbf{H}} = f(\mathbf{H}, \boldsymbol{\sigma}, \dot{\boldsymbol{\varepsilon}}^p, \dot{p})
 * \f]
 *
 * **Directional Hardening:**
 *
 * The model tracks directional hardening through internal variables that remember
 * the loading history orientation. This allows capturing:
 * - **Latent hardening**: Increased yield stress in inactive slip systems
 * - **Texture evolution effects**: Macroscopic manifestation of microstructural changes
 * - **Anisotropic damage accumulation**: Direction-dependent material degradation
 *
 * **Applications:**
 *
 * - Sheet metal forming with complex strain paths
 * - Materials exhibiting strong cross-hardening (aluminum alloys, TWIP steels)
 * - Multi-pass forming operations with changing strain directions
 * - Accurate springback prediction in automotive stamping
 *
 * **State Variables (statev):**
 *
 * | Index | Symbol | Description | Units |
 * |-------|--------|-------------|-------|
 * | statev[0] | \f$ T_{init} \f$ | Initial temperature | Temperature |
 * | statev[1] | \f$ p \f$ | Accumulated plastic strain | Strain |
 * | statev[2:7] | \f$ \boldsymbol{\varepsilon}^p \f$ | Plastic strain tensor (Voigt) | Strain |
 * | statev[8+6i:13+6i] | \f$ \mathbf{X}_i \f$ | Backstress i (Voigt) | Stress |
 * | statev[...] | \f$ H_{ij} \f$ | Evolving anisotropy components | - |
 *
 * @param Etot Total strain tensor at beginning of increment (Voigt notation: 6×1 vector)
 * @param DEtot Strain increment tensor (Voigt notation: 6×1 vector)
 * @param sigma Stress tensor (Voigt notation: 6×1 vector) [output]
 * @param Lt Consistent tangent modulus (6×6 matrix) [output]
 * @param L Elastic stiffness tensor (6×6 matrix) [output]
 * @param sigma_in Internal stress for explicit solvers (6×1 vector) [output]
 * @param DR Rotation increment matrix (3×3) for objective integration
 * @param nprops Number of material properties
 * @param props Material properties vector
 * @param nstatev Number of state variables
 * @param statev State variables vector [input/output]
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
 * @note Requires careful parameter identification from cross-hardening tests
 * @note Computationally more expensive than standard Chaboche due to evolving anisotropy
 * @note Useful for complex forming simulations with strain path changes
 *
 * @see umat_plasticity_chaboche_CCP() for standard Chaboche without distortional hardening
 * @see umat_ani_chaboche_CCP() for fixed anisotropy Chaboche
 *
 * **References:**
 * - Barlat, F., et al. (2011). "An alternative to kinematic hardening in classical plasticity." *Int. J. Plasticity*, 27(9), 1309-1327.
 * - Teodosiu, C., & Hu, Z. (1998). "Microstructure in the continuum modeling of plastic anisotropy." *Proc. 19th Riso Int. Symp.*, 149-168.
 * - Holmedal, B. (2019). "Bauschinger effect modelled by yield surface distortions." *Int. J. Plasticity*, 123, 86-100.
 */
void umat_dfa_chaboche_CCP(const arma::vec &Etot, const arma::vec &DEtot, arma::vec &sigma, arma::mat &Lt, arma::mat &L, arma::vec &sigma_in, const arma::mat &DR, const int &nprops, const arma::vec &props, const int &nstatev, arma::vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, const int &solver_type, double &tnew_dt);


/** @} */ // end of umat_mechanical group

} //namespace simcoon
