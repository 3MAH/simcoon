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
 * @file Generic_chaboche_ccp.hpp
 * @brief Generic elastic-plastic material framework with configurable yield criterion and Chaboche hardening
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
 * @brief Generic elastic-plastic framework with flexible yield criterion and Chaboche kinematic hardening
 *
 * @details This function implements a generic plasticity framework that allows for flexible
 * configuration of the yield criterion while using Chaboche-type hardening. It serves as a
 * template for implementing various plasticity models with different yield functions.
 *
 * **Key Features:**
 * - Configurable yield criterion through function pointers or template parameters
 * - Multiple Armstrong-Frederick backstresses for nonlinear kinematic hardening
 * - Optional Voce-type isotropic hardening
 * - Convex Cutting Plane algorithm for return mapping
 * - Consistent tangent modulus for implicit FE analysis
 *
 * **Generic Yield Function:**
 *
 * The yield function takes the general form:
 * \f[
 * \Phi(\boldsymbol{\sigma}, \mathbf{X}, p, \boldsymbol{\alpha}) = f(\boldsymbol{\sigma} - \mathbf{X}) - R(p) - \sigma_Y \leq 0
 * \f]
 * where:
 * - \f$ f(\boldsymbol{\eta}) \f$ is a general equivalent stress function (configurable)
 * - \f$ \boldsymbol{\alpha} \f$ represents additional internal variables
 *
 * **Supported Yield Criteria:**
 *
 * The generic framework can accommodate:
 * - **von Mises**: \f$ f = \sqrt{\frac{3}{2} \mathbf{s}:\mathbf{s}} \f$
 * - **Hill (1948)**: Quadratic orthotropic
 * - **Barlat Yld2000-2d**: Advanced anisotropic for sheet metals
 * - **Drucker-Prager**: Pressure-dependent for geomaterials
 * - **Hosford**: Non-quadratic isotropic
 * - **Custom user-defined**: Through function interface
 *
 * **Armstrong-Frederick Backstress Evolution:**
 *
 * Each backstress evolves according to:
 * \f[
 * \dot{\mathbf{X}}_i = \frac{2}{3} C_i \dot{\boldsymbol{\varepsilon}}^p - D_i \mathbf{X}_i \dot{p}
 * \f]
 *
 * **Isotropic Hardening (Voce Law):**
 *
 * \f[
 * R(p) = \sum_{j=1}^{N_{iso}} Q_j (1 - e^{-b_j p})
 * \f]
 *
 * **Framework Architecture:**
 *
 * The generic model separates:
 * 1. **Yield criterion evaluation**: \f$ f(\boldsymbol{\eta}) \f$ and \f$ \partial f / \partial \boldsymbol{\eta} \f$
 * 2. **Hardening laws**: Evolution of \f$ \mathbf{X}_i \f$ and \f$ R \f$
 * 3. **Return mapping**: CCP algorithm for stress update
 * 4. **Consistent tangent**: Algorithmic modulus computation
 *
 * This allows easy extension to new yield criteria while reusing the hardening and integration framework.
 *
 * **State Variables (statev):**
 *
 * | Index | Symbol | Description | Units |
 * |-------|--------|-------------|-------|
 * | statev[0] | \f$ T_{init} \f$ | Initial temperature | Temperature |
 * | statev[1] | \f$ p \f$ | Accumulated plastic strain | Strain |
 * | statev[2:7] | \f$ \boldsymbol{\varepsilon}^p \f$ | Plastic strain tensor (Voigt) | Strain |
 * | statev[8+6i:13+6i] | \f$ \mathbf{X}_i \f$ | Backstress i (Voigt) | Stress |
 * | statev[...] | \f$ \boldsymbol{\alpha} \f$ | Additional criterion-specific variables | - |
 *
 * @param Etot Total strain tensor at beginning of increment (Voigt notation: 6×1 vector)
 * @param DEtot Strain increment tensor (Voigt notation: 6×1 vector)
 * @param sigma Stress tensor (Voigt notation: 6×1 vector) [output]
 * @param Lt Consistent tangent modulus (6×6 matrix) [output]
 * @param L Elastic stiffness tensor (6×6 matrix) [output]
 * @param sigma_in Internal stress for explicit solvers (6×1 vector) [output]
 * @param DR Rotation increment matrix (3×3) for objective integration
 * @param nprops Number of material properties
 * @param props Material properties vector (criterion-dependent)
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
 * @note Use this as a base for implementing new yield criteria
 * @note Requires implementation of yield function and its derivatives
 * @note See specialized models (Hill, Ani) for specific applications
 *
 * @see umat_plasticity_chaboche_CCP() for standard von Mises Chaboche
 * @see umat_hill_chaboche_CCP() for Hill anisotropic Chaboche
 * @see umat_ani_chaboche_CCP() for general anisotropic Chaboche
 *
 * **References:**
 * - Chaboche, J. L. (2008). "A review of some plasticity and viscoplasticity constitutive theories." *Int. J. Plasticity*, 24(10), 1642-1693.
 * - Simo, J. C., & Hughes, T. J. R. (1998). *Computational Inelasticity*. Springer.
 */
void umat_generic_chaboche_CCP(const arma::vec &Etot, const arma::vec &DEtot, arma::vec &sigma, arma::mat &Lt, arma::mat &L, arma::vec &sigma_in, const arma::mat &DR, const int &nprops, const arma::vec &props, const int &nstatev, arma::vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, const int &solver_type, double &tnew_dt);


/** @} */ // end of umat_mechanical group

} //namespace simcoon
