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
 * @file plastic_kin_iso_ccp.hpp
 * @brief Elastic-plastic material with combined isotropic and kinematic hardening using CCP algorithm
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
 * @brief Elastic-plastic constitutive model with combined isotropic and linear kinematic hardening
 *
 * @details This function implements an elastic-plastic material model combining both isotropic
 * and linear kinematic hardening mechanisms. The model features:
 * - J2 (von Mises) plasticity with associative flow rule
 * - Power law isotropic hardening: \f$ R(p) = k \cdot p^m \f$
 * - Linear (Prager) kinematic hardening: \f$ \dot{\mathbf{X}} = \frac{2}{3} H_{kin} \dot{\boldsymbol{\varepsilon}}^p \f$
 * - Convex Cutting Plane algorithm for return mapping
 * - Consistent tangent modulus for implicit FE analysis
 *
 * **Constitutive Equations:**
 *
 * The yield function is defined as:
 * \f[
 * \Phi(\boldsymbol{\sigma}, \mathbf{X}, p) = \sigma_{eq}(\boldsymbol{\sigma} - \mathbf{X}) - R(p) - \sigma_Y \leq 0
 * \f]
 * where:
 * - \f$ \sigma_{eq}(\boldsymbol{\eta}) = \sqrt{\frac{3}{2} \boldsymbol{\eta}_{dev} : \boldsymbol{\eta}_{dev}} \f$ is the von Mises equivalent stress
 * - \f$ \boldsymbol{\eta} = \boldsymbol{\sigma} - \mathbf{X} \f$ is the shifted (effective) stress
 * - \f$ \mathbf{X} \f$ is the backstress (kinematic hardening variable)
 * - \f$ R(p) = k \cdot p^m \f$ is the isotropic hardening function
 * - \f$ \sigma_Y \f$ is the initial yield stress
 *
 * **Kinematic Hardening (Prager Linear):**
 *
 * The backstress evolves according to linear (Prager) kinematic hardening:
 * \f[
 * \dot{\mathbf{X}} = \frac{2}{3} H_{kin} \dot{\boldsymbol{\varepsilon}}^p = \frac{2}{3} H_{kin} \dot{p} \mathbf{n}
 * \f]
 * where:
 * - \f$ H_{kin} \f$ is the kinematic hardening modulus
 * - \f$ \mathbf{n} = \frac{3}{2} \frac{\boldsymbol{\eta}_{dev}}{\sigma_{eq}} \f$ is the flow direction
 *
 * **Combined Hardening Effect:**
 *
 * The combined hardening provides:
 * - **Isotropic hardening**: Expands the yield surface uniformly (captures strain hardening)
 * - **Kinematic hardening**: Translates the yield surface center (captures Bauschinger effect)
 *
 * **Plastic Flow Rule:**
 *
 * The plastic strain rate follows the associative flow rule:
 * \f[
 * \dot{\boldsymbol{\varepsilon}}^p = \dot{p} \frac{\partial \Phi}{\partial \boldsymbol{\sigma}} = \dot{p} \mathbf{n}
 * \f]
 *
 * **Incremental Update:**
 *
 * For an increment \f$ \Delta p \f$:
 * \f[
 * \boldsymbol{\sigma}_{n+1} = \boldsymbol{\sigma}^{trial} - \Delta p \cdot \mathbf{L} : \mathbf{n}
 * \f]
 * \f[
 * \mathbf{X}_{n+1} = \mathbf{X}_n + \frac{2}{3} H_{kin} \Delta p \cdot \mathbf{n}
 * \f]
 *
 * **Consistent Tangent Modulus:**
 *
 * The algorithmic tangent accounts for both hardening mechanisms:
 * \f[
 * \mathbf{L}_t = \mathbf{L} - \frac{(\mathbf{L}:\mathbf{n}) \otimes (\mathbf{n}:\mathbf{L})}{\mathbf{n}:\mathbf{L}:\mathbf{n} + H_{iso} + H_{kin}}
 * \f]
 * where \f$ H_{iso} = \frac{dR}{dp} = k \cdot m \cdot p^{m-1} \f$ is the isotropic hardening modulus.
 *
 * **Material Parameters (props):**
 *
 * | Index | Symbol | Description | Units |
 * |-------|--------|-------------|-------|
 * | props[0] | \f$ E \f$ | Young's modulus | Stress |
 * | props[1] | \f$ \nu \f$ | Poisson's ratio | - |
 * | props[2] | \f$ \alpha \f$ | Thermal expansion coefficient | 1/Temperature |
 * | props[3] | \f$ \sigma_Y \f$ | Initial yield stress | Stress |
 * | props[4] | \f$ k \f$ | Isotropic hardening coefficient | Stress |
 * | props[5] | \f$ m \f$ | Isotropic hardening exponent | - |
 * | props[6] | \f$ H_{kin} \f$ | Kinematic hardening modulus | Stress |
 *
 * **State Variables (statev):**
 *
 * Total state variables required: \f$ n_{statev} = 14 \f$
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
 * | statev[8] | \f$ X_{11} \f$ | Backstress component 11 | Stress |
 * | statev[9] | \f$ X_{22} \f$ | Backstress component 22 | Stress |
 * | statev[10] | \f$ X_{33} \f$ | Backstress component 33 | Stress |
 * | statev[11] | \f$ X_{12} \f$ | Backstress component 12 | Stress |
 * | statev[12] | \f$ X_{13} \f$ | Backstress component 13 | Stress |
 * | statev[13] | \f$ X_{23} \f$ | Backstress component 23 | Stress |
 *
 * @param Etot Total strain tensor at beginning of increment (Voigt notation: 6×1 vector)
 * @param DEtot Strain increment tensor (Voigt notation: 6×1 vector)
 * @param sigma Stress tensor (Voigt notation: 6×1 vector) [output]
 * @param Lt Consistent tangent modulus (6×6 matrix) [output]
 * @param L Elastic stiffness tensor (6×6 matrix) [output]
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
 * @param Wm_ir Irrecoverable work stored in hardening [output]
 * @param Wm_d Dissipated (plastic) work [output]
 * @param ndi Number of direct stress components (typically 3)
 * @param nshr Number of shear stress components (typically 3)
 * @param start Flag indicating first increment (true) or continuation (false)
 * @param solver_type Solver type: 0=implicit, 1=explicit, 2=dynamic implicit
 * @param tnew_dt Suggested new time step size for adaptive time stepping [output]
 *
 * @note Linear kinematic hardening cannot capture cyclic softening/hardening saturation
 * @note For saturating kinematic hardening, use the Chaboche model (Armstrong-Frederick)
 * @note The Bauschinger effect is proportional to the kinematic hardening modulus
 * @note Setting \f$ H_{kin} = 0 \f$ recovers pure isotropic hardening
 * @note Setting \f$ k = 0 \f$ (or \f$ m = 0 \f$) recovers pure kinematic hardening
 *
 * @see umat_plasticity_iso_CCP() for pure isotropic hardening
 * @see umat_plasticity_chaboche_CCP() for nonlinear kinematic hardening
 *
 * **References:**
 * - Prager, W. (1956). "A new method of analyzing stresses and strains in work-hardening plastic solids." *J. Appl. Mech.*, 23, 493-496.
 * - Ziegler, H. (1959). "A modification of Prager's hardening rule." *Q. Appl. Math.*, 17, 55-65.
 * - Lemaitre, J., & Chaboche, J. L. (1990). *Mechanics of Solid Materials*. Cambridge University Press.
 * - Simo, J. C., & Hughes, T. J. R. (1998). *Computational Inelasticity*. Springer.
 */
void umat_plasticity_kin_iso_CCP(const std::string &umat_name, const arma::vec &Etot, const arma::vec &DEtot, arma::vec &sigma, arma::mat &Lt, arma::mat &L, const arma::mat &DR, const int &nprops, const arma::vec &props, const int &nstatev, arma::vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, double &tnew_dt);


/** @} */ // end of umat_mechanical group

} //namespace simcoon
