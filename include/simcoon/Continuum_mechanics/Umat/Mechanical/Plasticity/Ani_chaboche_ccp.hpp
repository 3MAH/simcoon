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
 * @file Ani_chaboche_ccp.hpp
 * @brief Elastic-plastic material with general anisotropic yield and Chaboche kinematic hardening
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
 * @brief Elastic-plastic model with general anisotropic yield criterion and Chaboche kinematic hardening
 *
 * @details This function implements an advanced anisotropic plasticity model combining a general
 * quadratic anisotropic yield criterion with Chaboche nonlinear kinematic hardening. This model
 * provides maximum flexibility for capturing complex anisotropic cyclic plasticity behavior.
 *
 * **Key Features:**
 * - General quadratic anisotropic yield criterion (21 independent parameters)
 * - Multiple Armstrong-Frederick backstresses for nonlinear kinematic hardening
 * - Optional Voce-type isotropic hardening
 * - Convex Cutting Plane algorithm for return mapping
 * - Consistent tangent modulus for implicit FE analysis
 *
 * **General Anisotropic Yield Function:**
 *
 * The yield function is defined as:
 * \f[
 * \Phi(\boldsymbol{\sigma}, \mathbf{X}, p) = \sqrt{(\boldsymbol{\sigma} - \mathbf{X}) : \mathbf{H} : (\boldsymbol{\sigma} - \mathbf{X})} - R(p) - \sigma_Y \leq 0
 * \f]
 * where \f$ \mathbf{H} \f$ is a general fourth-order anisotropy tensor with up to 21 independent
 * components for a fully anisotropic material.
 *
 * **Anisotropy Tensor \f$ \mathbf{H} \f$:**
 *
 * In Voigt notation, the anisotropy tensor is represented as a 6×6 symmetric matrix:
 * \f[
 * \mathbf{H} = \begin{bmatrix}
 * H_{11} & H_{12} & H_{13} & H_{14} & H_{15} & H_{16} \\
 * H_{12} & H_{22} & H_{23} & H_{24} & H_{25} & H_{26} \\
 * H_{13} & H_{23} & H_{33} & H_{34} & H_{35} & H_{36} \\
 * H_{14} & H_{24} & H_{34} & H_{44} & H_{45} & H_{46} \\
 * H_{15} & H_{25} & H_{35} & H_{45} & H_{55} & H_{56} \\
 * H_{16} & H_{26} & H_{36} & H_{46} & H_{56} & H_{66}
 * \end{bmatrix}
 * \f]
 *
 * **Special Cases:**
 *
 * - **von Mises (isotropic)**: \f$ \mathbf{H} = \mathbf{I}_{dev} \f$ (deviatoric identity)
 * - **Hill orthotropic**: \f$ H_{ij} = 0 \f$ for \f$ i \leq 3 < j \f$ or \f$ j \leq 3 < i \f$
 * - **Transversely isotropic**: 5 independent parameters
 *
 * **Armstrong-Frederick Backstress Evolution:**
 *
 * Each backstress component evolves according to:
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
 * **Applications:**
 *
 * This model is suited for:
 * - Single crystals with complex slip system interactions
 * - Highly textured polycrystals requiring full anisotropy
 * - Materials with non-orthotropic symmetry (monoclinic, triclinic)
 * - Advanced composite materials with complex anisotropy
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
 * @note The anisotropy tensor must be positive semi-definite for convexity
 * @note For orthotropic materials, consider using Hill_chaboche_ccp instead
 * @note Parameter identification requires extensive multiaxial testing
 *
 * @see umat_hill_chaboche_CCP() for orthotropic Hill + Chaboche
 * @see umat_plasticity_chaboche_CCP() for isotropic von Mises + Chaboche
 *
 * **References:**
 * - Chaboche, J. L. (1986). "Time-independent constitutive theories for cyclic plasticity." *Int. J. Plasticity*, 2(2), 149-188.
 * - Barlat, F., et al. (2005). "Linear transformation-based anisotropic yield functions." *Int. J. Plasticity*, 21(5), 1009-1039.
 */
void umat_ani_chaboche_CCP(const arma::vec &Etot, const arma::vec &DEtot, arma::vec &sigma, arma::mat &Lt, arma::mat &L, arma::vec &sigma_in, const arma::mat &DR, const int &nprops, const arma::vec &props, const int &nstatev, arma::vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, const int &solver_type, double &tnew_dt);


/** @} */ // end of umat_mechanical group

} //namespace simcoon
