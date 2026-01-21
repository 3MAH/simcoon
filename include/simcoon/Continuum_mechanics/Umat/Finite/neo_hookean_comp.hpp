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

///@file neo_hookean_comp.hpp
///@brief User subroutine for Isotropic elastic materials in 3D case
///@version 1.0

#pragma once

#include <armadillo>

namespace simcoon{

/**
 * @file neo_hookean_comp.hpp
 * @brief Compressible Neo-Hookean hyperelastic material model for finite strain analysis
 * @author Yves Chemisky
 * @version 1.0
 */

/** @addtogroup umat_finite
 *  @{
 */

/**
 * @brief Compressible Neo-Hookean hyperelastic constitutive model for finite strain analysis
 *
 * @details This function implements the compressible Neo-Hookean hyperelastic model, one of the simplest
 * hyperelastic models suitable for large deformation analysis of rubber-like materials and soft tissues.
 *
 * **Strain Energy Function:**
 *
 * The total strain energy density is decomposed into isochoric (volume-preserving) and volumetric parts:
 * \f[
 * \Psi(\mathbf{F}) = \bar{W}(\bar{I}_1) + U(J)
 * \f]
 *
 * **Isochoric Part (Neo-Hookean):**
 * \f[
 * \bar{W}(\bar{I}_1) = \frac{\mu}{2} (\bar{I}_1 - 3)
 * \f]
 * where:
 * - \f$ \mu = \frac{E}{2(1+\nu)} \f$ is the shear modulus
 * - \f$ \bar{I}_1 = J^{-2/3} I_1 = J^{-2/3} \text{tr}(\mathbf{C}) \f$ is the first isochoric invariant
 * - \f$ I_1 = \lambda_1^2 + \lambda_2^2 + \lambda_3^2 \f$ is the first invariant of the right Cauchy-Green tensor
 * - \f$ J = \det(\mathbf{F}) \f$ is the volume ratio
 *
 * **Volumetric Part:**
 * \f[
 * U(J) = \frac{\kappa}{2} (J - 1)^2
 * \f]
 * where:
 * - \f$ \kappa = \frac{E}{3(1-2\nu)} \f$ is the bulk modulus
 *
 * **Cauchy Stress Tensor:**
 *
 * The Cauchy stress is computed as:
 * \f[
 * \boldsymbol{\sigma} = \boldsymbol{\sigma}_{\text{iso}} + \boldsymbol{\sigma}_{\text{vol}}
 * \f]
 *
 * Isochoric part:
 * \f[
 * \boldsymbol{\sigma}_{\text{iso}} = \frac{\mu}{J} J^{-2/3} \left( \mathbf{b} - \frac{I_1}{3} \mathbf{I} \right)
 * \f]
 *
 * Volumetric part:
 * \f[
 * \boldsymbol{\sigma}_{\text{vol}} = \kappa (J - 1) \mathbf{I}
 * \f]
 * where \f$ \mathbf{b} = \mathbf{F} \mathbf{F}^T \f$ is the left Cauchy-Green tensor.
 *
 * **Consistent Tangent Modulus:**
 *
 * The algorithmic tangent modulus is computed for implicit finite element analysis:
 * \f[
 * \mathbf{L}_t = \mathbf{L}_{\text{iso}} + \mathbf{L}_{\text{vol}}
 * \f]
 *
 * This ensures quadratic convergence in Newton-Raphson iterations.
 *
 * **Material Parameters (props):**
 *
 * | Index | Symbol | Description | Units |
 * |-------|--------|-------------|-------|
 * | props[0] | \f$ E \f$ | Young's modulus | Stress |
 * | props[1] | \f$ \nu \f$ | Poisson's ratio (should be < 0.5 for compressibility) | - |
 * | props[2] | \f$ \alpha \f$ | Thermal expansion coefficient | 1/Temperature |
 *
 * **State Variables (statev):**
 *
 * No internal state variables are required for this hyperelastic model (purely elastic response).
 *
 * @param Etot Total Green-Lagrange strain tensor at beginning of increment (Voigt notation: 6×1)
 * @param DEtot Green-Lagrange strain increment tensor (Voigt notation: 6×1)
 * @param F0 Deformation gradient at beginning of increment (3×3 matrix)
 * @param F1 Deformation gradient at end of increment (3×3 matrix)
 * @param sigma Cauchy stress tensor (Voigt notation: 6×1) [output]
 * @param Lt Consistent tangent modulus \f$ \mathbf{L}_t = \frac{\partial \boldsymbol{\sigma}}{\partial \boldsymbol{\varepsilon}} \f$ (6×6 matrix) [output]
 * @param L Elastic stiffness tensor (6×6 matrix) [output]
 * @param sigma_in Internal stress contribution for explicit solvers (6×1 vector) [output]
 * @param DR Rotation increment matrix (3×3) for objective integration
 * @param nprops Number of material properties
 * @param props Material properties vector (see table above)
 * @param nstatev Number of state variables (0 for this model)
 * @param statev State variables vector (unused for this model) [input/output]
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
 * @param solver_type Solver type: 0=implicit, 1=explicit, 2=dynamic implicit
 * @param tnew_dt Suggested new time step size for adaptive time stepping [output]
 *
 * @note The compressible Neo-Hookean model is suitable for moderate strains (< 100%)
 * @note For nearly incompressible materials (nu → 0.5), use the incompressible version instead
 * @note Thermal strains are handled through thermal expansion coefficient alpha
 * @note This is a purely hyperelastic model with no energy dissipation (Wm_ir = Wm_d = 0)
 *
 * @see isochoric_invariants() for computing isochoric strain invariants
 * @see sigma_iso_hyper_invariants() for isochoric Cauchy stress computation
 * @see sigma_vol_hyper() for volumetric Cauchy stress computation
 * @see L_iso_hyper_invariants() for isochoric tangent modulus
 * @see L_vol_hyper() for volumetric tangent modulus
 * @see L_Cauchy_Green() for left Cauchy-Green tensor computation
 *
 * @code
 * // Example usage:
 * mat F0 = eye(3,3);
 * mat F1 = {{1.1, 0.0, 0.0}, {0.0, 1.05, 0.0}, {0.0, 0.0, 0.95}};  // 10% stretch
 * vec Etot = Green_Lagrange(F1);
 * vec DEtot = Etot;  // Starting from reference
 * vec sigma = zeros(6);
 * mat Lt = zeros(6,6);
 * mat L = zeros(6,6);
 * vec sigma_in = zeros(6);
 * mat DR = eye(3,3);
 * vec props = {1000, 0.45, 1e-5};  // E=1000, nu=0.45 (nearly incompressible), alpha=1e-5
 * vec statev = zeros(1);  // Not used
 *
 * umat_neo_hookean_comp(Etot, DEtot, F0, F1, sigma, Lt, L, sigma_in, DR,
 *                       3, props, 0, statev, 20.0, 0.0, 0.0, 1.0,
 *                       Wm, Wm_r, Wm_ir, Wm_d, 3, 3, true, 0, tnew_dt);
 * @endcode
 *
 * **References:**
 * - Bonet, J., & Wood, R. D. (2008). *Nonlinear Continuum Mechanics for Finite Element Analysis*. Cambridge University Press.
 * - Holzapfel, G. A. (2000). *Nonlinear Solid Mechanics: A Continuum Approach for Engineering*. Wiley.
 * - Connolly, S. J., et al. (2019). "Automatic differentiation based formulation of computational models." *Computational Mechanics*, 64, 1273-1288.
 */
void umat_neo_hookean_comp(const arma::vec &Etot, const arma::vec &DEtot, const arma::mat &F0, const arma::mat &F1, arma::vec &sigma, arma::mat &Lt, arma::mat &L, arma::vec &sigma_in, const arma::mat &DR, const int &nprops, const arma::vec &props, const int &nstatev, arma::vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, const int &solver_type, double &tnew_dt);
                            

/** @} */ // end of umat_finite group

} //namespace simcoon
