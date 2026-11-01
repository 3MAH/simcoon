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
#include <simcoon/parameter.hpp>

namespace simcoon{

/** @addtogroup umat_mechanical
 *  @{
 */

/**
 * @brief Ladevèze-Le Dantec (LLD) damage model of the unidirectional ply, with coupled
 *        transverse-shear plasticity in effective stress (LLDM0)
 *
 * @details Matrix-dominated degradation of a unidirectional fibre-reinforced ply: two scalar
 * damage variables, \f$ d_{22} \f$ (transverse) and \f$ d_{12} \f$ (shear), and a plasticity of the
 * matrix written in the effective (undamaged) stress. No fibre failure. The equations below are in
 * the fibre frame, fibres along 1; @c props(0) puts the fibres along the global axis 1, 2 or 3, the
 * law then runs in that frame and returns its stress, tangent and plastic strain in the global one.
 *
 * **Damaged elasticity.** Transversely isotropic compliance \f$ \mathbf{S}_0 \f$ (moduli
 * \f$ E_L, E_T, \nu_{TL}, \nu_{TT}, G_{LT} \f$, \f$ G_{TT} = E_T/(2(1+\nu_{TT})) \f$), damaged as
 * \f[
 *   \tilde{\mathbf{S}}(d_{22}, d_{12}):\quad
 *   S_{22}, S_{33}, S_{23}, S_{66} \leftarrow \frac{S^0}{1-d_{22}}, \qquad
 *   S_{44}, S_{55} \leftarrow \frac{S^0}{1-d_{12}},
 * \f]
 * the fibre terms \f$ S_{11}, S_{12}, S_{13} \f$ unchanged, and
 * \f$ \boldsymbol{\sigma} = \tilde{\mathbf{L}} : (\boldsymbol{\varepsilon} - \boldsymbol{\alpha}\,\Delta\theta - \boldsymbol{\varepsilon}^p) \f$,
 * \f$ \tilde{\mathbf{L}} = \tilde{\mathbf{S}}^{-1} \f$. \f$ \nu_{TL} \f$ is the MINOR Poisson ratio
 * (\f$ \nu_{TL} = \nu_{LT} E_T/E_L \f$, about 0.02 for carbon/epoxy), the convention of
 * L_isotrans().
 *
 * **Plasticity in effective stress.** With \f$ \tilde{\boldsymbol{\sigma}} = \mathbf{L}_0 : (\boldsymbol{\varepsilon} - \boldsymbol{\alpha}\,\Delta\theta - \boldsymbol{\varepsilon}^p) \f$
 * and \f$ \boldsymbol{\Theta} = \mathrm{diag}(0, A_{ts}, A_{ts}, 1, 1, 0) \f$ (transverse and
 * fibre-plane shear stresses only),
 * \f[
 *   f = \sigma^{\mathrm{Mises}}\!\left(\boldsymbol{\Theta}\tilde{\boldsymbol{\sigma}}\right) - R(p) - \sigma_{ts,0} \leq 0,
 *   \qquad R(p) = \beta_{ts}\, p^{m_{ts}},
 * \f]
 * with the associated flow
 * \f$ \dot{\boldsymbol{\varepsilon}}^p = \dot p\, \partial f/\partial\tilde{\boldsymbol{\sigma}} = \dot p\, \boldsymbol{\Theta}\,\boldsymbol{\eta}(\boldsymbol{\Theta}\tilde{\boldsymbol{\sigma}}) \f$
 * (\f$ \boldsymbol{\eta} \f$ = eta_stress()). The effective stress does not depend on the damage:
 * the plasticity is solved first (convex cutting plane, Fischer-Burmeister), the damage then at
 * fixed plastic strain.
 *
 * **Damage.** The driving forces of the criterion are unilateral in the transverse direction,
 * \f[
 *   Y_{22} = \frac{\langle\sigma_{22}\rangle_+^2 + \langle\sigma_{33}\rangle_+^2}{2E_T(1-d_{22})^2}, \quad
 *   Y_{12} = \frac{\sigma_{12}^2 + \sigma_{13}^2}{2G_{LT}(1-d_{12})^2}, \quad
 *   Y_{ts} = \sqrt{Y_{12} + b\,Y_{22}}, \quad Y_t = \sqrt{Y_{22}},
 * \f]
 * and both damage variables follow \f$ Y_{ts} \f$ (Kuhn-Tucker, \f$ \dot d \geq 0 \f$):
 * \f[
 *   d_{22} + \lambda(d_{22}) = \frac{\langle Y_{ts} - Y_{22,0}\rangle_+}{Y_{22,c}}, \qquad
 *   d_{12} + \lambda(d_{12}) = \frac{\langle Y_{ts} - Y_{12,0}\rangle_+}{Y_{12,c}},
 * \f]
 * \f$ \lambda \f$ the penalty of lagrange_pow_1() keeping \f$ d < 1 \f$. The two are coupled through
 * \f$ Y_{ts} \f$ and through the stress they both relax (the damage directions
 * \f$ \boldsymbol{\Lambda}_j = \partial\tilde{\mathbf{S}}/\partial d_j : \boldsymbol{\sigma} \f$):
 * the 2x2 system is full and non-symmetric. Transverse failure: when \f$ Y_t > Y_{22,u} \f$,
 * \f$ d_{22} \f$ jumps to 0.99 and stays there; \f$ d_{12} \f$ keeps its own law.
 *
 * **Tangent.** The tangent of the two-stage update: \f$ \mathrm{d}\boldsymbol{\varepsilon}^e/\mathrm{d}\boldsymbol{\varepsilon} \f$
 * from the plasticity, the flow direction held fixed over the increment (continuum: its error
 * is first order in the increment), then the exact linearization of the active damage
 * conditions,
 * \f$ \mathbf{L}_t = \tilde{\mathbf{L}}\,\frac{\mathrm{d}\boldsymbol{\varepsilon}^e}{\mathrm{d}\boldsymbol{\varepsilon}}
 * - \sum_j \tilde{\mathbf{L}}\boldsymbol{\Lambda}_j \otimes \frac{\mathrm{d}d_j}{\mathrm{d}\boldsymbol{\varepsilon}} \f$.
 * @p tangent_mode is not used. @p L returns the damaged stiffness \f$ \tilde{\mathbf{L}} \f$.
 *
 * **Work split.** \f$ W_m^{ir} \f$ stores the hardening energy
 * (\f$ \tfrac12 (R_n + R_{n+1})\Delta p \f$); \f$ W_m^{d} \f$ the plastic dissipation
 * \f$ \bar{\boldsymbol{\sigma}}:\Delta\boldsymbol{\varepsilon}^p - \bar R\,\Delta p \f$ and the
 * damage dissipation \f$ \bar{Y}^{th}_j \Delta d_j \f$ with the thermodynamic forces
 * \f$ Y^{th}_j = \tfrac12\,\boldsymbol{\sigma}:\boldsymbol{\Lambda}_j \f$; \f$ W_m^{r} \f$ the rest
 * (overbars: increment averages). \f$ W_m = W_m^r + W_m^{ir} + W_m^d \f$ at every increment.
 *
 * **Material Parameters (props, 18):**
 * | Index | Symbol | Description |
 * |-------|--------|-------------|
 * | props[0] | axis | Fibre axis (1, 2 or 3) |
 * | props[1] | \f$ E_L \f$ | Longitudinal Young's modulus |
 * | props[2] | \f$ E_T \f$ | Transverse Young's modulus |
 * | props[3] | \f$ \nu_{TL} \f$ | Minor Poisson ratio \f$ \nu_{LT} E_T/E_L \f$ |
 * | props[4] | \f$ \nu_{TT} \f$ | Transverse Poisson ratio |
 * | props[5] | \f$ G_{LT} \f$ | In-plane shear modulus |
 * | props[6] | \f$ \alpha_L \f$ | Longitudinal CTE |
 * | props[7] | \f$ \alpha_T \f$ | Transverse CTE |
 * | props[8] | \f$ Y_{12,0} \f$ | Shear damage threshold |
 * | props[9] | \f$ Y_{12,c} \f$ | Shear damage characteristic force |
 * | props[10] | \f$ Y_{22,0} \f$ | Transverse damage threshold |
 * | props[11] | \f$ Y_{22,c} \f$ | Transverse damage characteristic force |
 * | props[12] | \f$ Y_{22,u} \f$ | Transverse failure force (brittle) |
 * | props[13] | \f$ b \f$ | Transverse-shear coupling of \f$ Y_{ts} \f$ |
 * | props[14] | \f$ A_{ts} \f$ | Transverse weight of the plasticity criterion |
 * | props[15] | \f$ \sigma_{ts,0} \f$ | Initial yield stress |
 * | props[16] | \f$ m_{ts} \f$ | Hardening exponent |
 * | props[17] | \f$ \beta_{ts} \f$ | Hardening coefficient |
 *
 * **State Variables (statev, 10):**
 * | Index | Symbol | Description |
 * |-------|--------|-------------|
 * | statev[0] | \f$ T_{init} \f$ | Reference temperature |
 * | statev[1] | \f$ d_{22} \f$ | Transverse damage |
 * | statev[2] | \f$ d_{12} \f$ | Shear damage |
 * | statev[3] | \f$ p \f$ | Accumulated plastic strain |
 * | statev[4-9] | \f$ \boldsymbol{\varepsilon}^p \f$ | Plastic strain (Voigt, engineering shears, global frame) |
 *
 * @param umat_name Name of the constitutive model
 * @param Etot Total strain at the beginning of the increment (6)
 * @param DEtot Strain increment (6)
 * @param stress Stress [in: start; out: end] (6)
 * @param Lt Tangent modulus [output] (6x6)
 * @param L Damaged elastic stiffness \f$ \tilde{\mathbf{L}} \f$ [output] (6x6)
 * @param DR Rotation increment (3x3)
 * @param nprops Number of material properties (18)
 * @param props Material properties (table above)
 * @param nstatev Number of state variables (10)
 * @param statev State variables [input/output] (table above)
 * @param T Temperature at the beginning of the increment
 * @param DT Temperature increment
 * @param Time Time at the beginning of the increment
 * @param DTime Time increment
 * @param Wm Total mechanical work [input/output]
 * @param Wm_r Recoverable work [input/output]
 * @param Wm_ir Stored (hardening) work [input/output]
 * @param Wm_d Dissipated work (plasticity and damage) [input/output]
 * @param ndi Number of direct stress components
 * @param nshr Number of shear stress components
 * @param start First increment flag
 * @param tnew_dt Suggested time step ratio [output] (0.2 when a local loop does not converge)
 * @param tangent_mode Not used: continuum tangent (see above)
 *
 * @code
 * // Carbon/epoxy (T300/914-like), fibres along 1
 * vec props = {1, 138000., 11000., 0.0223, 0.40, 5500., -0.3e-6, 30e-6,
 *              0.10, 0.30, 0.15, 0.40, 1.20, 2.5, 2.0, 50., 0.5, 300.};
 * vec statev = zeros(10);
 * @endcode
 *
 * **References:**
 * - Ladevèze, P., & Le Dantec, E. (1992). "Damage modelling of the elementary ply for laminated composites." *Composites Science and Technology*, 43(3), 257-267.
 * - Allix, O., & Ladevèze, P. (1992). "Interlaminar interface modelling for the prediction of delamination." *Composite Structures*, 22(4), 235-242.
 * - Ladevèze, P. (1992). "A damage computational method for composite structures." *Computers & Structures*, 44(1-2), 79-87.
 * - Matzenmiller, A., Lubliner, J., & Taylor, R. L. (1995). "A constitutive model for anisotropic damage in fiber-composites." *Mechanics of Materials*, 20(2), 125-152.
 * - Pinho, S. T., et al. (2012). "Material and structural response of polymer-matrix fibre-reinforced composites." *Journal of Composite Materials*, 46(19-20), 2313-2341.
 */
void umat_damage_LLD_0(const std::string &umat_name, const arma::vec &Etot, const arma::vec &DEtot, arma::vec &stress, arma::mat &Lt, arma::mat &L, const arma::mat &DR, const int &nprops, const arma::vec &props, const int &nstatev, arma::vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, double &tnew_dt, const int &tangent_mode = tangent_default);


/** @} */ // end of umat_mechanical group

} //namespace simcoon
