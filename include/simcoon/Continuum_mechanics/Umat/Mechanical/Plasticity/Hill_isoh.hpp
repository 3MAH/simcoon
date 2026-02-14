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
 * @file Hill_isoh.hpp
 * @brief Elastic-plastic material with Hill anisotropic yield criterion and isotropic hardening using CCP algorithm
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
 * @brief Elastic-plastic constitutive model with Hill anisotropic yield criterion and isotropic hardening solved by the Convex Cutting Plane (CCP) algorithm
 *
 * @details This function implements the Hill anisotropic plasticity model for small and finite strain analysis
 * of orthotropic materials. The model features:
 * - Hill (1948) quadratic anisotropic yield criterion
 * - Isotropic hardening with power law
 * - Associative flow rule
 * - Convex Cutting Plane algorithm for return mapping
 * - Thermal expansion effects
 * - Consistent tangent modulus for implicit FE analysis
 *
 * **Constitutive Equations:**
 *
 * The Hill anisotropic yield function is defined as:
 * \f[
 * \Phi(\boldsymbol{\sigma}, p) = \sigma_{eq}^{Hill} - H_p(p) - \sigma_Y \leq 0
 * \f]
 * where:
 * \f[
 * \sigma_{eq}^{Hill} = \sqrt{F(\sigma_{22} - \sigma_{33})^2 + G(\sigma_{33} - \sigma_{11})^2 + H(\sigma_{11} - \sigma_{22})^2 + 2L\sigma_{23}^2 + 2M\sigma_{13}^2 + 2N\sigma_{12}^2}
 * \f]
 *
 * **Hill Anisotropy Parameters:**
 *
 * The six Hill parameters \f$ (F, G, H, L, M, N) \f$ define the material's plastic anisotropy.
 * For orthotropic symmetry with principal axes aligned with material axes:
 * - \f$ F, G, H \f$ control yielding in normal stress components
 * - \f$ L, M, N \f$ control yielding in shear stress components
 *
 * **Relationship to Yield Stress Ratios:**
 *
 * The Hill parameters can be related to yield stresses in different directions:
 * \f[
 * F = \frac{1}{2}\left(\frac{1}{\sigma_{Y22}^2} + \frac{1}{\sigma_{Y33}^2} - \frac{1}{\sigma_{Y11}^2}\right)
 * \f]
 * \f[
 * G = \frac{1}{2}\left(\frac{1}{\sigma_{Y33}^2} + \frac{1}{\sigma_{Y11}^2} - \frac{1}{\sigma_{Y22}^2}\right)
 * \f]
 * \f[
 * H = \frac{1}{2}\left(\frac{1}{\sigma_{Y11}^2} + \frac{1}{\sigma_{Y22}^2} - \frac{1}{\sigma_{Y33}^2}\right)
 * \f]
 * \f[
 * L = \frac{3}{2\sigma_{Y23}^2}, \quad M = \frac{3}{2\sigma_{Y13}^2}, \quad N = \frac{3}{2\sigma_{Y12}^2}
 * \f]
 * where \f$ \sigma_{Yii} \f$ is the yield stress in uniaxial loading along direction \f$ i \f$, and
 * \f$ \sigma_{Yij} \f$ is the yield stress in pure shear in the \f$ ij \f$ plane.
 *
 * **Special Case - von Mises (Isotropic):**
 *
 * For isotropic materials, the Hill criterion reduces to von Mises when:
 * \f[
 * F = G = H = \frac{1}{2}, \quad L = M = N = \frac{3}{2}
 * \f]
 *
 * **Isotropic Hardening (Power Law):**
 *
 * The isotropic hardening follows a power law:
 * \f[
 * H_p(p) = k \cdot p^m
 * \f]
 * where:
 * - \f$ k \f$ is the hardening coefficient
 * - \f$ m \f$ is the hardening exponent
 * - \f$ p = \int_0^t \dot{p} \, dt \f$ is the accumulated plastic strain
 *
 * **Plastic Flow Rule:**
 *
 * The plastic strain rate follows the associative flow rule:
 * \f[
 * \dot{\boldsymbol{\varepsilon}}^p = \dot{p} \mathbf{n}
 * \f]
 * where:
 * \f[
 * \mathbf{n} = \frac{\partial \sigma_{eq}^{Hill}}{\partial \boldsymbol{\sigma}} = \frac{1}{\sigma_{eq}^{Hill}} \mathbf{P}^{Hill} : \boldsymbol{\sigma}
 * \f]
 * is the flow direction normal to the Hill yield surface, and \f$ \mathbf{P}^{Hill} \f$ is the Hill projection tensor.
 *
 * **Hill Projection Tensor (Voigt Notation):**
 *
 * \f[
 * \mathbf{P}^{Hill} = \begin{bmatrix}
 * G+H & -H & -G & 0 & 0 & 0 \\
 * -H & H+F & -F & 0 & 0 & 0 \\
 * -G & -F & F+G & 0 & 0 & 0 \\
 * 0 & 0 & 0 & 2N & 0 & 0 \\
 * 0 & 0 & 0 & 0 & 2M & 0 \\
 * 0 & 0 & 0 & 0 & 0 & 2L
 * \end{bmatrix}
 * \f]
 *
 * **Incremental Form (CCP Algorithm):**
 *
 * For an increment \f$ \Delta p \f$, the stress is updated as:
 * \f[
 * \boldsymbol{\sigma}_{n+1} = \boldsymbol{\sigma}^{trial} - \Delta p \left( \mathbf{L} : \mathbf{n} + \frac{\partial \mathbf{n}}{\partial \boldsymbol{\sigma}} : \boldsymbol{\sigma} \right)
 * \f]
 *
 * The CCP algorithm solves for \f$ \Delta p \f$ using the Fischer-Burmeister complementarity condition:
 * \f[
 * FB(\Delta p, \Phi) = \sqrt{(\Delta p)^2 + \Phi^2} - \Delta p - \Phi = 0
 * \f]
 *
 * **Consistent Tangent Modulus:**
 *
 * The algorithmic tangent modulus is:
 * \f[
 * \mathbf{L}_t = \mathbf{L} - \frac{(\mathbf{L}:\mathbf{n}) \otimes (\mathbf{n}:\mathbf{L})}{\mathbf{n}:\mathbf{L}:\mathbf{n} + H'}
 * \f]
 * where \f$ H' = \frac{dH_p}{dp} = k \cdot m \cdot p^{m-1} \f$ is the hardening modulus.
 *
 * **Material Parameters (props):**
 *
 * | Index | Symbol | Description | Units | Typical Range |
 * |-------|--------|-------------|-------|---------------|
 * | props[0] | \f$ E \f$ | Young's modulus | Stress | 50-500 GPa |
 * | props[1] | \f$ \nu \f$ | Poisson's ratio | - | 0.2-0.45 |
 * | props[2] | \f$ \alpha \f$ | Coefficient of thermal expansion | 1/Temperature | 1e-6 to 1e-4 /K |
 * | props[3] | \f$ \sigma_Y \f$ | Initial yield stress | Stress | 100-1000 MPa |
 * | props[4] | \f$ k \f$ | Hardening coefficient | Stress | 100-2000 MPa |
 * | props[5] | \f$ m \f$ | Hardening exponent | - | 0.05-0.5 |
 * | props[6] | \f$ F \f$ | Hill parameter F | 1/Stress² | 0-10 |
 * | props[7] | \f$ G \f$ | Hill parameter G | 1/Stress² | 0-10 |
 * | props[8] | \f$ H \f$ | Hill parameter H | 1/Stress² | 0-10 |
 * | props[9] | \f$ L \f$ | Hill parameter L | 1/Stress² | 0-10 |
 * | props[10] | \f$ M \f$ | Hill parameter M | 1/Stress² | 0-10 |
 * | props[11] | \f$ N \f$ | Hill parameter N | 1/Stress² | 0-10 |
 *
 * **Constraint on Hill Parameters:**
 * - For physical consistency: \f$ F + G + H > 0 \f$
 * - For von Mises isotropy: \f$ F = G = H = 0.5, L = M = N = 1.5 \f$
 * - For plane stress: Ensure appropriate ratios for the stress state
 *
 * **State Variables (statev):**
 *
 * Total state variables required: \f$ n_{statev} = 8 \f$
 *
 * | Index | Symbol | Description | Units |
 * |-------|--------|-------------|-------|
 * | statev[0] | \f$ T_{init} \f$ | Initial/reference temperature | Temperature |
 * | statev[1] | \f$ p \f$ | Accumulated plastic strain | Strain |
 * | statev[2] | \f$ \varepsilon^p_{11} \f$ | Plastic strain component 11 | Strain |
 * | statev[3] | \f$ \varepsilon^p_{22} \f$ | Plastic strain component 22 | Strain |
 * | statev[4] | \f$ \varepsilon^p_{33} \f$ | Plastic strain component 33 | Strain |
 * | statev[5] | \f$ \varepsilon^p_{12} \f$ | Plastic strain component 12 (engineering) | Strain |
 * | statev[6] | \f$ \varepsilon^p_{13} \f$ | Plastic strain component 13 (engineering) | Strain |
 * | statev[7] | \f$ \varepsilon^p_{23} \f$ | Plastic strain component 23 (engineering) | Strain |
 *
 * @param Etot Total strain tensor at beginning of increment (Voigt notation: 6×1 vector)
 * @param DEtot Strain increment tensor (Voigt notation: 6×1 vector)
 * @param sigma Stress tensor (Voigt notation: 6×1 vector) [output]
 * @param Lt Consistent tangent modulus \f$ \mathbf{L}_t \f$ (6×6 matrix) [output]
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
 * @note The Hill criterion is widely used for metals with crystallographic texture
 * @note Common applications: rolled sheets, drawn wires, forged components
 * @note Parameter identification requires yield stress measurements in multiple directions
 * @note For fiber composites, consider more advanced criteria (Tsai-Wu, Hoffman)
 * @note The model assumes plastic incompressibility (plastic volume change = 0)
 * @note Material axes must align with the global coordinate system or use rotation
 *
 * @see L_iso() for isotropic elastic stiffness construction
 * @see Ireal() for real identity tensor
 * @see Fischer_Burmeister() for complementarity solver
 * @see Eq_stress_Hill() for Hill equivalent stress computation
 * @see denom_FB_Hill() for CCP denominator with Hill criterion
 *
 * @code
 * // Example usage: Aluminum rolled sheet (transverse isotropy in 1-2 plane)
 * vec props(12);
 * props(0) = 70000;       // E = 70 GPa
 * props(1) = 0.33;        // nu = 0.33
 * props(2) = 2.3e-5;      // alpha = 23e-6 /K
 * props(3) = 200;         // sigma_Y = 200 MPa
 * props(4) = 400;         // k = 400 MPa
 * props(5) = 0.15;        // m = 0.15
 *
 * // Hill parameters for rolled sheet (stronger in rolling direction)
 * // Assume: sigma_Y11 = 200 MPa (rolling), sigma_Y22 = 220 MPa (transverse), sigma_Y33 = 210 MPa
 * // sigma_Y12 = sigma_Y13 = sigma_Y23 = 115.5 MPa (shear, sqrt(3)*sigma_Y/3)
 * props(6) = 0.00115;     // F
 * props(7) = 0.00116;     // G
 * props(8) = 0.00126;     // H
 * props(9) = 0.0130;      // L = 1.5 for isotropy
 * props(10) = 0.0130;     // M = 1.5 for isotropy
 * props(11) = 0.0130;     // N = 1.5 for isotropy
 *
 * vec statev = zeros(8);
 * statev(0) = 20.0;  // Reference temperature 20°C
 *
 * vec Etot = {0.003, 0.0, 0.0, 0.0, 0.0, 0.0};  // 0.3% uniaxial strain
 * vec DEtot = {0.0001, 0.0, 0.0, 0.0, 0.0, 0.0};
 * vec sigma = zeros(6);
 * mat Lt = zeros(6,6);
 * mat DR = eye(3,3);
 *
 * umat_plasticity_hill_isoh_CCP(Etot, DEtot, sigma, Lt, DR,
 *                               12, props, 8, statev, 20.0, 0.0, 0.0, 1.0,
 *                               Wm, Wm_r, Wm_ir, Wm_d, 3, 3, false, tnew_dt);
 *
 * // Expected: Different yield stress depending on loading direction
 * @endcode
 *
 * **References:**
 * - Hill, R. (1948). "A theory of the yielding and plastic flow of anisotropic metals." *Proceedings of the Royal Society of London A*, 193(1033), 281-297.
 * - Hill, R. (1950). *The Mathematical Theory of Plasticity*. Oxford University Press.
 * - Barlat, F., & Lian, J. (1989). "Plastic behavior and stretchability of sheet metals. Part I: A yield function for orthotropic sheets under plane stress conditions." *International Journal of Plasticity*, 5(1), 51-66.
 * - Banabic, D. (2010). *Sheet Metal Forming Processes: Constitutive Modelling and Numerical Simulation*. Springer.
 * - Ortiz, M., & Simo, J. C. (1986). "An analysis of a new class of integration algorithms for elastoplastic constitutive relations." *International Journal for Numerical Methods in Engineering*, 23(3), 353-366.
 */
void umat_plasticity_hill_isoh_CCP(const std::string &umat_name, const arma::vec &Etot, const arma::vec &DEtot, arma::vec &sigma, arma::mat &Lt, arma::mat &L, const arma::mat &DR, const int &nprops, const arma::vec &props, const int &nstatev, arma::vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, double &tnew_dt);

/** @} */ // end of umat_mechanical group
    
} //namespace smart
