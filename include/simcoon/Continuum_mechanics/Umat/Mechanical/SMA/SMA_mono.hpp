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
 * @file SMA_mono.hpp
 * @brief Unified thermomechanical constitutive model for shape memory alloys (SMAs)
 * @author Y. Chemisky, D. Chatziathanasiou, G. Chatzigeorgiou, F. Meraghni (LEM3-SMART group)
 * @version 1.0
 */

#pragma once
#include <armadillo>

namespace simcoon{

/** @addtogroup umat_mechanical
 *  @{
 */

/**
 * @brief Unified thermomechanical constitutive model for monocrystalline shape memory alloys (SMAs)
 *
 * @details This function implements a comprehensive SMA model capturing both pseudoelastic and
 * shape memory effects through stress- and temperature-induced martensitic phase transformation.
 * The model features:
 * - Thermomechanical coupling between stress, temperature, and phase transformation
 * - Martensitic volume fraction (MVF) as internal variable
 * - Transformation strain evolution with stress-dependent magnitude
 * - Forward (A→M) and reverse (M→A) transformations with distinct kinetics
 * - Smooth transformation surfaces with exponential hardening
 * - Penalty method for enforcing MVF bounds [0,1]
 * - Temperature-dependent elastic properties
 *
 * **Phase Transformation Framework:**
 *
 * SMAs undergo a diffusionless solid-state phase transformation between two crystal structures:
 * - **Austenite (A)**: High-temperature parent phase, cubic structure, higher stiffness
 * - **Martensite (M)**: Low-temperature product phase, lower symmetry, accommodates large strains
 *
 * **Martensitic Volume Fraction (MVF):**
 *
 * The material state is characterized by the martensitic volume fraction \f$ \xi \in [0,1] \f$:
 * - \f$ \xi = 0 \f$: Pure austenite
 * - \f$ 0 < \xi < 1 \f$: Mixed phase (austenite + martensite)
 * - \f$ \xi = 1 \f$: Pure martensite
 *
 * **Transformation Strain:**
 *
 * The total strain decomposes as:
 * \f[
 * \boldsymbol{\varepsilon} = \boldsymbol{\varepsilon}^{el} + \boldsymbol{\varepsilon}^{th} + \boldsymbol{\varepsilon}^{tr}
 * \f]
 * where:
 * - \f$ \boldsymbol{\varepsilon}^{el} \f$ is the elastic strain
 * - \f$ \boldsymbol{\varepsilon}^{th} \f$ is the thermal expansion strain
 * - \f$ \boldsymbol{\varepsilon}^{tr} \f$ is the transformation strain
 *
 * **Transformation Strain Magnitude:**
 *
 * The maximum transformation strain varies with stress level:
 * \f[
 * H(\sigma_{eq}) = H_{min} + (H_{max} - H_{min}) \left( 1 - e^{-k_1 (\sigma_{eq} - \sigma_{crit})} \right) \quad \text{for } \sigma_{eq} > \sigma_{crit}
 * \f]
 * where:
 * - \f$ H_{min} \f$ is the minimum transformation strain (low stress)
 * - \f$ H_{max} \f$ is the maximum transformation strain (high stress)
 * - \f$ k_1 \f$ controls the exponential transition rate
 * - \f$ \sigma_{crit} \f$ is the critical stress for strain magnitude evolution
 * - \f$ \sigma_{eq} \f$ is the von Mises equivalent stress
 *
 * **Transformation Surfaces:**
 *
 * **Forward Transformation (Austenite → Martensite):**
 * \f[
 * \Phi^{A \to M} = \sigma_{eq} - C_M (T - M_s) - Y_0^t - a_1 \xi^F \leq 0
 * \f]
 * where:
 * - \f$ C_M \f$ is the Clausius-Clapeyron slope for forward transformation
 * - \f$ M_s \f$ is the martensite start temperature at zero stress
 * - \f$ Y_0^t \f$ is the initial transformation threshold
 * - \f$ a_1 \f$ is the forward transformation hardening parameter
 * - \f$ \xi^F \f$ is the forward transformation MVF contribution
 *
 * **Reverse Transformation (Martensite → Austenite):**
 * \f[
 * \Phi^{M \to A} = -\sigma_{eq} + C_A (T - A_f) - Y_0^t - a_2 \xi^R \leq 0
 * \f]
 * where:
 * - \f$ C_A \f$ is the Clausius-Clapeyron slope for reverse transformation
 * - \f$ A_f \f$ is the austenite finish temperature at zero stress
 * - \f$ a_2 \f$ is the reverse transformation hardening parameter
 * - \f$ \xi^R \f$ is the reverse transformation MVF contribution
 *
 * **Critical Temperatures:**
 *
 * Four characteristic temperatures define the transformation behavior at zero stress:
 * - \f$ M_s \f$: Martensite start (cooling begins austenite → martensite)
 * - \f$ M_f \f$: Martensite finish (cooling completes austenite → martensite)
 * - \f$ A_s \f$: Austenite start (heating begins martensite → austenite)
 * - \f$ A_f \f$: Austenite finish (heating completes martensite → austenite)
 *
 * Typically: \f$ M_f < M_s < A_s < A_f \f$
 *
 * **Stress-Temperature Relationship:**
 *
 * The transformation temperatures shift linearly with applied stress:
 * \f[
 * M_s(\sigma) = M_s^0 + \frac{\sigma_{eq}}{C_M}, \quad A_f(\sigma) = A_f^0 + \frac{\sigma_{eq}}{C_A}
 * \f]
 *
 * **Mixed Elastic Properties:**
 *
 * Elastic properties vary linearly with MVF (rule of mixtures):
 * \f[
 * E(\xi) = E_A (1 - \xi) + E_M \xi
 * \f]
 * \f[
 * \nu(\xi) = \nu_A (1 - \xi) + \nu_M \xi
 * \f]
 * \f[
 * \alpha(\xi) = \alpha_A (1 - \xi) + \alpha_M \xi
 * \f]
 *
 * **Penalty Method for MVF Constraints:**
 *
 * To enforce \f$ 0 \leq \xi \leq 1 \f$, a penalty function is applied:
 * \f[
 * P(\xi) = p_0 \left[ \left( \frac{\xi - c}{1 - c} \right)^n + \left( \frac{c - \xi}{c} \right)^n \right]^\alpha
 * \f]
 * where \f$ c, p_0, n, \alpha \f$ are penalty parameters.
 *
 * **Material Parameters (props):**
 *
 * | Index | Symbol | Description | Units | Typical Range |
 * |-------|--------|-------------|-------|---------------|
 * | props[0] | flagT | Temperature smoothing flag (0=linear, 1=smooth) | - | 0 or 1 |
 * | props[1] | \f$ E_A \f$ | Young's modulus of austenite | Stress | 40-100 GPa |
 * | props[2] | \f$ E_M \f$ | Young's modulus of martensite | Stress | 20-60 GPa |
 * | props[3] | \f$ \nu_A \f$ | Poisson's ratio of austenite | - | 0.3-0.4 |
 * | props[4] | \f$ \nu_M \f$ | Poisson's ratio of martensite | - | 0.3-0.4 |
 * | props[5] | \f$ \alpha_A \f$ | CTE of austenite | 1/Temperature | 10-20 e-6 /K |
 * | props[6] | \f$ \alpha_M \f$ | CTE of martensite | 1/Temperature | 5-15 e-6 /K |
 * | props[7] | \f$ H_{min} \f$ | Minimum transformation strain magnitude | Strain | 0.01-0.03 |
 * | props[8] | \f$ H_{max} \f$ | Maximum transformation strain magnitude | Strain | 0.05-0.10 |
 * | props[9] | \f$ k_1 \f$ | Exponential rate for H evolution | 1/Stress | 0.001-0.01 |
 * | props[10] | \f$ \sigma_{crit} \f$ | Critical stress for H evolution | Stress | 50-200 MPa |
 * | props[11] | \f$ C_A \f$ | Clausius-Clapeyron slope (M→A) | Stress/Temperature | 5-15 MPa/K |
 * | props[12] | \f$ C_M \f$ | Clausius-Clapeyron slope (A→M) | Stress/Temperature | 5-15 MPa/K |
 * | props[13] | \f$ M_s \f$ | Martensite start temperature | Temperature | -50 to 100°C |
 * | props[14] | \f$ M_f \f$ | Martensite finish temperature | Temperature | -70 to 80°C |
 * | props[15] | \f$ A_s \f$ | Austenite start temperature | Temperature | -30 to 120°C |
 * | props[16] | \f$ A_f \f$ | Austenite finish temperature | Temperature | -10 to 140°C |
 * | props[17] | \f$ n_1 \f$ | Martensite start smooth exponent | - | 1-5 |
 * | props[18] | \f$ n_2 \f$ | Martensite finish smooth exponent | - | 1-5 |
 * | props[19] | \f$ n_3 \f$ | Austenite start smooth exponent | - | 1-5 |
 * | props[20] | \f$ n_4 \f$ | Austenite finish smooth exponent | - | 1-5 |
 * | props[21] | \f$ c_\lambda \f$ | Penalty function center point | - | 0.5 |
 * | props[22] | \f$ p_{0\lambda} \f$ | Penalty function magnitude | Stress | 1000-10000 |
 * | props[23] | \f$ n_\lambda \f$ | Penalty function power exponent | - | 2-10 |
 * | props[24] | \f$ \alpha_\lambda \f$ | Penalty function parameter | - | 0.5-2 |
 *
 * **State Variables (statev):**
 *
 * Total state variables required: \f$ n_{statev} = 17 \f$
 *
 * | Index | Symbol | Description | Units |
 * |-------|--------|-------------|-------|
 * | statev[0] | \f$ T_{init} \f$ | Initial/reference temperature | Temperature |
 * | statev[1] | \f$ \xi \f$ | Martensitic volume fraction | - |
 * | statev[2] | \f$ \varepsilon^{tr}_{11} \f$ | Transformation strain component 11 | Strain |
 * | statev[3] | \f$ \varepsilon^{tr}_{22} \f$ | Transformation strain component 22 | Strain |
 * | statev[4] | \f$ \varepsilon^{tr}_{33} \f$ | Transformation strain component 33 | Strain |
 * | statev[5] | \f$ \varepsilon^{tr}_{12} \f$ | Transformation strain component 12 (engineering) | Strain |
 * | statev[6] | \f$ \varepsilon^{tr}_{13} \f$ | Transformation strain component 13 (engineering) | Strain |
 * | statev[7] | \f$ \varepsilon^{tr}_{23} \f$ | Transformation strain component 23 (engineering) | Strain |
 * | statev[8] | \f$ \xi^F \f$ | Forward transformation MVF | - |
 * | statev[9] | \f$ \xi^R \f$ | Reverse transformation MVF | - |
 * | statev[10] | \f$ \Delta s_0 \f$ | Entropy difference (M - A) | Energy/Temperature |
 * | statev[11] | \f$ \Delta u_0 \f$ | Internal energy difference (M - A) | Energy |
 * | statev[12] | - | Stress dependence parameter | - |
 * | statev[13] | \f$ a_1 \f$ | Forward hardening parameter | Stress |
 * | statev[14] | \f$ a_2 \f$ | Reverse hardening parameter | Stress |
 * | statev[15] | \f$ a_3 \f$ | Equilibrium hardening parameter | Stress |
 * | statev[16] | \f$ Y_0^t \f$ | Initial transformation threshold | Stress |
 *
 * @param Etot Total strain tensor at beginning of increment (Voigt notation: 6×1 vector)
 * @param DEtot Strain increment tensor (Voigt notation: 6×1 vector)
 * @param sigma Cauchy stress tensor (Voigt notation: 6×1 vector) [output]
 * @param Lt Consistent tangent modulus (6×6 matrix) [output]
 * @param L Mixed elastic stiffness tensor (6×6 matrix) [output]
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
 * @param Wm_ir Irrecoverable work stored in transformation [output]
 * @param Wm_d Dissipated work (hysteresis) [output]
 * @param ndi Number of direct stress components (typically 3)
 * @param nshr Number of shear stress components (typically 3)
 * @param start Flag indicating first increment (true) or continuation (false)
 * @param tnew_dt Suggested new time step size for adaptive time stepping [output]
 *
 * @note The model captures the characteristic SMA behaviors:
 * @note - **Pseudoelasticity (superelasticity)**: Stress-induced transformation at T > Af
 * @note - **Shape memory effect**: Temperature-induced recovery of transformation strain
 * @note - **Transformation hysteresis**: Energy dissipation during cycling
 * @note - **Orientation effects**: Transformation strain aligns with stress direction
 * @note Typical SMA materials: NiTi (Nitinol), CuAlNi, CuZnAl, FeMnSi
 * @note The model requires experimental characterization of transformation temperatures
 * @note Parameter identification should include DSC, tensile tests at various temperatures
 * @note Convergence requires small load steps during phase transformation
 *
 * @see L_iso() for isotropic elastic stiffness construction
 * @see Lagrange_exp() for penalty function implementation
 * @see eta_stress() for effective stress computation
 * @see Mises_stress() for von Mises equivalent stress
 *
 * @code
 * // Example usage: NiTi SMA wire (pseudoelastic regime at room temperature)
 * vec props(25);
 * props(0) = 0;           // flagT = 0 (linear extrapolation)
 * props(1) = 70000;       // EA = 70 GPa (austenite)
 * props(2) = 40000;       // EM = 40 GPa (martensite)
 * props(3) = 0.33;        // nuA = 0.33
 * props(4) = 0.33;        // nuM = 0.33
 * props(5) = 1.1e-5;      // alphaA = 11e-6 /K
 * props(6) = 6.6e-6;      // alphaM = 6.6e-6 /K
 * props(7) = 0.02;        // Hmin = 2% (min transformation strain)
 * props(8) = 0.067;       // Hmax = 6.7% (max transformation strain)
 * props(9) = 0.005;       // k1 = 0.005
 * props(10) = 100;        // sigmacrit = 100 MPa
 * props(11) = 6.5;        // CA = 6.5 MPa/K (reverse transformation)
 * props(12) = 7.0;        // CM = 7.0 MPa/K (forward transformation)
 * props(13) = -20;        // Ms = -20°C
 * props(14) = -30;        // Mf = -30°C
 * props(15) = 10;         // As = 10°C
 * props(16) = 20;         // Af = 20°C
 * props(17) = 2;          // n1, n2, n3, n4 = 2 (smooth exponents)
 * props(18) = 2;
 * props(19) = 2;
 * props(20) = 2;
 * props(21) = 0.5;        // c_lambda = 0.5
 * props(22) = 5000;       // p0_lambda = 5000
 * props(23) = 5;          // n_lambda = 5
 * props(24) = 1.0;        // alpha_lambda = 1.0
 *
 * vec statev = zeros(17);
 * statev(0) = 25.0;  // Reference temperature 25°C (above Af, pseudoelastic regime)
 *
 * vec Etot = {0.05, -0.015, -0.015, 0.0, 0.0, 0.0};  // 5% axial strain (induces transformation)
 * vec DEtot = {0.001, -0.0003, -0.0003, 0.0, 0.0, 0.0};
 * vec sigma = zeros(6);
 * mat Lt = zeros(6,6);
 * mat L = zeros(6,6);
 * mat DR = eye(3,3);
 *
 * umat_sma_mono(Etot, DEtot, sigma, Lt, L, DR,
 *               25, props, 17, statev, 25.0, 0.0, 0.0, 1.0,
 *               Wm, Wm_r, Wm_ir, Wm_d, 3, 3, false, tnew_dt);
 *
 * // Check transformation state
 * double xi = statev(1);  // Martensitic volume fraction
 * cout << "MVF: " << xi << " (0=austenite, 1=martensite)" << endl;
 * @endcode
 *
 * **References:**
 * - Chemisky, Y., et al. (2011). "Constitutive model for shape memory alloys including phase transformation, martensitic reorientation and twins accommodation." *Mechanics of Materials*, 43(7), 361-376.
 * - Chatziathanasiou, D., et al. (2016). "Modeling of coupled phase transformation and reorientation in shape memory alloys under non-proportional thermomechanical loading." *International Journal of Plasticity*, 82, 192-224.
 * - Lagoudas, D. C. (Ed.). (2008). *Shape Memory Alloys: Modeling and Engineering Applications*. Springer.
 * - Boyd, J. G., & Lagoudas, D. C. (1996). "A thermodynamical constitutive model for shape memory materials. Part I: The monolithic shape memory alloy." *International Journal of Plasticity*, 12(6), 805-842.
 * - Auricchio, F., & Petrini, L. (2004). "A three-dimensional model describing stress-temperature induced solid phase transformations: solution algorithm and boundary value problems." *International Journal for Numerical Methods in Engineering*, 61(6), 807-836.
 */
void umat_sma_mono(const arma::vec &Etot, const arma::vec &DEtot, arma::vec &sigma, arma::mat &Lt, arma::mat &L, const arma::mat &DR, const int &nprops, const arma::vec &props, const int &nstatev, arma::vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, double &tnew_dt);

/** @} */ // end of umat_mechanical group
    
} //namespace simcoon
