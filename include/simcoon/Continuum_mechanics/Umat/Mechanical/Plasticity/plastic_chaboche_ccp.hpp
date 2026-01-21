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
 * @file plastic_chaboche_ccp.hpp
 * @brief Elastic-plastic material with Chaboche unified viscoplasticity using CCP algorithm
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
 * @brief Elastic-plastic constitutive model with Chaboche kinematic hardening solved by the Convex Cutting Plane (CCP) algorithm
 *
 * @details This function implements the Chaboche unified viscoplasticity model for small and finite strain analysis.
 * The model features:
 * - J2 (von Mises) plasticity with associative flow rule
 * - Multiple Armstrong-Frederick backstresses for nonlinear kinematic hardening
 * - Optional isotropic hardening with Voce law
 * - Dynamic recovery (recall terms) in backstress evolution
 * - Convex Cutting Plane algorithm for return mapping
 * - Thermal expansion effects
 * - Consistent tangent modulus for implicit FE analysis
 *
 * **Constitutive Equations:**
 *
 * The yield function is defined as:
 * \f[
 * \Phi(\boldsymbol{\sigma}, \mathbf{X}, R) = \sigma_{eq}(\boldsymbol{\sigma} - \mathbf{X}) - R - \sigma_Y \leq 0
 * \f]
 * where:
 * - \f$ \sigma_{eq}(\boldsymbol{\eta}) = \sqrt{\frac{3}{2} \boldsymbol{\eta}_{dev} : \boldsymbol{\eta}_{dev}} \f$ is the von Mises equivalent stress
 * - \f$ \boldsymbol{\eta} = \boldsymbol{\sigma} - \mathbf{X} \f$ is the shifted (effective) stress tensor
 * - \f$ \mathbf{X} = \sum_{i=1}^{N_{kin}} \mathbf{X}_i \f$ is the total backstress (sum of individual backstresses)
 * - \f$ R \f$ is the isotropic hardening stress
 * - \f$ \sigma_Y \f$ is the initial yield stress
 *
 * **Armstrong-Frederick Backstress Evolution:**
 *
 * Each backstress component evolves according to the Armstrong-Frederick law with dynamic recovery:
 * \f[
 * \dot{\mathbf{X}}_i = \frac{2}{3} C_i \dot{\boldsymbol{\varepsilon}}^p - D_i \mathbf{X}_i \dot{p}
 * \f]
 * where:
 * - \f$ C_i \f$ is the kinematic hardening modulus of the i-th backstress
 * - \f$ D_i \f$ is the dynamic recovery (recall) parameter of the i-th backstress
 * - \f$ \dot{\boldsymbol{\varepsilon}}^p \f$ is the plastic strain rate tensor
 * - \f$ \dot{p} = \sqrt{\frac{2}{3} \dot{\boldsymbol{\varepsilon}}^p : \dot{\boldsymbol{\varepsilon}}^p} \f$ is the accumulated plastic strain rate
 *
 * **Physical Interpretation:**
 * - The \f$ C_i \f$ term represents strain hardening (backstress growth)
 * - The \f$ D_i \mathbf{X}_i \dot{p} \f$ term represents dynamic recovery (backstress fading)
 * - Multiple backstresses capture different scales of cyclic hardening behavior
 * - Typically 2-3 backstress components are used
 *
 * **Isotropic Hardening (Optional):**
 *
 * The isotropic hardening follows the Voce law:
 * \f[
 * R(p) = \sum_{j=1}^{N_{iso}} Q_j \left( 1 - e^{-b_j p} \right)
 * \f]
 * where:
 * - \f$ Q_j \f$ is the saturation value of the j-th isotropic hardening component
 * - \f$ b_j \f$ is the hardening rate parameter
 * - \f$ p \f$ is the accumulated plastic strain
 *
 * **Plastic Flow Rule:**
 *
 * The plastic strain rate follows the associative flow rule:
 * \f[
 * \dot{\boldsymbol{\varepsilon}}^p = \dot{p} \mathbf{n}
 * \f]
 * where:
 * \f[
 * \mathbf{n} = \frac{3}{2} \frac{\boldsymbol{\eta}_{dev}}{\sigma_{eq}(\boldsymbol{\eta})}
 * \f]
 * is the flow direction (normal to the yield surface).
 *
 * **Incremental Form (CCP Algorithm):**
 *
 * For an increment \f$ \Delta p \f$, the stress and backstresses are updated as:
 * \f[
 * \boldsymbol{\sigma}_{n+1} = \boldsymbol{\sigma}^{trial} - \Delta p \left( \mathbf{L} : \mathbf{n} + \frac{\partial \mathbf{n}}{\partial \boldsymbol{\sigma}} : \boldsymbol{\sigma} \right)
 * \f]
 * \f[
 * \mathbf{X}_{i,n+1} = \mathbf{X}_{i,n} + \frac{2}{3} C_i \Delta \boldsymbol{\varepsilon}^p - D_i \mathbf{X}_{i,n} \Delta p
 * \f]
 *
 * The CCP algorithm solves for \f$ \Delta p \f$ using the Fischer-Burmeister complementarity condition:
 * \f[
 * FB(\Delta p, \Phi) = \sqrt{(\Delta p)^2 + \Phi^2} - \Delta p - \Phi = 0
 * \f]
 *
 * **Consistent Tangent Modulus:**
 *
 * The algorithmic tangent modulus accounts for both plastic flow and backstress evolution:
 * \f[
 * \mathbf{L}_t = \mathbf{L} - \frac{(\mathbf{L}:\mathbf{n}) \otimes (\mathbf{n}:\mathbf{L})}{\mathbf{n}:\mathbf{L}:\mathbf{n} + H_{tot}}
 * \f]
 * where:
 * \f[
 * H_{tot} = \sum_{i=1}^{N_{kin}} \frac{C_i}{1 + D_i \Delta p} + \sum_{j=1}^{N_{iso}} Q_j b_j e^{-b_j p}
 * \f]
 * is the total hardening modulus combining kinematic and isotropic contributions.
 *
 * **Material Parameters (props):**
 *
 * For \f$ N_{iso} \f$ isotropic hardening terms and \f$ N_{kin} \f$ kinematic hardening terms:
 *
 * | Index | Symbol | Description | Units | Typical Range |
 * |-------|--------|-------------|-------|---------------|
 * | props[0] | \f$ E \f$ | Young's modulus | Stress | 50-500 GPa |
 * | props[1] | \f$ \nu \f$ | Poisson's ratio | - | 0.2-0.45 |
 * | props[2] | \f$ \alpha \f$ | Coefficient of thermal expansion | 1/Temperature | 1e-6 to 1e-4 /K |
 * | props[3] | \f$ \sigma_Y \f$ | Initial yield stress | Stress | 100-1000 MPa |
 * | props[4] | \f$ N_{iso} \f$ | Number of isotropic hardening terms | - | 0-3 |
 * | props[5] | \f$ N_{kin} \f$ | Number of kinematic hardening terms | - | 1-3 |
 * | props[6+2j] | \f$ Q_j \f$ | Saturation stress of j-th isotropic term | Stress | 0-500 MPa |
 * | props[7+2j] | \f$ b_j \f$ | Hardening rate of j-th isotropic term | 1/Strain | 1-100 |
 * | props[6+2N_{iso}+2i] | \f$ C_i \f$ | Kinematic modulus of i-th backstress | Stress | 10-500 GPa |
 * | props[7+2N_{iso}+2i] | \f$ D_i \f$ | Dynamic recovery parameter of i-th backstress | - | 0-100 |
 *
 * **Notes on Parameter Selection:**
 * - For pure kinematic hardening, set \f$ N_{iso} = 0 \f$
 * - For pure isotropic hardening, set \f$ N_{kin} = 0 \f$ (not typical for Chaboche)
 * - First backstress (i=0) typically has high \f$ C_1 \f$, low \f$ D_1 \f$ (long-range backstress)
 * - Second backstress (i=1) typically has moderate \f$ C_2 \f$, moderate \f$ D_2 \f$
 * - Third backstress (i=2) typically has low \f$ C_3 \f$, high \f$ D_3 \f$ (short-range backstress)
 *
 * **State Variables (statev):**
 *
 * Total state variables: \f$ n_{statev} = 1 + 1 + 6 + 6 \times N_{kin} \f$
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
 * | statev[8+6i] | \f$ X_{i,11} \f$ | Backstress i, component 11 | Stress |
 * | statev[9+6i] | \f$ X_{i,22} \f$ | Backstress i, component 22 | Stress |
 * | statev[10+6i] | \f$ X_{i,33} \f$ | Backstress i, component 33 | Stress |
 * | statev[11+6i] | \f$ X_{i,12} \f$ | Backstress i, component 12 | Stress |
 * | statev[12+6i] | \f$ X_{i,13} \f$ | Backstress i, component 13 | Stress |
 * | statev[13+6i] | \f$ X_{i,23} \f$ | Backstress i, component 23 | Stress |
 *
 * For \f$ N_{kin} = 2 \f$ backstresses: \f$ n_{statev} = 1 + 1 + 6 + 12 = 20 \f$
 *
 * @param Etot Total strain tensor at beginning of increment (Voigt notation: 6×1 vector)
 * @param DEtot Strain increment tensor (Voigt notation: 6×1 vector)
 * @param sigma Stress tensor (Voigt notation: 6×1 vector) [output]
 * @param Lt Consistent tangent modulus \f$ \mathbf{L}_t \f$ (6×6 matrix) [output]
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
 * @param Wm_ir Irrecoverable work stored in plastic deformation [output]
 * @param Wm_d Dissipated (plastic) work [output]
 * @param ndi Number of direct stress components (typically 3)
 * @param nshr Number of shear stress components (typically 3)
 * @param start Flag indicating first increment (true) or continuation (false)
 * @param solver_type Solver type: 0=implicit, 1=explicit, 2=dynamic implicit
 * @param tnew_dt Suggested new time step size for adaptive time stepping [output]
 *
 * @note The Chaboche model excels at capturing ratcheting and cyclic plasticity
 * @note Multiple backstresses are essential for accurate multiaxial loading predictions
 * @note Parameter identification typically requires cyclic test data (tension-compression, torsion)
 * @note For monotonic loading, simpler isotropic hardening may suffice
 * @note The model assumes small strains (< 10%); use updated Lagrangian for larger strains
 * @note Dynamic recovery (D > 0) prevents unbounded backstress growth under cycling
 *
 * @see Ireal() for real identity tensor (Voigt 6×6)
 * @see Idev() for deviatoric projection tensor
 * @see Fischer_Burmeister() for complementarity solver
 * @see eta_stress() for shifted stress computation
 * @see denom_FB_N_Mises() for CCP denominator
 *
 * @code
 * // Example usage: 316 stainless steel with 2 backstresses, 1 isotropic hardening
 * int N_iso = 1;
 * int N_kin = 2;
 * vec props(6 + 2*N_iso + 2*N_kin);
 * props(0) = 200000;      // E = 200 GPa
 * props(1) = 0.3;         // nu = 0.3
 * props(2) = 1.7e-5;      // alpha = 17e-6 /K
 * props(3) = 150;         // sigma_Y = 150 MPa
 * props(4) = N_iso;       // 1 isotropic hardening term
 * props(5) = N_kin;       // 2 kinematic hardening terms
 *
 * // Isotropic hardening (Voce)
 * props(6) = 100;         // Q1 = 100 MPa (saturation stress)
 * props(7) = 10;          // b1 = 10 (hardening rate)
 *
 * // First backstress (long-range)
 * props(8) = 300000;      // C1 = 300 GPa
 * props(9) = 1000;        // D1 = 1000 (moderate recovery)
 *
 * // Second backstress (short-range)
 * props(10) = 50000;      // C2 = 50 GPa
 * props(11) = 100;        // D2 = 100 (strong recovery)
 *
 * vec statev = zeros(1 + 1 + 6 + 6*N_kin);  // 20 state variables
 * statev(0) = 20.0;  // Reference temperature 20°C
 *
 * vec Etot = {0.002, -0.0006, -0.0006, 0.0, 0.0, 0.0};  // 0.2% axial strain
 * vec DEtot = {0.0001, -0.00003, -0.00003, 0.0, 0.0, 0.0};
 * vec sigma = zeros(6);
 * mat Lt = zeros(6,6);
 * mat L = zeros(6,6);
 * vec sigma_in = zeros(6);
 * mat DR = eye(3,3);
 *
 * umat_plasticity_chaboche_CCP(Etot, DEtot, sigma, Lt, L, sigma_in, DR,
 *                              12, props, 20, statev, 20.0, 0.0, 0.0, 1.0,
 *                              Wm, Wm_r, Wm_ir, Wm_d, 3, 3, false, 0, tnew_dt);
 *
 * // Check backstress components
 * vec X1 = statev.subvec(8, 13);   // First backstress
 * vec X2 = statev.subvec(14, 19);  // Second backstress
 * @endcode
 *
 * **References:**
 * - Chaboche, J. L. (1986). "Time-independent constitutive theories for cyclic plasticity." *International Journal of Plasticity*, 2(2), 149-188.
 * - Chaboche, J. L. (1989). "Constitutive equations for cyclic plasticity and cyclic viscoplasticity." *International Journal of Plasticity*, 5(3), 247-302.
 * - Chaboche, J. L. (2008). "A review of some plasticity and viscoplasticity constitutive theories." *International Journal of Plasticity*, 24(10), 1642-1693.
 * - Armstrong, P. J., & Frederick, C. O. (1966). "A mathematical representation of the multiaxial Bauschinger effect." *CEGB Report RD/B/N731*.
 * - Lemaitre, J., & Chaboche, J. L. (1990). *Mechanics of Solid Materials*. Cambridge University Press.
 * - Ortiz, M., & Simo, J. C. (1986). "An analysis of a new class of integration algorithms for elastoplastic constitutive relations." *International Journal for Numerical Methods in Engineering*, 23(3), 353-366.
 */
void umat_plasticity_chaboche_CCP(const arma::vec &Etot, const arma::vec &DEtot, arma::vec &sigma, arma::mat &Lt, arma::mat &L, arma::vec &sigma_in, const arma::mat &DR, const int &nprops, const arma::vec &props, const int &nstatev, arma::vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, const int &solver_type, double &tnew_dt);
                                

/** @} */ // end of umat_mechanical group

} //namespace simcoon
