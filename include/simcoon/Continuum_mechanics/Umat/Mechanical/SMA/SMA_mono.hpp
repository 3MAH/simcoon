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
 * @file sma_mono.hpp
 * @brief Micromechanical monocrystal model for shape memory alloys based on Patoor et al. (1996)
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
 * @brief Micromechanical monocrystal model for SMA based on Patoor et al. (1996)
 *
 * @details This function implements the micromechanical monocrystal model of Patoor et al. (1996)
 * for shape memory alloys. The model describes the martensitic phase transformation in a single
 * crystal by tracking N individual martensite variants, where N is typically 12 or 24 depending
 * on the crystallographic system.
 *
 * **Key Features:**
 * - Explicit tracking of N martensite variant volume fractions \f$ f_n \f$ (n = 1, ..., N)
 * - Crystallographic transformation strains for each variant
 * - Variant selection driven by resolved thermodynamic driving force
 * - Inter-variant interactions through hardening matrix
 * - Applicable to single crystal SMA behavior
 *
 * **Physical Background:**
 *
 * In single crystal SMAs, the austenite-to-martensite transformation occurs through the formation
 * of distinct crystallographic variants. For cubic-to-monoclinic transformations (e.g., NiTi),
 * there are typically 12 or 24 habit plane variants, each characterized by:
 * - A specific transformation strain tensor \f$ \boldsymbol{\varepsilon}^{tr}_n \f$
 * - A habit plane normal and transformation direction
 *
 * **Variant Volume Fractions:**
 *
 * The microstructural state is characterized by N variant volume fractions:
 * \f[
 * f_n \geq 0, \quad \sum_{n=1}^{N} f_n \leq 1, \quad f_A = 1 - \sum_{n=1}^{N} f_n
 * \f]
 * where:
 * - \f$ f_n \f$ is the volume fraction of martensite variant n
 * - \f$ f_A \f$ is the austenite volume fraction
 * - The total martensite fraction is \f$ \xi = \sum_{n=1}^{N} f_n \f$
 *
 * **Transformation Strain:**
 *
 * The macroscopic transformation strain is the volume-weighted sum over all variants:
 * \f[
 * \boldsymbol{\varepsilon}^{tr} = \sum_{n=1}^{N} f_n \boldsymbol{\varepsilon}^{tr}_n
 * \f]
 * where \f$ \boldsymbol{\varepsilon}^{tr}_n \f$ is the crystallographic transformation strain
 * tensor of variant n, computed from the lattice correspondence.
 *
 * **Thermodynamic Driving Force:**
 *
 * For each variant n, the driving force for transformation is:
 * \f[
 * F_n = \boldsymbol{\sigma} : \boldsymbol{\varepsilon}^{tr}_n - \Delta G^{chem}(T) - \sum_{m=1}^{N} H_{nm} f_m
 * \f]
 * where:
 * - \f$ \boldsymbol{\sigma} : \boldsymbol{\varepsilon}^{tr}_n \f$ is the mechanical driving force
 * - \f$ \Delta G^{chem}(T) \f$ is the chemical free energy difference (temperature-dependent)
 * - \f$ H_{nm} \f$ is the interaction matrix describing variant-variant hardening
 *
 * **Transformation Criteria:**
 *
 * **Forward Transformation (A → M_n):**
 * \f[
 * \Phi_n^{fwd} = F_n - F_c^{fwd} \leq 0, \quad \dot{f}_n \geq 0
 * \f]
 *
 * **Reverse Transformation (M_n → A):**
 * \f[
 * \Phi_n^{rev} = -F_n - F_c^{rev} \leq 0, \quad \dot{f}_n \leq 0
 * \f]
 *
 * where \f$ F_c^{fwd} \f$ and \f$ F_c^{rev} \f$ are critical driving forces for forward
 * and reverse transformations.
 *
 * **Interaction Matrix:**
 *
 * The hardening matrix \f$ H_{nm} \f$ captures:
 * - Self-hardening (\f$ H_{nn} \f$): resistance to growth of variant n
 * - Latent hardening (\f$ H_{nm}, n \neq m \f$): interaction between different variants
 *
 * **Variant Selection:**
 *
 * Under applied stress, variants with favorable orientation (high resolved stress on
 * transformation system) are preferentially activated. This leads to:
 * - Single variant formation under uniaxial loading along specific orientations
 * - Multi-variant microstructures under complex loading
 * - Texture-dependent macroscopic response
 *
 * **Number of Variants:**
 *
 * Common crystallographic systems:
 * - Cubic → Orthorhombic: N = 6 variants
 * - Cubic → Monoclinic (NiTi): N = 12 or 24 variants
 * - Cubic → Tetragonal: N = 3 variants
 *
 * **State Variables (statev):**
 *
 * Total state variables required: \f$ n_{statev} = 1 + N + 6 \f$ (for N variants)
 *
 * | Index | Symbol | Description | Units |
 * |-------|--------|-------------|-------|
 * | statev[0] | \f$ T_{init} \f$ | Initial/reference temperature | Temperature |
 * | statev[1] | \f$ f_1 \f$ | Volume fraction of martensite variant 1 | - |
 * | statev[2] | \f$ f_2 \f$ | Volume fraction of martensite variant 2 | - |
 * | ... | ... | ... | ... |
 * | statev[N] | \f$ f_N \f$ | Volume fraction of martensite variant N | - |
 * | statev[N+1] | \f$ \varepsilon^{tr}_{11} \f$ | Macroscopic transformation strain component 11 | Strain |
 * | statev[N+2] | \f$ \varepsilon^{tr}_{22} \f$ | Macroscopic transformation strain component 22 | Strain |
 * | statev[N+3] | \f$ \varepsilon^{tr}_{33} \f$ | Macroscopic transformation strain component 33 | Strain |
 * | statev[N+4] | \f$ \varepsilon^{tr}_{12} \f$ | Macroscopic transformation strain component 12 | Strain |
 * | statev[N+5] | \f$ \varepsilon^{tr}_{13} \f$ | Macroscopic transformation strain component 13 | Strain |
 * | statev[N+6] | \f$ \varepsilon^{tr}_{23} \f$ | Macroscopic transformation strain component 23 | Strain |
 *
 * The macroscopic transformation strain is computed as:
 * \f$ \boldsymbol{\varepsilon}^{tr} = \sum_{n=1}^{N} f_n \boldsymbol{\varepsilon}^{tr}_n \f$
 *
 * @param Etot Total strain tensor at beginning of increment (Voigt notation: 6×1 vector)
 * @param DEtot Strain increment tensor (Voigt notation: 6×1 vector)
 * @param sigma Cauchy stress tensor (Voigt notation: 6×1 vector) [output]
 * @param Lt Consistent tangent modulus (6×6 matrix) [output]
 * @param L Elastic stiffness tensor (6×6 matrix) [output]
 * @param DR Rotation increment matrix (3×3) for objective integration
 * @param nprops Number of material properties
 * @param props Material properties vector
 * @param nstatev Number of state variables (includes N variant volume fractions)
 * @param statev State variables vector containing variant fractions [input/output]
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
 * @note This model is specifically for **single crystal** SMA behavior
 * @note For polycrystalline SMAs, use this model within a homogenization scheme
 * @note The number of variants N depends on the crystallographic transformation system
 * @note Variant transformation strains must be provided based on crystallographic data
 *
 * **References:**
 * - Patoor, E., Eberhardt, A., & Berveiller, M. (1996). "Micromechanical modelling of
 *   superelasticity in shape memory alloys." *Journal de Physique IV*, 6(C1), 277-292.
 * - Patoor, E., Lagoudas, D. C., Entchev, P. B., Brinson, L. C., & Gao, X. (2006).
 *   "Shape memory alloys, Part I: General properties and modeling of single crystals."
 *   *Mechanics of Materials*, 38(5-6), 391-429.
 * - Gall, K., & Sehitoglu, H. (1999). "The role of texture in tension-compression
 *   asymmetry in polycrystalline NiTi." *International Journal of Plasticity*, 15(1), 69-92.
 */
void umat_sma_mono(const arma::vec &Etot, const arma::vec &DEtot, arma::vec &sigma, arma::mat &Lt, arma::mat &L, const arma::mat &DR, const int &nprops, const arma::vec &props, const int &nstatev, arma::vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, double &tnew_dt);

/** @} */ // end of umat_mechanical group
    
} //namespace simcoon
