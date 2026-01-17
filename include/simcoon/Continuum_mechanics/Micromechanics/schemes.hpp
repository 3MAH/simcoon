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

///@file schemes.hpp
///@brief micromechanical schemes for non-linear N-phases heterogeneous materials:
///@brief Mori-Tanaka scheme 
///@version 1.0

#pragma once

#include <armadillo>
#include <simcoon/Simulation//Phase/phase_characteristics.hpp>

namespace simcoon{

/**
 * @file schemes.hpp
 * @brief Micromechanical homogenization schemes.
 */

/** @addtogroup micromechanics
 *  @{
 */


// The Mori-Tanaka function works like a UMAT with the following material properties
///@brief props[0] : Number of phases
///@brief props[1] : Number of the file NPhase[i].dat utilized
///@brief props[2] : Number of integration points in the 1 direction
///@brief props[3] : Number of integration points in the 2 direction

/**
 * @brief Compute tangent stiffness using the homogeneous (Voigt) scheme.
 * @param rve Phase characteristics of the RVE
 * 
 * @details This scheme assumes uniform strain in all phases (Voigt hypothesis):
 * \f[
 *     \mathbf{L}^{eff} = \sum_i f_i \mathbf{L}_i
 * \f]
 * where \f$ f_i \f$ is the volume fraction and \f$ \mathbf{L}_i \f$ the stiffness of phase \f$ i \f$.
 */
void Lt_Homogeneous_E(phase_characteristics &);

/**
 * @brief Compute strain increment using the homogeneous (Voigt) scheme.
 * @param rve Phase characteristics of the RVE
 */
void DE_Homogeneous_E(phase_characteristics &);

/**
 * @brief Compute tangent stiffness using the Mori-Tanaka scheme.
 * @param rve Phase characteristics of the RVE
 * @param n_matrix Index of the matrix phase
 * 
 * @details The Mori-Tanaka scheme accounts for inclusion interactions through the
 * matrix strain field. The effective stiffness is:
 * \f[
 *     \mathbf{L}^{MT} = \mathbf{L}_0 + \sum_i f_i (\mathbf{L}_i - \mathbf{L}_0) : \mathbf{A}_i^{MT}
 * \f]
 */
void Lt_Mori_Tanaka(phase_characteristics &, const int &);

/**
 * @brief Compute strain increment using the Mori-Tanaka scheme.
 * @param rve Phase characteristics of the RVE
 * @param n_matrix Index of the matrix phase
 */
void DE_Mori_Tanaka(phase_characteristics &, const int &);
    
/**
 * @brief Compute tangent stiffness using the isotropized Mori-Tanaka scheme.
 * @param rve Phase characteristics of the RVE
 * @param n_matrix Index of the matrix phase
 */
void Lt_Mori_Tanaka_iso(phase_characteristics &, const int &);

/**
 * @brief Compute strain increment using the isotropized Mori-Tanaka scheme.
 * @param rve Phase characteristics of the RVE
 * @param n_matrix Index of the matrix phase
 */
void DE_Mori_Tanaka_iso(phase_characteristics &, const int &);
    
/**
 * @brief Compute tangent stiffness using the Self-Consistent scheme.
 * @param rve Phase characteristics of the RVE
 * @param n_matrix Index of the matrix phase
 * @param start Flag indicating if this is the first iteration
 * @param max_iter Maximum number of iterations (default: 1)
 * 
 * @details The Self-Consistent scheme determines the effective medium iteratively:
 * \f[
 *     \mathbf{L}^{SC} = \sum_i f_i \mathbf{L}_i : \mathbf{A}_i^{SC}
 * \f]
 * where the localization tensors depend on the (unknown) effective medium.
 */
void Lt_Self_Consistent(phase_characteristics &, const int &, const bool &, const int & = 1);

/**
 * @brief Compute strain increment using the Self-Consistent scheme.
 * @param rve Phase characteristics of the RVE
 * @param n_matrix Index of the matrix phase
 * @param start Flag indicating if this is the first iteration
 * @param max_iter Maximum number of iterations (default: 1)
 */
void DE_Self_Consistent(phase_characteristics &, const int &, const bool &, const int & = 1);
    
/**
 * @brief Compute tangent stiffness using the Periodic Layer scheme.
 * @param rve Phase characteristics of the RVE
 * 
 * @details This scheme is used for periodic laminate composites where layers
 * are stacked in the direction 1.
 */
void Lt_Periodic_Layer(phase_characteristics &);

/**
 * @brief Compute strain increment using the Periodic Layer scheme.
 * @param rve Phase characteristics of the RVE
 * @param n_matrix Index of the matrix phase
 */
void dE_Periodic_Layer(phase_characteristics &, const int &);
    

/** @} */ // end of micromechanics group

} //namespace simcoon
