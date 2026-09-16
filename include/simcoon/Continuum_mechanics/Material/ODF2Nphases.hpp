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

///@file ODF2Nphases.hpp
///@brief ODF2Nphases discretization of ODFs
///@version 1.0

#pragma once

#include <iostream>
#include <string.h>
#include <armadillo>
#include <simcoon/Continuum_mechanics/Material/ODF.hpp>

namespace simcoon{

/**
 * @file ODF2Nphases.hpp
 * @brief ODF to N-phases discretization.
 */

/** @addtogroup material
 *  @{
 */


/**
 * @brief Densities of an ODF, summed over its peaks, at the angles x.
 * @param x the angles, in [0, pi] (radians) or [0, 180] (degrees, radian = false)
 * @param peaks the peaks; with radian = false their mean, s_dev and width are degrees
 * @param radian whether x and the peak angles are radians
 */
arma::vec get_densities_ODF(const arma::vec &x, const std::vector<peak> &peaks, const bool &radian);
    
// The discretisation of a phase into phases oriented along the ODF is done in Python
// (simcoon.solver.micromechanics.discretize_odf) from these densities.
    

/** @} */ // end of material group

} //namespace simcoon
