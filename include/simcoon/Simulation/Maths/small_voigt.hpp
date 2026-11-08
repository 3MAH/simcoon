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
///@file small_voigt.hpp
///@brief Voigt indexing of a symmetric 3x3 and the stack conversions between the 6 engineering-Voigt values and the 3x3
///@version 1.0

#pragma once
#include <armadillo>

namespace simcoon{

/** @addtogroup maths
 *  @{
 */

/// Voigt index of the pair (i,j): 11, 22, 33, 12, 13, 23 -> 0..5 (symmetric in i, j)
constexpr int voigt_index[3][3] = {{0, 3, 4}, {3, 1, 5}, {4, 5, 2}};
/// First index i of the Voigt entry I
constexpr int voigt_i[6] = {0, 1, 2, 0, 0, 1};
/// Second index j of the Voigt entry I
constexpr int voigt_j[6] = {0, 1, 2, 1, 2, 2};

/**
 * @brief Symmetric 3x3 from 6 engineering-Voigt values (11, 22, 33, 12, 13, 23).
 *
 * @p shear scales the three shear entries on the way in: 1 for a stress (the Voigt entries are
 * the tensor shears), 0.5 for a strain (the Voigt entries are \f$ \gamma_{ij} = 2\varepsilon_{ij} \f$).
 * One rounded product per shear entry, as tensor2::from_voigt.
 */
inline void mat_from_voigt(arma::mat::fixed<3,3> &m, const double *v, const double shear) {
    m(0,0) = v[0]; m(1,1) = v[1]; m(2,2) = v[2];
    m(0,1) = shear * v[3]; m(1,0) = shear * v[3];
    m(0,2) = shear * v[4]; m(2,0) = shear * v[4];
    m(1,2) = shear * v[5]; m(2,1) = shear * v[5];
}

/// Same, returned by value.
inline arma::mat::fixed<3,3> mat_from_voigt(const double *v, const double shear) {
    arma::mat::fixed<3,3> m;
    mat_from_voigt(m, v, shear);
    return m;
}

/**
 * @brief 6 engineering-Voigt values of a 3x3, symmetrised: \f$ v_I = w\,(m_{ij} + m_{ji}) \f$ for the shears.
 *
 * @p shear_sum is 0.5 for a stress (the average of the two off-diagonal entries) and 1 for a
 * strain (their sum, \f$ \gamma_{ij} \f$). One rounded product per shear entry, as tensor2::voigt.
 */
inline void voigt_from_mat(double *v, const arma::mat::fixed<3,3> &m, const double shear_sum) {
    v[0] = m(0,0); v[1] = m(1,1); v[2] = m(2,2);
    v[3] = shear_sum * (m(0,1) + m(1,0));
    v[4] = shear_sum * (m(0,2) + m(2,0));
    v[5] = shear_sum * (m(1,2) + m(2,1));
}

/** @}*/

} //namespace simcoon
