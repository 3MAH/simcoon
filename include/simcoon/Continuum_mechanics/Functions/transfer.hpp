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

#pragma once
#include <armadillo>

namespace simcoon{

/**
 * @file transfer.hpp
 * @author Yves Chemisky
 * @brief A set of functions that transform stress and strain tensors from 3x3 mat into vector with Voigt notation.
 */

/** @addtogroup transfer
 *  @{
 */
    
/**
 * @brief Provides the matrix (3x3 arma::mat) version of a strain tensor initially written with Voigt notation (arma::vec size=6)
 *
 * @param v (arma::vec size=6) the strain tensor written in with Voigt notation
 * @return (3x3 arma::mat) the strain tensor written as a 3x3 matrix
 *
 * @details Example:
 * @code
        vec v = randu(6);
        mat m = v2t_strain(v);
 * @endcode
*/
arma::mat v2t_strain(const arma::vec &v);

/**
 * @brief Provides the Voigt notation (arma::vec size=6) version of a strain tensor initially written using a matrix format (3x3 arma::mat)
 *
 * @param m (3x3 arma::mat) the strain tensor written as a 3x3 matrix 
 * @return (arma::vec size=6) the strain tensor written in with Voigt notation
 *
 * @details Example:
 * @code
        mat m = randu(3,3);
        vec v = t2v_strain(m);
 * @endcode
*/
arma::vec t2v_strain (const arma::mat &m);

/**
 * @brief Provides the matrix (3x3 arma::mat) version of a stress tensor initially written with Voigt notation (arma::vec size=6)
 *
 * @param v (arma::vec size=6) the stress tensor written in with Voigt notation
 * @return (3x3 arma::mat) the stress tensor written as a 3x3 matrix
 *
 * @details Example:
 * @code
        vec v = randu(6);
        mat m = v2t_stress(v);
 * @endcode
*/
arma::mat v2t_stress(const arma::vec &v);

/**
 * @brief Provides the Voigt notation (arma::vec size=6) version of a stress tensor initially written using a matrix format (3x3 arma::mat)
 *
 * @param m (3x3 arma::mat) the stress tensor written as a 3x3 matrix 
 * @return (arma::vec size=6) the stress tensor written in with Voigt notation
 *
 * @details Example:
 * @code
        mat m = randu(3,3);
        vec v = t2v_stress(m);
 * @endcode
*/
arma::vec t2v_stress (const arma::mat &m);

/**
 * @brief Provides the Voigt notation (arma::vec size=6) version of a symmetric tensor initially written using a matrix format (3x3 arma::mat)
 *
 * Note that the 6 components are organized as the following (it is a column vector) : \f$ \mathbf{v} \equiv \left( m_{11},m_{22},m_{33},m_{12},m_{13},m_{23} right) \f$
 * 
 * @param m (3x3 arma::mat) the symmetric tensor written as a 3x3 matrix 
 * @return (arma::vec size=6) the symmetric tensor written in with Voigt notation
 *
 * @details Example:
 * @code
        mat m = randu(3,3);
        vec v = t2v_sym(m);
 * @endcode
*/
arma::vec t2v_sym (const arma::mat &m);

/**
 * @brief Provides the matrix (3x3 arma::mat) version of a symmetric tensor initially written with Voigt notation (arma::vec size=6)
 *
 * Note that the 6 components are organized as the following (it is a column vector) : \f$ \mathbf{v} \equiv \left( m_{11},m_{22},m_{33},m_{12},m_{13},m_{23} right) \f$
 * 
 * @param v (arma::vec size=6) the symmetric tensor written in with Voigt notation
 * @return (3x3 arma::mat) the symmetric tensor written as a 3x3 matrix
 *
 * @details Example:
 * @code
        vec v = randu(6);
        mat m = v2t_sym(v);
 * @endcode
*/
arma::mat v2t_sym (const arma::vec &v);

//This function transforms a vector (6 components 11,22,33,12,13,23) into a skew-symmetric 3x3 stress matrix

/**
 * @brief Provides the matrix (3x3 arma::mat) version of a skew-symmetric tensor initially written with Voigt notation (arma::vec size=6)
 *
 * Note that the 6 components are organized as the following (it is a column vector) : \f$ \mathbf{v} \equiv \left( m_{11},m_{22},m_{33},m_{12},m_{13},m_{23} right) \f$
 * 
 * so that the components \f$ m_{21} = -m_{12}, \quad m_{31} = -m_{13}, \quad m_{32} = -m_{23} \f$ so that:
 * 
 * \f[  m = \left( \begin{array}{ccc}
        v_1 & v_4 & v_5 \\
        -v_4 & v_2 & v_6 \\
        v_5 & -v_6 & v_3 \end{array} \right)
 * \f]
 *  
 * @param v (arma::vec size=6) the skew-symmetric tensor written in with Voigt notation
 * @return (3x3 arma::mat) the skew-symmetric tensor written as a 3x3 matrix
 *
 * @details Example:
 * @code
        vec v = randu(6);
        mat m = v2t_skewsym(v);
 * @endcode
*/
arma::mat v2t_skewsym (const arma::vec &v);

/**
 * @brief Provides the matrix (3x3 arma::mat) version of a tensor initially written with Voigt notation (arma::vec size=9)
 *
 * Note that the 9 components are organized as the following (it is a column vector) : \f$ \mathbf{v} \equiv \left( m_{11},m_{12},m_{13},m_{21},m_{22},m_{23},m_{31},m_{32},m_{33} right) \f$
 * So that this operation is the opposite of a flatten (.as_col() for armadillo matrix to column vectors)
 *  
 * @param v (arma::vec size=9) the flatten tensor written in with Voigt notation
 * @return (3x3 arma::mat) the tensor written as a 3x3 matrix
 *
 * @details Example:
 * @code
        vec v = randu(9);
        mat m = v2t(v);
 * @endcode
*/
arma::mat v2t (const arma::vec &v);

/** @} */ // end of transfer group

} //namespace simcoon
