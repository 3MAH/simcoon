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
#include <simcoon/FTensor.hpp>

namespace simcoon{

/**
* @file transfer.hpp
* @author Yves Chemisky
* @brief A set of functions that transforms stress and strain tensors from 3x3 mat into vector with Voigt notation 
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

//This function transforms an armadillo 3 colvec to a FTensor Tensor of the 1st rank

/**
 * @brief Provides the FTensor vector (Tensor1<double,3>) of a dimension 3 colvec (arma::vec size=3)
 *
 * @param v (arma::vec size=3) the column tensor 
 * @return (Tensor1<double,3>) the Ftensor of dimension 1, size=3
 *
 * @details Example:
 * @code
        vec v = randu(3);
        FTensor::Tensor1<double,3> F1 = vec_FTensor1(v);
 * @endcode
*/
FTensor::Tensor1<double,3> vec_FTensor1(const arma::vec &v);

/**
 * @brief Provides the FTensor matrix (Tensor2<double,3,3>) of a 3x3 armadillo matrix (arma::mat 3x3)
 *
 * @param m (arma::mat 3x3) the armadillo matrix 
 * @return (Tensor2<double,3,3>) the Ftensor of dimension 2, size=3x3
 *
 * @details Example:
 * @code
        mat m = randu(3,3);
        FTensor::Tensor2<double,3,3> FT2 = mat_FTensor2(m);
 * @endcode
*/
FTensor::Tensor2<double,3,3> mat_FTensor2(const arma::mat &m);

/**
 * @brief Provides the FTensor matrix (Tensor2<double,3,3>) of armadillo column vector that correspond to a strain vector in Voigt notation (arma::vec size=6)
 *
 * @param v (arma::vec size=6) the strain tensor written in with Voigt notation
 * @return (Tensor2<double,3,3>) the Ftensor of dimension 2, size=3x3
 *
 * @details Example:
 * @code
        mat m = randu(3,3);
        vec v = t2v_strain(m);
        FTensor::Tensor2<double,3,3> FT2_strain = v_FTensor2_strain(v);
 * @endcode
*/
FTensor::Tensor2<double,3,3> v_FTensor2_strain(const arma::vec &v);

/**
 * @brief Provides the FTensor matrix (Tensor2<double,3,3>) of armadillo column vector that correspond to a stress vector in Voigt notation (arma::vec size=6)
 *
 * @param v (arma::vec size=6) the strain tensor written in with Voigt notation
 * @return (Tensor2<double,3,3>) the Ftensor of dimension 2, size=3x3
 *
 * @details Example:
 * @code
        mat m = randu(3,3);
        vec v = t2v_strain(m);
        FTensor::Tensor2<double,3,3> FT2_strain = v_FTensor2_strain(v);
 * @endcode
*/
FTensor::Tensor2<double,3,3> v_FTensor2_stress(const arma::vec &v);

/**
 * @brief Provides an armadillo vector (arma::vec size=3) from a FTensor Tensor of the 1st rank (Tensor1<double,3>)
 *
 * @param FT1 (FTensor::Tensor1<double,3>) the FTensor of rank 1, and dimension 3
 * @return (arma::vec size=3) the corresponding armadillo vector
 *
 * @details Example:
 * @code
        FTensor::Tensor1<double,3> FT1 {0,1,2};
        arma::vec v = FTensor1_vec(FT1);
 * @endcode
*/
arma::vec FTensor1_vec(const FTensor::Tensor1<double,3> &FT1);

/**
 * @brief Provides an armadillo matrix (arma::mat 3x3) to a FTensor Tensor of the 2nd rank (Tensor2<double,3,3>)
 *
 * @param FT2 (Tensor2<double,3,3>) the Ftensor of dimension 2, size=3x3
 * @return (arma::mat 3x3) the corresponding armadillo matrix
 *
 * @details Example:
 * @code
        FTensor::Tensor2<double,3,3> FT2 {00, 01, 02
                                          10, 11, 12
                                          20, 21, 22};
        arma::mat m = FTensor2_mat(FT2);
 * @endcode
*/
arma::mat FTensor2_mat(const FTensor::Tensor2<double,3,3> &FT2);

/**
 * @brief Provides an armadillo vec in Voigt notation (arma::vec size=6) from a strain tensor expressed as a FTensor Tensor of the 2nd rank (Tensor2<double,3,3>)
 *
 * @param FT2 (Tensor2<double,3,3>) the Ftensor of dimension 2, size=3x3
 * @return (arma::vec size=6) the corresponding armadillo strain vector in Voigt notation
 *
 * @details Example:
 * @code
        FTensor::Tensor2<double,3,3> FT2 {00, 01, 02
                                          10, 11, 12
                                          20, 21, 22};
        arma::vec E = FTensor2_v_strain(FT2);
 * @endcode
*/
arma::vec FTensor2_v_strain(const FTensor::Tensor2<double,3,3> &FT2);
    
/**
 * @brief Provides an armadillo vec in Voigt notation (arma::vec size=6) from a stress tensor expressed as a FTensor Tensor of the 2nd rank (Tensor2<double,3,3>)
 *
 * @param FT2 (Tensor2<double,3,3>) the Ftensor of dimension 2, size=3x3
 * @return (arma::vec size=6) the corresponding armadillo strain vector in Voigt notation
 *
 * @details Example:
 * @code
        FTensor::Tensor2<double,3,3> FT2 {00, 01, 02
                                          10, 11, 12
                                          20, 21, 22};
        arma::vec sigma = FTensor2_v_stress(FT2);
 * @endcode
*/
arma::vec FTensor2_v_stress(const FTensor::Tensor2<double,3,3> &FT2);

/**
 * @brief Provides an armadillo (stiffness) mat in Voigt notation (arma::mat 6x6) from a  4th order tensor expressed as a FTensor (Tensor2<double,3,3,3,3>)
 *
 * @param FT4 (Tensor4<double,3,3,3,3>) the (stifness) tensor Ftensor of dimension 4, size=3x3x3x3
 * @return (arma::mat 6x6) the corresponding armadillo (stiffness) matrix in Voigt notation
 *
 * @details Example:
 * @code
       Tensor2<double,3,3> A_ = mat_FTensor2(A);
       Tensor4<double,3,3,3,3> C_;

       Index<'i', 3> i;
       Index<'j', 3> j;
       Index<'k', 3> k;
       Index<'l', 3> l;
              
       C_(i,j,k,l) = A_(i,j)*A_(k,l);
       return FTensor4_mat(C_);
 * @endcode
*/
arma::mat FTensor4_mat(const FTensor::Tensor4<double,3,3,3,3> &FT4);

//This function transforms a stiffness symmetric matrix 6x6 to a FTensor Tensor of the 4nd rank


/**
 * @brief Provides a 4th order (stiffness) tensor expressed as a FTensor (Tensor2<double,3,3,3,3>) from an armadillo (stiffness) mat in Voigt notation (arma::mat 6x6) 
 *
 * @param L (arma::mat 6x6) the corresponding armadillo (stiffness) matrix in Voigt notation
 * @return (Tensor4<double,3,3,3,3>) the (stifness) tensor Ftensor of dimension 4, size=3x3x3x3
 *
 * @details Example:
 * @code
	mat L = randu(6,6);
       Tensor4<double,3,3,3,3> C;	
	C = mat_FTensor4(L);
 * @endcode
*/
FTensor::Tensor4<double,3,3,3,3> mat_FTensor4(const arma::mat &L);
    
} //namespace simcoon
