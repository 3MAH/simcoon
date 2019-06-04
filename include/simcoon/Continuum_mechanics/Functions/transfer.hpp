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

///@file objective_rate.hpp
///@brief A set of function that help to define different quantities, depending on a selected objective rate
///@version 1.0

#pragma once
#include <armadillo>
#include <FTensor.hpp>

namespace simcoon{
    
//This function transforms the strain Voigt vector into a 3x3 strain matrix
arma::mat v2t_strain(const arma::vec &v);

//This function transforms a 3x3 strain matrix into a strain Voigt vector
arma::vec t2v_strain (const arma::mat &);

//This function transforms the stress Voigt vector into a 3x3 stress matrix
arma::mat v2t_stress(const arma::vec &);

//This function transforms a 3x3 stress matrix into a stress Voigt vector
arma::vec t2v_stress (const arma::mat &);

//This function transforms a symmetric 3x3 stress matrix into vector (6 components 11,22,33,12,13,23)
arma::vec t2v_sym (const arma::mat &);

//This function transforms a vector (6 components 11,22,33,12,13,23) into a symmetric 3x3 stress matrix
arma::mat v2t_sym (const arma::vec &);

//This function transforms a vector (6 components 11,22,33,12,13,23) into a skew-symmetric 3x3 stress matrix
arma::mat v2t_skewsym (const arma::vec &);

//This function transforms a 9 component vector into a 3x3 matrix
arma::mat v2t (const arma::vec &);

//This function transforms an armadillo 3 colvec to a FTensor Tensor of the 1st rank
FTensor::Tensor1<double,3> vec_FTensor1(const arma::vec &);

//This function transforms an armadillo 3x3 matrix to a FTensor Tensor of the 2nd rank
FTensor::Tensor2<double,3,3> mat_FTensor2(const arma::mat &);

//This function transforms a strain Voigt vector to a FTensor Tensor of the 2nd rank
FTensor::Tensor2<double,3,3> v_FTensor2_strain(const arma::vec &);

//This function transforms a stress Voigt vector to a FTensor Tensor of the 2nd rank
FTensor::Tensor2<double,3,3> v_FTensor2_stress(const arma::vec &);

//This function transforms an armadillo 3x3 matrix to a FTensor Tensor of the 2nd rank
arma::vec FTensor1_vec(const FTensor::Tensor1<double,3> &);

//This function transforms an armadillo 3x3 matrix to a FTensor Tensor of the 2nd rank
arma::mat FTensor2_mat(const FTensor::Tensor2<double,3,3> &);

//This function transforms a FTensor Tensor of the 2nd rank to a strain Voigt vector
arma::vec FTensor2_v_strain(const FTensor::Tensor2<double,3,3> &);
    
//This function transforms a FTensor Tensor of the 2nd rank to a stress Voigt vector
arma::vec FTensor2_v_stress(const FTensor::Tensor2<double,3,3> &);

//This function transforms a FTensor Tensor of the 4nd rank to a stiffness symmetric matrix 6x6
arma::mat FTensor4_mat(const FTensor::Tensor4<double,3,3,3,3> &);

//This function transforms a stiffness symmetric matrix 6x6 to a FTensor Tensor of the 4nd rank
FTensor::Tensor4<double,3,3,3,3> mat_FTensor4(const arma::mat &);

arma::mat B_klmn(const arma::vec &b_i, const arma::vec &b_j);
    
} //namespace simcoon
