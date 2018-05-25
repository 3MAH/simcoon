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

///@file rotation.hpp
///@brief rotation of a Voigt tensor
///@version 1.0

#pragma once

#include <armadillo>

namespace simcoon{

//To rotate a vector, given a rotation matrix
arma::vec rotate_vec(const arma::vec &, const arma::mat &);
arma::vec rotate_vec(const arma::vec &, const double &, const int &);
//To rotate a matrix, given a rotation matrix
arma::mat rotate_mat(const arma::mat &, const arma::mat &);
arma::mat rotate_mat(const arma::mat &, const double &, const int &);

//To generate a 3x3 rotation matrix (with 1 angle)
arma::mat fillR(const double &, const int &, const bool & = true);
//To generate a 3x3 rotation matrix (with 3 angles and a convention, by default null string means the generic convention indicated in the parameters.hpp). By default active , alibi rotation is considered
arma::mat fillR(const double &, const double &, const double &, const bool & = true, const std::string & = "");
arma::mat fillR(const double &, const double &, const double phi, const int & , const int & , const int & , const bool &);
    
//To generate a 6x6 rotation matrix for stress tensors
arma::mat fillQS(const double &, const int &, const bool & = true);
arma::mat fillQS(const arma::mat &, const bool & = true);
    
//To generate a 6x6 rotation matrix for strain tensors
arma::mat fillQE(const double &, const int &, const bool & = true);
arma::mat fillQE(const arma::mat &, const bool & = true);

//To rotate a stiffness matrix (6,6)
arma::mat rotateL(const arma::mat &, const double &, const int &, const bool & = true);
arma::mat rotateL(const arma::mat &, const arma::mat &, const bool & = true);

//To rotate a compliance matrix (6,6)
arma::mat rotateM(const arma::mat &, const double &, const int &, const bool & = true);
arma::mat rotateM(const arma::mat &, const arma::mat &, const bool & = true);
    
//To rotate an interaction matrix A (6,6)
arma::mat rotateA(const arma::mat &, const double &, const int &, const bool & = true);
arma::mat rotateA(const arma::mat &, const arma::mat &, const bool & = true);

//To rotate a interaction matrix B (6,6)
arma::mat rotateB(const arma::mat &, const double &, const int &, const bool & = true);
arma::mat rotateB(const arma::mat &, const arma::mat &, const bool & = true);
    
//To rotate a stress vector (6)
arma::vec rotate_stress(const arma::vec &, const double &, const int &, const bool & = true);
arma::vec rotate_stress(const arma::vec &, const arma::mat &, const bool & = true);
    
//To rotate a strain vector (6)
arma::vec rotate_strain(const arma::vec &, const double &, const int &, const bool & = true);
arma::vec rotate_strain(const arma::vec &, const arma::mat &, const bool & = true);

//To rotate from local to global a strain tensor (6)
arma::mat rotate_l2g_strain(const arma::vec &, const double &, const double &, const double &);

//To rotate from global to local a strain tensor (6)
arma::mat rotate_g2l_strain(const arma::vec &, const double &, const double &, const double &);

//To rotate from local to global a stress tensor (6)
arma::mat rotate_l2g_stress(const arma::vec &, const double &, const double &, const double &);

//To rotate from global to local a stress tensor (6)
arma::mat rotate_g2l_stress(const arma::vec &, const double &, const double &, const double &);
    
//To rotate from local to global a stiffness matrix (6,6)
arma::mat rotate_l2g_L(const arma::mat &, const double &, const double &, const double &);

//To rotate from global to local a stiffness matrix (6,6)
arma::mat rotate_g2l_L(const arma::mat &, const double &, const double &, const double &);

//To rotate from local to global a strain localisation matrix (6,6)
arma::mat rotate_l2g_A(const arma::mat &, const double &, const double &, const double &);

//To rotate from global to local a strain localisation matrix (6,6)
arma::mat rotate_g2l_A(const arma::mat &, const double &, const double &, const double &);

//To rotate from local to global a stress localisation matrix (6,6)
arma::mat rotate_l2g_B(const arma::mat &, const double &, const double &, const double &);

//To rotate from global to local a stress localisation matrix (6,6)
arma::mat rotate_g2l_B(const arma::mat &, const double &, const double &, const double &);

//To rotate from local to global a compliance matrix (6,6)
arma::mat rotate_l2g_M(const arma::mat &, const double &, const double &, const double &);

//To rotate from global to local a compliance matrix (6,6)
arma::mat rotate_g2l_M(const arma::mat &, const double &, const double &, const double &);
    
    
} //namespace simcoon
