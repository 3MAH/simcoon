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

namespace simcoon{

//This function computes the Rotation matrix R, the increment of the rotation, the rate of deformation D and the spin W depending on F_0 and F_1 (F at the beginning and end of an increment) using the Jauman corotational framework
void Jaumann(arma::mat &, arma::mat &, arma::mat &, arma::mat &, const double &, const arma::mat &, const arma::mat &);
    
//This function computes the Rotation matrix R, the increment of the rotation, the rate of deformation D and the spin W depending on F_0 and F_1 (F at the beginning and end of an increment) using the Green-Naghdi corotational framework
void Green_Naghdi(arma::mat &, arma::mat &, arma::mat &, arma::mat &, const double &, const arma::mat &, const arma::mat &);

//This function computes the tangent modulus that links the Lie derivative of the Kirchoff stress tau to the rate of deformation D, from the Saint-Venant Kirchoff elastic tensor (that links the Piola-Kirchoff II stress S to the Green-Lagrange stress E) and the transformation gradient F
arma::mat DSDE_2_Dtau_LieDD(const arma::mat &, const arma::mat &);

//This function computes the tangent modulus that links the Jaumann rate of the Kirchoff stress tau to the rate of deformation D, from the Saint-Venant Kirchoff elastic tensor (that links the Piola-Kirchoff II stress S to the Green-Lagrange stress E), the transformation gradient F and the Kirchoff stress tau
arma::mat DSDE_2_Dtau_JaumannDD(const arma::mat &, const arma::mat &, const arma::mat &);

//This function computes the tangent modulus that links the Jaumann rate of the Kirchoff stress tau to the rate of deformation D, from the tangent modulus that links the Jaumann rate of the Kirchoff stress tau to the rate of deformation D and the Kirchoff stress tau
arma::mat Dtau_LieDD_Dtau_JaumannDD(const arma::mat &, const arma::mat &);

//This function computes the tangent modulus that computes the Jaumann rate of the Kirchoff stress tau to the rate of deformation D, weighted by J (the determinant of F), from the tangent modulus that links the Lie derivative of the Kirchoff stress tau to the rate of deformation D
//Note : This is the tangent modulus utilized by Abaqus
arma::mat Dtau_LieDD_Lt_Aba(const arma::mat &Dtau_JaumannDD, const double &J);

//This function computes the tangent modulus that computes the Jaumann rate of the Kirchoff stress tau to the rate of deformation D, weighted by J (the determinant of F), from the tangent modulus that links the Jaumann rate of the Kirchoff stress tau to the rate of deformation D
//Note : This is the tangent modulus utilized by Abaqus
arma::mat Dtau_JaumannDD_Lt_Aba(const arma::mat &Dtau_JaumannDD, const double &J);
 
} //namespace simcoon
