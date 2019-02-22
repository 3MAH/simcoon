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

///@file kinematics.hpp
///@brief A set of function that allows various strain transformation (in Finite strains)
///@version 1.0

#pragma once
#include <armadillo>

namespace simcoon{

//This function returns the deviatoric part of m
arma::mat dev(const arma::mat &m);

//This function returns the spherical part of m
arma::mat sph(const arma::mat &m);

//This function returns F (in a vectorized), from E (Green-Lagrange strain) and R (Rotation matrix), according to a RU decomposition
arma::vec ER_to_F(const mat&E, const mat&R)

//This function computes the gradient of displacement (Lagrangian) from the deformation gradient tensor
arma::mat G_UdX(const arma::mat &);

//This function computes the gradient of displacement (Eulerian) from the deformation gradient tensor
arma::mat G_Udx(const arma::mat &);

//This function computes the Right Cauchy-Green C
arma::mat R_Cauchy_Green(const arma::mat &);

//This function computes the Left Cauchy-Green B
arma::mat L_Cauchy_Green(const arma::mat &);

//This funtion computes the RU decomposition of the deformation gradient F (order of arguments: R, U, F)
void RU_decomposition(arma::mat &, arma::mat &, const arma::mat &);

//This funtion computes the VR decomposition of the deformation gradient F (order of arguments: R, V, F)
void VR_decomposition(arma::mat &, arma::mat &, const arma::mat &);
    
//This function computes the common Right (or Left) Cauchy-Green invariants
arma::vec Inv_X(const arma::mat &);

//This function computes the Cauchy deformation tensor c
arma::mat Cauchy(const arma::mat &);
    
//This function computes the Green-Lagrange finite strain tensor E
arma::mat Green_Lagrange(const arma::mat &);

//This function computes the Euler-Almansi finite strain tensor A
arma::mat Euler_Almansi(const arma::mat &);
    
//This function computes the velocity difference (F,DF,DTime)
arma::mat finite_L(const arma::mat &, const arma::mat &, const double &);

//This function computes the spin tensor W (correspond to Jaumann rate) (F,DF,DTime)
arma::mat finite_W(const arma::mat &, const arma::mat &, const double &);
    
//This function computes the spin tensor Omega (corrspond to Green-Naghdi rate)
// Note : here R is the is the rigid body rotation in the RU or VR polar decomposition of the deformation gradient F (R,DR,DTime)
arma::mat finite_Omega(const arma::mat &, const arma::mat &, const double &);

//This function computes the deformation rate D (F,DF,DTime)
arma::mat finite_D(const arma::mat &, const arma::mat &, const double &);
    
//This function computes the increment of finite rotation (Omega0, Omega1, DTime)
arma::mat finite_DQ(const arma::mat &, const arma::mat &, const double &);
    
} //namespace simcoon
