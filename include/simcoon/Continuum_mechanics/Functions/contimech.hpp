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

///@file contimech.hpp
///@brief Functions that computes Mises stress/strains, directions, etc
///@version 1.0

#pragma once
#include <armadillo>

namespace simcoon{

//This function returns the deviatoric part of m
arma::mat dev(const arma::mat &);

//This function returns the spherical part of m
arma::mat sph(const arma::mat &);

//This function returns the trace of the tensor v
double tr(const arma::vec &);

//This function returns the deviatoric part of v
arma::vec dev(const arma::vec &);

//This function determines the Mises equivalent of a stress tensor, according to the Voigt convention for stress 
double Mises_stress(const arma::vec &);

//This function determines the strain flow (direction) from a stress tensor (Mises convention), according to the Voigt convention for strains
arma::vec eta_stress(const arma::vec &);
    
//This function determines the strain flow (direction) from a stress tensor, according to the Voigt convention for strains
arma::vec eta_norm_stress(const arma::vec &);

//This function determines the strain flow (direction) from a strain tensor, according to the Voigt convention for strains
arma::vec eta_norm_strain(const arma::vec &);
    
//This function determines the norm of a stress tensor
double norm_stress(const arma::vec &);

//This function determines the norm of a strain tensor
double norm_strain(const arma::vec &);
    
//This function determines the Mises equivalent of a strain tensor, according to the Voigt convention for strains 
double Mises_strain(const arma::vec &);

//This function determines the strain flow (direction) from a strain tensor, according to the Voigt convention for strains
arma::vec eta_strain(const arma::vec &);

//Returns the second invariant of the deviatoric part of a second order stress tensor written as a Voigt vector
double J2_stress(const arma::vec &);

//Returns the second invariant of the deviatoric part of a second order strain tensor written as a Voigt vector
double J2_strain(const arma::vec &);

//Returns the third invariant of the deviatoric part of a second order stress tensor written as a Voigt vector
double J3_stress(const arma::vec &);

//Returns the third invariant of the deviatoric part of a second order stress tensor written as a Voigt vector
double J3_strain(const arma::vec &);

//This function returns the value if it's positive, zero if it's negative (Macaulay brackets <>+)
double Macaulay_p(const double &);

//This function returns the value if it's negative, zero if it's positive (Macaulay brackets <>-)
double Macaulay_n(const double &);

//This function returns the value if it's negative, zero if it's positive (Macaulay brackets <>-)
double sign(const double &);
    
//Returns the normalized vector normal to an ellipsoid with semi-principal axes of length a1, a2, a3. The direction of the normalized vector is set by angles u
arma::vec normal_ellipsoid(const double &, const double &, const double &, const double &, const double &);

//Returns the curvature of an ellipsoid with semi-principal axes of length a1, a2, a3 at the angle u,v.
double curvature_ellipsoid(const double &, const double &, const double &, const double &, const double &);

//Returns the normal and tangent components of the stress vector in the normal direction n to an ellipsoid with axes a1, a2, a3. The direction of the normalized vector is set by angles u
arma::vec sigma_int(const arma::vec &, const double &, const double &, const double &, const double &, const double &);

///This computes the Hill interfacial operator according to a normal a (see papers of Siredey and Entemeyer phD dissertation)
arma::mat p_ikjl(const arma::vec &);

//This function returns the dyadic product of a symmetric tensor A in a 6*6 matrix format
arma::mat auto_sym_dyadic(const arma::mat &);

//This function returns the dyadic product of 2 symmetric tensors A and B, in a 6*6 matrix format
arma::mat sym_dyadic(const arma::mat &, const arma::mat &);

//This function returns the dyadic product of a symmetric tensors A and B, in a 6*6 matrix format
arma::mat auto_dyadic(const arma::mat &);

//This function returns the dyadic product of 2 symmetric tensors A and B, in a 6*6 matrix format
arma::mat dyadic(const arma::mat &, const arma::mat &);


} //namespace simcoon
