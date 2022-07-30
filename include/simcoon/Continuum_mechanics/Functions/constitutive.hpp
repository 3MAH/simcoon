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

///@file constitutive.hpp
///@brief Constitutive tensors in Voigt notation
///@version 1.0

#pragma once
#include <string.h>
#include <armadillo>

namespace simcoon{

/**
* @file constitutive.hpp
* @author Yves Chemisky 
* @section DESCRIPTION *
* The time class represents a moment of time. 
*/

///@brief Returns the fourth order identity tensor Ireal written in Voigt notation 
///@param None
///@return The following 6x6 mat (arma::mat) 
/**
\f[
    I_{real} = \left( \begin{array}{cccccc}
    1 & 0 & 0 & 0 & 0 & 0 \\
    0 & 1 & 0 & 0 & 0 & 0 \\
    0 & 0 & 1 & 0 & 0 & 0 \\
    0 & 0 & 0 & 0.5 & 0 & 0 \\
    0 & 0 & 0 & 0 & 0.5 & 0 \\
    0 & 0 & 0 & 0 & 0 & 0.5 \end{array} \right)
\f]
*/
arma::mat Ireal();

//Returns the volumic of the identity tensor Ireal written in Voigt notation
///@brief Returns the fourth order identity tensor written in Voigt notation
///@param None
///@return The following 6x6 mat (arma::mat) 
/**
\f[
    I_{vol} = \left( \begin{array}{cccccc}
    1/3 & 1/3 & 1/3 & 0 & 0 & 0 \\
    1/3 & 1/3 & 1/3 & 0 & 0 & 0 \\
    1/3 & 1/3 & 1/3 & 0 & 0 & 0 \\
    0 & 0 & 0 & 0 & 0 & 0 \\
    0 & 0 & 0 & 0 & 0 & 0 \\
    0 & 0 & 0 & 0 & 0 & 0 \end{array} \right)
\f]
*/
arma::mat Ivol();

///@brief Returns the deviatoric of the identity tensor Ireal written in Voigt notation
///@param None
///@return The following 6x6 mat (arma::mat) 
/**
\f[ 
    I_{dev} = I_{real} - I_{vol} = \left( \begin{array}{cccccc}
    2/3 & -1/3 & -1/3 & 0 & 0 & 0 \\
    -1/3 & 2/3 & -1/3 & 0 & 0 & 0 \\
    -1/3 & -1/3 & 2/3 & 0 & 0 & 0 \\
    0 & 0 & 0 & 0.5 & 0 & 0 \\
    0 & 0 & 0 & 0 & 0.5 & 0 \\
    0 & 0 & 0 & 0 & 0 & 0.5 \end{array} \right)
\f]
*/
arma::mat Idev();

///@brief Returns the fourth order identity tensor \f$ \widehat{I} \f$ written in Voigt notation
///@param None
///@return The following 6x6 mat (arma::mat) 
/**
\f[ 
        \widehat{I} = \left( \begin{array}{cccccc}
        1 & 0 & 0 & 0 & 0 & 0 \\
        0 & 1 & 0 & 0 & 0 & 0 \\
        0 & 0 & 1 & 0 & 0 & 0 \\
        0 & 0 & 0 & 2 & 0 & 0 \\
        0 & 0 & 0 & 0 & 2 & 0 \\
        0 & 0 & 0 & 0 & 0 & 2 \end{array} \right)
\f]
*/
arma::mat Ireal2();

///@brief Returns the deviatoric part of the identity tensor, in the form of \f$ \widehat{I} \f$ considering the Voigt notation
///@param None
///@return The following 6x6 mat (arma::mat) 
/**
\f[ 
    I_{dev2} = \left( \begin{array}{cccccc}
    2/3 & -1/3 & -1/3 & 0 & 0 & 0 \\
    -1/3 & 2/3 & -1/3 & 0 & 0 & 0 \\
    -1/3 & -1/3 & 2/3 & 0 & 0 & 0 \\
    0 & 0 & 0 & 2 & 0 & 0 \\
    0 & 0 & 0 & 0 & 2 & 0 \\
    0 & 0 & 0 & 0 & 0 & 2 \end{array} \right)
\f]
*/
arma::mat Idev2();

///@brief Returns the expansion vector
///@param None
///@return The following 6 vec (arma::vec) 
/**
\f[ 
    I_{th} = \left( \begin{array}{c}
    1 \\
    1 \\
    1 \\
    0 \\
    0 \\
    0 \end{array} \right)
\f]
*/
arma::vec Ith();

///@brief Returns the stress to strain operator
///@param None
///@return The following 6 vec (arma::vec) 
/**
\f[ 
    I_{r2} = \left( \begin{array}{ccc}
    1 \\
    1 \\
    1 \\
    2 \\
    2 \\
    2 \end{array} \right)
\f]
*/
arma::vec Ir2();

///@brief Returns the stress 2 strain operator
///@param None
///@return The following 6 vec (arma::vec) 
/**
\f[ 
    I_{r05} = \left( \begin{array}{ccc}
    1 \\
    1 \\
    1 \\
    0.5 \\
    0.5 \\
    0.5 \end{array} \right)
\f]
*/
arma::vec Ir05();

//Provides the elastic stiffness tensor for an isotropic material.
//The two first arguments are a couple of Lamé coefficients. The third argument specify which couple has been provided and the order of coefficients.
//Exhaustive list of possible third argument :
// ‘Enu’,’nuE,’Kmu’,’muK’, ‘KG’, ‘GK’, ‘lambdamu’, ‘mulambda’, ‘lambdaG’, ‘Glambda’.

///@brief Provides the elastic stiffness tensor for an isotropic material.
/*
* The two first arguments are a couple of Lamé coefficients. The third argument specify which couple has been provided and the order of coefficients.
* Exhaustive list of possible third argument :
* ‘Enu’,’nuE,’Kmu’,’muK’, ‘KG’, ‘GK’, ‘lambdamu’, ‘mulambda’, ‘lambdaG’, ‘Glambda’.
*/
///@param C1, C2, conv
///@return The 6x6 stiffness matrix considering a Voigt notation (arma::mat) 

arma::mat L_iso(const double &, const double &, const std::string& = "Enu");

//Provides the elastic compliance tensor for an isotropic material.
//The two first arguments are a couple of Lamé coefficients. The third argument specify which couple has been provided and the order of coefficients.
//Exhaustive list of possible third argument :
//‘Enu’,’nuE,’Kmu’,’muK’, ‘KG’, ‘GK’, ‘lambdamu’, ‘mulambda’, ‘lambdaG’, ‘Glambda’.
arma::mat M_iso(const double &, const double &, const std::string& = "Enu");

//Returns the elastic stiffness tensor for a cubic material.
//Arguments are the stiffness coefficients C11, C12 and C44 or E, nu and G.
//‘EnuG’,’Cii’.
arma::mat L_cubic(const double &, const double &, const double &, const std::string& = "EnuG");

//Returns the elastic compliance tensor for an isotropic material.
//Arguments are the stiffness coefficients C11, C12 and C44, or E, nu and G.
arma::mat M_cubic(const double &, const double &, const double &, const std::string& = "EnuG");

//Returns the elastic stiffness tensor for an orthotropic material.
//Arguments are the stiffness coefficients Cii or E and nu's
arma::mat L_ortho(const double &, const double &, const double &, const double &, const double &, const double &, const double &, const double &, const double &, const std::string& = "EnuG");

//Returns the elastic compliance tensor for an orthotropic material.
//Arguments are the stiffness coefficients Cii or E and nu's
arma::mat M_ortho(const double &, const double &, const double &, const double &, const double &, const double &, const double &, const double &, const double &, const std::string& = "EnuG");

//Returns the elastic stiffness tensor for an isotropic transverse material.
//Arguments are longitudinal Young modulus EL, transverse young modulus, Poisson’s ratio for loading along the longitudinal axis nuTL, Poisson’s ratio for loading along the transverse axis nuTT, shear modulus GLT and the axis of symmetry.
arma::mat L_isotrans(const double &, const double &, const double &, const double &, const double &, const int &);

//Returns the elastic compliance tensor for an isotropic transverse material.
//Arguments are longitudinal Young modulus EL, transverse young modulus, Poisson’s ratio for loading along the longitudinal axis nuTL, Poisson’s ratio for loading along the transverse axis nuTT, shear modulus GLT and the axis of symmetry.
arma::mat M_isotrans(const double &, const double &, const double &, const double &, const double &, const int &);

//Provides the viscous tensor H an isotropic material.
//The two first arguments are a couple of viscous coefficients (the first is bulk, the second is shear).
arma::mat H_iso(const double &, const double &);

//Update the elastic prediction, providing the stiffness tensor and the trial elastic strain
arma::vec el_pred(const arma::vec &, const arma::mat &, const arma::vec &, const int & = 3);
    
//Return the elastic prediction stress, providing the stiffness tensor and the trial elastic strain
arma::vec el_pred(const arma::mat &, const arma::vec &, const int &ndi = 3);

//Return the isotropized tangent modulus from the spectral decomposition of Bornert.etal (2001)
arma::mat Isotropize(const arma::mat &);

} //namespace simcoon
