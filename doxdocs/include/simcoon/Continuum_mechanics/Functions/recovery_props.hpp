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

///@file recovery_props.cpp
///@brief A set of function that allow to get the Elastic properties from stiffness/compliance tensors
///@version 1.0

#pragma once
#include <iostream>
#include <armadillo>

namespace simcoon{

//This function check the the properties of an unknown stiffness tensor
void check_symetries(const arma::mat &, std::string &, int &, arma::vec &, int &);

//This function recovers the properties of an isotropic stiffness tensor
arma::vec L_iso_props(const arma::mat &);

//This function recovers the properties of an isotropic compliance tensor
arma::vec M_iso_props(const arma::mat &);
    
//This function recovers the properties of a transversely isotropic stiffness tensor
arma::vec L_isotrans_props(const arma::mat &, const int &);

//This function recovers the properties of a transversely isotropic compliance tensor
arma::vec M_isotrans_props(const arma::mat &, const int &);

//This function recovers the properties of a cubic isotropic stiffness tensor
arma::vec L_cubic_props(const arma::mat &);

//This function recovers the properties of a cubic isotropic compliance tensor
arma::vec M_cubic_props(const arma::mat &);
    
//This function recovers the properties of an orthotropic compliance tensor
arma::vec L_ortho_props(const arma::mat &);
    
//This function recovers the properties of an orthotropic compliance tensor
arma::vec M_ortho_props(const arma::mat &);
    
//This function recovers the properties of a fully anisotropic compliance tensor
arma::vec M_aniso_props(const arma::mat &);

} //namespace simcoon
