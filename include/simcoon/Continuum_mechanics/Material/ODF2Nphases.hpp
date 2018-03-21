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

///@file ODF2Nphases.hpp
///@brief ODF2Nphases discretization of ODFs
///@version 1.0

#pragma once

#include <iostream>
#include <string.h>
#include <armadillo>
#include <simcoon/Simulation/Phase/phase_characteristics.hpp>
#include <simcoon/Continuum_mechanics/Material/ODF.hpp>

namespace simcoon{

//Fill the ODF from a vector of angles, providing a file with the peak informations
arma::vec get_densities_ODF(const arma::vec &, const std::string &, const std::string &, const bool &);
    
//Fill the angles of the geom and material (if indicated 1 in angles_mat)
void fill_angles(const double &, phase_characteristics &, const ODF &, const int & = 1);
    
//This function computes the ODF of the selected angle, according to different methods (Lorentzian, Pearson...)
phase_characteristics discretize_ODF(const phase_characteristics &, ODF &, const int &, const int &, const int & = 1);

//Writes the Nphases.dat file for multiphase modeling, according to specific ODFs
//void ODF2Nphases(const arma::Col<int> &, const arma::Col<int> &, const arma::Col<int> &, const std::vector<std::string> &, const arma::mat &, const bool& = false, const double& = 0.);

} //namespace simcoon
