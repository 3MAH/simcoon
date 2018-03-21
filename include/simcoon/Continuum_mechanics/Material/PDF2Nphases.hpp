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
#include <simcoon/Continuum_mechanics/Material/PDF.hpp>

namespace simcoon{

//Fill the PDF from a vector of parameter, providing a file with the peak informations
arma::vec get_densities_PDF(const arma::vec &, const std::string &, const std::string &);
    
//Fill the parameters of the geom and material
void fill_parameters(const double &, phase_characteristics &, const PDF &);
    
//This function computes the PDF of the selected angle, according to different methods (Lorentzian, Pearson...)
phase_characteristics discretize_PDF(const phase_characteristics &, PDF &, const int &, const int &);

//Writes the Nphases.dat file for multiphase modeling, according to specific ODFs
//void PDF2Nphases(const arma::Col<int> &, const arma::Col<int> &, const arma::Col<int> &, const std::vector<std::string> &, const arma::mat &, const bool& = false, const double& = 0.);

} //namespace simcoon
