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

///@file identification.hpp
///@brief Function that run solver identification algorithms based on Smart+ solver
///@version 1.0

#pragma once

#include <iostream>
#include <armadillo>
#include "parameters.hpp"
#include "opti_data.hpp"
#include "generation.hpp"

namespace simcoon{
    
void run_identification(const std::string &, const int &, const int &, const int &, const int &, const int &, int &, int &, const int &, const int &, const int & = 6, const double & = 1.E-12, const std::string & = "data/", const std::string & = "keys/", const std::string & = "results/", const std::string & = "material.dat", const std::string & = "id_params.txt", const std::string & = "simul.txt", const double & = 5, const double & = 0.01, const double & = 0.001, const double & = 10, const double & = 0.01);

} //namespace simcoon
