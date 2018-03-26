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

///@file interpolate.hpp
///@brief Interpolation tools for geometry used in Finite Element Analysis
///@version 1.0

#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <assert.h>
#include <armadillo>
#include <simcoon/Continuum_mechanics/Unit_cell/node.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/equation.hpp>

namespace simcoon{
    
void find_neighbours_set(std::vector<std::vector<Node> > &, std::vector<arma::mat> &, const std::vector<Node> &, const std::vector<Node> &, const unsigned int &);
    
void set_weights(std::vector<arma::mat> &, const std::vector<arma::mat> &, const unsigned int &, const double &);
    
equation set_equation(const std::vector<Node> &neigh, const arma::mat &weight, const unsigned int &, const unsigned int &, const unsigned int &);

} //namespace simcoon