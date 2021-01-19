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

///@file solver.hpp
///@brief To solver an homogeneous thermomechanical problem
///@version 1.0

#pragma once
#include <armadillo>
#include <string>

namespace simcoon{

//function that solves a
void solver(const std::string &, const arma::vec &, const unsigned int &, const double &, const double &, const double &, const int &, const int &, const double & = 0.5, const double & = 2., const int & = 10, const int & = 100, const int & = 1, const double & = 1.E-6, const double & = 10000., const std::string& = "data", const std::string& = "results", const std::string& = "path.txt", const std::string& = "result_job.txt");

} //namespace simcoon
