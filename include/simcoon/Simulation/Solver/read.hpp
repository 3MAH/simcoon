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

///@file read.hpp
///@brief To read from material.dat and path.dat
///@version 1.0

#pragma once
#include <armadillo>
#include <string>
#include "block.hpp"
#include "output.hpp"

namespace simcoon{

arma::Col<int> subdiag2vec();

/// Function that fills the matrix Tdsde for mix strain/stress conditions
void Lt_2_K(const arma::mat &, arma::mat &, const arma::Col<int> &, const double &);

/// Function that fills the matrix Tdsde for mix strain/stress conditions
void Lth_2_K(const arma::mat &, arma::mat &, arma::mat &, arma::mat &, arma::mat &, const arma::Col<int> &, const int &, const double &);

/// Function that reads the material properties
void solver_essentials(int &, int &, const std::string & = "data", const std::string & = "solver_essentials.inp");

/// Function that reads the material properties
void solver_control(double &, double &, int &, int &, int &, double &, double &, const std::string & = "data", const std::string & = "solver_control.inp");
    
/// Function that reads the material properties
void read_matprops(std::string &, unsigned int &, arma::vec &, unsigned int &, double &, double &, double &, const std::string & = "data", const std::string & = "material.dat");
    
/// Function that reads the output parameters
void read_output(solver_output &, const int &, const int &, const std::string & = "data", const std::string & = "output.dat");

/// Function that checks the coherency between the path and the step increments provided
void check_path_output(const std::vector<block> &, const solver_output &);
    
/// Function that reads the loading path
void read_path(std::vector<block> &, double &, const std::string & = "data", const std::string & = "path.txt");

} //namespace simcoon
