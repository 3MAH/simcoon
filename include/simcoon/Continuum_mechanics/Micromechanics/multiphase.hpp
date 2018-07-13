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

///@file multiphase.hpp
///@brief User subroutine for non-linear N-phases heterogeneous materials using the method
///@version 1.0

#pragma once

#include <armadillo>
#include <simcoon/Simulation/Phase/phase_characteristics.hpp>

namespace simcoon{

// The multiphase function works with the following material properties
///@brief props[0] : Number of phases
///@brief props[1] : Number of the file NPhase[i].dat utilized
///@brief props[2] : Number of integration points in the 1 direction
///@brief props[3] : Number of integration points in the 2 direction

void umat_multi(phase_characteristics &, const arma::mat &, const double &,const double &, const int &, const int &, const bool &, const unsigned int &, double &, const int &);

} //namespace simcoon
