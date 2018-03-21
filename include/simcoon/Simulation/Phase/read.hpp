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
///@brief To read from NphasesX.dat and NlayerX.dat
///@version 1.0

#pragma once
#include <armadillo>
#include <string>
#include <simcoon/Simulation/Phase/phase_characteristics.hpp>

namespace simcoon{
    
/// Function that reads the characteristics of a phase
void read_phase(phase_characteristics &, const std::string & = "data", const std::string & = "Nphases0.dat");

/// Function that reads the characteristics of a layer
void read_layer(phase_characteristics &, const std::string & = "data", const std::string & = "Nlayers0.dat");
    
/// Function that reads the characteristics of an ellipsoid
void read_ellipsoid(phase_characteristics &, const std::string & = "data", const std::string & = "Nellipsoids0.dat");

/// Function that reads the characteristics of a cylinder
void read_cylinder(phase_characteristics &, const std::string & = "data", const std::string & = "Ncylinders0.dat");

} //namespace simcoon