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

///@file random.hpp
///@brief random number generators
///@version 1.0

#pragma once

namespace simcoon{

//This function returns a random in number between 0 and a
int alea(const int &);

//This function returns a random in number between a and b
int aleab(const int &, const int &);

//This function returns a random double number between a and b
double alead(const double &, const double &);

} //namespace simcoon