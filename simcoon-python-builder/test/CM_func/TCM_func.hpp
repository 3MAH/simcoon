/* This file is part of arma2numpy.
 
 arma2numpy is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 arma2numpy is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with arma2numpy.  If not, see <http://www.gnu.org/licenses/>.
 
 */

///@file Tarma2numpy.hpp
///@brief Test for the arma2numpy library
///@version 1.0

#pragma once
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

namespace arma2numpy {

boost::python::numpy::ndarray L_iso(const double &, const double &, const std::string &);
    
} //end of namespace arma2numpy
