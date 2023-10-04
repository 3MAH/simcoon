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
 along with arma2numpy.  If not, see <http://www.gnu.org/licenses/>.
 
 */

///@file Troundtrip.hpp
///@brief Test for the Carma library
///@version 1.0

#pragma once
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <armadillo>
#include <carma>

namespace py = pybind11;

namespace testcarma {

py::array_t<double> test_mat_roundtrip(py::array_t<double>& arr);
py::array_t<double> test_row_roundtrip(py::array_t<double>& arr);
py::array_t<double> test_col_roundtrip(py::array_t<double>& arr);
py::array_t<double> test_cube_roundtrip(py::array_t<double>& arr);
    
} //end of namespace testcarma
