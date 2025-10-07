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

///@file Tarr_to_mat.hpp
///@brief Test for the Carma library
///@version 1.0

#pragma once
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <armadillo>
#include <carma>

namespace testcarma {

int test_arr_to_mat_double_copy(const py::array_t<double>& arr);
int test_arr_to_mat_1d(py::array_t<double>& arr, bool copy);
int test_arr_to_mat_long(py::array_t<int64_t>& arr, bool copy);
int test_arr_to_mat_double(py::array_t<double>& arr, bool copy);
int test_arr_to_col(py::array_t<double>& arr, bool copy);
int test_arr_to_row(py::array_t<double>& arr, bool copy);
int test_arr_to_cube(py::array_t<double>& arr, bool copy);
int test_to_arma_mat(py::array_t<double>& arr, bool copy);
int test_to_arma_col(py::array_t<double>& arr, bool copy);
int test_to_arma_row(py::array_t<double>& arr, bool copy);
int test_to_arma_cube(py::array_t<double>& arr, bool copy);

}  // namespace testcarma