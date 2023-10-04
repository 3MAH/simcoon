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

///@file Troundtrip.cpp
///@brief Test for the carma library
///@version 1.0

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <armadillo>
#include <carma>

namespace py = pybind11;

namespace testcarma {

py::array_t<double> test_mat_roundtrip(py::array_t<double>& arr) {
    // call function to be tested
    arma::Mat<double> M = carma::arr_to_mat<double>(std::move(arr));
    return carma::mat_to_arr(M);
}

py::array_t<double> test_row_roundtrip(py::array_t<double>& arr) {
    // call function to be tested
    arma::Row<double> M = carma::arr_to_row<double>(std::move(arr));
    return carma::row_to_arr(M);
}

py::array_t<double> test_col_roundtrip(py::array_t<double>& arr) {
    // call function to be tested
    arma::Col<double> M = carma::arr_to_col<double>(std::move(arr));
    return carma::col_to_arr(M);
}

py::array_t<double> test_cube_roundtrip(py::array_t<double>& arr) {
    // call function to be tested
    arma::Cube<double> M = carma::arr_to_cube<double>(std::move(arr));
    return carma::cube_to_arr(M);
}    
    
} //end of namespace arma2numpy

