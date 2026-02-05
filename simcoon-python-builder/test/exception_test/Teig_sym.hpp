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

///@file Teig_sym.hpp
///@brief Test for the exceptions raised with simcoon
///@version 1.0

#pragma once
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <armadillo>
#include <carma>

namespace testexception {

pybind11::array_t<double> test_eig_sym_val_affect(const py::array_t<double> &matrix_fail_eigen, const bool &copy=true);

pybind11::array_t<double> test_eig_sym_val_modify(const py::array_t<double> &matrix_fail_eigen, const bool &copy=true);
    
pybind11::tuple test_eig_sym_val_vec(const py::array_t<double> &matrix_fail_eigen, const bool &copy=true);

//py::tuple test_eig_sym(py::array_t<double>& arr, bool copy);

}  // namespace testexception