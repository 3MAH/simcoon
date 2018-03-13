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

///@file Tarma2numpy.cpp
///@brief Test for the arma2numpy library
///@version 1.0

#include <armadillo>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <simcoon/arma2numpy/numpy_arma.hpp>
#include <simcoon/arma2numpy/list_vector.hpp>
#include <simcoon/Simulation/Identification/constants.hpp>

namespace bp = boost::python;
namespace bn = boost::python::numpy;
using namespace std;
using namespace arma;

namespace arma2numpy {

bn::ndarray test_vec_int(bn::ndarray const &y) {
    Col<int> v = array2Col_int(y);
    return Col_int2array(v);
}

bn::ndarray test_mat_int(bn::ndarray const &y) {
    Mat<int> m = array2Mat_int(y);
    return Mat_int2array(m);
}
    
bn::ndarray test_vec_double(bn::ndarray const &y) {
    vec v = array2vec(y);
    return vec2array(v);
}

bn::ndarray test_mat_double(bn::ndarray const &y) {
    mat m = array2mat(y);
    return mat2array(m);
}
    
bp::list test_vector_list_double(bp::list const &y) {
 std::vector<double> v = py_list_to_std_vector_double(y);
 return std_vector_to_py_list_double(v);
}
    
bp::list test_vector_list_int(bp::list const &y) {
    std::vector<int> v = py_list_to_std_vector_int(y);
    return std_vector_to_py_list_int(v);
}

bp::list test_vector_list_string(bp::list const &y) {
    std::vector<std::string> v = py_list_to_std_vector_string(y);
    return std_vector_to_py_list_string(v);
}

bp::list test_vector_list_constants(bp::list const &y) {
    std::vector<simcoon::constants> v = py_list_to_std_vector_constants(y);
    return std_vector_to_py_list_constants(v);
}
    
bp::list test_vector_list_parameters(bp::list const &y) {
    std::vector<simcoon::constants> v = py_list_to_std_vector_constants(y);
    return std_vector_to_py_list_constants(v);
}
    
    
} //end of namespace arma2numpy

