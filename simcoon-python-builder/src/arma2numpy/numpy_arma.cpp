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

#include <string.h>
#include <assert.h>
#include <memory>
#include <armadillo>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <boost/python/extract.hpp>

namespace bp = boost::python;
namespace bn = boost::python::numpy;
using namespace arma;

namespace arma2numpy {

vec array2vec(bn::ndarray const &array) {
    
    assert(array.get_nd() == 1);
    Py_intptr_t const *shape = array.get_shape();
    int n_rows = shape[0];
    vec v(reinterpret_cast<double*>(array.get_data()), n_rows, false);
    return v;
}

bn::ndarray vec2array(const vec &v) {
    
    //create a tuple with the size of v
    bp::tuple shape = bp::make_tuple(v.n_elem);
    bp::tuple stride = bp::make_tuple(sizeof(double));
    bp::object own;
    //as well as a type for C++ double
    bn::dtype dtype = bn::dtype::get_builtin<double>();
    //Construct an array with the above shape and type
    bn::ndarray py_array = bn::from_data(v.memptr(),dtype,shape,stride,own);
    return py_array.copy();
}


mat array2mat(bn::ndarray const &array) {
    
    assert(array.get_nd() == 2);
    Py_intptr_t const *shape = array.get_shape();
    int n_rows = shape[0];
    int n_cols = shape[1];
    mat m(reinterpret_cast<double*>(array.get_data()), n_rows, n_cols, false);
    return m.t();
}

bn::ndarray mat2array(const mat &m) {
        
    //create a tuple with the size of m
    bp::tuple shape = bp::make_tuple(m.n_rows, m.n_cols);
    bp::tuple stride = bp::make_tuple(sizeof(double)*m.n_rows, sizeof(double));
    bp::object own;
    //as well as a type for C++ double
    bn::dtype dtype = bn::dtype::get_builtin<double>();
    mat m_t = m.t();
    bn::ndarray py_array = bn::from_data(m_t.memptr(),dtype,shape,stride,own);
    return py_array.copy();
}
    
Col<int> array2Col_int(bn::ndarray const &array) {
    
    assert(array.get_nd() == 1);
    Py_intptr_t const *shape = array.get_shape();
    int n_rows = shape[0];
    Col<int> v(n_rows);
    for (unsigned int i = 0; i < n_rows; ++i) {
        v(i) = atoi(bp::extract<char const *>(bp::str(array[i])));
    }
    return v;
}

bn::ndarray Col_int2array(const Col<int> &v) {
    
    //create a tuple with the size of v
    bp::tuple shape = bp::make_tuple(v.n_elem);
    bp::tuple stride = bp::make_tuple(sizeof(int));
    bp::object own;
    //as well as a type for C++ int
    bn::dtype dtype = bn::dtype::get_builtin<int>();
    //Construct an array with the above shape and type
    bn::ndarray py_array = bn::from_data(v.memptr(),dtype,shape,stride,own);
    return py_array.copy();
}


Mat<int> array2Mat_int(bn::ndarray const &array) {
    
    assert(array.get_nd() == 2);
    Py_intptr_t const *shape = array.get_shape();
    int n_rows = shape[0];
    int n_cols = shape[1];
    Mat<int> m(n_rows, n_cols);
    
    //    Py_intptr_t const * strides = array.get_strides();
    for (int i = 0; i < n_rows; ++i) {
        for (int j = 0; j < n_cols; ++j) {
            m(i,j) = atoi(bp::extract<char const *>(bp::str(array[i][j])));
        }
    }
    return m;
}

bn::ndarray Mat_int2array(const Mat<int> &m) {
    
    //create a tuple with the size of m
    bp::tuple shape = bp::make_tuple(m.n_rows, m.n_cols);
    bp::tuple stride = bp::make_tuple(sizeof(int)*m.n_cols, sizeof(int));
    bp::object own;
    //as well as a type for C++ double
    bn::dtype dtype = bn::dtype::get_builtin<int>();
    //Construct an array with the above shape and type
    Mat<int> m_t = m.t();
    bn::ndarray py_array = bn::from_data(m_t.memptr(),dtype,shape,stride,own);
    return py_array.copy();
}
    

} //end of namespace arma2numpy
