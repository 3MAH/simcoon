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
#include <CGAL/Simple_cartesian.h>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <simcoon/arma2numpy/numpy_arma.hpp>

namespace bp = boost::python;
namespace bn = boost::python::numpy;

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;

namespace cgal2numpy {

Point array2Point(bn::ndarray const &array) {
    
    assert(array.get_nd() == 1);
    Py_intptr_t const *shape = array.get_shape();
    assert(shape[0]==3);
    arma::vec coords = arma2numpy::array2vec(array);
    return Point(coords(0),coords(1),coords(2));
}

bn::ndarray Point2array(const Point &p) {
    
    //create a tuple with the size of v
    bp::tuple shape = bp::make_tuple(3);
    //as well as a type for C++ double
    bn::dtype dtype = bn::dtype::get_builtin<double>();
    //Construct an array with the above shape and type
    bn::ndarray pyp = bn::zeros(shape, dtype);
    
    pyp[0] = p.x();
    pyp[1] = p.y();
    pyp[2] = p.z();
    return pyp;
}


} //end of namespace cgal2numpy