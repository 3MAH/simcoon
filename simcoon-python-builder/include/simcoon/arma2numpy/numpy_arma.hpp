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

#pragma once
#include <armadillo>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

namespace arma2numpy {

arma::vec array2vec(const boost::python::numpy::ndarray&);

arma::mat array2mat(const boost::python::numpy::ndarray&);

boost::python::numpy::ndarray vec2array(const arma::vec&);

boost::python::numpy::ndarray mat2array(const arma::mat&);

arma::Col<int> array2Col_int(const boost::python::numpy::ndarray&);

boost::python::numpy::ndarray Col_int2array(const arma::Col<int>&);

arma::Mat<int> array2Mat_int(const boost::python::numpy::ndarray&);

boost::python::numpy::ndarray Mat_int2array(const arma::Mat<int>&);
    
} //end of namespace arma2numpy