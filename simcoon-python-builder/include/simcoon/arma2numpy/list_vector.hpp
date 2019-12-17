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
#include <vector>
#include <boost/python.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/node.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/section_characteristics.hpp>
#include <simcoon/Simulation/Identification/constants.hpp>
#include <simcoon/Simulation/Identification/parameters.hpp>
#include <simcoon/Simulation/Solver/block.hpp>

namespace arma2numpy {

template<typename T> std::vector<T> py_list_to_std_vector(const boost::python::object &);
    
std::vector<double> py_list_to_std_vector_double(const boost::python::object &);
    
std::vector<int> py_list_to_std_vector_int(const boost::python::object &);
    
std::vector<std::string> py_list_to_std_vector_string(const boost::python::object &);
    
std::vector<simcoon::constants> py_list_to_std_vector_constants(const boost::python::object &);
    
std::vector<simcoon::parameters> py_list_to_std_vector_parameters(const boost::python::object &);

std::vector<simcoon::Node> py_list_to_std_vector_Node(const boost::python::object &);

std::vector<simcoon::block> py_list_to_std_vector_block(const boost::python::object &);

template <typename T> boost::python::object transfer_to_python(T* t);    
    
template <class T> boost::python::list std_vector_to_py_list(const std::vector<T> &);
    
template <class T> boost::python::list std_vector_to_py_list_class(const std::vector<T> &);
    
boost::python::list std_vector_to_py_list_double(const std::vector<double> &);
    
boost::python::list std_vector_to_py_list_int(const std::vector<int> &);
    
boost::python::list std_vector_to_py_list_string(const std::vector<std::string> &);
    
boost::python::list std_vector_to_py_list_constants(const std::vector<simcoon::constants> &);
    
boost::python::list std_vector_to_py_list_parameters(const std::vector<simcoon::parameters> &);

boost::python::list std_vector_to_py_list_Node(const std::vector<simcoon::Node> &);

boost::python::list std_vector_to_py_list_section_characteristics(const std::vector<simcoon::section_characteristics> &);

boost::python::list std_vector_to_py_list_block(const std::vector<simcoon::block> &);
    
} //end of namespace arma2numpy