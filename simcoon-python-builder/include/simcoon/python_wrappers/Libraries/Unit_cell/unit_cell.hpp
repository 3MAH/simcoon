#pragma once
#include <boost/python.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/node.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/cubic_mesh.hpp>

namespace simpy{
    
simcoon::Node build_node(const int &, const boost::python::numpy::ndarray &);

boost::python::numpy::ndarray Node_get_input_coords(simcoon::Node &);

void Node_set_input_coords(simcoon::Node &, const boost::python::numpy::ndarray &);
    
simcoon::cubic_mesh build_cubic_mesh(const boost::python::str &, const boost::python::list &);

void get_domain(simcoon::cubic_mesh &);

void construct_lists(simcoon::cubic_mesh &);
    
boost::python::list read_mesh(const boost::python::str &, const boost::python::str &);

boost::python::list read_sections(const int &, const boost::python::str &, const boost::python::str &);

//boost::python::list read_path(double &T, const boost::python::str &, const boost::python::str &);
    
} //namespace simpy