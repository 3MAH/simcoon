#pragma once
#include <boost/python.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/node.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/cubic_mesh.hpp>

namespace simpy{

boost::python::tuple test_mesh(const boost::python::numpy::ndarray &nodes_coords_py, const double & = 1.E-6, const double & = 1.E-4);

boost::python::list build_MPC_from_cubic_mesh(const boost::python::numpy::ndarray &, const boost::python::list &, const int & = 4, const double & = 1.0);

simcoon::Node build_node(const int &, const boost::python::numpy::ndarray &);

boost::python::numpy::ndarray Node_get_input_coords(simcoon::Node &);

void Node_set_input_coords(simcoon::Node &, const boost::python::numpy::ndarray &);
    
simcoon::cubic_mesh build_cubic_mesh(const std::string &, const boost::python::list &);

void get_domain(simcoon::cubic_mesh &);

void construct_lists(simcoon::cubic_mesh &);
    
boost::python::list read_nodes_file(const std::string &, const std::string &);

boost::python::list read_sections(const int &, const std::string &, const std::string &s);

//boost::python::list read_path(double &T, const boost::python::str &, const boost::python::str &);
    
} //namespace simpy
