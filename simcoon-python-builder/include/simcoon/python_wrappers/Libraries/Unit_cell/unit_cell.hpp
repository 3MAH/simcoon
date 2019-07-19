#pragma once
#include <boost/python.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/node.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/cubic_mesh.hpp>

namespace simpy{
    
simcoon::Node build_node(const int &, const boost::python::numpy::ndarray &);

boost::python::numpy::ndarray readwrite_coords(simcoon::Node &);
    
void get_domain(simcoon::cubic_mesh);

void construct_lists(simcoon::cubic_mesh &);
    
//simcoon::cubic_mesh build_cubic_mesh(const boost::python::list &, const boost::python::list &);
    
} //namespace simpy