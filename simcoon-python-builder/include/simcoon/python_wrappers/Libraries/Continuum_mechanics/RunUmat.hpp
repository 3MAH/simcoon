#pragma once
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

namespace bn = boost::python::numpy;
namespace bp = boost::python;

namespace simpy {

	bp::tuple get_Detot_fedoo(const bn::ndarray&, const bn::ndarray&, const double, const int);

	bp::tuple RunUmat_fedoo(const std::string&, const bn::ndarray&, const bn::ndarray&, const bn::ndarray&, const bn::ndarray&, bn::ndarray&, const double, const double, const double, const double, const int, const int, bn::ndarray&);

	bp::tuple Log_strain_fedoo(const bn::ndarray&);

}