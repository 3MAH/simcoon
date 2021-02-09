#pragma once
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

namespace bn = boost::python::numpy;
namespace bp = boost::python;

namespace simpy {

	bp::tuple RunUmat_fedoo(const std::string&, bn::ndarray&, const bn::ndarray&, const bn::ndarray&, bn::ndarray&, bn::ndarray&, bn::ndarray&, const bn::ndarray&, bn::ndarray&, const double, const double, const double, const double, const int, const int);

	//bp::tuple RunUmat(const std::string& umat_name_py, const bn::ndarray& Etot_py, const bn::ndarray& DEtot_py, bn::ndarray& sigma_py, bn::ndarray& Lt_py, bn::ndarray& L_py, const bn::ndarray& DR_py, const bn::ndarray& props_py, bn::ndarray& statev_py, const double T, const double DT, const double Time, const double DTime, const int ndi);
}