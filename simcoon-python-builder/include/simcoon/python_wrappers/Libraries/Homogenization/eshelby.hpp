#pragma once
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

namespace simpy{

//Eshelby tensor for a sphere
boost::python::numpy::ndarray Eshelby_sphere(const double &);

//	Eshelby tensor determination. The cylinder is oriented in such a way that the axis direction is the 1 direction. a2=a3 here
boost::python::numpy::ndarray Eshelby_cylinder(const double &);

//	Eshelby tensor determination. The prolate shape is oriented in such a way that the axis direction is the 1 direction. a1>a2=a3 here
boost::python::numpy::ndarray Eshelby_prolate(const double &, const double &);

//	Eshelby tensor determination. The oblate shape is oriented in such a way that the axis direction is the 1 direction. a1<a2=a3 here
boost::python::numpy::ndarray Eshelby_oblate(const double &, const double &);

//Numerical Eshelby tensor determination
boost::python::numpy::ndarray Eshelby(const boost::python::numpy::ndarray &, const double &, const double &, const double &, const int &, const int &);

//Numerical Hill Interaction tensor determination
boost::python::numpy::ndarray T_II(const boost::python::numpy::ndarray &, const double &, const double &, const double &, const int &, const int &);

} //namespace simcoon
