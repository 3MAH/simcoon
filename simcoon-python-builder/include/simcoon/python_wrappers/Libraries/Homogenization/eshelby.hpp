#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace simpy{

//Eshelby tensor for a sphere
pybind11::array_t<double> Eshelby_sphere(const double &nu, const bool &copy=true);

//	Eshelby tensor determination. The cylinder is oriented in such a way that the axis direction is the 1 direction. a2=a3 here
pybind11::array_t<double> Eshelby_cylinder(const double &nu, const bool &copy=true);

//	Eshelby tensor determination. The prolate shape is oriented in such a way that the axis direction is the 1 direction. a1>a2=a3 here
pybind11::array_t<double> Eshelby_prolate(const double &nu, const double &aspect_ratio, const bool &copy=true);

//	Eshelby tensor determination. The oblate shape is oriented in such a way that the axis direction is the 1 direction. a1<a2=a3 here
pybind11::array_t<double> Eshelby_oblate(const double &nu, const double &aspect_ratio, const bool &copy=true);

//	Eshelby tensor determination for a penny-shaped crack (limit of oblate as ar->0). The crack normal is the 1 direction.
pybind11::array_t<double> Eshelby_penny(const double &nu, const bool &copy=true);

//Numerical Eshelby tensor determination
pybind11::array_t<double> Eshelby(const pybind11::array_t<double> &L, const double &a1=1., const double &a2=1., const double &a3=1., const int &mp=50, const int &np=50, const bool &copy=true);

//Numerical Hill Interaction tensor determination
pybind11::array_t<double> T_II(const pybind11::array_t<double> &L, const double &a1=1., const double &a2=1., const double &a3=1., const int &mp=50, const int &np=50, const bool &copy=true);

} //namespace simcoon
