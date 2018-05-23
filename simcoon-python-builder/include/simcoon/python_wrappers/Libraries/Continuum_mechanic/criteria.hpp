#pragma once
#include <string>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

namespace simpy{
    
//This function returns the Prager equivalent stress.
double Prager_stress(const boost::python::numpy::ndarray &, const double &, const double &);

//This function returns the derivative of the Prager equivalent stress.
boost::python::numpy::ndarray dPrager_stress(const boost::python::numpy::ndarray &, const double &, const double &);

//This function returns the Prager equivalent stress.
double Tresca_stress(const boost::python::numpy::ndarray &);

//This function returns the derivative of the Prager equivalent stress.
boost::python::numpy::ndarray dTresca_stress(const boost::python::numpy::ndarray &);

//This function returns the Hill equivalent stress.
double Hill_stress(const boost::python::numpy::ndarray &, const boost::python::numpy::ndarray &);

//This function returns the derivative of the Hill equivalent stress.
boost::python::numpy::ndarray dHill_stress(const boost::python::numpy::ndarray &, const boost::python::numpy::ndarray &);
    
//This function returns the Ani equivalent stress.
double Ani_stress(const boost::python::numpy::ndarray &, const boost::python::numpy::ndarray &);

//This function returns the derivative of the Ani equivalent stress.
boost::python::numpy::ndarray dAni_stress(const boost::python::numpy::ndarray &, const boost::python::numpy::ndarray &);
    
//This function computes the selected equivalent stress function
double Eq_stress(const boost::python::numpy::ndarray &, const std::string &, const boost::python::numpy::ndarray &);

//This function computes the deriavtive of the selected equivalent stress function
boost::python::numpy::ndarray dEq_stress(const boost::python::numpy::ndarray &, const std::string &, const boost::python::numpy::ndarray &);
    
} //namespace simpy