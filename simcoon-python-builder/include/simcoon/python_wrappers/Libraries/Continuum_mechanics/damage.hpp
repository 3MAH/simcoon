#pragma once
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

namespace simpy{
    
//This function returns damage evolution (/dt) considering a Weibull damage law
double damage_weibull(const boost::python::numpy::ndarray &, const double &, const double &, const double &, const double &, const std::string& = "vonmises");

//This function returns damage evolution (/dt) considering Kachanov's creep damage law
double damage_kachanov(const boost::python::numpy::ndarray &, const boost::python::numpy::ndarray &, const double &, const double &, const double &, const std::string &);

//This function returns the constant damage evolution (/dN) considering Woehler- Miner's damage law
double damage_miner(const double &, const double &, const double &, const double &, const double &, const double &, const double & =0.);

//This function returns the constant damage evolution (/dN) considering Coffin-Manson's damage law
double damage_manson(const double &, const double &, const double &);
    
} //namespace simpy
