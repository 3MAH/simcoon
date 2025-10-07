#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace simpy{
    
//This function returns damage evolution (/dt) considering a Weibull damage law
double damage_weibull(const pybind11::array_t<double> &stress, const double &damage, const double &alpha, const double &beta, const double &DTime, const std::string &criterion = "vonmises");

//This function returns damage evolution (/dt) considering Kachanov's creep damage law
double damage_kachanov(const pybind11::array_t<double> &stress, const pybind11::array_t<double> &strain, const double &damage, const double &A0, const double &r, const std::string &criterion);

//This function returns the constant damage evolution (/dN) considering Woehler- Miner's damage law
double damage_miner(const double &S_max, const double &S_mean, const double &S_ult, const double &b, const double &B0, const double &beta, const double &Sl_0 = 0.);

//This function returns the constant damage evolution (/dN) considering Coffin-Manson's damage law
double damage_manson(const double &S_amp, const double &C2, const double &gamma2);
    
} //namespace simpy
