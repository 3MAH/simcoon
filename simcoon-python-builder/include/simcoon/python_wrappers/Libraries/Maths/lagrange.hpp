#pragma once
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

namespace simpy{

//This function is used to determine an exponential Lagrange Multiplier (like contact in Abaqus)
double lagrange_exp(const double &, const double &, const double &);

//This function is used to determine the first derivative of an exponential Lagrange Multiplier
double dlagrange_exp(const double &, const double &, const double &);

//This function is used to determine a power-law Lagrange Multiplier for problem such x >= 0
double lagrange_pow_0(const double &, const double &, const double &, const double &, const double &);

//This function is used to determine the first derivative of a power-law Lagrange Multiplier for problem such x >= 0
double dlagrange_pow_0(const double &, const double &, const double &, const double &, const double &);

//This function is used to determine a power-law Lagrange Multiplier for problem such x <= 1
double lagrange_pow_1(const double &, const double &, const double &, const double &, const double &);

//This function is used to determine the first derivative of a power-law Lagrange Multiplier for problem such x <= 1
double dlagrange_pow_1(const double &, const double &, const double &, const double &, const double &);

//This function is used to determine the SECOND derivative of a power-law Lagrange Multiplier for problem such x <= 1
double d2lagrange_pow_1(const double &, const double &, const double &, const double &, const double &);
    
} //namespace simpy