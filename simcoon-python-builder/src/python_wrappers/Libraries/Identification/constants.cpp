///@file constants.cpp
///@brief Handle of input constants exposed in python
///@version 0.9


#include <iostream>
#include <fstream>
#include <assert.h>
#include <math.h>
#include <armadillo>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

#include <simcoon/arma2numpy/list_vector.hpp>
#include <simcoon/arma2numpy/numpy_arma.hpp>
#include <simcoon/Simulation/Identification/constants.hpp>
#include <simcoon/python_wrappers/Libraries/Identification/constants.hpp>

using namespace std;
using namespace arma;
using namespace arma2numpy;
namespace bp = boost::python;
namespace bn = boost::python::numpy;

namespace simpy{
    
//-------------------------------------------------------------
simcoon::constants build_constants_full(const int &mnumber, const double &mvalue, const bn::ndarray &minput_values, const std::string &mkey, const int &mninput_files, const bp::list &minput_files)
//-------------------------------------------------------------
{
    simcoon::constants a;
    a.number = mnumber;
    a.value = mvalue;
    a.input_values = array2vec(minput_values);
    a.key = mkey;
    a.ninput_files = mninput_files;
    a.input_files = py_list_to_std_vector_string(minput_files);
    return a;
}

//------------------------------------------------------
bn::ndarray constants_get_input_values(simcoon::constants &self) {
    return vec2array(self.input_values);
}
//------------------------------------------------------

//------------------------------------------------------
bp::list constants_get_input_files(simcoon::constants &self) {
    return std_vector_to_py_list_string(self.input_files);
}
//------------------------------------------------------

//------------------------------------------------------
void constants_set_input_values(simcoon::constants &self, const bn::ndarray &minput_values) {
    self.input_values = array2vec(minput_values);
}
//------------------------------------------------------

//------------------------------------------------------
void constants_set_input_files(simcoon::constants &self, const bp::list &minput_values) {
    self.input_files = py_list_to_std_vector_string(minput_values);
    self.ninput_files = self.input_files.size();
}
//------------------------------------------------------
    
    
} //namespace simpy