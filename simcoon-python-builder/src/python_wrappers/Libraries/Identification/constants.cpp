#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <string>
#include <armadillo>
#include <simcoon/python_wrappers/conversion_helpers.hpp>
#include <assert.h>

#include <simcoon/Simulation/Identification/constants.hpp>
#include <simcoon/python_wrappers/Libraries/Identification/constants.hpp>

using namespace std;
using namespace arma;
namespace py=pybind11;

namespace simpy{
    
//-------------------------------------------------------------
simcoon::constants build_constants_full(const int &mnumber, const double &mvalue, const py::array_t<double> &minput_values, const std::string &mkey, const int &mninput_files, const py::list &minput_files)
//-------------------------------------------------------------
{
    simcoon::constants a;
    a.number = mnumber;
    a.value = mvalue;
    a.input_values = simpy::arr_to_col(minput_values);
    a.key = mkey;
    a.ninput_files = mninput_files;
    a.input_files = minput_files.cast<std::vector<std::string>>();
    return a;
}

//------------------------------------------------------
py::array_t<double> constants_get_input_values(simcoon::constants &self) {
    return simpy::col_to_arr(self.input_values);
}
//------------------------------------------------------

//------------------------------------------------------
py::list constants_get_input_files(simcoon::constants &self) {
    py::list list_to_return = py::cast(self.input_files);
    return list_to_return;
}
//------------------------------------------------------

//------------------------------------------------------
void constants_set_input_values(simcoon::constants &self, const py::array_t<double> &minput_values) {
    self.input_values = simpy::arr_to_col(minput_values);
}
//------------------------------------------------------

//------------------------------------------------------
void constants_set_input_files(simcoon::constants &self, const py::list &minput_files) {
    self.input_files = minput_files.cast<std::vector<std::string>>();
    self.ninput_files = self.input_files.size();
}
//------------------------------------------------------
    
    
} //namespace simpy