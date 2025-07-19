#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <string>
#include <carma>
#include <armadillo>
#include <assert.h>

#include <simcoon/Simulation/Identification/parameters.hpp>
#include <simcoon/python_wrappers/Libraries/Identification/parameters.hpp>

using namespace std;
using namespace arma;
namespace py=pybind11;

namespace simpy{
    
//-------------------------------------------------------------
simcoon::parameters build_parameters_full(const int &mnumber, const double &mmin_value, const double &mmax_value, const std::string &mkey, const int &mninput_files, const py::list &minput_files)
//-------------------------------------------------------------
{
    simcoon::parameters a;
    a.number = mnumber;
    a.min_value = mmin_value;
    a.max_value = mmax_value;
    a.value = (a.min_value+a.max_value)/2.;
    a.key = mkey;
    a.ninput_files = mninput_files;
    a.input_files = minput_files.cast<std::vector<std::string>>();
    return a;
}

//------------------------------------------------------
py::list parameters_get_input_files(simcoon::parameters &self) {
    py::list list_to_return = py::cast(self.input_files);
    return list_to_return;
}
//------------------------------------------------------

//------------------------------------------------------
void parameters_set_input_files(simcoon::parameters &self, const py::list &minput_values) {
    self.input_files = minput_values.cast<std::vector<std::string>>();
    self.ninput_files = self.input_files.size();
}
//------------------------------------------------------
    
} //namespace simpy