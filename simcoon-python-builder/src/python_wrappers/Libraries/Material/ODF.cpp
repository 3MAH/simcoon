#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <stdexcept>
#include <string>
#include <vector>
#include <carma>
#include <armadillo>

#include <simcoon/Continuum_mechanics/Material/peak.hpp>
#include <simcoon/Continuum_mechanics/Material/ODF.hpp>
#include <simcoon/Continuum_mechanics/Material/ODF2Nphases.hpp>

#include <simcoon/python_wrappers/dict_get.hpp>
#include <simcoon/python_wrappers/Libraries/Material/ODF.hpp>

using namespace std;
using namespace arma;
namespace py=pybind11;

namespace simpy{

namespace {

//The peaks of an ODF/PDF from their Python description (a sequence of dicts). Angles are
//left as given: the callers say whether they are degrees.
std::vector<simcoon::peak> make_peaks(const py::object &peaks) {
    if (!static_cast<bool>(peaks) || peaks.is_none()) {
        throw std::invalid_argument("peaks: give the peaks of the distribution (a sequence of dicts)");
    }
    py::sequence seq = peaks.cast<py::sequence>();
    std::vector<simcoon::peak> out(py::len(seq));
    for (size_t i = 0; i < out.size(); i++) {
        py::dict d = seq[i].cast<py::dict>();
        out[i].number = dget(d, "number", static_cast<int>(i));
        out[i].method = dget(d, "method", 0);
        out[i].mean = dget(d, "mean", 0.);
        out[i].s_dev = dget(d, "s_dev", 0.);
        out[i].width = dget(d, "width", 0.);
        out[i].ampl = dget(d, "ampl", 1.);
        if (d.contains("params") && !d["params"].is_none()) {
            auto arr = py::array_t<double, py::array::c_style | py::array::forcecast>::ensure(d["params"]);
            out[i].params = vec(static_cast<const double *>(arr.data()), arr.size());
        }
    }
    return out;
}

} //anonymous namespace

py::array_t<double> get_densities_ODF(const py::array_t<double> &x_py, const py::object &peaks, const bool &radian) {
    vec x = carma::arr_to_col(x_py);
    vec y = simcoon::get_densities_ODF(x, make_peaks(peaks), radian);
    return carma::col_to_arr(y);
}

} //namespace simpy
