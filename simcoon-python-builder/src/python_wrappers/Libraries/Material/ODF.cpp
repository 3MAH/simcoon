#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <stdexcept>
#include <string>
#include <vector>
#include <carma>
#include <armadillo>

#include <simcoon/Simulation/Phase/phase_characteristics.hpp>
#include <simcoon/Continuum_mechanics/Material/peak.hpp>
#include <simcoon/Continuum_mechanics/Material/ODF.hpp>
#include <simcoon/Continuum_mechanics/Material/ODF2Nphases.hpp>
#include <simcoon/Continuum_mechanics/Micromechanics/multiphase.hpp>

#include <simcoon/python_wrappers/dict_get.hpp>
#include <simcoon/python_wrappers/Libraries/Phase/phases.hpp>
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

py::list ODF_discretization(const py::object &phases, const py::object &peaks,
                            const std::string &umat_name, const py::array_t<double> &props_py,
                            const int &num_phase_disc, const int &nphases_disc,
                            const double &angle_min, const double &angle_max,
                            const bool &angles_mat, const int &angle) {

    vec props = carma::arr_to_col(props_py);
    const int shape_type = simcoon::sub_phase_shape(umat_name);
    if (shape_type == 0) {
        throw std::invalid_argument(umat_name + " is a homogeneous model: only a mean-field RVE has sub-phases to discretise");
    }

    simcoon::phase_characteristics rve_init;
    rve_init.sptr_matprops->update(0, umat_name, 1, 0., 0., 0., props.n_elem, props);
    rve_init.construct(shape_type, 1);   //mechanical only
    rve_init.sub_phases = make_sub_phases(phases, umat_name, 0.);
    if (num_phase_disc < 0 || num_phase_disc >= static_cast<int>(rve_init.sub_phases.size())) {
        throw std::invalid_argument("num_phase_disc = " + to_string(num_phase_disc) + " is outside the "
                                    + to_string(rve_init.sub_phases.size()) + " phases given");
    }

    //The peaks and the limits are degrees; the ODF stores radians.
    simcoon::ODF odf_rve(angle, false, angle_min, angle_max);
    odf_rve.peaks = make_peaks(peaks);
    for (auto &p : odf_rve.peaks) {
        p.mean = simcoon::deg2rad(p.mean);
        p.s_dev = simcoon::deg2rad(p.s_dev);
        p.width = simcoon::deg2rad(p.width);
    }

    simcoon::phase_characteristics rve = simcoon::discretize_ODF(rve_init, odf_rve, num_phase_disc, nphases_disc, angles_mat ? 1 : 0);
    return phases_to_list(rve.sub_phases);
}

} //namespace simpy
