#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>
//<carma> first: ARMA_ALIEN_MEM routes armadillo through numpy's allocator, in every _core TU.
#include <carma>
#include <armadillo>

#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Micromechanics/multiphase.hpp>
#include <simcoon/Simulation/Geometry/ellipsoid.hpp>
#include <simcoon/Simulation/Geometry/layer.hpp>
#include <simcoon/Simulation/Phase/phase_characteristics.hpp>
#include <simcoon/Simulation/Phase/state_variables_M.hpp>

#include <simcoon/python_wrappers/dict_get.hpp>
#include <simcoon/python_wrappers/Libraries/Phase/phases.hpp>

using namespace std;
using namespace arma;
namespace py=pybind11;

namespace simpy{

namespace {

void dict_angles(const py::dict &d, const char *key, double &psi, double &theta, double &phi) {
    psi = theta = phi = 0.;
    if (!d.contains(key) || d[key].is_none()) {
        return;
    }
    py::dict angles = d[key].cast<py::dict>();
    psi = simcoon::deg2rad(dget(angles, "psi", 0.));
    theta = simcoon::deg2rad(dget(angles, "theta", 0.));
    phi = simcoon::deg2rad(dget(angles, "phi", 0.));
}

std::string dict_string(const py::dict &d, const char *key, const int &number) {
    if (!d.contains(key) || d[key].is_none()) {
        throw std::invalid_argument("phase " + to_string(number) + ": no '" + key + "' entry");
    }
    return d[key].cast<std::string>();
}

vec dict_props(const py::dict &d, const int &number) {
    if (!d.contains("props")) {
        throw std::invalid_argument("phase " + to_string(number) + ": no 'props' entry");
    }
    auto arr = py::array_t<double, py::array::c_style | py::array::forcecast>::ensure(d["props"]);
    if (!arr) {
        throw std::invalid_argument("phase " + to_string(number) + ": 'props' is not a sequence of numbers");
    }
    if (arr.size() == 0) {
        //update() only asserts, and asserts are compiled out of the release wheels: props(0)
        //would then read out of bounds.
        throw std::invalid_argument("phase " + to_string(number) + ": 'props' is empty");
    }
    //Copying constructor: the buffer belongs to Python and must not be aliased into arma.
    return vec(static_cast<const double *>(arr.data()), arr.size());
}

} //anonymous namespace

std::vector<simcoon::phase_characteristics> make_sub_phases(const py::object &phases,
                                                            const std::string &umat_name,
                                                            const double &T_init) {

    const int shape_type = simcoon::sub_phase_shape(umat_name);
    const bool given = static_cast<bool>(phases) && !phases.is_none();

    if (shape_type == 0) {
        if (given) {
            throw std::invalid_argument(umat_name + " is a homogeneous model and takes no sub-phases");
        }
        return {};
    }
    if (!given) {
        throw std::invalid_argument(umat_name + " is a mean-field model: give its sub-phases through "
                                    "`phases` (the Nellipsoids/Nlayers files are no longer read)");
    }

    py::sequence seq = phases.cast<py::sequence>();
    const int nphases = static_cast<int>(py::len(seq));

    //A holder yields the sub_phases fully constructed, without the RVE they attach to.
    simcoon::phase_characteristics holder;
    holder.sub_phases_construct(nphases, shape_type, 1);

    const char *expected_kind = (shape_type == 2) ? "ellipsoid" : "layer";
    double total_concentration = 0.;

    for (int i = 0; i < nphases; i++) {
        py::dict p = seq[i].cast<py::dict>();
        //By reference: the phase is filled in place.
        simcoon::phase_characteristics &sub = holder.sub_phases[i];

        const int number = dget(p, "number", i);
        if (number != i) {
            throw std::invalid_argument("phase at position " + to_string(i) + " carries number "
                                        + to_string(number) + ": the phases must be numbered by "
                                        "their position in the list (0, 1, ... n-1)");
        }
        //The geometry the dict describes must be the one the model builds: an ellipsoidal
        //scheme fed layers or cylinders would otherwise homogenise unit spheres, silently.
        if (p.contains("kind") && !p["kind"].is_none()) {
            const std::string kind = p["kind"].cast<std::string>();
            if (kind != expected_kind) {
                throw std::invalid_argument("phase " + to_string(number) + " is a " + kind + ", but "
                                            + umat_name + " builds " + expected_kind + "s");
            }
        }
        const std::string umat_i = dict_string(p, "umat_name", number);
        const vec props_i = dict_props(p, number);
        const int nstatev_i = dget(p, "nstatev", 1);

        double psi_mat, theta_mat, phi_mat;
        dict_angles(p, "material_orientation", psi_mat, theta_mat, phi_mat);

        sub.sptr_matprops->update(number, umat_i, dget(p, "save", 1),
                                  psi_mat, theta_mat, phi_mat, props_i.n_elem, props_i);

        //A freshly constructed state is already the zero state at DT = 0: only the state
        //variables' size and the temperature are to set.
        for (auto &sv : {sub.sptr_sv_global, sub.sptr_sv_local}) {
            sv->resize(nstatev_i);
            sv->T = T_init;
        }

        if (!p.contains("concentration") || p["concentration"].is_none()) {
            throw std::invalid_argument("phase " + to_string(number) + ": no 'concentration' entry");
        }
        sub.sptr_shape->concentration = p["concentration"].cast<double>();
        total_concentration += sub.sptr_shape->concentration;

        //A mean-field sub-phase carries its own sub-phases; a homogeneous one must not.
        const py::object nested = p.contains("phases") ? p["phases"].cast<py::object>()
                                                       : py::object(py::none());
        sub.sub_phases = make_sub_phases(nested, umat_i, T_init);

        double psi_geom, theta_geom, phi_geom;
        dict_angles(p, "geometry_orientation", psi_geom, theta_geom, phi_geom);

        if (shape_type == 2) {
            auto sptr_ellipsoid = std::dynamic_pointer_cast<simcoon::ellipsoid>(sub.sptr_shape);
            const int coatingof = dget(p, "coatingof", 0);
            if (coatingof < 0 || coatingof >= nphases) {
                throw std::invalid_argument("phase " + to_string(number) + ": 'coatingof' is "
                                            + to_string(coatingof) + ", outside the "
                                            + to_string(nphases) + " phases given");
            }
            sptr_ellipsoid->coatingof = coatingof;
            py::dict axes = (p.contains("semi_axes") && !p["semi_axes"].is_none())
                          ? p["semi_axes"].cast<py::dict>() : p;
            sptr_ellipsoid->a1 = dget(axes, "a1", 1.);
            sptr_ellipsoid->a2 = dget(axes, "a2", 1.);
            sptr_ellipsoid->a3 = dget(axes, "a3", 1.);
            sptr_ellipsoid->psi_geom = psi_geom;
            sptr_ellipsoid->theta_geom = theta_geom;
            sptr_ellipsoid->phi_geom = phi_geom;
        }
        else if (shape_type == 1) {
            auto sptr_layer = std::dynamic_pointer_cast<simcoon::layer>(sub.sptr_shape);
            //0/0, as layer::layer() leaves them: a -1 default would change the stacking state.
            sptr_layer->layerup = dget(p, "layerup", 0);
            sptr_layer->layerdown = dget(p, "layerdown", 0);
            sptr_layer->psi_geom = psi_geom;
            sptr_layer->theta_geom = theta_geom;
            sptr_layer->phi_geom = phi_geom;
        }
    }

    if (std::abs(total_concentration - 1.) > 1.e-6) {
        throw std::invalid_argument(umat_name + ": the concentrations of the " + to_string(nphases)
                                    + " phases sum to " + to_string(total_concentration)
                                    + ", not 1");
    }

    //coatedby is only resolvable once every phase is known.
    if (shape_type == 2) {
        for (int i = 0; i < nphases; i++) {
            auto sptr_ellipsoid = std::dynamic_pointer_cast<simcoon::ellipsoid>(holder.sub_phases[i].sptr_shape);
            if (sptr_ellipsoid->coatingof != 0) {
                auto sptr_coated = std::dynamic_pointer_cast<simcoon::ellipsoid>(holder.sub_phases[sptr_ellipsoid->coatingof].sptr_shape);
                sptr_coated->coatedby = i;
            }
        }
    }

    return holder.sub_phases;
}

py::list phases_to_list(const std::vector<simcoon::phase_characteristics> &phases) {
    py::list out;
    for (size_t i = 0; i < phases.size(); i++) {
        const simcoon::phase_characteristics &ph = phases[i];
        const auto &mp = *ph.sptr_matprops;
        py::dict d;
        d["number"] = mp.number;
        d["umat_name"] = mp.umat_name;
        d["save"] = mp.save;
        d["concentration"] = ph.sptr_shape->concentration;
        py::dict mat_angles;
        mat_angles["psi"] = simcoon::rad2deg(mp.psi_mat);
        mat_angles["theta"] = simcoon::rad2deg(mp.theta_mat);
        mat_angles["phi"] = simcoon::rad2deg(mp.phi_mat);
        d["material_orientation"] = mat_angles;
        d["nstatev"] = ph.sptr_sv_global->nstatev;
        vec props_copy = mp.props;   //a copy numpy owns
        d["props"] = carma::col_to_arr(props_copy, true);
        if (auto ell = std::dynamic_pointer_cast<simcoon::ellipsoid>(ph.sptr_shape)) {
            d["kind"] = "ellipsoid";
            d["coatingof"] = ell->coatingof;
            py::dict axes;
            axes["a1"] = ell->a1; axes["a2"] = ell->a2; axes["a3"] = ell->a3;
            d["semi_axes"] = axes;
            py::dict geo;
            geo["psi"] = simcoon::rad2deg(ell->psi_geom);
            geo["theta"] = simcoon::rad2deg(ell->theta_geom);
            geo["phi"] = simcoon::rad2deg(ell->phi_geom);
            d["geometry_orientation"] = geo;
        }
        else if (auto lay = std::dynamic_pointer_cast<simcoon::layer>(ph.sptr_shape)) {
            d["kind"] = "layer";
            d["layerup"] = lay->layerup;
            d["layerdown"] = lay->layerdown;
            py::dict geo;
            geo["psi"] = simcoon::rad2deg(lay->psi_geom);
            geo["theta"] = simcoon::rad2deg(lay->theta_geom);
            geo["phi"] = simcoon::rad2deg(lay->phi_geom);
            d["geometry_orientation"] = geo;
        }
        else {
            d["kind"] = "phase";
        }
        if (!ph.sub_phases.empty()) {
            d["phases"] = phases_to_list(ph.sub_phases);
        }
        out.append(d);
    }
    return out;
}

} //namespace simpy
