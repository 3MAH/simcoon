#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

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
        const std::string umat_i = dict_string(p, "umat_name", number);
        if (simcoon::sub_phase_shape(umat_i) != 0) {
            throw std::invalid_argument("phase " + to_string(number) + " is itself a mean-field "
                                        "model (" + umat_i + "): nested composites are not "
                                        "supported through `phases`");
        }
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

        sub.sptr_shape->concentration = dget(p, "concentration", 0.);

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

} //namespace simpy
