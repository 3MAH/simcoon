#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <stdexcept>
#include <string>
#include <vector>
// Before <armadillo>, as in every other _core translation unit: carma defines the
// ARMA_ALIEN_MEM macros so armadillo allocates through numpy. A TU that omits it
// instantiates the same header-only armadillo memory symbols WITHOUT them, and the
// linker keeps a single definition module-wide (ODR) — mixing allocators.
#include <carma>
#include <armadillo>

#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Functions/natural_basis.hpp>
#include <simcoon/Simulation/Geometry/ellipsoid.hpp>
#include <simcoon/Simulation/Geometry/layer.hpp>
#include <simcoon/Simulation/Phase/phase_characteristics.hpp>
#include <simcoon/Simulation/Phase/state_variables_M.hpp>

#include <simcoon/python_wrappers/Libraries/Phase/phases.hpp>

using namespace std;
using namespace arma;
namespace py=pybind11;

namespace simpy{

namespace {

double dict_double(const py::dict &d, const char *key, const double &fallback) {
    return (d.contains(key) && !d[key].is_none()) ? d[key].cast<double>() : fallback;
}

int dict_int(const py::dict &d, const char *key, const int &fallback) {
    return (d.contains(key) && !d[key].is_none()) ? d[key].cast<int>() : fallback;
}

//Euler angles arrive in degrees, as they were written in the .dat files; everything below
//material_characteristics and the geometry classes is in radians.
void dict_angles(const py::dict &d, const char *key, double &psi, double &theta, double &phi) {
    psi = theta = phi = 0.;
    if (!d.contains(key) || d[key].is_none()) {
        return;
    }
    py::dict angles = d[key].cast<py::dict>();
    psi = simcoon::deg2rad(dict_double(angles, "psi", 0.));
    theta = simcoon::deg2rad(dict_double(angles, "theta", 0.));
    phi = simcoon::deg2rad(dict_double(angles, "phi", 0.));
}

vec dict_props(const py::dict &d, const int &number) {
    if (!d.contains("props")) {
        throw std::invalid_argument("phase " + to_string(number) + ": no 'props' entry");
    }
    auto arr = py::array_t<double, py::array::c_style | py::array::forcecast>::ensure(d["props"]);
    if (!arr) {
        throw std::invalid_argument("phase " + to_string(number) + ": 'props' is not a sequence of numbers");
    }
    //Copying constructor: the buffer belongs to Python and must not be aliased into arma.
    return vec(static_cast<const double *>(arr.data()), arr.size());
}

} //anonymous namespace

int shape_type_of(const std::string &umat_name) {
    if (umat_name == "MIHEN" || umat_name == "MIMTN" || umat_name == "MISCN") {
        return 2;
    }
    if (umat_name == "MIPLN") {
        return 1;
    }
    return 0;
}

std::vector<simcoon::phase_characteristics> make_sub_phases(const py::object &phases,
                                                            const std::string &umat_name,
                                                            const int &announced_nphases,
                                                            const double &T_init) {

    const int shape_type = shape_type_of(umat_name);
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
    if (nphases == 0) {
        throw std::invalid_argument(umat_name + ": `phases` is empty");
    }
    if (announced_nphases >= 0 && announced_nphases != nphases) {
        throw std::invalid_argument("props[0] announces " + to_string(announced_nphases)
                                    + " phases but `phases` holds " + to_string(nphases));
    }

    //A local holder gives us the sub_phases vector fully constructed (geometry, multi and
    //state variables), without needing the RVE it will later be attached to.
    simcoon::phase_characteristics holder;
    holder.sub_phases_construct(nphases, shape_type, 1);

    for (int i = 0; i < nphases; i++) {
        py::dict p = seq[i].cast<py::dict>();
        //By reference: read.cpp iterated its sub_phases BY VALUE and only got away with it
        //because every member is a shared_ptr. Do not copy the phase to fill it.
        simcoon::phase_characteristics &sub = holder.sub_phases[i];

        const int number = dict_int(p, "number", i);
        const vec props_i = dict_props(p, number);
        const int nstatev_i = dict_int(p, "nstatev", 1);

        double psi_mat, theta_mat, phi_mat;
        dict_angles(p, "material_orientation", psi_mat, theta_mat, phi_mat);

        sub.sptr_matprops->resize(props_i.n_elem);
        sub.sptr_matprops->update(number, p["umat_name"].cast<std::string>(), dict_int(p, "save", 1),
                                  psi_mat, theta_mat, phi_mat, props_i.n_elem, props_i);

        simcoon::natural_basis nb;
        sub.sptr_sv_global->update(zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), eye(3,3), eye(3,3), eye(3,3), eye(3,3), eye(3,3), eye(3,3), T_init, 0., nstatev_i, zeros(nstatev_i), zeros(nstatev_i), nb);
        sub.sptr_sv_local->update(zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), eye(3,3), eye(3,3), eye(3,3), eye(3,3), eye(3,3), eye(3,3), T_init, 0., nstatev_i, zeros(nstatev_i), zeros(nstatev_i), nb);

        sub.sptr_shape->concentration = dict_double(p, "concentration", 0.);

        double psi_geom, theta_geom, phi_geom;
        dict_angles(p, "geometry_orientation", psi_geom, theta_geom, phi_geom);

        if (shape_type == 2) {
            auto sptr_ellipsoid = std::dynamic_pointer_cast<simcoon::ellipsoid>(sub.sptr_shape);
            sptr_ellipsoid->coatingof = dict_int(p, "coatingof", 0);
            py::dict axes = (p.contains("semi_axes") && !p["semi_axes"].is_none())
                          ? p["semi_axes"].cast<py::dict>() : p;
            sptr_ellipsoid->a1 = dict_double(axes, "a1", 1.);
            sptr_ellipsoid->a2 = dict_double(axes, "a2", 1.);
            sptr_ellipsoid->a3 = dict_double(axes, "a3", 1.);
            sptr_ellipsoid->psi_geom = psi_geom;
            sptr_ellipsoid->theta_geom = theta_geom;
            sptr_ellipsoid->phi_geom = phi_geom;
        }
        else if (shape_type == 1) {
            auto sptr_layer = std::dynamic_pointer_cast<simcoon::layer>(sub.sptr_shape);
            sptr_layer->layerup = dict_int(p, "layerup", -1);
            sptr_layer->layerdown = dict_int(p, "layerdown", -1);
            sptr_layer->psi_geom = psi_geom;
            sptr_layer->theta_geom = theta_geom;
            sptr_layer->phi_geom = phi_geom;
        }
    }

    //Fill the coatedby parameter, as read_ellipsoid did once every phase was known.
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
