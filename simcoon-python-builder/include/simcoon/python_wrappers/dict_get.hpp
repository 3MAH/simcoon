#pragma once
#include <initializer_list>
#include <stdexcept>
#include <string>
#include <armadillo>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <simcoon/parameter.hpp>

namespace simpy {

//A dict entry with a default: absent or None both mean "not given".
template <typename T>
T dget(const pybind11::dict &d, const char *key, const T &dflt) {
    return (d.contains(key) && !d[key].is_none()) ? d[key].cast<T>() : dflt;
}

//A misspelt key would otherwise fall back on its default, silently: a 'semi_axis' entry
//turned an inclusion into a sphere, a 'psy' angle into 0.
inline void check_keys(const pybind11::dict &d, std::initializer_list<const char *> allowed, const std::string &what) {
    for (auto item : d) {
        const std::string key = pybind11::str(item.first);
        bool known = false;
        for (const char *a : allowed) {
            if (key == a) { known = true; break; }
        }
        if (!known) {
            std::string list;
            for (const char *a : allowed) { list += std::string(list.empty() ? "" : ", ") + a; }
            throw std::invalid_argument(what + ": unknown entry '" + key + "' (known: " + list + ")");
        }
    }
}

//Euler angles {psi, theta, phi} in DEGREES (the convention of every dict of the bindings),
//returned in radians; None is the identity.
inline void angles_of(const pybind11::object &angles, const std::string &what, double &psi, double &theta, double &phi) {
    psi = theta = phi = 0.;
    if (!static_cast<bool>(angles) || angles.is_none()) {
        return;
    }
    const pybind11::dict d = angles.cast<pybind11::dict>();
    check_keys(d, {"psi", "theta", "phi"}, what);
    psi = simcoon::deg2rad(dget(d, "psi", 0.));
    theta = simcoon::deg2rad(dget(d, "theta", 0.));
    phi = simcoon::deg2rad(dget(d, "phi", 0.));
}

//A sequence of numbers as an arma vec (a copy: the buffer belongs to Python). ensure()
//returns a null handle for a ragged list or a string, which must not be dereferenced.
inline arma::vec vec_of(const pybind11::handle &values, const std::string &what) {
    auto arr = pybind11::array_t<double, pybind11::array::c_style | pybind11::array::forcecast>::ensure(values);
    if (!arr) {
        PyErr_Clear();
        throw std::invalid_argument(what + " is not a sequence of numbers");
    }
    return arma::vec(static_cast<const double *>(arr.data()), arr.size());
}

} //namespace simpy
