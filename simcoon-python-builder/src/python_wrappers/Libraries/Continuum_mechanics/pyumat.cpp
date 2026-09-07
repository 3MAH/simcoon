/* This file is part of simcoon.

 simcoon is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 simcoon is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with simcoon.  If not, see <http://www.gnu.org/licenses/>.

 */

///@file pyumat.cpp
///@brief Bridge between the C++ "PYEXT" callback slot (umat_callback.hpp) and a Python
///       callable: numpy marshalling, GIL handling, StepCut and error propagation.

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/gil_safe_call_once.h>

#include <algorithm>
#include <cstring>
#include <sstream>
#include <string>
#include <vector>

#include <carma>
#include <armadillo>

#include <simcoon/exception.hpp>
#include <simcoon/Continuum_mechanics/Umat/umat_callback.hpp>
#include <simcoon/python_wrappers/Libraries/Continuum_mechanics/pyumat.hpp>

namespace py = pybind11;
using namespace pybind11::literals;

namespace simpy {

namespace {

struct step_cut_tag {};

// Module-lifetime Python objects. Neither is ever destroyed: a static py::object destructor
// running at interpreter teardown (after Py_Finalize) is a documented pybind11 crash. The
// exception type uses the gil_safe_call_once_and_store idiom; the callable lives in a
// heap slot that is intentionally leaked (unregister assigns None to release the user object).
PYBIND11_CONSTINIT static py::gil_safe_call_once_and_store<py::object> g_step_cut;
static py::object *g_fn = nullptr;   // touched only with the GIL held

using farr = py::array_t<double, py::array::f_style | py::array::forcecast>;

// (rows, cols) column-major copy — cheaper than carma::mat_to_arr(arma::mat(m)), which
// heap-allocates an arma copy plus a capsule for the 9 doubles of DR on every call
farr np2d(const arma::mat &m) {
    farr a({static_cast<py::ssize_t>(m.n_rows), static_cast<py::ssize_t>(m.n_cols)});
    if (m.n_elem > 0) {
        std::memcpy(a.mutable_data(), m.memptr(), m.n_elem * sizeof(double));
    }
    return a;
}

// (n,) copy — carma::col_to_arr would give (n,1)
py::array_t<double> np1d(const arma::vec &v) {
    py::array_t<double> a(static_cast<py::ssize_t>(v.n_elem));
    if (v.n_elem > 0) {
        std::memcpy(a.mutable_data(), v.memptr(), v.n_elem * sizeof(double));
    }
    return a;
}

std::string shape_str(const std::vector<py::ssize_t> &shape) {
    std::ostringstream os;
    os << "(";
    for (size_t i = 0; i < shape.size(); i++) {
        if (i > 0) os << ", ";
        os << shape[i];
    }
    os << (shape.size() == 1 ? ",)" : ")");
    return os.str();
}

std::string shape_str(const farr &a) {
    std::vector<py::ssize_t> s(a.ndim());
    for (py::ssize_t k = 0; k < a.ndim(); k++) s[k] = a.shape(k);
    return shape_str(s);
}

// Normalise one returned value to a float64 Fortran-ordered array of the expected shape.
// Accepts anything numpy can convert (list, float32, CPU torch tensor via __array__).
farr take(py::handle obj, const char *what, const std::vector<py::ssize_t> &shape) {
    farr a = farr::ensure(obj);
    if (!a) {
        throw py::type_error(std::string("PYEXT: returned '") + what +
                             "' is not convertible to a float64 array (return numpy arrays; "
                             "for torch tensors use .detach().cpu().numpy())");
    }
    bool ok = (a.ndim() == static_cast<py::ssize_t>(shape.size()));
    for (size_t k = 0; ok && k < shape.size(); k++) {
        ok = (a.shape(static_cast<py::ssize_t>(k)) == shape[k]);
    }
    if (!ok) {
        throw py::value_error(std::string("PYEXT: returned '") + what + "' must have shape " +
                              shape_str(shape) + ", got " + shape_str(a));
    }
    return a;
}

void python_umat_bridge(const std::string &umat_name, const arma::vec &Etot, const arma::vec &DEtot,
                        arma::vec &sigma, arma::mat &Lt, arma::mat &L, const arma::mat &DR,
                        const int &nprops, const arma::vec &props, const int &nstatev, arma::vec &statev,
                        const double &T, const double &DT, const double &Time, const double &DTime,
                        double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d,
                        const int &ndi, const int &nshr, const bool &start, double &tnew_dt,
                        const int &tangent_mode)
{
    (void)umat_name;
    (void)nprops;
    // Re-entrant on the calling thread: a no-op when the GIL is already held (batch path),
    // an acquire when the solver released it (solver_run).
    py::gil_scoped_acquire gil;

    if (!g_fn || g_fn->is_none()) {
        throw simcoon::exception_solver("PYEXT: no Python UMAT registered");
    }

    py::array_t<double> Wm4(4);
    {
        double *w = Wm4.mutable_data();
        w[0] = Wm; w[1] = Wm_r; w[2] = Wm_ir; w[3] = Wm_d;
    }
    py::object ret;
    try {
        ret = (*g_fn)("Etot"_a = np1d(Etot), "DEtot"_a = np1d(DEtot), "sigma"_a = np1d(sigma),
                      "DR"_a = np2d(DR), "props"_a = np1d(props),
                      "statev"_a = np1d(statev),
                      "T"_a = T, "DT"_a = DT, "Time"_a = Time, "DTime"_a = DTime,
                      "Wm"_a = Wm4, "ndi"_a = ndi, "nshr"_a = nshr, "start"_a = start,
                      "tangent_mode"_a = tangent_mode);
    } catch (py::error_already_set &e) {
        if (e.matches(g_step_cut.get_stored())) {
            // Requested step cut: leave sigma/statev/Wm untouched (the solver restores the
            // start-of-increment state before retrying); it only tests tnew_dt < 1.
            double ratio = 0.5;
            try {
                ratio = e.value().attr("ratio").cast<double>();
            } catch (...) {
            }
            tnew_dt = std::min(std::max(ratio, 1.e-3), 0.99);
            return;
        }
        // Propagate the ORIGINAL Python exception unchanged: never inspect it here
        // (what() needs the GIL); the pybind11 dispatcher restores it on the way out.
        throw;
    }

    if (!(py::isinstance<py::tuple>(ret) || py::isinstance<py::list>(ret))) {
        throw py::type_error("PYEXT: integrate() must return a tuple "
                             "(sigma, Lt, statev, Wm[, L])");
    }
    py::sequence out = py::reinterpret_borrow<py::sequence>(ret);
    const py::ssize_t n = static_cast<py::ssize_t>(py::len(out));
    if (n != 4 && n != 5) {
        throw py::value_error("PYEXT: expected 4 or 5 return values (sigma, Lt, statev, Wm[, L]), got "
                              + std::to_string(n));
    }
    farr s = take(out[0], "sigma", {6});
    farr lt = take(out[1], "Lt", {6, 6});
    farr sv = take(out[2], "statev", {static_cast<py::ssize_t>(nstatev)});
    farr wm = take(out[3], "Wm", {4});

    // Build locals first (f_style => column-major, so mat(ptr, 6, 6) reads the returned layout
    // as is), validate, and only then commit: sigma/Lt/statev/Wm* alias the caller's live state,
    // so a throw after a partial assignment would leave it corrupted.
    const arma::vec sigma_new(s.data(), 6);
    const arma::mat Lt_new(lt.data(), 6, 6);
    const arma::vec statev_new(sv.data(), static_cast<arma::uword>(std::max(nstatev, 0)));
    const arma::vec Wm_new(wm.data(), 4);
    arma::mat L_new;
    if (n == 5) {
        farr l = take(out[4], "L", {6, 6});
        L_new = arma::mat(l.data(), 6, 6);
    } else {
        L_new = Lt_new;
    }
    if (!sigma_new.is_finite() || !Lt_new.is_finite() || !statev_new.is_finite()
            || !Wm_new.is_finite() || !L_new.is_finite()) {
        throw py::value_error(
            "PYEXT: non-finite value in the returned sigma/Lt/statev/Wm/L (the state was left "
            "untouched). Raise simcoon.StepCut instead to ask the material-point solver for a "
            "smaller increment.");
    }
    sigma = sigma_new;
    Lt = Lt_new;
    if (nstatev > 0) {
        statev = statev_new;
    }
    Wm = Wm_new(0);
    Wm_r = Wm_new(1);
    Wm_ir = Wm_new(2);
    Wm_d = Wm_new(3);
    L = L_new;
}

} // namespace

void init_pyumat(py::module_ &m)
{
    py::object step_cut = py::exception<step_cut_tag>(m, "StepCut", PyExc_RuntimeError);
    step_cut.attr("__doc__") = "Raised by a Python UMAT to request a smaller increment "
                               "(see simcoon.pyumat.StepCut).";
    g_step_cut.call_once_and_store_result([&]() { return step_cut; });

    m.def("register_python_umat", [](py::object fn) {
        if (!PyCallable_Check(fn.ptr())) {
            throw py::type_error("register_python_umat: object is not callable");
        }
        if (!g_fn) g_fn = new py::object();
        *g_fn = std::move(fn);
        // Plain function pointer: no Python object is captured in the C++ std::function.
        simcoon::set_umat_callback(&python_umat_bridge);
    }, "fn"_a,
    "Register the Python callable served under the 'PYEXT' UMAT name (keyword call: "
    "Etot, DEtot, sigma, DR, props, statev, T, DT, Time, DTime, Wm, ndi, nshr, start, "
    "tangent_mode -> (sigma, Lt, statev, Wm[, L])). Prefer simcoon.registered().");

    m.def("unregister_python_umat", []() {
        simcoon::clear_umat_callback();
        if (g_fn) *g_fn = py::none();   // release the user object now, keep the slot
    }, "Remove the Python UMAT served under 'PYEXT'.");

    m.def("has_python_umat", []() {
        return g_fn && !g_fn->is_none() && simcoon::has_umat_callback();
    }, "True when a Python UMAT is registered under 'PYEXT'.");

    // Module teardown (GIL held, interpreter alive): drop the user object (e.g. a torch
    // model) deterministically. Correctness does not depend on this.
    m.add_object("_pyumat_cleanup", py::capsule([]() {
        simcoon::clear_umat_callback();
        if (g_fn) *g_fn = py::none();
    }));
}

} // namespace simpy
