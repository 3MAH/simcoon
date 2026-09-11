#pragma once
#include <pybind11/pybind11.h>

namespace simpy {

/// @brief Bind the Python-callback UMAT ("PYEXT") machinery on module `m`:
///        the `StepCut` exception type, `register_python_umat(fn)`,
///        `unregister_python_umat()` and `has_python_umat()`.
///
/// The registered Python callable is invoked, with the GIL acquired, by the
/// C++ dispatchers through simcoon::umat_callback_M. Call once from PYBIND11_MODULE.
void init_pyumat(pybind11::module_ &m);

} // namespace simpy
