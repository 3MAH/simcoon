// Armadillo memory through numpy's allocator, for every translation unit of a module.
//
// Force-included (MSVC /FI, GCC/Clang -include) into every translation unit of _core and,
// on Windows, of libsimcoon as well: armadillo buffers are stolen into numpy arrays
// (carma, zero-copy), so whoever allocates them must use the allocator numpy frees with,
// and on Windows the two DLLs must agree or the heap is corrupted (carma issue #89).
//
// This replaces carma's cnalloc.h, which did the same with a hidden flaw. Numpy's C-API
// table (PyArray_API) was a STATIC variable of each translation unit, imported lazily by
// the first allocation made in that unit through _import_array(), a Python call that needs
// the GIL. The solver engine runs with the GIL released (solver_run binding), so the first
// allocation of a not-yet-touched unit of libsimcoon inside the solver was an access
// violation (PyImport with a NULL thread state). Whether it hit depended on which units the
// bindings called earlier in the process had already warmed with the GIL held.
//
// Here the table is ONE variable per module (PY_ARRAY_UNIQUE_SYMBOL), defined by that
// module's owner translation unit (compiled with SIMCOON_NUMPY_API_OWNER: python_module.cpp
// for _core, numpy_alloc.cpp for libsimcoon) and imported once, at _core import time, with
// the GIL held. After that no allocation ever calls into Python: PyDataMem_NEW/FREE are
// plain malloc/free plus tracemalloc bookkeeping, and tracemalloc takes the GIL itself.
#pragma once

// pybind11 before numpy: numpy's headers would otherwise pull pythonXX_d.lib on a
// Py_DEBUG build (pybind11 issue #1295).
#include <pybind11/pybind11.h>
#ifndef NPY_NO_DEPRECATED_API
#define NPY_NO_DEPRECATED_API NPY_1_14_API_VERSION
#endif
#define PY_ARRAY_UNIQUE_SYMBOL simcoon_numpy_ARRAY_API
#ifndef SIMCOON_NUMPY_API_OWNER
#define NO_IMPORT_ARRAY
#endif
#include <numpy/arrayobject.h>

#include <cstddef>

namespace simcoon {
namespace numpy_alloc {

// Import numpy's C-API table into this module. The GIL must be held. Throws
// pybind11::error_already_set when numpy cannot be imported. Defined once per module, in
// its owner translation unit (below).
void import_api();

namespace detail {

// PyGILState_Ensure/Release as RAII: valid from a thread that released the GIL and from a
// thread Python has never seen.
struct gil_state {
    PyGILState_STATE state;
    gil_state() : state(PyGILState_Ensure()) {}
    ~gil_state() { PyGILState_Release(state); }
    gil_state(const gil_state &) = delete;
    gil_state &operator=(const gil_state &) = delete;
};

// Safety net only: _core imports the table before anything allocates. It stays for an
// allocation made before that (a C++ program using a bindings-enabled libsimcoon inside an
// embedded interpreter): the import then takes the GIL itself instead of crashing.
inline void ensure_api() {
    if (PyArray_API == nullptr) {
        gil_state gil;
        import_api();
    }
}

} // namespace detail

inline void *npy_malloc(std::size_t bytes) {
    detail::ensure_api();
    return PyDataMem_NEW(bytes);
}

inline void npy_free(void *ptr) {
    detail::ensure_api();
    PyDataMem_FREE(ptr);
}

#ifdef SIMCOON_NUMPY_API_OWNER
// The one definition per module: this header is included once per translation unit and
// SIMCOON_NUMPY_API_OWNER is set on a single translation unit of each module. Only the
// owner has _import_array() (NO_IMPORT_ARRAY strips it everywhere else).
void import_api() {
    if (_import_array() < 0) {
        throw pybind11::error_already_set();
    }
}
#endif

} // namespace numpy_alloc
} // namespace simcoon

#define ARMA_ALIEN_MEM_ALLOC_FUNCTION simcoon::numpy_alloc::npy_malloc
#define ARMA_ALIEN_MEM_FREE_FUNCTION simcoon::numpy_alloc::npy_free
// Tells <carma> the allocator is set, so it does not include its own cnalloc.h.
#define CARMA_ARMA_ALIEN_MEM_FUNCTIONS_SET
