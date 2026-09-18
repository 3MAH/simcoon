// Armadillo memory through numpy's allocator, for every translation unit of a module.
//
// Force-included by CMake (simcoon_use_numpy_allocator) into every translation unit of
// _core and, on Windows, of libsimcoon: armadillo buffers are stolen into numpy arrays
// (carma, zero-copy), so whoever allocates them must use the allocator numpy frees with,
// and on Windows the two DLLs must agree or the heap is corrupted (carma issue #89).
//
// The contract, and why this is not carma's cnalloc.h: numpy's C-API table is ONE
// variable per module (PY_ARRAY_UNIQUE_SYMBOL), defined by the module's owner translation
// unit (compiled with SIMCOON_NUMPY_API_OWNER) and imported once, with the GIL held, when
// _core is imported. Allocations never import it: cnalloc.h kept a static table per
// translation unit and imported it lazily from the first allocation, a Python call that
// needs the GIL, which the solver had released (see test_solver_run.py,
// test_solver_is_safe_as_the_first_call_of_a_process).
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
#include <stdexcept>

namespace simcoon {
namespace numpy_alloc {

// Import numpy's C-API table into this module. The GIL must be held. Throws
// pybind11::error_already_set when numpy cannot be imported. Defined once per module, in
// its owner translation unit (below).
void import_api();

inline void *npy_malloc(std::size_t bytes) {
    if (PyArray_API == nullptr) {
        // Loud on purpose: the lazy import this replaces was the bug. Every supported
        // build imports the table at `import simcoon`; an allocation before that is a
        // sequencing error, not a case to paper over with a GIL-taking import.
        throw std::logic_error("simcoon: numpy's C-API table is not imported in this module; "
                               "import simcoon._core before allocating armadillo memory");
    }
    return PyDataMem_NEW(bytes);
}

// A buffer reaching here came from npy_malloc of the same module: the table is imported.
inline void npy_free(void *ptr) {
    PyDataMem_FREE(ptr);
}

#ifdef SIMCOON_NUMPY_API_OWNER
// The one definition per module: SIMCOON_NUMPY_API_OWNER is set on a single translation
// unit of each module, the only one where _import_array() exists (NO_IMPORT_ARRAY strips
// it everywhere else).
void import_api() {
    if (_import_array() < 0) {
        throw pybind11::error_already_set();
    }
}
#endif

} // namespace numpy_alloc
} // namespace simcoon

// libsimcoon's own table, when libsimcoon uses this allocator (Windows): exported from the
// DLL by numpy_alloc.cpp, called by _core at import. Returns 0 on success, -1 with the
// Python error set (a C-linkage frame must not carry a C++ exception across DLLs).
extern "C" int simcoon_numpy_alloc_import(void);

#define ARMA_ALIEN_MEM_ALLOC_FUNCTION simcoon::numpy_alloc::npy_malloc
#define ARMA_ALIEN_MEM_FREE_FUNCTION simcoon::numpy_alloc::npy_free
// Tells <carma> the allocator is set, so it does not include its own cnalloc.h.
#define CARMA_ARMA_ALIEN_MEM_FUNCTIONS_SET
