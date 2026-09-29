// Armadillo memory through the C runtime allocator, for every translation unit of a module.
//
// Force-included by CMake (simcoon_use_numpy_allocator) into every translation unit of
// _core and, on Windows, of libsimcoon: carma returns zero-copy numpy arrays whose
// capsules own Armadillo objects. Their destructors must use the same allocator as the
// code that created the buffers across DLL boundaries (carma issue #89).
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
#include <cstdlib>
#include <stdexcept>

namespace simcoon {
namespace numpy_alloc {

// Import numpy's C-API table into this module. The GIL must be held. Throws
// pybind11::error_already_set when numpy cannot be imported. Defined once per module, in
// its owner translation unit (below).
void import_api();

// The C runtime's malloc/free (Armadillo's own MSVC allocator uses _aligned_malloc,
// which free cannot release; the builds share the C runtime heap, /MD). Not PyDataMem_NEW/FREE: those
// add a tracemalloc notification that, since CPython 3.13.2 (gh-129185), takes the GIL on
// every call, traced or not -- every allocation of a parallel loop would queue on it. The
// cost is that simcoon's buffers are not counted by tracemalloc. Needs no numpy C-API table.
inline void *npy_malloc(std::size_t bytes) {
    return std::malloc(bytes);
}

inline void npy_free(void *ptr) {
    std::free(ptr);
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
