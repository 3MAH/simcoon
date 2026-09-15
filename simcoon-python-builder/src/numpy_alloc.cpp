// Owner translation unit of libsimcoon's numpy C-API table (Windows, bindings on).
//
// Compiled into libsimcoon by simcoon_use_numpy_allocator with SIMCOON_NUMPY_API_OWNER set
// on this file alone: the force-included numpy_alloc.hpp then defines the table (and
// import_api) here and declares it extern in every other unit of the DLL.
#ifndef SIMCOON_NUMPY_API_OWNER
#error "numpy_alloc.cpp must be compiled with SIMCOON_NUMPY_API_OWNER (see CMakeLists.txt)"
#endif
#ifndef CARMA_ARMA_ALIEN_MEM_FUNCTIONS_SET
#error "numpy_alloc.hpp must be force-included into every translation unit of libsimcoon"
#endif

extern "C" int simcoon_numpy_alloc_import(void) {
    try {
        simcoon::numpy_alloc::import_api();
    } catch (pybind11::error_already_set &e) {
        e.restore();
        return -1;
    }
    return 0;
}
