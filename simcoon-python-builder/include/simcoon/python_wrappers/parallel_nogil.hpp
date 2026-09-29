#pragma once

// simcoon_parallel_for_safe for the bindings, with the GIL released around the loop: the
// kernels allocate Armadillo memory through malloc/free in _core (and on Windows in
// libsimcoon). Releasing the GIL lets other Python threads run while the caller waits for
// workers and protects against future kernels that may acquire it. The loop body must not
// touch Python objects; exceptions are rethrown once the GIL is back.

#include <utility>
#include <pybind11/pybind11.h>
#include <simcoon/parallel.hpp>

namespace simpy {

template<typename F>
void parallel_for_nogil(int N, F&& func, unsigned int n_threads = 0) {
    pybind11::gil_scoped_release release;
    simcoon_parallel_for_safe(N, std::forward<F>(func), simcoon_parallel_cutoff, n_threads);
}

}  // namespace simpy
