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

///@file parallel.hpp
///@brief Cross-platform parallel-for: GCD on macOS, OpenMP on Linux, serial on Windows.

#pragma once

#if defined(__APPLE__)
  #include <dispatch/dispatch.h>

  template<typename F>
  void simcoon_parallel_for(int N, F&& func) {
      dispatch_apply(static_cast<size_t>(N), DISPATCH_APPLY_AUTO, ^(size_t i) {
          func(static_cast<int>(i));
      });
  }

#elif defined(_OPENMP)
  #include <omp.h>

  template<typename F>
  void simcoon_parallel_for(int N, F&& func) {
      #pragma omp parallel for
      for (int i = 0; i < N; i++) {
          func(i);
      }
  }

#else

  template<typename F>
  void simcoon_parallel_for(int N, F&& func) {
      for (int i = 0; i < N; i++) {
          func(i);
      }
  }

#endif

#include <exception>

/// @brief Exception-safe parallel loop over [0,N) for batch kernels that can throw.
///
/// OpenMP when available (parallel only past @p cutoff items; the pragma is inert otherwise),
/// serial fallback elsewhere. The first exception thrown by @p func is captured and rethrown
/// AFTER the loop: an exception escaping an active OpenMP parallel region is undefined
/// behavior (std::terminate), which would kill e.g. a Python session on the first singular
/// slice of a batch. Deliberately NOT routed through the GCD backend above: these loops are
/// reached from Python bindings, where GCD blocks interacting with the GIL have a deadlock
/// history (see the parallel-UMAT fix).
template<typename F>
void simcoon_parallel_for_safe(int N, F&& func, int cutoff = 100) {
    std::exception_ptr eptr = nullptr;
    #pragma omp parallel for schedule(static) if(N > cutoff)
    for (int i = 0; i < N; i++) {
        try {
            func(i);
        } catch (...) {
            #pragma omp critical (simcoon_parallel_for_safe_eptr)
            { if (!eptr) eptr = std::current_exception(); }
        }
    }
    if (eptr) std::rethrow_exception(eptr);
}
