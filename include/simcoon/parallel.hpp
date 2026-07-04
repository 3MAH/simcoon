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
#include <mutex>
#include <thread>
#include <algorithm>

/// @brief Exception-safe parallel loop over [0,N) for batch kernels that can throw.
///
/// GCD on macOS, OpenMP on Linux (both parallel only past @p cutoff items), serial fallback
/// elsewhere. The first exception thrown by @p func is captured and rethrown AFTER the loop:
/// an exception escaping an active OpenMP parallel region (or a GCD block) is undefined
/// behavior (std::terminate), which would kill e.g. a Python session on the first singular
/// slice of a batch.
///
/// @warning CONTRACT: @p func must not touch Python memory (numpy allocation, refcounts,
/// anything needing the GIL). Worker threads acquiring the GIL while the calling thread
/// blocks on the loop is the lock cycle behind the 1.11.2 macOS parallel-UMAT deadlock
/// (carma copy -> PyDataMem_NEW -> GIL inside GCD). Do all numpy<->arma conversion BEFORE
/// the loop, and release the GIL around the C++ call in the Python bindings.
#if defined(__APPLE__)

template<typename F>
void simcoon_parallel_for_safe(int N, F&& func, int cutoff = 100) {
    if (N <= cutoff) {
        for (int i = 0; i < N; i++) func(i);
        return;
    }
    // Manual chunking (contiguous ranges, like OpenMP schedule(static)): dispatch_apply
    // hands out its iterations through a contended atomic counter, which at ~1 microsecond
    // per item costs more than the parallelism gains. ~8 chunks per core balances load
    // without paying that per-item toll.
    const int hw = std::max(1u, std::thread::hardware_concurrency());
    const int chunk = std::max(N / (8 * hw), 32);
    const size_t nblocks = static_cast<size_t>((N + chunk - 1) / chunk);
    std::exception_ptr eptr = nullptr;
    std::mutex mtx;
    // dispatch_apply is synchronous, so pointers to these stack locals stay valid; the
    // block captures the pointers by value (no __block C++-object machinery needed).
    std::exception_ptr *pe = &eptr;
    std::mutex *pm = &mtx;
    dispatch_apply(nblocks, DISPATCH_APPLY_AUTO, ^(size_t b) {
        const int start = static_cast<int>(b) * chunk;
        const int end = std::min(N, start + chunk);
        try {
            for (int i = start; i < end; i++) func(i);
        } catch (...) {
            // Remaining items of THIS chunk are skipped; other chunks still run
            // (same semantics as the OpenMP branch: no cancellation, first
            // exception rethrown after the join).
            std::lock_guard<std::mutex> lk(*pm);
            if (!*pe) *pe = std::current_exception();
        }
    });
    if (eptr) std::rethrow_exception(eptr);
}

#else

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

#endif
