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
///@brief Exception-safe parallel-for: GCD on macOS, OpenMP where available, std::thread otherwise.

#pragma once

#include <algorithm>
#include <atomic>
#include <exception>
#include <mutex>
#include <thread>
#include <vector>

#if defined(__APPLE__)
  #include <dispatch/dispatch.h>
#endif

/// Batches of at most this many items run serially: the parallel set-up would cost more.
constexpr int simcoon_parallel_cutoff = 100;

namespace simcoon_parallel_detail {

/// Run func over [begin, end), keeping the first exception: the remaining items of this range
/// are skipped, the other ranges still run, and the caller rethrows after the join.
template<typename F>
void run_range(int begin, int end, F& func, std::exception_ptr& eptr, std::mutex& mtx) {
    try {
        for (int i = begin; i < end; i++) func(i);
    } catch (...) {
        std::lock_guard<std::mutex> lk(mtx);
        if (!eptr) eptr = std::current_exception();
    }
}

/// Contiguous ranges of about 8 per hardware thread, at least 32 items (like OpenMP
/// schedule(static) but load-balanced): a per-item hand-out through a shared counter costs more
/// than ~1 microsecond kernels gain, one range per thread leaves the plastic hot spots of a mesh
/// to a single worker.
inline int chunk_size(int N) {
    const int hw = static_cast<int>(std::max(1u, std::thread::hardware_concurrency()));
    return std::max(N / (8 * hw), 32);
}

}  // namespace simcoon_parallel_detail

/// @brief Exception-safe parallel loop over [0,N) for batch kernels that can throw.
///
/// GCD on macOS, OpenMP where the build has it (Linux), std::thread otherwise (Windows, whose
/// Python stack would carry a second OpenMP runtime); serial up to @p cutoff items.
/// @p n_threads = 1 runs serially on every platform; otherwise it caps the std::thread workers
/// (0 = one per hardware thread), while GCD and OpenMP keep their own runtime sizing. The first
/// exception thrown by @p func is rethrown AFTER the loop: an exception escaping an active
/// OpenMP parallel region (or a GCD block, or a thread) would std::terminate, killing e.g. a
/// Python session on the first singular slice of a batch.
///
/// @warning CONTRACT: @p func must not touch Python memory (numpy allocation, refcounts,
/// anything needing the GIL). Worker threads acquiring the GIL while the calling thread
/// blocks on the loop is the lock cycle behind the 1.11.2 macOS parallel-UMAT deadlock
/// (carma copy -> PyDataMem_NEW -> GIL inside GCD). Do all numpy<->arma conversion BEFORE
/// the loop, and release the GIL around the C++ call in the Python bindings.
template<typename F>
void simcoon_parallel_for_safe(int N, F&& func, int cutoff = simcoon_parallel_cutoff,
                               unsigned int n_threads = 0) {
    if (N <= cutoff || n_threads == 1) {
        for (int i = 0; i < N; i++) func(i);
        return;
    }
    std::exception_ptr eptr = nullptr;
    std::mutex mtx;

#if defined(__APPLE__)
    const int chunk = simcoon_parallel_detail::chunk_size(N);
    const size_t nblocks = static_cast<size_t>((N + chunk - 1) / chunk);
    // dispatch_apply is synchronous, so pointers to these stack locals stay valid; the block
    // captures the pointers by value (no __block C++-object machinery needed).
    std::exception_ptr *pe = &eptr;
    std::mutex *pm = &mtx;
    auto *pf = &func;
    dispatch_apply(nblocks, DISPATCH_APPLY_AUTO, ^(size_t b) {
        const int start = static_cast<int>(b) * chunk;
        simcoon_parallel_detail::run_range(start, std::min(N, start + chunk), *pf, *pe, *pm);
    });

#elif defined(_OPENMP)
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < N; i++) {
        simcoon_parallel_detail::run_range(i, i + 1, func, eptr, mtx);
    }

#else
    const int chunk = simcoon_parallel_detail::chunk_size(N);
    const int nchunks = (N + chunk - 1) / chunk;
    const unsigned int hw = std::max(1u, std::thread::hardware_concurrency());
    const int workers = std::min(nchunks, static_cast<int>(n_threads ? std::min(n_threads, hw) : hw));
    std::atomic<int> next{0};
    auto worker = [&]() {
        for (int c = next++; c < nchunks; c = next++) {
            simcoon_parallel_detail::run_range(c * chunk, std::min(N, (c + 1) * chunk), func, eptr, mtx);
        }
    };
    std::vector<std::thread> threads;
    threads.reserve(workers - 1);
    try {
        for (int w = 1; w < workers; w++) threads.emplace_back(worker);
    } catch (...) {
        // thread creation failed: finish on the threads we have
    }
    worker();
    for (auto &t : threads) t.join();
#endif

    if (eptr) std::rethrow_exception(eptr);
}
