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

///@file Tparallel.cpp
///@brief Test of simcoon_parallel_for_safe: coverage, n_threads cap, exceptions
///@version 1.0

#include <gtest/gtest.h>
#include <atomic>
#include <chrono>
#include <mutex>
#include <set>
#include <stdexcept>
#include <thread>
#include <vector>
#include <simcoon/parallel.hpp>

// Every item runs once, and no more than n_threads threads ever run the loop body.
TEST(Tparallel, n_threads_caps_the_threads)
{
    const int N = 4000;
    for (unsigned int n_threads : {1u, 2u, 3u, 0u}) {
        std::vector<int> count(N, 0);
        std::set<std::thread::id> ids;
        std::mutex mtx;
        std::atomic<int> running{0};
        std::atomic<int> peak{0};
        simcoon_parallel_for_safe(N, [&](int i) {
            const int now = ++running;
            int p = peak.load();
            while (now > p && !peak.compare_exchange_weak(p, now)) {}
            count[i]++;
            {
                std::lock_guard<std::mutex> lk(mtx);
                ids.insert(std::this_thread::get_id());
            }
            std::this_thread::sleep_for(std::chrono::microseconds(20));
            --running;
        }, simcoon_parallel_cutoff, n_threads);
        for (int i = 0; i < N; i++) ASSERT_EQ(count[i], 1) << "item " << i << ", n_threads " << n_threads;
        if (n_threads > 0) {
            EXPECT_LE(ids.size(), n_threads) << "n_threads " << n_threads;
            EXPECT_LE(peak.load(), static_cast<int>(n_threads)) << "n_threads " << n_threads;
        }
    }
}

// The first exception of the loop body is rethrown once the loop is over, whatever n_threads.
TEST(Tparallel, exception_is_rethrown_after_the_loop)
{
    for (unsigned int n_threads : {1u, 2u, 0u}) {
        EXPECT_THROW(simcoon_parallel_for_safe(1000, [](int i) {
            if (i == 517) throw std::runtime_error("item 517");
        }, simcoon_parallel_cutoff, n_threads), std::runtime_error);
    }
}
