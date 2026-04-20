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
