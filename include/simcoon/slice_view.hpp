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

///@file slice_view.hpp
///@brief Lock-free aliases of cube slices, for the parallel batch loops.

#pragma once

#include <stdexcept>
#include <armadillo>

namespace simcoon {

/// @brief Alias of slice @p s of @p c: the memory and fixed size of Cube::slice(s).
///
/// Cube::slice(s) builds its Mat accessor lazily, on the first access of each slice, behind a
/// per-cube std::mutex (an OpenMP critical section in OpenMP builds) and a heap allocation
/// (Cube_meat.hpp, get_mat_ptr). In a parallel loop over the slices of a fresh cube that is one
/// lock and one malloc per item and per cube, which serializes the loop. This alias is a plain
/// pointer: use it for every slice access inside parallel loops.
///
/// Initialize a named local from it (guaranteed copy elision keeps the alias); a copy of an
/// aliasing Mat owns fresh memory, and writes to it would not reach the cube. Bounds-checked
/// like Cube::slice: an out-of-range index throws std::out_of_range (the parallel loops rethrow
/// it after the join) instead of aliasing past the end of the buffer.
inline arma::mat slice_view(arma::cube &c, const arma::uword s) {
    if (s >= c.n_slices) {
        throw std::out_of_range("simcoon::slice_view(): slice index out of bounds");
    }
    return arma::mat(c.slice_memptr(s), c.n_rows, c.n_cols, false, true);
}

}  // namespace simcoon
