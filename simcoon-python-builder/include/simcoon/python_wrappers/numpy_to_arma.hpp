// NumPy -> Armadillo input conversions for the python wrappers.
//
// Two rules, the reason this header exists instead of carma's arr_to_*:
//   1. Armadillo never owns NumPy memory. Armadillo allocates and frees with malloc/free
//      (numpy_alloc.hpp); NumPy frees with the PyDataMem_Handler of the array, which may
//      not be malloc/free. carma's input converters always end with Armadillo owning a
//      NumPy-allocated buffer (stolen, or a copy made with PyArray_NewLikeArray).
//   2. The NumPy array is never modified: no shape, stride or flag change, no ownership
//      transfer.
//
// Two entry points per shape, decided at the call site, not by the constness of the
// argument expression:
//   arr_to_{mat,col,cube}(src)       always an Armadillo-owned copy (any layout, order,
//                                    alignment or writeability).
//   arr_to_{mat,col,cube}_view(src)  a strict alias of the NumPy buffer when it is
//                                    F-contiguous, aligned and writeable; a copy otherwise.
//                                    The caller keeps src alive and unmodified until the
//                                    C++ work is done. A strict alias cannot be resized:
//                                    initialise the Armadillo object from it once and
//                                    never re-assign it with another size.
// A 0-d array is accepted as a single element (mat 1x1, col of one).
// Call with the GIL held.
#pragma once

#include <pybind11/numpy.h>
#include <armadillo>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <type_traits>

namespace simpy::numpy_to_arma {

namespace detail {

template<class A, typename T, int Flags>
A convert(const pybind11::array_t<T, Flags>& src, bool view) {
    const auto info = src.request();
    constexpr bool cube = std::is_same_v<A, arma::Cube<T>>;
    constexpr bool col = std::is_same_v<A, arma::Col<T>>;
    if constexpr (cube) {
        if (info.ndim != 3)
            throw std::invalid_argument("Expected a 3-D array for an Armadillo cube");
    } else {
        if (info.ndim > 2)
            throw std::invalid_argument("Expected a 0-D, 1-D or 2-D array");
        if constexpr (col) {
            if (info.ndim == 2 && info.shape[1] != 1)
                throw std::invalid_argument("Expected a column vector");
        }
    }
    const auto nrows = static_cast<arma::uword>(info.ndim > 0 ? info.shape[0] : 1);
    const auto ncols = static_cast<arma::uword>(info.ndim > 1 ? info.shape[1] : 1);
    const auto nslices = static_cast<arma::uword>(cube ? info.shape[2] : 1);
    const auto n_elem = static_cast<arma::uword>(info.size);
    const bool contiguous = info.ndim == 0 || (src.flags() & pybind11::array::f_style) != 0;
    const bool aligned = reinterpret_cast<std::uintptr_t>(info.ptr) % alignof(T) == 0;
    auto* data = static_cast<T*>(info.ptr);
    // No writable Armadillo alias of a read-only NumPy buffer.
    if (view && n_elem != 0 && contiguous && aligned && src.writeable()) {
        if constexpr (cube)
            return A(data, nrows, ncols, nslices, false, true);
        else if constexpr (col)
            return A(data, n_elem, false, true);
        else
            return A(data, nrows, ncols, false, true);
    }
    A result;
    if constexpr (cube)
        result.set_size(nrows, ncols, nslices);
    else if constexpr (col)
        result.set_size(n_elem);
    else
        result.set_size(nrows, ncols);
    if (n_elem == 0) return result;
    if (contiguous) {
        std::memcpy(result.memptr(), info.ptr, static_cast<std::size_t>(n_elem) * sizeof(T));
    } else {
        // Element gather in Fortran order. Signed offsets support reversed slices; memcpy
        // supports unaligned inputs without dereferencing an unaligned T*.
        const auto* base = static_cast<const char*>(info.ptr);
        T* dest = result.memptr();
        const auto stride1 = info.ndim > 1 ? info.strides[1] : 0;
        const auto stride2 = cube ? info.strides[2] : 0;
        for (arma::uword k = 0; k < nslices; ++k)
            for (arma::uword j = 0; j < ncols; ++j)
                for (arma::uword i = 0; i < nrows; ++i) {
                    const auto offset = static_cast<pybind11::ssize_t>(i) * info.strides[0]
                                      + static_cast<pybind11::ssize_t>(j) * stride1
                                      + static_cast<pybind11::ssize_t>(k) * stride2;
                    std::memcpy(dest++, base + offset, sizeof(T));
                }
    }
    return result;
}

} // namespace detail

template<typename T, int Flags>
arma::Mat<T> arr_to_mat(const pybind11::array_t<T, Flags>& src) {
    return detail::convert<arma::Mat<T>>(src, false);
}
template<typename T, int Flags>
arma::Mat<T> arr_to_mat_view(const pybind11::array_t<T, Flags>& src) {
    return detail::convert<arma::Mat<T>>(src, true);
}

template<typename T, int Flags>
arma::Col<T> arr_to_col(const pybind11::array_t<T, Flags>& src) {
    return detail::convert<arma::Col<T>>(src, false);
}
template<typename T, int Flags>
arma::Col<T> arr_to_col_view(const pybind11::array_t<T, Flags>& src) {
    return detail::convert<arma::Col<T>>(src, true);
}

template<typename T, int Flags>
arma::Cube<T> arr_to_cube(const pybind11::array_t<T, Flags>& src) {
    return detail::convert<arma::Cube<T>>(src, false);
}
template<typename T, int Flags>
arma::Cube<T> arr_to_cube_view(const pybind11::array_t<T, Flags>& src) {
    return detail::convert<arma::Cube<T>>(src, true);
}

} // namespace simpy::numpy_to_arma
