#pragma once

#include <pybind11/numpy.h>
#include <armadillo>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <type_traits>

namespace simpy::numpy_to_arma {

// Input conversion must never give Armadillo ownership of NumPy memory: NumPy
// may use a custom PyDataMem_Handler, whereas Armadillo uses malloc/free.
// Copies below are allocated by Armadillo, including the fallback for strided
// views. Call with the GIL held. For a borrowed view, the caller must keep src
// alive and must not resize/mutate it concurrently until the C++ work finishes.
namespace detail {

template<class A, typename T>
A convert(const pybind11::array_t<T>& src, bool copy) {
    const auto info = src.request();
    constexpr bool cube = std::is_same_v<A, arma::Cube<T>>;
    constexpr bool col = std::is_same_v<A, arma::Col<T>>;
    constexpr bool row = std::is_same_v<A, arma::Row<T>>;
    if constexpr (cube) {
        if (info.ndim != 3)
            throw std::invalid_argument("Expected a 3-D array for an Armadillo cube");
    } else {
        if (info.ndim < 1 || info.ndim > 2)
            throw std::invalid_argument("Expected a 1-D or 2-D array");
        if constexpr (col) {
            if (info.ndim == 2 && info.shape[1] != 1)
                throw std::invalid_argument("Expected a column vector");
        }
        if constexpr (row) {
            if (info.ndim == 2 && info.shape[0] != 1)
                throw std::invalid_argument("Expected a row vector");
        }
    }
    const auto nrows = static_cast<arma::uword>(info.shape[0]);
    const auto ncols = static_cast<arma::uword>(info.ndim > 1 ? info.shape[1] : 1);
    const auto nslices = static_cast<arma::uword>(cube ? info.shape[2] : 1);
    const bool contiguous = (src.flags() & pybind11::array::f_style) != 0;
    const bool aligned = reinterpret_cast<std::uintptr_t>(info.ptr) % alignof(T) == 0;
    auto* data = static_cast<T*>(info.ptr);
    // Do not expose a writable Armadillo view of a read-only NumPy buffer.
    if (!copy && info.size != 0 && contiguous && aligned && src.writeable()) {
        if constexpr (cube)
            return A(data, nrows, ncols, nslices, false, true);
        else if constexpr (col || row)
            return A(data, static_cast<arma::uword>(info.size), false, true);
        else
            return A(data, nrows, ncols, false, true);
    }
    A result;
    if constexpr (cube)
        result.set_size(nrows, ncols, nslices);
    else if constexpr (col || row)
        result.set_size(static_cast<arma::uword>(info.size));
    else
        result.set_size(nrows, ncols);
    if (info.size == 0) return result;
    if (contiguous) {
        std::memcpy(result.memptr(), info.ptr, static_cast<std::size_t>(info.size) * sizeof(T));
    } else {
        const auto* base = static_cast<const char*>(info.ptr);
        T* dest = result.memptr();
        const auto stride1 = info.ndim > 1 ? info.strides[1] : 0;
        const auto stride2 = cube ? info.strides[2] : 0;
        for (arma::uword k = 0; k < nslices; ++k)
            for (arma::uword j = 0; j < ncols; ++j)
                for (arma::uword i = 0; i < nrows; ++i) {
                    // Signed offsets support reversed slices; memcpy supports
                    // unaligned inputs without dereferencing an unaligned T*.
                    const auto offset = static_cast<pybind11::ssize_t>(i) * info.strides[0]
                                      + static_cast<pybind11::ssize_t>(j) * stride1
                                      + static_cast<pybind11::ssize_t>(k) * stride2;
                    std::memcpy(dest++, base + offset, sizeof(T));
                }
    }
    return result;
}

} // namespace detail

// Match the call sites' existing copy/view distinction, without any ownership
// stealing or modification of the NumPy array's shape, strides or flags.
template<typename T>
arma::Mat<T> arr_to_mat(const pybind11::array_t<T>& src) {
    return detail::convert<arma::Mat<T>>(src, true);
}
template<typename T>
arma::Mat<T> arr_to_mat(pybind11::array_t<T>& src, bool copy = false) {
    return detail::convert<arma::Mat<T>>(src, copy);
}
template<typename T>
arma::Mat<T> arr_to_mat(pybind11::array_t<T>&& src) {
    return detail::convert<arma::Mat<T>>(src, true);
}
template<typename T>
arma::Mat<T> arr_to_mat_view(const pybind11::array_t<T>& src) {
    return detail::convert<arma::Mat<T>>(src, false);
}

template<typename T>
arma::Col<T> arr_to_col(const pybind11::array_t<T>& src) {
    return detail::convert<arma::Col<T>>(src, true);
}
template<typename T>
arma::Col<T> arr_to_col(pybind11::array_t<T>& src, bool copy = false) {
    return detail::convert<arma::Col<T>>(src, copy);
}
template<typename T>
arma::Col<T> arr_to_col(pybind11::array_t<T>&& src) {
    return detail::convert<arma::Col<T>>(src, true);
}
template<typename T>
arma::Col<T> arr_to_col_view(const pybind11::array_t<T>& src) {
    return detail::convert<arma::Col<T>>(src, false);
}

template<typename T>
arma::Row<T> arr_to_row(const pybind11::array_t<T>& src) {
    return detail::convert<arma::Row<T>>(src, true);
}
template<typename T>
arma::Row<T> arr_to_row(pybind11::array_t<T>& src, bool copy = false) {
    return detail::convert<arma::Row<T>>(src, copy);
}
template<typename T>
arma::Row<T> arr_to_row(pybind11::array_t<T>&& src) {
    return detail::convert<arma::Row<T>>(src, true);
}
template<typename T>
arma::Row<T> arr_to_row_view(const pybind11::array_t<T>& src) {
    return detail::convert<arma::Row<T>>(src, false);
}

template<typename T>
arma::Cube<T> arr_to_cube(const pybind11::array_t<T>& src) {
    return detail::convert<arma::Cube<T>>(src, true);
}
template<typename T>
arma::Cube<T> arr_to_cube(pybind11::array_t<T>& src, bool copy = false) {
    return detail::convert<arma::Cube<T>>(src, copy);
}
template<typename T>
arma::Cube<T> arr_to_cube(pybind11::array_t<T>&& src) {
    return detail::convert<arma::Cube<T>>(src, true);
}
template<typename T>
arma::Cube<T> arr_to_cube_view(const pybind11::array_t<T>& src) {
    return detail::convert<arma::Cube<T>>(src, false);
}

} // namespace simpy::numpy_to_arma
