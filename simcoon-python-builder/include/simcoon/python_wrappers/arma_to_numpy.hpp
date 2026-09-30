// Armadillo -> NumPy output conversions for the python wrappers.
//
// Zero-copy by default: the Armadillo object is moved to the heap and owned by a capsule
// set as the array's base. NumPy never frees the buffer (OWNDATA is false); the capsule
// deletes the object, so allocation and release both go through Armadillo's own
// allocator, identical in libsimcoon and _core. Do not override it (ARMA_ALIEN_MEM_*) on
// one side only: a buffer allocated by one allocator would be freed by the other.
// Layout: Fortran order, writeable; Mat (r,c), Col (n,1), Row (1,n), Cube (r,c,s).
// Call with the GIL held.
#pragma once

#include <pybind11/numpy.h>
#include <armadillo>
#include <memory>
#include <type_traits>
#include <utility>

#ifdef ARMA_ALIEN_MEM_ALLOC_FUNCTION
#error "arma_to_numpy.hpp: libsimcoon and _core must share Armadillo's default allocator"
#endif

namespace simpy::arma_to_numpy {

namespace detail {

namespace py = pybind11;

// Heap object owning its memory. Moving an object on auxiliary memory (mem_state 1/2,
// e.g. a numpy_to_arma view) would keep aliasing that memory, so it is copied instead;
// fixed-size objects (mem_state 3) are copied by Armadillo's move anyway.
template<class A>
std::unique_ptr<A> take(A& src, bool copy) {
    if (copy || src.mem_state != 0)
        return std::make_unique<A>(src);
    return std::make_unique<A>(std::move(src));
}

template<class A>
py::array_t<typename A::elem_type> wrap(std::unique_ptr<A> owned) {
    using T = typename A::elem_type;
    constexpr auto s = static_cast<py::ssize_t>(sizeof(T));
    const auto r = static_cast<py::ssize_t>(owned->n_rows);
    const auto c = static_cast<py::ssize_t>(owned->n_cols);
    py::ssize_t sl = 1;
    if constexpr (std::is_same_v<A, arma::Cube<T>>)
        sl = static_cast<py::ssize_t>(owned->n_slices);
    T* data = owned->memptr();
    // Ownership moves to the capsule only once it exists: no leak if it throws.
    py::capsule base(owned.get(), [](void* p) { delete static_cast<A*>(p); });
    owned.release();
    if constexpr (std::is_same_v<A, arma::Cube<T>>) {
        return py::array_t<T>({r, c, sl}, {s, r * s, r * c * s}, data, base);
    } else if constexpr (std::is_same_v<A, arma::Row<T>>) {
        return py::array_t<T>({r, c}, {s, s}, data, base);
    } else {
        return py::array_t<T>({r, c}, {s, r * s}, data, base);
    }
}

} // namespace detail

// A const lvalue or `copy = true` copies, an rvalue or a non-const lvalue is moved (left empty).
template<typename T>
pybind11::array_t<T> mat_to_arr(const arma::Mat<T>& src) {
    return detail::wrap(std::make_unique<arma::Mat<T>>(src));
}
template<typename T>
pybind11::array_t<T> mat_to_arr(arma::Mat<T>&& src) {
    return detail::wrap(detail::take(src, false));
}
template<typename T>
pybind11::array_t<T> mat_to_arr(arma::Mat<T>& src, bool copy = false) {
    return detail::wrap(detail::take(src, copy));
}

template<typename T>
pybind11::array_t<T> col_to_arr(const arma::Col<T>& src) {
    return detail::wrap(std::make_unique<arma::Col<T>>(src));
}
template<typename T>
pybind11::array_t<T> col_to_arr(arma::Col<T>&& src) {
    return detail::wrap(detail::take(src, false));
}
template<typename T>
pybind11::array_t<T> col_to_arr(arma::Col<T>& src, bool copy = false) {
    return detail::wrap(detail::take(src, copy));
}

template<typename T>
pybind11::array_t<T> row_to_arr(const arma::Row<T>& src) {
    return detail::wrap(std::make_unique<arma::Row<T>>(src));
}
template<typename T>
pybind11::array_t<T> row_to_arr(arma::Row<T>&& src) {
    return detail::wrap(detail::take(src, false));
}
template<typename T>
pybind11::array_t<T> row_to_arr(arma::Row<T>& src, bool copy = false) {
    return detail::wrap(detail::take(src, copy));
}

template<typename T>
pybind11::array_t<T> cube_to_arr(const arma::Cube<T>& src) {
    return detail::wrap(std::make_unique<arma::Cube<T>>(src));
}
template<typename T>
pybind11::array_t<T> cube_to_arr(arma::Cube<T>&& src) {
    return detail::wrap(detail::take(src, false));
}
template<typename T>
pybind11::array_t<T> cube_to_arr(arma::Cube<T>& src, bool copy = false) {
    return detail::wrap(detail::take(src, copy));
}

} // namespace simpy::arma_to_numpy
