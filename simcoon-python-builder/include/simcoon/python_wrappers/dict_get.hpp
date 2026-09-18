#pragma once
#include <pybind11/pybind11.h>

namespace simpy {

//A dict entry with a default: absent or None both mean "not given".
template <typename T>
T dget(const pybind11::dict &d, const char *key, const T &dflt) {
    return (d.contains(key) && !d[key].is_none()) ? d[key].cast<T>() : dflt;
}

} //namespace simpy
