#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace simpy {

// Function to register the Rotation class with pybind11
void register_rotation(pybind11::module_& m);

} // namespace simpy
