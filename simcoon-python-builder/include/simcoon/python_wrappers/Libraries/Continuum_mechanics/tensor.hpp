#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace simpy {

// Function to register tensor2 and tensor4 classes with pybind11
void register_tensor(pybind11::module_& m);

} // namespace simpy
