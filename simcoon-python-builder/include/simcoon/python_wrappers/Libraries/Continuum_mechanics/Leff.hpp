#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace simpy{

//Return the elastic stiffness tensor of a composite material.
//`phases` describes the sub-phases in memory: a sequence of dicts, one per phase, as
//produced by simcoon.solver.micromechanics (ellipsoids, layers, ...). It replaces the
//Nellipsoids<N>.dat / Nlayers<N>.dat files the C++ side used to read from a "data"
//directory; pass None for the homogeneous models (ELISO, ELIST, ELORT), which have none.
//`phases` comes last on purpose: it is a new argument, and the RVE angles keep the positions
//callers already pass them in.
pybind11::array_t<double> L_eff(const std::string &umat_name, const pybind11::array_t<double> &props, const int &nstatev, const pybind11::object &orientation=pybind11::none(), const pybind11::object &phases=pybind11::none());
    
} //namespace simpy
