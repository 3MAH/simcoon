#pragma once
#include <pybind11/pybind11.h>

#include <string>
#include <vector>

#include <simcoon/Simulation/Phase/phase_characteristics.hpp>

namespace simpy{

//Geometry of the sub-phases a mean-field model builds, with the codes of
//phase_characteristics::construct: 1 = layer, 2 = ellipsoid. 0 means the model is
//homogeneous and has no sub-phases at all.
int shape_type_of(const std::string &umat_name);

//Build the sub-phases of a mean-field model from their Python description: a sequence of
//dicts, one per phase, in the shape simcoon.solver.micromechanics produces. This replaces
//the Nellipsoids<N>.dat / Nlayers<N>.dat files the C++ side used to read from a "data"
//directory, and is shared by L_eff and by the solver.
//
//Angles arrive in degrees, as they did in those files, and are stored in radians.
//
//Returns an empty vector for a homogeneous model. Throws std::invalid_argument when the
//phases are missing, superfluous, or do not match the count props[0] announces; pass
//announced_nphases < 0 to skip that last check.
std::vector<simcoon::phase_characteristics> make_sub_phases(const pybind11::object &phases,
                                                            const std::string &umat_name,
                                                            const int &announced_nphases,
                                                            const double &T_init);

} //namespace simpy
