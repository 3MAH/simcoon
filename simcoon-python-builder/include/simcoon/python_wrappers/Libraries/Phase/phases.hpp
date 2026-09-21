#pragma once
#include <pybind11/pybind11.h>

#include <string>
#include <vector>

#include <simcoon/Simulation/Phase/phase_characteristics.hpp>

namespace simpy{

//Build the sub-phases of a mean-field model from their Python description: a sequence of
//dicts, one per phase, in the shape simcoon.solver.micromechanics produces. This replaces
//the Nellipsoids<N>.dat / Nlayers<N>.dat files the C++ side used to read from a "data"
//directory, and is shared by L_eff and by the solver.
//
//Angles arrive in degrees, as they did in those files, and are stored in radians.
//
//The phases must be numbered by their position in the list (0, 1, ... n-1): the schemes
//select the matrix by `number` and then index `sub_phases[n_matrix]`, so the two agree only
//when the number IS the position. Phases sharing a number would make every concentration
//tensor the identity, and L_eff would quietly return the Voigt average.
//
//Recursive: a sub-phase that is itself a mean-field model carries its own sub-phases under
//
//Each dict may carry a `kind` ('ellipsoid', 'layer', 'cylinder', 'phase', what
//to_phase_dict writes): it must be the geometry the model builds. `concentration` is
//required, and the concentrations of one level must sum to 1.
//
//Returns an empty vector for a homogeneous model. Throws std::invalid_argument when the
//phases are missing, superfluous, misnumbered or of the wrong kind. Their consistency with
//the props of the model is the library's check (check_sub_phases, multiphase.hpp).
std::vector<simcoon::phase_characteristics> make_sub_phases(const pybind11::object &phases,
                                                            const std::string &umat_name,
                                                            const double &T_init, const int &depth = 0);


} //namespace simpy
