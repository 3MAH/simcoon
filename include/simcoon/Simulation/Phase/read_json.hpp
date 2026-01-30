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

///@file read_json.hpp
///@brief JSON-based I/O for phase configurations (ellipsoids, layers, cylinders)
///@version 2.0
///
///@note JSON files can be created/edited using the Python interface:
///      @code{.py}
///      from simcoon.solver.micromechanics import Ellipsoid, save_ellipsoids_json
///      @endcode

#pragma once
#include <armadillo>
#include <string>
#include <simcoon/Simulation/Phase/phase_characteristics.hpp>

namespace simcoon{

/**
 * @file read_json.hpp
 * @brief JSON-based phase configuration I/O.
 *
 * This module provides functions to read phase configurations from JSON files.
 * The JSON format is:
 * - Self-documenting with named properties
 * - Easy to edit programmatically
 * - Compatible with the Python `simcoon.solver.micromechanics` module
 */

/** @addtogroup phase
 *  @{
 */

/**
 * @brief Read ellipsoid phase characteristics from a JSON file.
 *
 * The JSON file should have the following structure:
 * @code{.json}
 * {
 *   "ellipsoids": [
 *     {
 *       "number": 0,
 *       "coatingof": 0,
 *       "umat_name": "ELISO",
 *       "save": 1,
 *       "concentration": 0.7,
 *       "material_orientation": {"psi": 0, "theta": 0, "phi": 0},
 *       "semi_axes": {"a1": 1, "a2": 1, "a3": 1},
 *       "geometry_orientation": {"psi": 0, "theta": 0, "phi": 0},
 *       "nstatev": 1,
 *       "props": {"E": 70000, "nu": 0.3, "alpha": 1e-5}
 *     }
 *   ]
 * }
 * @endcode
 *
 * @param rve Reference to phase_characteristics to populate
 * @param path_data Directory containing the JSON file (default: "data")
 * @param inputfile JSON filename (default: "ellipsoids.json")
 */
void read_ellipsoid_json(phase_characteristics &rve,
                         const std::string &path_data = "data",
                         const std::string &inputfile = "ellipsoids.json");

/**
 * @brief Read layer phase characteristics from a JSON file.
 *
 * The JSON file should have the following structure:
 * @code{.json}
 * {
 *   "layers": [
 *     {
 *       "number": 0,
 *       "umat_name": "ELISO",
 *       "save": 1,
 *       "concentration": 0.5,
 *       "material_orientation": {"psi": 0, "theta": 0, "phi": 0},
 *       "geometry_orientation": {"psi": 0, "theta": 90, "phi": -90},
 *       "nstatev": 1,
 *       "props": {"E": 70000, "nu": 0.3, "alpha": 1e-5}
 *     }
 *   ]
 * }
 * @endcode
 *
 * @param rve Reference to phase_characteristics to populate
 * @param path_data Directory containing the JSON file (default: "data")
 * @param inputfile JSON filename (default: "layers.json")
 */
void read_layer_json(phase_characteristics &rve,
                     const std::string &path_data = "data",
                     const std::string &inputfile = "layers.json");

/**
 * @brief Read cylinder phase characteristics from a JSON file.
 *
 * @param rve Reference to phase_characteristics to populate
 * @param path_data Directory containing the JSON file (default: "data")
 * @param inputfile JSON filename (default: "cylinders.json")
 */
void read_cylinder_json(phase_characteristics &rve,
                        const std::string &path_data = "data",
                        const std::string &inputfile = "cylinders.json");

/**
 * @brief Read generic phase characteristics from a JSON file.
 *
 * @param rve Reference to phase_characteristics to populate
 * @param path_data Directory containing the JSON file (default: "data")
 * @param inputfile JSON filename (default: "phases.json")
 */
void read_phase_json(phase_characteristics &rve,
                     const std::string &path_data = "data",
                     const std::string &inputfile = "phases.json");

/**
 * @brief Check if a JSON file exists for the given configuration.
 *
 * @param path_data Directory to check
 * @param inputfile JSON filename to check
 * @return true if the JSON file exists, false otherwise
 */
bool json_file_exists(const std::string &path_data, const std::string &inputfile);

/**
 * @brief Write phase characteristics to a JSON file.
 *
 * @param rve Reference to phase_characteristics to write
 * @param path_data Directory to write the JSON file (default: "data")
 * @param outputfile JSON filename (default: "phases.json")
 */
void write_phase_json(phase_characteristics &rve,
                      const std::string &path_data = "data",
                      const std::string &outputfile = "phases.json");

/**
 * @brief Write ellipsoid phase characteristics to a JSON file.
 *
 * @param rve Reference to phase_characteristics to write
 * @param path_data Directory to write the JSON file (default: "data")
 * @param outputfile JSON filename (default: "ellipsoids.json")
 */
void write_ellipsoid_json(phase_characteristics &rve,
                          const std::string &path_data = "data",
                          const std::string &outputfile = "ellipsoids.json");

/**
 * @brief Write layer phase characteristics to a JSON file.
 *
 * @param rve Reference to phase_characteristics to write
 * @param path_data Directory to write the JSON file (default: "data")
 * @param outputfile JSON filename (default: "layers.json")
 */
void write_layer_json(phase_characteristics &rve,
                      const std::string &path_data = "data",
                      const std::string &outputfile = "layers.json");

/**
 * @brief Write cylinder phase characteristics to a JSON file.
 *
 * @param rve Reference to phase_characteristics to write
 * @param path_data Directory to write the JSON file (default: "data")
 * @param outputfile JSON filename (default: "cylinders.json")
 */
void write_cylinder_json(phase_characteristics &rve,
                         const std::string &path_data = "data",
                         const std::string &outputfile = "cylinders.json");

/** @} */ // end of phase group

} //namespace simcoon
