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

///@file read_json.cpp
///@brief JSON-based I/O for phase configurations
///@version 1.0

// Define _USE_MATH_DEFINES before cmath for M_PI on Windows MSVC
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#include <cmath>

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <filesystem>
#include <nlohmann/json.hpp>
#include <simcoon/Simulation/Phase/read_json.hpp>
#include <simcoon/Simulation/Geometry/ellipsoid.hpp>
#include <simcoon/Simulation/Geometry/layer.hpp>
#include <simcoon/Simulation/Geometry/cylinder.hpp>

using namespace arma;
using json = nlohmann::json;

namespace simcoon {

// Helper to convert degrees to radians
constexpr double deg2rad = M_PI / 180.0;

// Helper to get props array from JSON (handles both dict and array formats)
static vec get_props_from_json(const json& j) {
    if (j.contains("props")) {
        const auto& props_json = j["props"];
        if (props_json.is_object()) {
            // Props as dict: {"E": 70000, "nu": 0.3, ...}
            std::vector<double> props_vec;
            for (auto& [key, val] : props_json.items()) {
                props_vec.push_back(val.get<double>());
            }
            return vec(props_vec);
        } else if (props_json.is_array()) {
            // Props as array: [70000, 0.3, ...]
            std::vector<double> props_vec = props_json.get<std::vector<double>>();
            return vec(props_vec);
        }
    }
    return vec();
}

// Helper to get orientation from JSON
static void get_orientation(const json& j, const std::string& key,
                           double& psi, double& theta, double& phi,
                           bool to_radians = true) {
    if (j.contains(key)) {
        const auto& orient = j[key];
        psi = orient.value("psi", 0.0);
        theta = orient.value("theta", 0.0);
        phi = orient.value("phi", 0.0);
        if (to_radians) {
            psi *= deg2rad;
            theta *= deg2rad;
            phi *= deg2rad;
        }
    }
}

bool json_file_exists(const std::string &path_data, const std::string &inputfile) {
    std::filesystem::path filepath = std::filesystem::path(path_data) / inputfile;
    return std::filesystem::exists(filepath);
}

void read_ellipsoid_json(phase_characteristics &rve,
                         const std::string &path_data,
                         const std::string &inputfile) {

    std::filesystem::path filepath = std::filesystem::path(path_data) / inputfile;

    if (!std::filesystem::exists(filepath)) {
        throw std::runtime_error("JSON file not found: " + filepath.string());
    }

    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open JSON file: " + filepath.string());
    }

    json j;
    try {
        file >> j;
    } catch (const json::parse_error& e) {
        throw std::runtime_error("JSON parse error in " + filepath.string() + ": " + e.what());
    }

    if (!j.contains("ellipsoids") || !j["ellipsoids"].is_array()) {
        throw std::runtime_error("JSON file must contain 'ellipsoids' array");
    }

    const auto& ellipsoids = j["ellipsoids"];
    unsigned int nphases = ellipsoids.size();

    // Construct sub-phases with ellipsoid geometry (type 2)
    rve.sub_phases_construct(nphases, 2, 1);

    for (unsigned int i = 0; i < nphases; i++) {
        const auto& ell_json = ellipsoids[i];
        auto& phase = rve.sub_phases[i];

        // Get ellipsoid geometry pointer
        auto sptr_ellipsoid = std::dynamic_pointer_cast<ellipsoid>(phase.sptr_shape);
        if (!sptr_ellipsoid) {
            throw std::runtime_error("Failed to cast shape to ellipsoid");
        }

        // Read material properties
        phase.sptr_matprops->number = ell_json.value("number", static_cast<int>(i));
        phase.sptr_matprops->umat_name = ell_json.value("umat_name", std::string("ELISO"));
        phase.sptr_matprops->save = ell_json.value("save", 1);

        // Read concentration
        sptr_ellipsoid->concentration = ell_json.value("concentration", 1.0);

        // Read coating info
        sptr_ellipsoid->coatingof = ell_json.value("coatingof", 0);
        sptr_ellipsoid->coatedby = ell_json.value("coatedby", 0);

        // Read semi-axes (can be in "semi_axes" dict or directly)
        if (ell_json.contains("semi_axes")) {
            const auto& axes = ell_json["semi_axes"];
            sptr_ellipsoid->a1 = axes.value("a1", 1.0);
            sptr_ellipsoid->a2 = axes.value("a2", 1.0);
            sptr_ellipsoid->a3 = axes.value("a3", 1.0);
        } else {
            sptr_ellipsoid->a1 = ell_json.value("a1", 1.0);
            sptr_ellipsoid->a2 = ell_json.value("a2", 1.0);
            sptr_ellipsoid->a3 = ell_json.value("a3", 1.0);
        }

        // Read material orientation (degrees in JSON, stored as radians)
        get_orientation(ell_json, "material_orientation",
                       phase.sptr_matprops->psi_mat,
                       phase.sptr_matprops->theta_mat,
                       phase.sptr_matprops->phi_mat, true);

        // Read geometry orientation (degrees in JSON, stored as radians)
        get_orientation(ell_json, "geometry_orientation",
                       sptr_ellipsoid->psi_geom,
                       sptr_ellipsoid->theta_geom,
                       sptr_ellipsoid->phi_geom, true);

        // Read props array
        phase.sptr_matprops->props = get_props_from_json(ell_json);
        phase.sptr_matprops->nprops = phase.sptr_matprops->props.n_elem;
    }
}

void read_layer_json(phase_characteristics &rve,
                     const std::string &path_data,
                     const std::string &inputfile) {

    std::filesystem::path filepath = std::filesystem::path(path_data) / inputfile;

    if (!std::filesystem::exists(filepath)) {
        throw std::runtime_error("JSON file not found: " + filepath.string());
    }

    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open JSON file: " + filepath.string());
    }

    json j;
    try {
        file >> j;
    } catch (const json::parse_error& e) {
        throw std::runtime_error("JSON parse error in " + filepath.string() + ": " + e.what());
    }

    if (!j.contains("layers") || !j["layers"].is_array()) {
        throw std::runtime_error("JSON file must contain 'layers' array");
    }

    const auto& layers = j["layers"];
    unsigned int nphases = layers.size();

    // Construct sub-phases with layer geometry (type 1)
    rve.sub_phases_construct(nphases, 1, 1);

    for (unsigned int i = 0; i < nphases; i++) {
        const auto& layer_json = layers[i];
        auto& phase = rve.sub_phases[i];

        // Get layer geometry pointer
        auto sptr_layer = std::dynamic_pointer_cast<layer>(phase.sptr_shape);
        if (!sptr_layer) {
            throw std::runtime_error("Failed to cast shape to layer");
        }

        // Read material properties
        phase.sptr_matprops->number = layer_json.value("number", static_cast<int>(i));
        phase.sptr_matprops->umat_name = layer_json.value("umat_name", std::string("ELISO"));
        phase.sptr_matprops->save = layer_json.value("save", 1);

        // Read concentration
        sptr_layer->concentration = layer_json.value("concentration", 1.0);

        // Read layer connectivity
        sptr_layer->layerup = layer_json.value("layerup", -1);
        sptr_layer->layerdown = layer_json.value("layerdown", -1);

        // Read material orientation (degrees in JSON, stored as radians)
        get_orientation(layer_json, "material_orientation",
                       phase.sptr_matprops->psi_mat,
                       phase.sptr_matprops->theta_mat,
                       phase.sptr_matprops->phi_mat, true);

        // Read geometry orientation (degrees in JSON, stored as radians)
        get_orientation(layer_json, "geometry_orientation",
                       sptr_layer->psi_geom,
                       sptr_layer->theta_geom,
                       sptr_layer->phi_geom, true);

        // Read props array
        phase.sptr_matprops->props = get_props_from_json(layer_json);
        phase.sptr_matprops->nprops = phase.sptr_matprops->props.n_elem;
    }
}

void read_cylinder_json(phase_characteristics &rve,
                        const std::string &path_data,
                        const std::string &inputfile) {

    std::filesystem::path filepath = std::filesystem::path(path_data) / inputfile;

    if (!std::filesystem::exists(filepath)) {
        throw std::runtime_error("JSON file not found: " + filepath.string());
    }

    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open JSON file: " + filepath.string());
    }

    json j;
    try {
        file >> j;
    } catch (const json::parse_error& e) {
        throw std::runtime_error("JSON parse error in " + filepath.string() + ": " + e.what());
    }

    if (!j.contains("cylinders") || !j["cylinders"].is_array()) {
        throw std::runtime_error("JSON file must contain 'cylinders' array");
    }

    const auto& cylinders = j["cylinders"];
    unsigned int nphases = cylinders.size();

    // Construct sub-phases with cylinder geometry (type 3)
    rve.sub_phases_construct(nphases, 3, 1);

    for (unsigned int i = 0; i < nphases; i++) {
        const auto& cyl_json = cylinders[i];
        auto& phase = rve.sub_phases[i];

        // Get cylinder geometry pointer
        auto sptr_cylinder = std::dynamic_pointer_cast<cylinder>(phase.sptr_shape);
        if (!sptr_cylinder) {
            throw std::runtime_error("Failed to cast shape to cylinder");
        }

        // Read material properties
        phase.sptr_matprops->number = cyl_json.value("number", static_cast<int>(i));
        phase.sptr_matprops->umat_name = cyl_json.value("umat_name", std::string("ELISO"));
        phase.sptr_matprops->save = cyl_json.value("save", 1);

        // Read concentration
        sptr_cylinder->concentration = cyl_json.value("concentration", 1.0);

        // Read coating info
        sptr_cylinder->coatingof = cyl_json.value("coatingof", 0);
        sptr_cylinder->coatedby = cyl_json.value("coatedby", 0);

        // Read geometry (can be in "geometry" dict or directly)
        if (cyl_json.contains("geometry")) {
            const auto& geom = cyl_json["geometry"];
            sptr_cylinder->L = geom.value("L", 1.0);
            sptr_cylinder->R = geom.value("R", 1.0);
        } else {
            sptr_cylinder->L = cyl_json.value("L", 1.0);
            sptr_cylinder->R = cyl_json.value("R", 1.0);
        }

        // Read material orientation (degrees in JSON, stored as radians)
        get_orientation(cyl_json, "material_orientation",
                       phase.sptr_matprops->psi_mat,
                       phase.sptr_matprops->theta_mat,
                       phase.sptr_matprops->phi_mat, true);

        // Read geometry orientation (degrees in JSON, stored as radians)
        get_orientation(cyl_json, "geometry_orientation",
                       sptr_cylinder->psi_geom,
                       sptr_cylinder->theta_geom,
                       sptr_cylinder->phi_geom, true);

        // Read props array
        phase.sptr_matprops->props = get_props_from_json(cyl_json);
        phase.sptr_matprops->nprops = phase.sptr_matprops->props.n_elem;
    }
}

void read_phase_json(phase_characteristics &rve,
                     const std::string &path_data,
                     const std::string &inputfile) {

    std::filesystem::path filepath = std::filesystem::path(path_data) / inputfile;

    if (!std::filesystem::exists(filepath)) {
        throw std::runtime_error("JSON file not found: " + filepath.string());
    }

    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open JSON file: " + filepath.string());
    }

    json j;
    try {
        file >> j;
    } catch (const json::parse_error& e) {
        throw std::runtime_error("JSON parse error in " + filepath.string() + ": " + e.what());
    }

    if (!j.contains("phases") || !j["phases"].is_array()) {
        throw std::runtime_error("JSON file must contain 'phases' array");
    }

    const auto& phases = j["phases"];
    unsigned int nphases = phases.size();

    // Construct sub-phases with generic geometry (type 0)
    rve.sub_phases_construct(nphases, 0, 1);

    for (unsigned int i = 0; i < nphases; i++) {
        const auto& phase_json = phases[i];
        auto& phase = rve.sub_phases[i];

        // Read material properties
        phase.sptr_matprops->number = phase_json.value("number", static_cast<int>(i));
        phase.sptr_matprops->umat_name = phase_json.value("umat_name", std::string("ELISO"));
        phase.sptr_matprops->save = phase_json.value("save", 1);

        // Read concentration
        phase.sptr_shape->concentration = phase_json.value("concentration", 1.0);

        // Read material orientation (degrees in JSON, stored as radians)
        get_orientation(phase_json, "material_orientation",
                       phase.sptr_matprops->psi_mat,
                       phase.sptr_matprops->theta_mat,
                       phase.sptr_matprops->phi_mat, true);

        // Read props array
        phase.sptr_matprops->props = get_props_from_json(phase_json);
        phase.sptr_matprops->nprops = phase.sptr_matprops->props.n_elem;
    }
}

// Helper to convert radians to degrees
constexpr double rad2deg = 180.0 / M_PI;

// Helper to write props array to JSON
static json props_to_json(const vec& props) {
    json j = json::array();
    for (unsigned int i = 0; i < props.n_elem; i++) {
        j.push_back(props(i));
    }
    return j;
}

// Helper to write orientation to JSON
static json orientation_to_json(double psi, double theta, double phi) {
    return {
        {"psi", psi * rad2deg},
        {"theta", theta * rad2deg},
        {"phi", phi * rad2deg}
    };
}

void write_phase_json(phase_characteristics &rve,
                      const std::string &path_data,
                      const std::string &outputfile) {

    std::filesystem::path filepath = std::filesystem::path(path_data) / outputfile;

    json j;
    j["phases"] = json::array();

    for (auto& r : rve.sub_phases) {
        json phase_json;
        phase_json["number"] = r.sptr_matprops->number;
        phase_json["umat_name"] = r.sptr_matprops->umat_name;
        phase_json["save"] = r.sptr_matprops->save;
        phase_json["concentration"] = r.sptr_shape->concentration;
        phase_json["material_orientation"] = orientation_to_json(
            r.sptr_matprops->psi_mat,
            r.sptr_matprops->theta_mat,
            r.sptr_matprops->phi_mat
        );
        phase_json["nprops"] = r.sptr_matprops->nprops;
        phase_json["nstatev"] = r.sptr_sv_global->nstatev;
        phase_json["props"] = props_to_json(r.sptr_matprops->props);
        j["phases"].push_back(phase_json);
    }

    std::ofstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open JSON file for writing: " + filepath.string());
    }
    file << j.dump(2);
}

void write_ellipsoid_json(phase_characteristics &rve,
                          const std::string &path_data,
                          const std::string &outputfile) {

    std::filesystem::path filepath = std::filesystem::path(path_data) / outputfile;

    json j;
    j["ellipsoids"] = json::array();

    for (auto& r : rve.sub_phases) {
        auto sptr_ellipsoid = std::dynamic_pointer_cast<ellipsoid>(r.sptr_shape);
        if (!sptr_ellipsoid) {
            throw std::runtime_error("Failed to cast shape to ellipsoid");
        }

        json ell_json;
        ell_json["number"] = r.sptr_matprops->number;
        ell_json["coatingof"] = sptr_ellipsoid->coatingof;
        ell_json["umat_name"] = r.sptr_matprops->umat_name;
        ell_json["save"] = r.sptr_matprops->save;
        ell_json["concentration"] = sptr_ellipsoid->concentration;
        ell_json["material_orientation"] = orientation_to_json(
            r.sptr_matprops->psi_mat,
            r.sptr_matprops->theta_mat,
            r.sptr_matprops->phi_mat
        );
        ell_json["semi_axes"] = {
            {"a1", sptr_ellipsoid->a1},
            {"a2", sptr_ellipsoid->a2},
            {"a3", sptr_ellipsoid->a3}
        };
        ell_json["geometry_orientation"] = orientation_to_json(
            sptr_ellipsoid->psi_geom,
            sptr_ellipsoid->theta_geom,
            sptr_ellipsoid->phi_geom
        );
        ell_json["nprops"] = r.sptr_matprops->nprops;
        ell_json["nstatev"] = r.sptr_sv_global->nstatev;
        ell_json["props"] = props_to_json(r.sptr_matprops->props);
        j["ellipsoids"].push_back(ell_json);
    }

    std::ofstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open JSON file for writing: " + filepath.string());
    }
    file << j.dump(2);
}

void write_layer_json(phase_characteristics &rve,
                      const std::string &path_data,
                      const std::string &outputfile) {

    std::filesystem::path filepath = std::filesystem::path(path_data) / outputfile;

    json j;
    j["layers"] = json::array();

    for (auto& r : rve.sub_phases) {
        auto sptr_layer = std::dynamic_pointer_cast<layer>(r.sptr_shape);
        if (!sptr_layer) {
            throw std::runtime_error("Failed to cast shape to layer");
        }

        json layer_json;
        layer_json["number"] = r.sptr_matprops->number;
        layer_json["umat_name"] = r.sptr_matprops->umat_name;
        layer_json["save"] = r.sptr_matprops->save;
        layer_json["concentration"] = sptr_layer->concentration;
        layer_json["material_orientation"] = orientation_to_json(
            r.sptr_matprops->psi_mat,
            r.sptr_matprops->theta_mat,
            r.sptr_matprops->phi_mat
        );
        layer_json["geometry_orientation"] = orientation_to_json(
            sptr_layer->psi_geom,
            sptr_layer->theta_geom,
            sptr_layer->phi_geom
        );
        layer_json["nprops"] = r.sptr_matprops->nprops;
        layer_json["nstatev"] = r.sptr_sv_global->nstatev;
        layer_json["props"] = props_to_json(r.sptr_matprops->props);
        j["layers"].push_back(layer_json);
    }

    std::ofstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open JSON file for writing: " + filepath.string());
    }
    file << j.dump(2);
}

void write_cylinder_json(phase_characteristics &rve,
                         const std::string &path_data,
                         const std::string &outputfile) {

    std::filesystem::path filepath = std::filesystem::path(path_data) / outputfile;

    json j;
    j["cylinders"] = json::array();

    for (auto& r : rve.sub_phases) {
        auto sptr_cylinder = std::dynamic_pointer_cast<cylinder>(r.sptr_shape);
        if (!sptr_cylinder) {
            throw std::runtime_error("Failed to cast shape to cylinder");
        }

        json cyl_json;
        cyl_json["number"] = r.sptr_matprops->number;
        cyl_json["coatingof"] = sptr_cylinder->coatingof;
        cyl_json["umat_name"] = r.sptr_matprops->umat_name;
        cyl_json["save"] = r.sptr_matprops->save;
        cyl_json["concentration"] = sptr_cylinder->concentration;
        cyl_json["material_orientation"] = orientation_to_json(
            r.sptr_matprops->psi_mat,
            r.sptr_matprops->theta_mat,
            r.sptr_matprops->phi_mat
        );
        cyl_json["geometry"] = {
            {"L", sptr_cylinder->L},
            {"R", sptr_cylinder->R}
        };
        cyl_json["geometry_orientation"] = orientation_to_json(
            sptr_cylinder->psi_geom,
            sptr_cylinder->theta_geom,
            sptr_cylinder->phi_geom
        );
        cyl_json["nprops"] = r.sptr_matprops->nprops;
        cyl_json["nstatev"] = r.sptr_sv_global->nstatev;
        cyl_json["props"] = props_to_json(r.sptr_matprops->props);
        j["cylinders"].push_back(cyl_json);
    }

    std::ofstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open JSON file for writing: " + filepath.string());
    }
    file << j.dump(2);
}

} // namespace simcoon
