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

///@file umat_dispatch.hpp
///@brief Static UMAT dispatch singleton for optimized solver
///@version 1.0

#pragma once

#include <string>
#include <unordered_map>
#include <armadillo>

namespace simcoon {

/**
 * @brief Singleton class for efficient UMAT dispatch.
 *
 * This class replaces the pattern of rebuilding std::map<string,int> on every
 * UMAT call. The map is built once at first access and reused for all
 * subsequent calls, providing significant performance improvement for
 * simulations with many increments.
 *
 * Thread-safe initialization via C++11 magic statics.
 */
class UmatDispatch {
public:
    /**
     * @brief Get the singleton instance.
     * @return Reference to the singleton UmatDispatch instance
     */
    static UmatDispatch& instance();

    // Prevent copying
    UmatDispatch(const UmatDispatch&) = delete;
    UmatDispatch& operator=(const UmatDispatch&) = delete;

    /**
     * @brief Call the appropriate UMAT function for mechanical problems.
     *
     * @param umat_name 5-character UMAT identifier
     * @param Etot Total strain
     * @param DEtot Strain increment
     * @param sigma Stress (output)
     * @param Lt Algorithmic tangent (output)
     * @param L Elastic stiffness (output)
     * @param sigma_in Internal stress (output)
     * @param DR Rotation increment
     * @param nprops Number of material properties
     * @param props Material properties
     * @param nstatev Number of state variables
     * @param statev State variables (input/output)
     * @param T Temperature
     * @param DT Temperature increment
     * @param Time Current time
     * @param DTime Time increment
     * @param Wm Mechanical work (output)
     * @param Wm_r Recoverable work (output)
     * @param Wm_ir Irrecoverable work (output)
     * @param Wm_d Dissipated work (output)
     * @param ndi Number of direct stress components
     * @param nshr Number of shear stress components
     * @param start Flag for first increment
     * @param solver_type Solver type identifier
     * @param tnew_dt New time step (output)
     * @return true if UMAT was found and called, false otherwise
     */
    bool call_umat_M(
        const std::string& umat_name,
        const arma::vec& Etot, const arma::vec& DEtot,
        arma::vec& sigma, arma::mat& Lt,
        arma::mat& L, arma::vec& sigma_in,
        const arma::mat& DR, int nprops, const arma::vec& props,
        int nstatev, arma::vec& statev,
        double T, double DT, double Time, double DTime,
        double& Wm, double& Wm_r, double& Wm_ir, double& Wm_d,
        int ndi, int nshr, bool start, int solver_type, double& tnew_dt);

    /**
     * @brief Check if a UMAT name is registered.
     * @param umat_name UMAT identifier to check
     * @return true if UMAT exists, false otherwise
     */
    bool has_umat(const std::string& umat_name) const;

private:
    UmatDispatch();
    ~UmatDispatch() = default;

    std::unordered_map<std::string, int> dispatch_map_;
};

} // namespace simcoon
