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
///@brief Static UMAT dispatch singleton - single source of truth for all UMAT selection
///@version 2.0

#pragma once

#include <string>
#include <unordered_map>
#include <armadillo>

namespace simcoon {

/**
 * @brief Singleton class for efficient UMAT dispatch.
 *
 * This class is the single source of truth for UMAT name-to-function mapping.
 * The maps are built once at first access and reused for all subsequent calls,
 * providing significant performance improvement for simulations with many increments.
 *
 * Supports three UMAT categories:
 * - Mechanical small strain (call_umat_M)
 * - Mechanical finite strain (call_umat_M_finite)
 * - Thermomechanical (call_umat_T)
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
     * @brief Check if a UMAT name is registered for mechanical small strain.
     * @param umat_name UMAT identifier to check
     * @return true if UMAT exists, false otherwise
     */
    bool has_umat_M(const std::string& umat_name) const;

    /**
     * @brief Check if a UMAT name is registered for mechanical finite strain.
     * @param umat_name UMAT identifier to check
     * @return true if UMAT exists, false otherwise
     */
    bool has_umat_M_finite(const std::string& umat_name) const;

    /**
     * @brief Check if a UMAT name is registered for thermomechanical.
     * @param umat_name UMAT identifier to check
     * @return true if UMAT exists, false otherwise
     */
    bool has_umat_T(const std::string& umat_name) const;

    /**
     * @brief Get the dispatch ID for a mechanical small strain UMAT.
     * @param umat_name UMAT identifier
     * @return dispatch ID, or -1 if not found
     */
    int get_umat_M_id(const std::string& umat_name) const;

    /**
     * @brief Get the dispatch ID for a mechanical finite strain UMAT.
     * @param umat_name UMAT identifier
     * @return dispatch ID, or -1 if not found
     */
    int get_umat_M_finite_id(const std::string& umat_name) const;

    /**
     * @brief Get the dispatch ID for a thermomechanical UMAT.
     * @param umat_name UMAT identifier
     * @return dispatch ID, or -1 if not found
     */
    int get_umat_T_id(const std::string& umat_name) const;

    /**
     * @brief Call the appropriate UMAT function for mechanical small strain problems.
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
     * @brief Call the appropriate UMAT function for mechanical finite strain problems.
     *
     * @param umat_name 5-character UMAT identifier
     * @param etot Total logarithmic strain
     * @param Detot Strain increment
     * @param F0 Deformation gradient at start of increment
     * @param F1 Deformation gradient at end of increment
     * @param sigma Cauchy stress (output)
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
    bool call_umat_M_finite(
        const std::string& umat_name,
        const arma::vec& etot, const arma::vec& Detot,
        const arma::mat& F0, const arma::mat& F1,
        arma::vec& sigma, arma::mat& Lt,
        arma::mat& L, arma::vec& sigma_in,
        const arma::mat& DR, int nprops, const arma::vec& props,
        int nstatev, arma::vec& statev,
        double T, double DT, double Time, double DTime,
        double& Wm, double& Wm_r, double& Wm_ir, double& Wm_d,
        int ndi, int nshr, bool start, int solver_type, double& tnew_dt);

    /**
     * @brief Call the appropriate UMAT function for thermomechanical problems.
     *
     * @param umat_name 5-character UMAT identifier
     * @param Etot Total strain
     * @param DEtot Strain increment
     * @param sigma Stress (output)
     * @param r Heat source (output)
     * @param dSdE Mechanical tangent (output)
     * @param dSdT Thermal-mechanical tangent (output)
     * @param drdE Heat source/strain tangent (output)
     * @param drdT Heat source/temperature tangent (output)
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
     * @param Wt0 Thermal work component 0 (output)
     * @param Wt1 Thermal work component 1 (output)
     * @param Wt2 Thermal work component 2 (output)
     * @param ndi Number of direct stress components
     * @param nshr Number of shear stress components
     * @param start Flag for first increment
     * @param tnew_dt New time step (output)
     * @return true if UMAT was found and called, false otherwise
     */
    bool call_umat_T(
        const std::string& umat_name,
        const arma::vec& Etot, const arma::vec& DEtot,
        arma::vec& sigma, double& r,
        arma::mat& dSdE, arma::mat& dSdT,
        arma::mat& drdE, arma::mat& drdT,
        const arma::mat& DR, int nprops, const arma::vec& props,
        int nstatev, arma::vec& statev,
        double T, double DT, double Time, double DTime,
        double& Wm, double& Wm_r, double& Wm_ir, double& Wm_d,
        double& Wt0, double& Wt1, double& Wt2,
        int ndi, int nshr, bool start, double& tnew_dt);

    // Legacy alias for backward compatibility
    bool has_umat(const std::string& umat_name) const { return has_umat_M(umat_name); }

private:
    UmatDispatch();
    ~UmatDispatch() = default;

    std::unordered_map<std::string, int> dispatch_map_M_;
    std::unordered_map<std::string, int> dispatch_map_M_finite_;
    std::unordered_map<std::string, int> dispatch_map_T_;
};

} // namespace simcoon
