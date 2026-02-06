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

/**
 * @file umat_plugin_api.hpp
 * @brief Plugin API for external user material subroutines
 * @author Yves Chemisky
 * @version 1.0
 *
 * This file defines abstract base classes for implementing custom user material
 * subroutines (UMATs) as dynamically loadable plugins. Three plugin interfaces
 * are provided for different use cases:
 *
 * - **umat_plugin_ext_api**: Generic external UMAT interface
 * - **umat_plugin_aba_api**: Abaqus-compatible UMAT interface
 * - **umat_plugin_ans_api**: ANSYS-compatible UMAT interface
 *
 * To create a custom plugin, derive from the appropriate base class and implement
 * all pure virtual methods. The plugin can then be loaded at runtime using the
 * simcoon plugin loader.
 *
 * @see phase_characteristics for material state storage
 */

#pragma once

#include <string>
#include <armadillo>
#include <simcoon/Simulation/Phase/phase_characteristics.hpp>

/** @addtogroup umat_plugins
 *  @{
 */

/**
 * @brief Abstract base class for external mechanical UMAT plugins
 *
 * @details This interface allows users to implement custom constitutive models
 * that can be loaded as shared libraries at runtime. The interface follows
 * the standard simcoon UMAT calling convention.
 *
 * **Implementation Example:**
 * @code
 * class MyCustomUmat : public umat_plugin_ext_api {
 * public:
 *     std::string name() const override { return "my_custom_umat"; }
 *
 *     void umat_external_M(
 *         const arma::vec &Etot, const arma::vec &DEtot, arma::vec &sigma,
 *         arma::mat &Lt, arma::mat &L, arma::vec &sigma_in,
 *         const arma::mat &DR, const int &nprops,
 *         const arma::vec &props, const int &nstatev, arma::vec &statev,
 *         const double &T, const double &DT, const double &Time, const double &DTime,
 *         double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d,
 *         const int &ndi, const int &nshr, const bool &start, const int &solver_type, double &tnew_dt
 *     ) override {
 *         // Your constitutive model implementation here
 *     }
 * };
 *
 * // Export function for plugin loading
 * extern "C" umat_plugin_ext_api* create_plugin() {
 *     return new MyCustomUmat();
 * }
 * @endcode
 */
class umat_plugin_ext_api {
    public:
        /**
         * @brief Returns the unique name identifier for this plugin
         * @return String containing the plugin name
         */
        virtual std::string name() const = 0;

        /**
         * @brief External mechanical UMAT subroutine
         *
         * @details Implements a constitutive model following the simcoon UMAT convention.
         * This function is called at each integration point to update stress and
         * compute the consistent tangent modulus.
         *
         * @param umat_name Name of the constitutive model
         * @param Etot Total strain tensor at beginning of increment (Voigt notation: 6×1)
         * @param DEtot Strain increment tensor (Voigt notation: 6×1)
         * @param sigma Stress tensor [output] (Voigt notation: 6×1)
         * @param Lt Consistent tangent modulus [output] (6×6)
         * @param L Elastic stiffness tensor [output] (6×6)
         * @param DR Rotation increment matrix (3×3)
         * @param nprops Number of material properties
         * @param props Material properties vector
         * @param nstatev Number of state variables
         * @param statev State variables vector [input/output]
         * @param T Temperature at beginning of increment
         * @param DT Temperature increment
         * @param Time Time at beginning of increment
         * @param DTime Time increment
         * @param Wm Total mechanical work [output]
         * @param Wm_r Recoverable work [output]
         * @param Wm_ir Irrecoverable work [output]
         * @param Wm_d Dissipated work [output]
         * @param ndi Number of direct stress components
         * @param nshr Number of shear stress components
         * @param start Flag indicating first increment
         * @param tnew_dt Suggested new time step [output]
         */
        virtual void umat_external_M(
            const std::string &umat_name, const arma::vec &Etot, const arma::vec &DEtot, arma::vec &sigma,
            arma::mat &Lt, arma::mat &L,
            const arma::mat &DR, const int &nprops,
            const arma::vec &props, const int &nstatev, arma::vec &statev,
            const double &T, const double &DT, const double &Time, const double &DTime,
            double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d,
            const int &ndi, const int &nshr, const bool &start, double &tnew_dt
        ) = 0;

        /** @brief Virtual destructor for proper cleanup */
        virtual ~umat_plugin_ext_api() = default;
};

/**
 * @brief Abstract base class for Abaqus-compatible UMAT plugins
 *
 * @details This interface provides compatibility with Abaqus user material
 * subroutines. It wraps the phase_characteristics structure to provide
 * a familiar interface for users migrating from Abaqus.
 *
 * **Abaqus Compatibility:**
 * - Uses the same state variable convention as Abaqus
 * - Supports the same solver types
 * - Handles Abaqus-specific rotation conventions
 */
class umat_plugin_aba_api {
    public:
        /**
         * @brief Returns the unique name identifier for this plugin
         * @return String containing the plugin name
         */
        virtual std::string name() const = 0;

        /**
         * @brief Abaqus-compatible UMAT subroutine
         *
         * @details Implements a constitutive model compatible with Abaqus UMAT conventions.
         * The phase_characteristics object contains all material state information.
         *
         * @param phase Material phase containing state, properties, strain, stress [input/output]
         * @param DR Rotation increment matrix (3×3) for objective integration
         * @param Time Time at beginning of increment
         * @param DTime Time increment
         * @param ndi Number of direct stress components
         * @param nshr Number of shear stress components
         * @param start Flag indicating first increment [input/output]
         * @param solver_type Solver type (0=implicit, 1=explicit)
         * @param tnew_dt Suggested new time step [output]
         */
        virtual void umat_abaqus(
            simcoon::phase_characteristics &phase,
            const arma::mat &DR,
            const double &Time, const double &DTime,
            const int &ndi, const int &nshr,
            bool &start, const int &solver_type, double &tnew_dt
        ) = 0;

        /** @brief Virtual destructor for proper cleanup */
        virtual ~umat_plugin_aba_api() = default;
};

/**
 * @brief Abstract base class for ANSYS-compatible UMAT plugins
 *
 * @details This interface provides compatibility with ANSYS user material
 * subroutines (USERMAT). It wraps the phase_characteristics structure to
 * provide a familiar interface for users migrating from ANSYS.
 *
 * **ANSYS Compatibility:**
 * - Uses conventions compatible with ANSYS USERMAT
 * - Supports ANSYS solver types
 * - Handles ANSYS-specific state variable conventions
 */
class umat_plugin_ans_api {
    public:
        /**
         * @brief Returns the unique name identifier for this plugin
         * @return String containing the plugin name
         */
        virtual std::string name() const = 0;

        /**
         * @brief ANSYS-compatible UMAT subroutine
         *
         * @details Implements a constitutive model compatible with ANSYS USERMAT conventions.
         * The phase_characteristics object contains all material state information.
         *
         * @param phase Material phase containing state, properties, strain, stress [input/output]
         * @param DR Rotation increment matrix (3×3) for objective integration
         * @param Time Time at beginning of increment
         * @param DTime Time increment
         * @param ndi Number of direct stress components
         * @param nshr Number of shear stress components
         * @param start Flag indicating first increment [input/output]
         * @param solver_type Solver type
         * @param tnew_dt Suggested new time step [output]
         */
        virtual void umat_ansys(
            simcoon::phase_characteristics &phase,
            const arma::mat &DR,
            const double &Time, const double &DTime,
            const int &ndi, const int &nshr,
            bool &start, const int &solver_type, double &tnew_dt
        ) = 0;

        /** @brief Virtual destructor for proper cleanup */
        virtual ~umat_plugin_ans_api() = default;
};

/** @} */ // end of umat_plugins group
