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

///@file external_umat.hpp
///@brief External user material subroutine interface for custom mechanical models
///@version 1.0

#pragma once

#include <armadillo>

/**
 * @file external_umat.hpp
 * @brief External user material subroutine interface for custom mechanical models
 * @author Yves Chemisky
 * @version 1.0
 *
 * This file declares the external interface for user-defined mechanical constitutive
 * models that can be loaded as plugins at runtime. Users implement their material
 * models through the plugin API system.
 */

/** @addtogroup umat_mechanical
 *  @{
 */

/**
 * @brief External mechanical UMAT subroutine interface
 *
 * @details This function declaration provides the interface for external user-defined
 * material subroutines (UMATs). The actual implementation is provided by the user
 * through the plugin system, allowing custom constitutive models to be loaded
 * dynamically at runtime.
 *
 * **Plugin System:**
 *
 * To implement a custom material model, users should:
 * 1. Create a class inheriting from `umat_plugin_ext_api`
 * 2. Implement the `umat_external_M()` method with the constitutive equations
 * 3. Compile as a shared library (`.so`, `.dylib`, or `.dll`)
 * 4. Configure the material name in `material.dat` to match the plugin's `name()`
 *
 * **Available Plugin Formats:**
 *
 * | Format | Class | Description |
 * |--------|-------|-------------|
 * | UMEXT | `umat_plugin_ext_api` | Native simcoon format using Armadillo |
 * | UMABA | `umat_plugin_aba_api` | Abaqus UMAT-compatible format |
 * | UMANS | `umat_plugin_ans_api` | ANSYS USERMAT-compatible format |
 *
 * **Implementation Example:**
 *
 * @code
 * #include <simcoon/Continuum_mechanics/Umat/umat_plugin_api.hpp>
 *
 * class MyMaterial : public umat_plugin_ext_api {
 * public:
 *     std::string name() const override { return "my_material"; }
 *
 *     void umat_external_M(
 *         const arma::vec &Etot, const arma::vec &DEtot, arma::vec &sigma,
 *         arma::mat &Lt, arma::mat &L, arma::vec &sigma_in,
 *         const arma::mat &DR, const int &nprops,
 *         const arma::vec &props, const int &nstatev, arma::vec &statev,
 *         const double &T, const double &DT, const double &Time, const double &DTime,
 *         double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d,
 *         const int &ndi, const int &nshr, const bool &start,
 *         const int &solver_type, double &tnew_dt
 *     ) override {
 *         // Implement your constitutive model here
 *         // 1. Extract material properties from props
 *         // 2. Update state variables in statev
 *         // 3. Compute stress tensor sigma
 *         // 4. Compute tangent modulus Lt (for implicit solvers)
 *         // 5. Update work quantities Wm, Wm_r, Wm_ir, Wm_d
 *     }
 * };
 *
 * // Required factory functions for plugin loading
 * extern "C" umat_plugin_ext_api* create_api() {
 *     return new MyMaterial();
 * }
 *
 * extern "C" void destroy_api(umat_plugin_ext_api* p) {
 *     delete p;
 * }
 * @endcode
 *
 * **Compiling the Plugin:**
 *
 * @code{.sh}
 * # On macOS
 * clang++ -c -fPIC -std=c++14 my_material.cpp -I/path/to/simcoon/include
 * clang++ -std=c++14 -shared -lsimcoon -larmadillo -o libmy_material.dylib my_material.o
 *
 * # On Linux
 * g++ -c -fPIC -std=c++14 my_material.cpp -I/path/to/simcoon/include
 * g++ -std=c++14 -shared -lsimcoon -larmadillo -o libmy_material.so my_material.o
 * @endcode
 *
 * **Material Configuration (material.dat):**
 *
 * @code{.txt}
 * Material
 * Name    MY_MATERIAL
 * Number_of_material_parameters   3
 * Number_of_internal_variables    1
 *
 * #Mechanical
 * E       210000
 * nu      0.3
 * alpha   1.2e-5
 * @endcode
 *
 * **Material Parameters (props):**
 *
 * User-defined. The number and meaning of material parameters depends on
 * the specific constitutive model implemented in the plugin.
 *
 * **State Variables (statev):**
 *
 * User-defined. The number and meaning of state variables depends on
 * the specific constitutive model implemented in the plugin.
 *
 * @param Etot Total strain tensor at beginning of increment (Voigt notation: 6×1 vector)
 * @param DEtot Strain increment tensor (Voigt notation: 6×1 vector)
 * @param sigma Stress tensor (Voigt notation: 6×1 vector) [output]
 * @param Lt Consistent tangent modulus (6×6 matrix) [output]
 * @param L Elastic stiffness tensor (6×6 matrix) [output]
 * @param sigma_in Internal stress contribution for explicit solvers (6×1 vector) [output]
 * @param DR Rotation increment matrix (3×3) for objective integration
 * @param nprops Number of material properties
 * @param props Material properties vector
 * @param nstatev Number of state variables
 * @param statev State variables vector [input/output]
 * @param T Temperature at beginning of increment
 * @param DT Temperature increment
 * @param Time Time at beginning of increment
 * @param DTime Time increment
 * @param Wm Total mechanical work [output]
 * @param Wm_r Recoverable (elastic) work [output]
 * @param Wm_ir Irrecoverable work [output]
 * @param Wm_d Dissipated work [output]
 * @param ndi Number of direct stress components (typically 3)
 * @param nshr Number of shear stress components (typically 3)
 * @param start Flag indicating first increment (true) or continuation (false)
 * @param solver_type Solver type: 0=implicit (Newton), 1=explicit (RNL), 2=dynamic implicit
 * @param tnew_dt Suggested new time step ratio for adaptive time stepping [output]
 *
 * @note This is an external declaration; the function must be provided by a loaded plugin
 * @note The plugin shared library must be in the runtime library path
 * @note Use `start` flag to initialize state variables on first call
 * @note For implicit solvers (solver_type=0,2), compute consistent tangent Lt
 * @note For explicit solvers (solver_type=1), compute internal stress sigma_in
 *
 * @see umat_plugin_ext_api for the plugin base class interface
 * @see umat_plugin_aba_api for Abaqus-compatible plugins
 * @see umat_plugin_ans_api for ANSYS-compatible plugins
 */
extern void umat_external(const arma::vec &Etot, const arma::vec &DEtot, arma::vec &sigma, arma::mat &Lt, arma::mat &L, arma::vec &sigma_in, const arma::mat &DR, const int &nprops, const arma::vec &props, const int &nstatev, arma::vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, const int &solver_type, double &tnew_dt);

/** @} */ // end of umat_mechanical group
