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
* @file fea_transfer.hpp
* @author Yves Chemisky 
* @brief Functions to transfer data between FEA software formats and simcoon internal format
* @section FEA Transfer Functions
*
* This file contains functions for converting between:
* - Abaqus UMAT format ↔ simcoon Armadillo format
* - Ansys USERMAT format ↔ simcoon Armadillo format
*
* Voigt notation conventions:
* - simcoon/Abaqus: (11, 22, 33, 12, 13, 23) → indices (0, 1, 2, 3, 4, 5)
* - Ansys:          (11, 22, 33, 12, 23, 13) → indices (0, 1, 2, 3, 4, 5)
*   Note: Ansys swaps components 4 and 5 compared to simcoon/Abaqus
*/

#pragma once
#include <armadillo>

namespace simcoon {

/** @addtogroup umat_mechanical
 *  @{
 */


//=============================================================================
// ABAQUS UMAT TRANSFER FUNCTIONS
//=============================================================================

/**
 * @brief Light transfer from Abaqus output arrays to simcoon format (mechanical)
 * @details Converts stress, tangent, and state variables from Abaqus format.
 *          Use this after calling an external Abaqus UMAT to retrieve results.
 * 
 * @param stress Abaqus stress array (ntens)
 * @param ddsdde Abaqus tangent operator array (ntens×ntens, column-major)
 * @param nstatev Number of state variables
 * @param statev Abaqus state variables array
 * @param ndi Number of direct stress components
 * @param nshr Number of shear stress components
 * @param sigma [out] simcoon stress vector (6)
 * @param Lt [out] simcoon tangent matrix (6×6)
 * @param Wm [out] simcoon work quantities vector (4): Wm, Wm_r, Wm_ir, Wm_d
 * @param statev_smart [out] simcoon state variables vector
 */
void abaqus2smart_M_light(const double *stress, const double *ddsdde, 
    const int &nstatev, double *statev, const int &ndi, const int &nshr, 
    arma::vec &sigma, arma::mat &Lt, arma::vec &Wm, arma::vec &statev_smart);

/**
 * @brief Full transfer from Abaqus input arrays to simcoon format (mechanical)
 * @details Converts all UMAT inputs from Abaqus format to simcoon Armadillo format.
 *          Use this before calling a simcoon constitutive model from an Abaqus UMAT wrapper.
 * 
 * @param stress Abaqus stress array
 * @param ddsdde Abaqus tangent operator array
 * @param stran Abaqus total strain array
 * @param dstran Abaqus strain increment array
 * @param time Abaqus time array [step_time, total_time]
 * @param dtime Time increment
 * @param temperature Temperature at start of increment
 * @param Dtemperature Temperature increment
 * @param nprops Number of material properties
 * @param props Abaqus material properties array
 * @param nstatev Number of state variables
 * @param statev Abaqus state variables array
 * @param ndi Number of direct stress components
 * @param nshr Number of shear stress components
 * @param drot Abaqus rotation increment matrix (3×3, column-major)
 * @param sigma [out] simcoon stress vector
 * @param Lt [out] simcoon tangent matrix
 * @param Etot [out] simcoon total strain vector
 * @param DEtot [out] simcoon strain increment vector
 * @param T [out] simcoon temperature
 * @param DT [out] simcoon temperature increment
 * @param Time [out] simcoon total time
 * @param DTime [out] simcoon time increment
 * @param props_smart [out] simcoon properties vector
 * @param Wm [out] simcoon work quantities vector
 * @param statev_smart [out] simcoon state variables vector
 * @param DR [out] simcoon rotation matrix
 * @param start [out] true if this is the first increment
 */
void abaqus2smart_M(const double *stress, const double *ddsdde, 
    const double *stran, const double *dstran, const double *time, 
    const double &dtime, const double &temperature, const double &Dtemperature, 
    const int &nprops, const double *props, const int &nstatev, double *statev, 
    const int &ndi, const int &nshr, const double *drot, 
    arma::vec &sigma, arma::mat &Lt, arma::vec &Etot, arma::vec &DEtot, 
    double &T, double &DT, double &Time, double &DTime, 
    arma::vec &props_smart, arma::vec &Wm, arma::vec &statev_smart, 
    arma::mat &DR, bool &start);

/**
 * @brief Full transfer from Abaqus input arrays to simcoon format (thermomechanical)
 * @details Converts all thermomechanical UMAT inputs from Abaqus format.
 *          Includes thermal tangent operators and heat generation terms.
 * 
 * @param stress Abaqus stress array
 * @param ddsdde Abaqus mechanical tangent array
 * @param ddsddt Abaqus stress-temperature tangent array
 * @param drplde Abaqus heat generation-strain tangent array
 * @param drpldt Abaqus heat generation-temperature tangent
 * @param stran Abaqus total strain array
 * @param dstran Abaqus strain increment array
 * @param time Abaqus time array
 * @param dtime Time increment
 * @param temperature Temperature at start of increment
 * @param Dtemperature Temperature increment
 * @param nprops Number of material properties
 * @param props Abaqus material properties array
 * @param nstatev Number of state variables
 * @param statev Abaqus state variables array
 * @param ndi Number of direct stress components
 * @param nshr Number of shear stress components
 * @param drot Abaqus rotation increment matrix
 * @param sigma [out] simcoon stress vector
 * @param dSdE [out] simcoon mechanical tangent matrix
 * @param dSdT [out] simcoon stress-temperature tangent matrix
 * @param drpldE [out] simcoon heat generation-strain tangent matrix
 * @param drpldT [out] simcoon heat generation-temperature tangent matrix
 * @param Etot [out] simcoon total strain vector
 * @param DEtot [out] simcoon strain increment vector
 * @param T [out] simcoon temperature
 * @param DT [out] simcoon temperature increment
 * @param Time [out] simcoon total time
 * @param DTime [out] simcoon time increment
 * @param props_smart [out] simcoon properties vector
 * @param Wm [out] simcoon mechanical work quantities vector
 * @param Wt [out] simcoon thermal work quantities vector
 * @param statev_smart [out] simcoon state variables vector
 * @param DR [out] simcoon rotation matrix
 * @param start [out] true if this is the first increment
 */
void abaqus2smart_T(const double *stress, const double *ddsdde, 
    const double *ddsddt, const double *drplde, const double &drpldt, 
    const double *stran, const double *dstran, const double *time, 
    const double &dtime, const double &temperature, const double &Dtemperature, 
    const int &nprops, const double *props, const int &nstatev, double *statev, 
    const int &ndi, const int &nshr, const double *drot, 
    arma::vec &sigma, arma::mat &dSdE, arma::mat &dSdT, 
    arma::mat &drpldE, arma::mat &drpldT, 
    arma::vec &Etot, arma::vec &DEtot, double &T, double &DT, 
    double &Time, double &DTime, arma::vec &props_smart, 
    arma::vec &Wm, arma::vec &Wt, arma::vec &statev_smart, 
    arma::mat &DR, bool &start);

/**
 * @brief Transfer simcoon output to Abaqus format (mechanical, simple)
 * @details Converts stress, tangent, and state variables to Abaqus format.
 *          Use this to return results from a simcoon model to an Abaqus UMAT.
 * 
 * @param stress [out] Abaqus stress array
 * @param ddsdde [out] Abaqus tangent operator array
 * @param statev [out] Abaqus state variables array
 * @param ndi Number of direct stress components
 * @param nshr Number of shear stress components
 * @param sigma simcoon stress vector
 * @param statev_smart simcoon state variables vector
 * @param Wm simcoon work quantities vector
 * @param Lt simcoon tangent matrix
 */
void smart2abaqus_M(double *stress, double *ddsdde, double *statev, 
    const int &ndi, const int &nshr, 
    const arma::vec &sigma, const arma::vec &statev_smart, 
    const arma::vec &Wm, const arma::mat &Lt);

/**
 * @brief Full transfer from simcoon to Abaqus format (mechanical)
 * @details Converts all simcoon outputs to Abaqus UMAT format.
 *          Use this for complete state transfer in external UMAT wrappers.
 * 
 * @param stress [out] Abaqus stress array
 * @param ddsdde [out] Abaqus tangent operator array
 * @param stran [out] Abaqus total strain array
 * @param dstran [out] Abaqus strain increment array
 * @param time [out] Abaqus time array
 * @param dtime [out] Time increment
 * @param temperature [out] Temperature
 * @param Dtemperature [out] Temperature increment
 * @param nprops Number of material properties
 * @param props [out] Abaqus material properties array
 * @param nstatev Number of state variables
 * @param statev [out] Abaqus state variables array
 * @param ndi Number of direct stress components
 * @param nshr Number of shear stress components
 * @param drot [out] Abaqus rotation increment matrix
 * @param sigma simcoon stress vector
 * @param Lt simcoon tangent matrix
 * @param Etot simcoon total strain vector
 * @param DEtot simcoon strain increment vector
 * @param T simcoon temperature
 * @param DT simcoon temperature increment
 * @param Time simcoon total time
 * @param DTime simcoon time increment
 * @param props_smart simcoon properties vector
 * @param Wm simcoon work quantities vector
 * @param statev_smart simcoon state variables vector
 * @param DR simcoon rotation matrix
 * @param start simcoon start flag
 */
void smart2abaqus_M_full(double *stress, double *ddsdde, double *stran, 
    double *dstran, double *time, double &dtime, 
    double &temperature, double &Dtemperature, 
    int &nprops, double *props, int &nstatev, double *statev, 
    const int &ndi, const int &nshr, double *drot, 
    const arma::vec &sigma, const arma::mat &Lt, 
    const arma::vec &Etot, const arma::vec &DEtot, 
    const double &T, const double &DT, const double &Time, const double &DTime, 
    const arma::vec &props_smart, const arma::vec &Wm, 
    const arma::vec &statev_smart, const arma::mat &DR, bool &start);

/**
 * @brief Transfer simcoon output to Abaqus format (thermomechanical)
 * @details Converts thermomechanical outputs to Abaqus UMAT format.
 *          Includes thermal tangent operators and heat generation.
 * 
 * @param stress [out] Abaqus stress array
 * @param ddsdde [out] Abaqus mechanical tangent array
 * @param ddsddt [out] Abaqus stress-temperature tangent array
 * @param drplde [out] Abaqus heat generation-strain tangent array
 * @param drpldt [out] Abaqus heat generation-temperature tangent
 * @param rpl [out] Abaqus heat generation rate
 * @param statev [out] Abaqus state variables array
 * @param ndi Number of direct stress components
 * @param nshr Number of shear stress components
 * @param sigma simcoon stress vector
 * @param statev_smart simcoon state variables vector
 * @param r simcoon heat generation rate
 * @param Wm simcoon mechanical work quantities vector
 * @param Wt simcoon thermal work quantities vector
 * @param dSdE simcoon mechanical tangent matrix
 * @param dSdT simcoon stress-temperature tangent matrix
 * @param drpldE simcoon heat generation-strain tangent matrix
 * @param drpldT simcoon heat generation-temperature tangent matrix
 */
void smart2abaqus_T(double *stress, double *ddsdde, double *ddsddt, 
    double *drplde, double &drpldt, double &rpl, double *statev, 
    const int &ndi, const int &nshr, 
    const arma::vec &sigma, const arma::vec &statev_smart, const double &r, 
    const arma::vec &Wm, const arma::vec &Wt, 
    const arma::mat &dSdE, const arma::mat &dSdT, 
    const arma::mat &drpldE, const arma::mat &drpldT);

//=============================================================================
// ANSYS USERMAT TRANSFER FUNCTIONS
//=============================================================================

/**
 * @brief Transfer Ansys USERMAT input arrays to simcoon format (mechanical)
 * @details Converts all USERMAT inputs from Ansys format to simcoon format.
 *          Handles Ansys Voigt notation: (11,22,33,12,23,13) by swapping 
 *          components 4 and 5.
 * 
 * @param stress Ansys stress array (ncomp)
 * @param dstran Ansys strain increment array (ncomp)
 * @param sedEl Specific elastic strain energy
 * @param sedPl Specific plastic strain energy
 * @param epseq Equivalent plastic strain
 * @param statev Ansys state variables array
 * @param props Ansys material properties array
 * @param Time Total time
 * @param DTime Time increment
 * @param temperature Temperature
 * @param Dtemperature Temperature increment
 * @param ncomp Number of stress/strain components
 * @param nprops Number of material properties
 * @param nstatev Number of state variables
 * @param sigma [out] simcoon stress vector
 * @param DEtot [out] simcoon strain increment vector
 * @param Wm [out] simcoon work quantities vector
 * @param statev_smart [out] simcoon state variables vector
 * @param props_smart [out] simcoon properties vector
 * @param T [out] simcoon temperature
 * @param DT [out] simcoon temperature increment
 */
void ansys2smart_M(const double *stress, const double *dstran,
    const double &sedEl, const double &sedPl, const double &epseq,
    const double *statev, const double *props,
    const double &Time, const double &DTime,
    const double &temperature, const double &Dtemperature,
    const int &ncomp, const int &nprops, const int &nstatev,
    arma::vec &sigma, arma::vec &DEtot, arma::vec &Wm, 
    arma::vec &statev_smart, arma::vec &props_smart,
    double &T, double &DT);

/**
 * @brief Transfer simcoon output to Ansys USERMAT format (mechanical)
 * @details Converts stress, tangent, and state variables to Ansys format.
 *          Handles Ansys Voigt notation by swapping components 4 and 5.
 * 
 * @param stress [out] Ansys stress array
 * @param ddsdde [out] Ansys tangent operator array
 * @param sedEl [out] Specific elastic strain energy
 * @param sedPl [out] Specific plastic strain energy
 * @param statev [out] Ansys state variables array
 * @param ncomp Number of stress/strain components
 * @param sigma simcoon stress vector
 * @param Lt simcoon tangent matrix
 * @param Wm simcoon work quantities vector
 * @param statev_smart simcoon state variables vector
 */
void smart2ansys_M(double *stress, double *ddsdde,
    double &sedEl, double &sedPl, double *statev,
    const int &ncomp,
    const arma::vec &sigma, const arma::mat &Lt,
    const arma::vec &Wm, const arma::vec &statev_smart);


/** @} */ // end of umat_mechanical group

} // namespace simcoon
