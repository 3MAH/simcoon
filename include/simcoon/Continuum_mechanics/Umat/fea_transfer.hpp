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
* - Abaqus UMAT format \f$\leftrightarrow\f$ simcoon Armadillo format
* - Ansys USERMAT format \f$\leftrightarrow\f$ simcoon Armadillo format
*
* Voigt notation conventions:
* - simcoon/Abaqus: (11, 22, 33, 12, 13, 23) \f$\to\f$ indices (0, 1, 2, 3, 4, 5)
* - Ansys:          (11, 22, 33, 12, 23, 13) \f$\to\f$ indices (0, 1, 2, 3, 4, 5)
*   Note: Ansys swaps components 4 and 5 compared to simcoon/Abaqus
*/

#pragma once
#include <exception>
#include <string>
#include <armadillo>

namespace simcoon {

class state_variables;
class state_variables_T;

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
 * @param ddsdde Abaqus tangent operator array (ntens \f$\times\f$ ntens, column-major)
 * @param nstatev Number of state variables
 * @param statev Abaqus state variables array
 * @param ndi Number of direct stress components
 * @param nshr Number of shear stress components
 * @param sigma [out] simcoon stress vector (6)
 * @param Lt [out] simcoon tangent matrix (\f$6 \times 6\f$)
 * @param Wm [out] simcoon work quantities vector (4): Wm, Wm_r, Wm_ir, Wm_d
 * @param statev_smart [out] simcoon state variables vector
 */
void abaqus2smart_M_light(const double *stress, const double *ddsdde, 
    const int &nstatev, double *statev, const int &ndi, const int &nshr, 
    arma::vec &sigma, arma::mat &Lt, arma::vec &Wm, arma::vec &statev_smart);

/**
 * @brief Full transfer from Abaqus input arrays to simcoon format (mechanical)
 * @details Converts all UMAT inputs from Abaqus format to simcoon Armadillo format.
 *          STATEV layout: Wm (4), then, when ndi < 3, the 6 - ntens stress-free strain
 *          components the material point solves (condensed_state.hpp, simcoon order, rotated
 *          with DROT and returned in Etot), then the law's state variables (nstatev - 4 - (6 - ntens)).
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
 * @param drot Abaqus rotation increment matrix (\f$3 \times 3\f$, column-major)
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
 *          STATEV layout: Wm (4), Wt (3), then the stress-free strain components when
 *          ndi < 3 (as abaqus2smart_M), then the law's state variables.
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
 * @param Lt simcoon tangent matrix (dsigma/deps of the corotated kernel), condensed on the
 *        element's components when ndi < 3 (select_umat_M): ddsdde is its block on them
 * @param nlgeom whether the step is geometrically nonlinear (abaqus_nlgeom(dfgrd1)): ddsdde
 *        then gets abaqus_jacobian of Lt for the 3D, plane-strain and axisymmetric cases;
 *        plane stress and 1D take the condensed tangent as it is
 * @param Etot_end total strain at the end of the increment (Etot + DEtot): when ndi < 3 its
 *        stress-free components are stored in STATEV after Wm (layout in abaqus2smart_M);
 *        empty leaves those slots untouched
 */
void smart2abaqus_M(double *stress, double *ddsdde, double *statev, 
    const int &ndi, const int &nshr, 
    const arma::vec &sigma, const arma::vec &statev_smart, 
    const arma::vec &Wm, const arma::mat &Lt, const bool &nlgeom = true,
    const arma::vec &Etot_end = arma::vec());

/**
 * @brief The Abaqus material Jacobian of a corotated kernel tangent.
 * @details Abaqus defines DDSDDE through the Jaumann rate of the Kirchhoff stress
 *          divided by J: \f$ \frac{1}{J}\frac{\partial (J\boldsymbol{\sigma})}{\partial \boldsymbol{\varepsilon}}
 *          = \frac{\partial \boldsymbol{\sigma}}{\partial \boldsymbol{\varepsilon}} + \boldsymbol{\sigma} \otimes \mathbf{I} \f$.
 *          The kernels return the first term in Abaqus's corotated frame; this adds the
 *          symmetric part of the second (the default symmetric solver keeps that part).
 *          Used by smart2abaqus_M and smart2abaqus_T for NLGEOM steps. Without NLGEOM the
 *          Jacobian is the kernel tangent (J = 1 by assumption), which also keeps
 *          linear-perturbation procedures (*FREQUENCY, *BUCKLE) on the right stiffness.
 * @param Lt kernel tangent (6x6, Voigt, engineering shears)
 * @param sigma Cauchy stress at the end of the increment (6)
 * @return DDSDDE as a 6x6 matrix
 */
arma::mat abaqus_jacobian(const arma::mat &Lt, const arma::vec &sigma);

/**
 * @brief Whether Abaqus runs the step with NLGEOM: DFGRD1 is the identity otherwise.
 * @param dfgrd1 the DFGRD1 array of the UMAT interface (9 values)
 */
bool abaqus_nlgeom(const double *dfgrd1);

/**
 * @brief The kinematic form an Abaqus material name asks for under NLGEOM.
 * @details The law is the first five characters of CMNAME. The total form (strain space,
 *          \f$ \ln\mathbf{U} \f$ and its dual stress, read from DFGRD0/DFGRD1, no history) is
 *          the default; the suffix `_R` (e.g. `EPICP_R`) asks for the rate form (log_R,
 *          corate 3), whose accumulated strain and frame rotation take the last
 *          abaqus_kinematic_statev() STATEV entries. `_R` must be the whole suffix (case
 *          insensitive): any other tail after the five characters (`EPICP_REF`, `EPICP_R2`)
 *          is ignored and gives the total form, as names were read before the suffix existed.
 *          A NUL-terminated name (a C caller) is read up to its terminator only.
 * @param cmname the CMNAME argument of the UMAT (space padded, 80 characters)
 * @return corate_type_total, or 3 for `_R`
 */
int abaqus_kinematic_form(const char *cmname);

/**
 * @brief STATEV entries the kinematic form keeps at the end of STATEV: 15 for the rate form
 *        (the accumulated box strain, 6, and the log_R frame rotation, 9, column-major), 0 for
 *        the total form.
 */
int abaqus_kinematic_statev(const int &corate_type);

/**
 * @brief STATEV entries left to the law behind an FE bridge (Abaqus UMAT, Ansys USERMAT).
 * @details STATEV holds, in order, the bridge's header (@p nheader: 4 for Wm, 7 for Wm and Wt),
 *          the stress-free strain of a condensed state (ndi < 3, make_condensed_state()), the
 *          law's variables, then the abaqus_kinematic_statev() entries of the kinematic form
 *          (Abaqus only). The STATEV size (Abaqus DEPVAR, Ansys TB,STATE) must cover the
 *          bridge's own entries: switching an Abaqus material to `_R` asks for 15 more. The
 *          law's own count is not known here; a size short of it shifts the trailing entries
 *          into the law's variables.
 * @param nstatev the STATEV size (NSTATV, nStatev)
 * @param nheader 4 (umat_singleM, usermat_singleM) or 7 (umat_singleT)
 * @param ndi number of direct stress components
 * @param nshr number of shear stress components
 * @param corate_type abaqus_kinematic_form(), corate_type_total for Ansys
 * @return the law's share of STATEV, \f$ \geq 0 \f$
 * @throw std::invalid_argument when the size is below the bridge's own entries
 */
int fea_law_statev(const int &nstatev, const int &nheader, const int &ndi, const int &nshr, const int &corate_type);

/**
 * @brief Whether an exception out of a material-point call is a numerical failure of the trial
 *        state, which a smaller increment may cure.
 * @details True for a std::runtime_error (a simcoon_error of a matrix function or decomposition,
 *          an Armadillo failure) other than simcoon::exception_solver, which the dispatch raises
 *          for an unsupported configuration (a refused corate or formulation); false for those and
 *          for an input error (std::logic_error: an unknown law, a short STATEV, an unsupported
 *          element layout). The FE entries (umat_singleM, umat_singleT, usermat_singleM) request a
 *          smaller increment on true and stop the analysis on false.
 */
bool fea_numerical_failure(const std::exception &e);

/**
 * @brief Whether an Abaqus call runs the finite-strain route (the deformation gradient driving
 *        select_umat_M_finite / the thermomechanical finite dispatch).
 * @details True for an NLGEOM step (abaqus_nlgeom()) with three direct components (3D, plane
 *          strain, axisymmetric) and a law the finite dispatch serves in the requested form
 *          (every thermomechanical law; the mechanical ones of finite_umat_names(), the total
 *          form only where total_formulation_serves()). Otherwise the bridge keeps the
 *          corotational small-strain route on STRAN/DSTRAN (plane stress and 1D, NLGEOM=NO,
 *          the small-strain-only laws).
 */
bool abaqus_finite_route(const std::string &umat_name, const int &ndi, const int &corate_type, const bool &thermomechanical, const double *dfgrd1);

/**
 * @brief Abaqus -> the finite-strain state of a material point, from the deformation gradients.
 * @details Called after abaqus2smart_M / abaqus2smart_T, whose STRESS read (rotated by DROT by
 *          Abaqus) it turns back into the start Kirchhoff stress in the global basis,
 *          \f$ \boldsymbol\tau_n = J_0\,\Delta\mathbf{R}^{T}\boldsymbol\sigma\,\Delta\mathbf{R} \f$.
 *          F0, F1 = DFGRD0, DFGRD1; the box strain is \f$ \ln\mathbf{U}_0 \f$ (total form) or the
 *          accumulated log_R strain and frame read from @p kin_statev (rate form; zero entries,
 *          Abaqus's initial STATEV, read as an undeformed start); DR and Detot follow from
 *          corate_kinematics, the Green-Lagrange pair from F. STRAN, DSTRAN and DROT play no
 *          further part.
 * @param sv [in/out] the global state variables of the point
 * @param drot DROT (9, column-major), the rotation Abaqus applied to STRESS
 * @param dfgrd0 DFGRD0 (9, column-major)
 * @param dfgrd1 DFGRD1 (9, column-major)
 * @param kin_statev the last abaqus_kinematic_statev() entries of STATEV (rate form), unused otherwise
 * @param corate_type corate_type_total or 3 (abaqus_kinematic_form())
 * @param DTime time increment
 */
void abaqus2smart_finite(state_variables &sv, const double *drot, const double *dfgrd0, const double *dfgrd1, const double *kin_statev, const int &corate_type, const double &DTime);

/**
 * @brief The finite-strain point -> the Abaqus stress and material Jacobian (mechanical).
 * @details STRESS is the Cauchy stress \f$ \boldsymbol\tau/J \f$ in the global basis, and
 *          DDSDDE the Jacobian Abaqus defines under NLGEOM,
 *          \f$ \frac{1}{J}\,\partial\Delta(J\boldsymbol\sigma)/\partial\Delta\boldsymbol\varepsilon \f$
 *          with the Jaumann rate of the Kirchhoff stress: the box tangent mapped to
 *          \f$ \partial\mathbf{S}/\partial\mathbf{E} \f$ (DtauDe_corate_2_DSDE) then to that rate
 *          (DSDE_2_Dsigma_JaumannDD). It replaces abaqus_jacobian, which serves the
 *          corotational route. The rate form writes its end strain and frame to @p kin_statev.
 * @param sv the global state variables after the finite dispatch
 * @param Lt the box tangent of the dispatch (state_variables_M::Lt or state_variables_T::dSdE)
 * @param corate_type corate_type_total or 3
 * @param sigma [out] Cauchy stress (6)
 * @param ddsdde [out] Abaqus Jacobian (6x6)
 * @param kin_statev [out] the last abaqus_kinematic_statev() entries of STATEV (rate form)
 */
void smart2abaqus_finite_M(const state_variables &sv, const arma::mat &Lt, const int &corate_type, arma::vec &sigma, arma::mat &ddsdde, double *kin_statev);

/**
 * @brief The finite-strain point -> the Abaqus outputs of a coupled temperature-displacement step.
 * @details smart2abaqus_finite_M for STRESS and DDSDDE, then the heat per unit CURRENT volume
 *          and time, as Abaqus integrates it over the deformed body. With \f$ r \f$ the heat
 *          of the increment per reference volume and \f$ J = \det\mathbf{F}_1 \f$:
 *          \f$ \mathrm{RPL} = r/(J\Delta t) \f$, \f$ \mathrm{DRPLDT} = (\partial r/\partial\theta)/(J\Delta t) \f$,
 *          \f$ \mathrm{DDSDDT} = (\partial\boldsymbol\tau/\partial\theta)/J \f$ (\f$ J \f$ fixed
 *          at fixed \f$ \mathbf{F} \f$) and, from \f$ \delta\mathbf{E} = \mathbf{F}^{T}\delta\boldsymbol\varepsilon\,\mathbf{F} \f$
 *          and \f$ \delta J = J\,\mathrm{tr}\,\delta\boldsymbol\varepsilon \f$,
 *          \f$ \mathrm{DRPLDE} = \left(\mathbf{F}\,\frac{\partial r}{\partial\mathbf{E}}\,\mathbf{F}^{T} - r\,\mathbf{I}\right)/(J\Delta t) \f$,
 *          \f$ \partial r/\partial\mathbf{E} \f$ from box_covector_2_dE(). All zero rates at
 *          \f$ \Delta t = 0 \f$.
 * @param sv the global state variables after the finite dispatch
 * @param corate_type corate_type_total or 3
 * @param DTime time increment
 * @param sigma [out] Cauchy stress (6)
 * @param ddsdde [out] Abaqus Jacobian (6x6)
 * @param ddsddt [out] DDSDDT (6)
 * @param rpl [out] RPL
 * @param drplde [out] DRPLDE (6)
 * @param drpldt [out] DRPLDT
 * @param kin_statev [out] the last abaqus_kinematic_statev() entries of STATEV (rate form)
 */
void smart2abaqus_finite_T(const state_variables_T &sv, const int &corate_type, const double &DTime, arma::vec &sigma, arma::mat &ddsdde, arma::vec &ddsddt, double &rpl, arma::vec &drplde, double &drpldt, double *kin_statev);

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
 * @param dSdE simcoon mechanical tangent matrix (condensed when ndi < 3, select_umat_T)
 * @param dSdT simcoon stress-temperature tangent matrix (condensed when ndi < 3)
 * @param drpldE simcoon heat generation-strain tangent matrix (condensed when ndi < 3)
 * @param drpldT simcoon heat generation-temperature tangent matrix
 * @param nlgeom as smart2abaqus_M
 * @param Etot_end as smart2abaqus_M, stored after Wm and Wt (layout in abaqus2smart_T)
 */
void smart2abaqus_T(double *stress, double *ddsdde, double *ddsddt, 
    double *drplde, double &drpldt, double &rpl, double *statev, 
    const int &ndi, const int &nshr, 
    const arma::vec &sigma, const arma::vec &statev_smart, const double &r, 
    const arma::vec &Wm, const arma::vec &Wt, 
    const arma::mat &dSdE, const arma::mat &dSdT, 
    const arma::mat &drpldE, const arma::mat &drpldT, const bool &nlgeom = true,
    const arma::vec &Etot_end = arma::vec());

//=============================================================================
// ANSYS USERMAT TRANSFER FUNCTIONS
//=============================================================================

/**
 * @brief Ansys USERMAT -> simcoon (mechanical).
 * @details Ansys passes the element's ncomp components (11, 22, 33, 12, 23, 13 in 3D: 13 and 23
 *          swapped against simcoon; 11, 22, 33, 12 in plane strain and axisymmetry; 11, 22, 12 in
 *          plane stress; 11 in 1D). USTATEV holds the four works Wm, then the stress-free strain
 *          of a condensed state (plane stress, 1D: the components select_umat_M solves, in simcoon
 *          order), then the law's variables, as the Abaqus bridge (abaqus2smart_M). rotateM is
 *          read row by row as the rotation increment DR of the state variables; an array that is
 *          not a rotation (zero, as left by an Ansys release that does not fill it) reads as the
 *          identity.
 * @param stress Ansys stress (ncomp)
 * @param dsdePl Ansys tangent on input (ncomp x ncomp, column-major), as abaqus2smart_M reads DDSDDE
 * @param strain Ansys total strain at the start of the increment (ncomp)
 * @param dstrain Ansys strain increment (ncomp)
 * @param Time time at the start of the increment
 * @param Temp temperature at the start of the increment
 * @param dTemp temperature increment
 * @param nprops number of law properties
 * @param props the law properties (the bridge's model code excluded)
 * @param nstatev USTATEV size
 * @param statev USTATEV
 * @param nDirect number of direct components
 * @param nShear number of shear components
 * @param rotateM Ansys rotateM (9)
 * @param sigma [out] stress (6)
 * @param Lt [out] the input tangent (6x6, zero outside the element's components)
 * @param Etot [out] total strain at the start (6), the stored stress-free strain included (rotated by DR)
 * @param DEtot [out] strain increment (6)
 * @param T [out] temperature
 * @param DT [out] temperature increment
 * @param props_smart [out] law properties
 * @param Wm [out] works (4)
 * @param statev_smart [out] law state variables
 * @param DR [out] rotation increment
 * @param start [out] first increment (Time = 0)
 * @throw std::invalid_argument for a component layout the bridge does not serve
 */
void ansys2smart_M(const double *stress, const double *dsdePl, const double *strain, const double *dstrain,
    const double &Time, const double &Temp, const double &dTemp,
    const int &nprops, const double *props, const int &nstatev, const double *statev,
    const int &nDirect, const int &nShear, const double *rotateM,
    arma::vec &sigma, arma::mat &Lt, arma::vec &Etot, arma::vec &DEtot, double &T, double &DT,
    arma::vec &props_smart, arma::vec &Wm, arma::vec &statev_smart, arma::mat &DR, bool &start);

/**
 * @brief simcoon -> Ansys USERMAT (mechanical).
 * @details The stress and dsdePl (ncomp x ncomp, column-major) on the element's components; the
 *          tangent arrives condensed from select_umat_M for plane stress and 1D and goes out as it
 *          is (the kernel tangent, no further term). USTATEV gets Wm, the solved stress-free strain
 *          and the law's variables (ansys2smart_M()).
 * @param stress [out] Ansys stress (ncomp)
 * @param dsdePl [out] Ansys tangent (ncomp x ncomp, column-major)
 * @param statev [out] USTATEV
 * @param nDirect number of direct components
 * @param nShear number of shear components
 * @param sigma stress (6)
 * @param Lt tangent (6x6)
 * @param Wm works (4)
 * @param statev_smart law state variables
 * @param Etot_end total strain at the end of the increment (6)
 */
void smart2ansys_M(double *stress, double *dsdePl, double *statev,
    const int &nDirect, const int &nShear,
    const arma::vec &sigma, const arma::mat &Lt,
    const arma::vec &Wm, const arma::vec &statev_smart, const arma::vec &Etot_end);


/** @} */ // end of umat_mechanical group

} // namespace simcoon
