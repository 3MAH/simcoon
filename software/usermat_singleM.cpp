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

///@file usermat_singleM.cpp
///@brief USERMAT template to run simcoon subroutines using Ansys MAPDL
///@brief Implemented in 1D-2D-3D
///@author Chemisky
///@version 1.0
///@date 01/2026

///@details This file provides a bridge between Ansys MAPDL's USERMAT interface
///         and simcoon's constitutive models. The model is selected by the
///         first material property (props[0]) which should be set to the
///         model code (e.g., 2 for ELISO, 5 for EPICP, etc.).
///
///         Voigt notation difference:
///         - simcoon/Abaqus: (11, 22, 33, 12, 13, 23) -> indices (0,1,2,3,4,5)
///         - Ansys:          (11, 22, 33, 12, 23, 13) -> indices (0,1,2,3,4,5)
///         Indices 4 and 5 are swapped between conventions.
///
///         Key limitations of Ansys USERMAT (from MFront/TFEL documentation):
///         - External state variables are not supported
///         - Internal state variables cannot be initialized to non-zero values
///         - No standard way of defining orthotropic axes
///         - No way of controlling time step increase/decrease
///
///         Material properties in Ansys:
///         - props[0]: Model code (see select_umat_M for available models)
///         - props[1:nProp-1]: Model-specific material properties
///
///         State variables layout (same as Abaqus):
///         - ustatev[0:3]: Work quantities (Wm, Wm_r, Wm_ir, Wm_d)
///         - ustatev[4:nStatev-1]: Model-specific state variables

#include <iostream>
#include <fstream>
#include <assert.h>
#include <string.h>
#include <armadillo>
#include <map>

#include <simcoon/parameter.hpp>
#include <simcoon/Simulation/Phase/phase_characteristics.hpp>
#include <simcoon/Simulation/Phase/state_variables_M.hpp>
#include <simcoon/Continuum_mechanics/Umat/umat_smart.hpp>

using namespace std;
using namespace arma;
using namespace simcoon;

///@brief Map from model code to 5-character UMAT name
static const std::map<int, std::string> model_code_to_name = {
    {0, "UMEXT"},   // External plugin
    {1, "UMABA"},   // Abaqus plugin
    {2, "ELISO"},   // Isotropic elasticity
    {3, "ELIST"},   // Transversely isotropic elasticity
    {4, "ELORT"},   // Orthotropic elasticity
    {5, "EPICP"},   // Isotropic plasticity (isotropic hardening)
    {6, "EPKCP"},   // Kinematic + isotropic hardening
    {7, "EPCHA"},   // Chaboche cyclic plasticity
    {8, "SMAUT"},   // SMA unified model
    {9, "SMANI"},   // SMA anisotropic model
    {10, "LLDM0"},  // Lemaitre-Chaboche damage
    {11, "ZENER"},  // Zener viscoelastic (single branch)
    {12, "ZENNK"},  // Zener viscoelastic (N branches)
    {13, "PRONK"},  // Prony series viscoelastic
    {17, "EPHIL"},  // Hill anisotropic plasticity (iso hardening)
    {18, "EPHAC"},  // Hill + Chaboche
    {19, "EPANI"},  // Anisotropic plasticity
    {20, "EPDFA"},  // 
    {21, "EPCHG"},  // 
    {22, "EPHIN"},  // 
    {23, "SMAMO"},  // SMA Morin
    {24, "SMAMC"},  // SMA Morin Casciato
    {100, "MIHEN"}, // Mori-Tanaka (Eshelby)
    {101, "MIMTN"}, // Mori-Tanaka N phases
    {103, "MISCN"}, // Self-consistent N phases
    {104, "MIPLN"}  // Periodic layered
};

///@brief Voigt index mapping arrays
/// Ansys uses (11,22,33,12,23,13) while simcoon uses (11,22,33,12,13,23)
/// So for 3D: ansys index 4 (23) -> simcoon index 5
///            ansys index 5 (13) -> simcoon index 4
static const int ansys2smart_idx_3D[6] = {0, 1, 2, 3, 5, 4};
static const int smart2ansys_idx_3D[6] = {0, 1, 2, 3, 5, 4};

///@brief Convert Ansys arrays to simcoon format
///@param stress Ansys stress array (input)
///@param dsdePl Ansys tangent modulus array (input, column-major)
///@param strain Ansys total strain array (input)
///@param dstrain Ansys strain increment array (input)
///@param Time Current time
///@param dTime Time increment
///@param Temp Temperature at start of increment
///@param dTemp Temperature increment
///@param nProp Number of material properties
///@param prop Material properties array
///@param nStatev Number of state variables
///@param ustatev State variables array
///@param nDirect Number of direct stress components
///@param nShear Number of shear stress components
///@param rotateM Rotation matrix (3x3, row-major)
///@param sigma_out simcoon stress vector (output)
///@param Lt_out simcoon tangent modulus matrix (output)
///@param Etot_out simcoon total strain vector (output)
///@param DEtot_out simcoon strain increment vector (output)
///@param T_out Temperature (output)
///@param DT_out Temperature increment (output)
///@param Time_out Time (output)
///@param DTime_out Time increment (output)
///@param props_smart Material properties vector (output)
///@param Wm_out Work quantities vector (output)
///@param statev_smart State variables vector (output)
///@param DR_out Rotation matrix (output)
///@param start First increment flag (output)
void ansys2smart_M(
    const double *stress, const double *dsdePl,
    const double *strain, const double *dstrain,
    const double &Time, const double &dTime,
    const double &Temp, const double &dTemp,
    const int &nProp, const double *prop,
    const int &nStatev, double *ustatev,
    const int &nDirect, const int &nShear,
    const double *rotateM,
    vec &sigma_out, mat &Lt_out,
    vec &Etot_out, vec &DEtot_out,
    double &T_out, double &DT_out,
    double &Time_out, double &DTime_out,
    vec &props_smart, vec &Wm_out, vec &statev_smart,
    mat &DR_out, bool &start)
{
    int ncomp = nDirect + nShear;
    
    // Handle different dimensionalities
    if (nDirect == 1) {
        // 1D case
        sigma_out(0) = stress[0];
        Etot_out(0) = strain[0];
        DEtot_out(0) = dstrain[0];
        Lt_out(0,0) = dsdePl[0];
    }
    else if (nDirect == 2) {
        // 2D Plane Stress (nDirect=2, nShear=1, ncomp=3)
        sigma_out(0) = stress[0];
        sigma_out(1) = stress[1];
        sigma_out(3) = stress[2];  // shear 12
        
        Etot_out(0) = strain[0];
        Etot_out(1) = strain[1];
        Etot_out(3) = strain[2];
        
        DEtot_out(0) = dstrain[0];
        DEtot_out(1) = dstrain[1];
        DEtot_out(3) = dstrain[2];
        
        // Tangent modulus (column-major in Ansys)
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                Lt_out(j, i) = dsdePl[i * 3 + j];
            }
        }
    }
    else if (nDirect == 3) {
        if (nShear == 1) {
            // 2D Generalized Plane Strain or Axisymmetric (nDirect=3, nShear=1, ncomp=4)
            sigma_out(0) = stress[0];
            sigma_out(1) = stress[1];
            sigma_out(2) = stress[2];
            sigma_out(3) = stress[3];  // shear 12
            
            Etot_out(0) = strain[0];
            Etot_out(1) = strain[1];
            Etot_out(2) = strain[2];
            Etot_out(3) = strain[3];
            
            DEtot_out(0) = dstrain[0];
            DEtot_out(1) = dstrain[1];
            DEtot_out(2) = dstrain[2];
            DEtot_out(3) = dstrain[3];
            
            // Tangent modulus (column-major, 4x4)
            Lt_out(0,0) = dsdePl[0];
            Lt_out(0,1) = dsdePl[4];
            Lt_out(0,2) = dsdePl[8];
            Lt_out(0,3) = dsdePl[12];
            Lt_out(1,0) = dsdePl[1];
            Lt_out(1,1) = dsdePl[5];
            Lt_out(1,2) = dsdePl[9];
            Lt_out(1,3) = dsdePl[13];
            Lt_out(2,0) = dsdePl[2];
            Lt_out(2,1) = dsdePl[6];
            Lt_out(2,2) = dsdePl[10];
            Lt_out(2,3) = dsdePl[14];
            Lt_out(3,0) = dsdePl[3];
            Lt_out(3,1) = dsdePl[7];
            Lt_out(3,2) = dsdePl[11];
            Lt_out(3,3) = dsdePl[15];
        }
        else {
            // 3D case (nDirect=3, nShear=3, ncomp=6)
            // Apply Voigt index remapping for shear components
            for (int i = 0; i < 6; i++) {
                int i_smart = ansys2smart_idx_3D[i];
                sigma_out(i_smart) = stress[i];
                Etot_out(i_smart) = strain[i];
                DEtot_out(i_smart) = dstrain[i];
            }
            
            // Tangent modulus with index remapping (column-major)
            for (int i = 0; i < 6; i++) {
                for (int j = 0; j < 6; j++) {
                    int i_smart = ansys2smart_idx_3D[i];
                    int j_smart = ansys2smart_idx_3D[j];
                    Lt_out(i_smart, j_smart) = dsdePl[j * 6 + i];
                }
            }
        }
    }
    
    // Rotation matrix (row-major in Ansys, stored as 3x3)
    DR_out(0,0) = rotateM[0];
    DR_out(0,1) = rotateM[1];
    DR_out(0,2) = rotateM[2];
    DR_out(1,0) = rotateM[3];
    DR_out(1,1) = rotateM[4];
    DR_out(1,2) = rotateM[5];
    DR_out(2,0) = rotateM[6];
    DR_out(2,1) = rotateM[7];
    DR_out(2,2) = rotateM[8];
    
    // Temperature
    T_out = Temp;
    DT_out = dTemp;
    
    // Time
    Time_out = Time;
    DTime_out = dTime;
    
    // First increment detection
    if (Time < 1E-12) {
        start = true;
    }
    else {
        start = false;
    }
    
    // Material properties (skip first one which is model code)
    for (int i = 1; i < nProp; i++) {
        props_smart(i-1) = prop[i];
    }
    
    // State variables: first 4 are work quantities
    for (int i = 0; i < 4; i++) {
        Wm_out(i) = ustatev[i];
    }
    // Remaining are model-specific state variables
    for (int i = 0; i < nStatev - 4; i++) {
        statev_smart(i) = ustatev[i + 4];
    }
}

///@brief Convert simcoon results to Ansys format
///@param stress Ansys stress array (output)
///@param dsdePl Ansys tangent modulus array (output, column-major)
///@param ustatev Ansys state variables array (output)
///@param nDirect Number of direct stress components
///@param nShear Number of shear stress components
///@param sigma simcoon stress vector (input)
///@param statev_smart simcoon state variables vector (input)
///@param Wm simcoon work quantities vector (input)
///@param Lt simcoon tangent modulus matrix (input)
void smart2ansys_M(
    double *stress, double *dsdePl, double *ustatev,
    const int &nDirect, const int &nShear,
    const vec &sigma, const vec &statev_smart,
    const vec &Wm, const mat &Lt)
{
    int ncomp = nDirect + nShear;
    
    if (nDirect == 1) {
        // 1D case
        stress[0] = sigma(0);
        dsdePl[0] = Lt(0,0);
    }
    else if (nDirect == 2) {
        // 2D Plane Stress
        stress[0] = sigma(0);
        stress[1] = sigma(1);
        stress[2] = sigma(3);  // shear 12
        
        // Plane stress condensation for tangent modulus
        assert(Lt(2,2) > 0.);
        dsdePl[0] = Lt(0,0) - Lt(0,2)*Lt(2,0)/Lt(2,2);
        dsdePl[4] = Lt(1,1) - Lt(1,2)*Lt(2,1)/Lt(2,2);
        dsdePl[3] = Lt(0,1) - Lt(0,2)*Lt(2,1)/Lt(2,2);
        dsdePl[1] = Lt(1,0) - Lt(1,2)*Lt(2,0)/Lt(2,2);
        dsdePl[6] = Lt(0,3) - Lt(0,2)*Lt(2,3)/Lt(2,2);
        dsdePl[7] = Lt(1,3) - Lt(1,2)*Lt(2,3)/Lt(2,2);
        dsdePl[2] = Lt(3,0) - Lt(3,2)*Lt(2,0)/Lt(2,2);
        dsdePl[5] = Lt(3,1) - Lt(3,2)*Lt(2,1)/Lt(2,2);
        dsdePl[8] = Lt(3,3) - Lt(3,2)*Lt(2,3)/Lt(2,2);
    }
    else if (nDirect == 3) {
        if (nShear == 1) {
            // 2D Generalized Plane Strain or Axisymmetric
            stress[0] = sigma(0);
            stress[1] = sigma(1);
            stress[2] = sigma(2);
            stress[3] = sigma(3);
            
            // Tangent modulus (column-major, 4x4)
            dsdePl[0] = Lt(0,0);
            dsdePl[4] = Lt(0,1);
            dsdePl[8] = Lt(0,2);
            dsdePl[12] = Lt(0,3);
            dsdePl[1] = Lt(1,0);
            dsdePl[5] = Lt(1,1);
            dsdePl[9] = Lt(1,2);
            dsdePl[13] = Lt(1,3);
            dsdePl[2] = Lt(2,0);
            dsdePl[6] = Lt(2,1);
            dsdePl[10] = Lt(2,2);
            dsdePl[14] = Lt(2,3);
            dsdePl[3] = Lt(3,0);
            dsdePl[7] = Lt(3,1);
            dsdePl[11] = Lt(3,2);
            dsdePl[15] = Lt(3,3);
        }
        else {
            // 3D case - apply Voigt index remapping
            for (int i = 0; i < 6; i++) {
                int i_smart = smart2ansys_idx_3D[i];
                stress[i] = sigma(i_smart);
            }
            
            // Tangent modulus with index remapping (column-major)
            for (int i = 0; i < 6; i++) {
                for (int j = 0; j < 6; j++) {
                    int i_smart = smart2ansys_idx_3D[i];
                    int j_smart = smart2ansys_idx_3D[j];
                    dsdePl[j * 6 + i] = Lt(i_smart, j_smart);
                }
            }
        }
    }
    
    // State variables: first 4 are work quantities
    for (int i = 0; i < 4; i++) {
        ustatev[i] = Wm(i);
    }
    // Remaining are model-specific state variables
    for (unsigned int i = 0; i < statev_smart.n_elem; i++) {
        ustatev[i + 4] = statev_smart(i);
    }
}

///@brief Ansys USERMAT interface
///
///@param matId Material ID number
///@param elemId Element number
///@param kDomIntPt Domain integration point number
///@param kLayer Layer number
///@param kSectPt Section point number
///@param ldStep Load step number
///@param iSubst Substep number
///@param keycut Cutback flag (output: set to 1 to request cutback - not supported)
///@param nDirect Number of direct stress components (1, 2, or 3)
///@param nShear Number of shear stress components (0, 1, or 3)
///@param ncomp Total stress components (nDirect + nShear)
///@param nStatev Number of state variables
///@param nProp Number of material properties (first one is model code)
///@param nPt Total number of integration points
///@param Time Current time
///@param dTime Time increment
///@param Temp Temperature at start of increment
///@param dTemp Temperature increment
///@param stress Stress array (input/output, dimension ncomp)
///@param ustatev State variables array (input/output, dimension nStatev)
///@param dsdePl Tangent modulus array (output, dimension ncomp*ncomp, column-major)
///@param sedEl Elastic strain energy density (output)
///@param sedPl Plastic strain energy density (output)
///@param epseq Equivalent plastic strain (output)
///@param strain Total strain array (input, dimension ncomp)
///@param dstrain Strain increment array (input, dimension ncomp)
///@param epsPl Plastic strain array (input, dimension ncomp)
///@param prop Material properties array (input, dimension nProp)
///@param coords Integration point coordinates (input, dimension 3)
///@param rotateM Rotation matrix (input, dimension 3x3, row-major)
///@param defGrad_t Deformation gradient at time t (input, dimension 3x3)
///@param defGrad Deformation gradient at time t+dt (input, dimension 3x3)
///@param tsstif Transverse shear stiffness (for shells)

extern "C" void usermat_(
    int *matId, int *elemId, int *kDomIntPt, int *kLayer, int *kSectPt,
    int *ldStep, int *iSubst, int *keycut,
    int *nDirect, int *nShear, int *ncomp, int *nStatev, int *nProp,
    int *nPt, double *Time, double *dTime, double *Temp, double *dTemp,
    double *stress, double *ustatev, double *dsdePl,
    double *sedEl, double *sedPl, double *epseq,
    double *strain, double *dstrain, double *epsPl,
    double *prop, double *coords, double *rotateM,
    double *defGrad_t, double *defGrad, double *tsstif)
{
    // Mark unused parameters
    UNUSED(matId);
    UNUSED(elemId);
    UNUSED(kDomIntPt);
    UNUSED(kLayer);
    UNUSED(kSectPt);
    UNUSED(ldStep);
    UNUSED(iSubst);
    UNUSED(nPt);
    UNUSED(sedEl);
    UNUSED(sedPl);
    UNUSED(epseq);
    UNUSED(epsPl);
    UNUSED(coords);
    UNUSED(defGrad_t);
    UNUSED(defGrad);
    UNUSED(tsstif);
    
    // keycut: set to 1 to request time step cutback (not supported in simcoon bridge)
    *keycut = 0;
    
    bool start = false;
    double Time_smart = 0.;
    double DTime_smart = 0.;
    int solver_type = 0;
    double pnewdt = 1.0;  // Not used by Ansys but required by select_umat_M
    
    // Get model code from first material property
    int model_code = static_cast<int>(prop[0]);
    
    // Look up model name
    auto it = model_code_to_name.find(model_code);
    if (it == model_code_to_name.end()) {
        std::cerr << "ERROR: Unknown model code " << model_code << " in USERMAT" << std::endl;
        *keycut = 1;  // Request cutback (will likely cause convergence failure)
        return;
    }
    std::string umat_name = it->second;
    
    // Number of model-specific state variables (total - 4 work quantities)
    int nstatev_smart = *nStatev - 4;
    
    // Number of model-specific properties (total - 1 for model code)
    int nprops_smart = *nProp - 1;
    
    // Allocate simcoon data structures
    vec props_smart = zeros(nprops_smart);
    vec statev_smart = zeros(nstatev_smart);
    vec sigma = zeros(6);
    vec Etot = zeros(6);
    vec DEtot = zeros(6);
    mat Lt = zeros(6, 6);
    mat DR = zeros(3, 3);
    vec Wm = zeros(4);
    
    // Set up phase characteristics
    phase_characteristics rve;
    rve.construct(0, 1);
    auto rve_sv_M = std::dynamic_pointer_cast<state_variables_M>(rve.sptr_sv_global);
    rve_sv_M->resize(nstatev_smart);
    rve.sptr_matprops->resize(nprops_smart);
    
    // Convert Ansys data to simcoon format
    ansys2smart_M(
        stress, dsdePl, strain, dstrain,
        *Time, *dTime, *Temp, *dTemp,
        *nProp, prop, *nStatev, ustatev,
        *nDirect, *nShear, rotateM,
        rve_sv_M->sigma, rve_sv_M->Lt,
        rve_sv_M->Etot, rve_sv_M->DEtot,
        rve_sv_M->T, rve_sv_M->DT,
        Time_smart, DTime_smart,
        props_smart, rve_sv_M->Wm, rve_sv_M->statev,
        DR, start
    );
    
    // Update material properties in the phase
    rve.sptr_matprops->update(0, umat_name, 1., 0., 0., 0., nprops_smart, props_smart);
    
    // Call the simcoon constitutive model
    select_umat_M(rve, DR, Time_smart, DTime_smart, *nDirect, *nShear, start, solver_type, pnewdt);
    
    // Convert results back to Ansys format
    smart2ansys_M(
        stress, dsdePl, ustatev,
        *nDirect, *nShear,
        rve_sv_M->sigma, rve_sv_M->statev, rve_sv_M->Wm, rve_sv_M->Lt
    );
    
    // Note: pnewdt from simcoon is ignored since Ansys USERMAT cannot control time step
    // If simcoon requested a smaller time step (pnewdt < 1), we cannot honor it
    // The user should ensure sufficiently small time increments are used
}
