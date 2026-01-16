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

///@file umat_plugin_ans.cpp
///@brief Ansys USERMAT compatibility plugin for simcoon
///@author Y. Chemisky
///@version 1.0
///@date 2026-01-13

#include <iostream>
#include <fstream>
#include <assert.h>
#include <string.h>
#include <armadillo>

#include <simcoon/parameter.hpp>
#include <simcoon/Simulation/Phase/phase_characteristics.hpp>
#include <simcoon/Simulation/Phase/state_variables_M.hpp>
#include <simcoon/Continuum_mechanics/Umat/umat_smart.hpp>
#include <simcoon/Continuum_mechanics/Umat/umat_plugin_api.hpp>

using namespace std;
using namespace arma;
using namespace simcoon;

#if defined(_WIN32) || defined(_WIN64)
    #define LIB_EXPORT __declspec(dllexport)
#elif defined(__GNUC__) || defined(__clang__)
    #if __GNUC__ >= 4
        #define LIB_EXPORT __attribute__((visibility("default")))
    #else
        #define LIB_EXPORT
    #endif
#else
    #define LIB_EXPORT
#endif

///@brief Declaration of the external Ansys USERMAT subroutine
/// This is the standard Ansys USERMAT interface (Fortran)
/// Reference: Ansys Mechanical APDL Programmer's Reference
extern "C" {
    void usermat_(
        int *matId,         // Material ID number
        int *elemId,        // Element number
        int *kDomIntPt,     // Current integration point in element
        int *kLayer,        // Current layer number
        int *kSectPt,       // Current section point within layer
        int *ldstep,        // Current load step
        int *isubst,        // Current substep
        int *keycut,        // Cutback flag (output: set to 1 to request cutback)
        int *nDirect,       // Number of direct stress components
        int *nShear,        // Number of shear stress components
        int *ncomp,         // Number of stress/strain components (nDirect + nShear)
        int *nStatev,       // Number of state variables
        int *nProp,         // Number of material properties
        double *Time,       // Current time
        double *dTime,      // Time increment
        double *Temp,       // Current temperature
        double *dTemp,      // Temperature increment
        double *stress,     // Stress tensor (in: at t, out: at t+dt)
        double *ustatev,    // State variables (in/out)
        double *dsdde,      // Material tangent modulus (output)
        double *sedEl,      // Elastic strain energy density (in/out)
        double *sedPl,      // Plastic strain energy density (in/out)
        double *epseq,      // Equivalent plastic strain (in/out)
        double *Strain,     // Total strain at t
        double *dStrain,    // Strain increment
        double *epsPl,      // Plastic strain components
        double *prop,       // Material property array
        double *coords,     // Integration point coordinates (x, y, z)
        double *rotateM,    // Rotation matrix (3x3)
        double *defGrad_t,  // Deformation gradient at t (3x3)
        double *defGrad,    // Deformation gradient at t+dt (3x3)
        double *tsstif,     // Transverse shear stiffness (shells)
        double *epsZZ,      // Out-of-plane strain for plane stress
        int *var1,          // Reserved
        int *var2,          // Reserved
        int *var3,          // Reserved
        int *var4,          // Reserved
        int *var5           // Reserved
    );
}

///@brief Convert simcoon state to Ansys USERMAT format
///@details Transforms simcoon Armadillo arrays to C-style arrays for Ansys interface
void smart2ansys_M(
    double *stress, double *dsdde, double *ustatev,
    double *Strain, double *dStrain, double *Time, double *dTime,
    double *Temp, double *dTemp, double *prop, double *rotateM,
    int &nDirect, int &nShear, int &ncomp, int &nStatev, int &nProp,
    const vec &sigma_smart, const mat &Lt_smart,
    const vec &Etot_smart, const vec &DEtot_smart,
    const double &T_smart, const double &DT_smart,
    const double &Time_smart, const double &DTime_smart,
    const vec &props_smart, const vec &Wm_smart, const vec &statev_smart,
    const mat &DR_smart, const bool &start
) {
    // Dimensions
    nDirect = 3;
    nShear = 3;
    ncomp = 6;
    nStatev = statev_smart.n_elem;
    nProp = props_smart.n_elem;
    
    // Time
    *Time = Time_smart;
    *dTime = DTime_smart;
    
    // Temperature
    *Temp = T_smart;
    *dTemp = DT_smart;
    
    // Material properties
    for (int i = 0; i < nProp; i++) {
        prop[i] = props_smart(i);
    }
    
    // Ansys uses Voigt notation: 11, 22, 33, 12, 23, 13
    // simcoon uses: 11, 22, 33, 12, 13, 23
    // Need to swap indices 4<->5 (13 and 23)
    int ansys_map[6] = {0, 1, 2, 3, 5, 4};
    
    // Stress (initialize if start)
    if (start) {
        for (int i = 0; i < 6; i++) {
            stress[i] = 0.0;
        }
    } else {
        for (int i = 0; i < 6; i++) {
            stress[i] = sigma_smart(ansys_map[i]);
        }
    }
    
    // Strain and strain increment
    for (int i = 0; i < 6; i++) {
        Strain[i] = Etot_smart(ansys_map[i]);
        dStrain[i] = DEtot_smart(ansys_map[i]);
    }
    
    // Rotation matrix (row-major for Ansys)
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            rotateM[i * 3 + j] = DR_smart(i, j);
        }
    }
    
    // State variables (first 4 are Wm components, then user statev)
    ustatev[0] = Wm_smart(0);  // Wm
    ustatev[1] = Wm_smart(1);  // Wm_r
    ustatev[2] = Wm_smart(2);  // Wm_ir
    ustatev[3] = Wm_smart(3);  // Wm_d
    for (int i = 0; i < nStatev; i++) {
        ustatev[i + 4] = statev_smart(i);
    }
    
    // Initialize tangent modulus
    for (int i = 0; i < 36; i++) {
        dsdde[i] = 0.0;
    }
}

///@brief Convert Ansys USERMAT output back to simcoon format
void ansys2smart_M(
    const double *stress, const double *dsdde,
    const int &nStatev, const double *ustatev,
    const int &nDirect, const int &nShear,
    vec &sigma_smart, mat &Lt_smart, vec &Wm_smart, vec &statev_smart
) {
    // Ansys to simcoon index mapping
    int ansys_map[6] = {0, 1, 2, 3, 5, 4};
    
    // Stress
    for (int i = 0; i < 6; i++) {
        sigma_smart(ansys_map[i]) = stress[i];
    }
    
    // Tangent modulus with index reordering
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            Lt_smart(ansys_map[i], ansys_map[j]) = dsdde[i * 6 + j];
        }
    }
    
    // State variables (extract Wm and user statev)
    Wm_smart(0) = ustatev[0];  // Wm
    Wm_smart(1) = ustatev[1];  // Wm_r
    Wm_smart(2) = ustatev[2];  // Wm_ir
    Wm_smart(3) = ustatev[3];  // Wm_d
    for (int i = 0; i < nStatev; i++) {
        statev_smart(i) = ustatev[i + 4];
    }
}

class LIB_EXPORT umat_plugin_ans : public umat_plugin_ans_api {
public:

    std::string name() const {
        return "umans";
    }
    
    void umat_ansys(simcoon::phase_characteristics &rve, const arma::mat &DR, const double &Time, const double &DTime, const int &ndi, const int &nshr, bool &start, const int &solver_type, double &tnew_dt) {
        
        // Get state variables pointer
        auto rve_sv_M = std::dynamic_pointer_cast<state_variables_M>(rve.sptr_sv_local);
        
        // Ansys USERMAT variables
        int matId = 1;
        int elemId = 1;
        int kDomIntPt = 1;
        int kLayer = 1;
        int kSectPt = 1;
        int ldstep = 1;
        int isubst = 1;
        int keycut = 0;
        int nDirect = 3;
        int nShear = 3;
        int ncomp = 6;
        int nStatev = rve.sptr_sv_global->nstatev;
        int nProp = rve.sptr_matprops->nprops;
        
        // Allocate arrays
        double stress[6];
        double ustatev[nStatev + 4];  // +4 for Wm components
        double dsdde[36];
        double sedEl = 0.0;
        double sedPl = 0.0;
        double epseq = 0.0;
        double Strain[6];
        double dStrain[6];
        double epsPl[6] = {0.0};
        double prop[nProp];
        double coords[3] = {0.0, 0.0, 0.0};
        double rotateM[9];
        double defGrad_t[9] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
        double defGrad[9] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
        double tsstif[2] = {0.0, 0.0};
        double epsZZ = 0.0;
        double dTime_loc = DTime;
        double Time_loc = Time;
        double Temp_loc = rve_sv_M->T;
        double dTemp_loc = rve_sv_M->DT;
        int var1 = 0, var2 = 0, var3 = 0, var4 = 0, var5 = 0;
        
        // Convert simcoon to Ansys format
        smart2ansys_M(stress, dsdde, ustatev, Strain, dStrain, &Time_loc, &dTime_loc,
                      &Temp_loc, &dTemp_loc, prop, rotateM,
                      nDirect, nShear, ncomp, nStatev, nProp,
                      rve_sv_M->sigma, rve_sv_M->Lt,
                      rve_sv_M->Etot, rve_sv_M->DEtot,
                      rve_sv_M->T, rve_sv_M->DT,
                      Time, DTime,
                      rve.sptr_matprops->props, rve_sv_M->Wm, rve_sv_M->statev,
                      DR, start);
        
        // Call Ansys USERMAT
        usermat_(&matId, &elemId, &kDomIntPt, &kLayer, &kSectPt,
                 &ldstep, &isubst, &keycut,
                 &nDirect, &nShear, &ncomp, &nStatev, &nProp,
                 &Time_loc, &dTime_loc, &Temp_loc, &dTemp_loc,
                 stress, ustatev, dsdde,
                 &sedEl, &sedPl, &epseq,
                 Strain, dStrain, epsPl,
                 prop, coords, rotateM,
                 defGrad_t, defGrad,
                 tsstif, &epsZZ,
                 &var1, &var2, &var3, &var4, &var5);
        
        // Convert Ansys output back to simcoon format
        ansys2smart_M(stress, dsdde, nStatev, ustatev, nDirect, nShear,
                      rve_sv_M->sigma, rve_sv_M->Lt, rve_sv_M->Wm, rve_sv_M->statev);
        
        // Handle cutback request
        if (keycut == 1) {
            tnew_dt = 0.5;  // Request smaller time step
        }
    }

    ~umat_plugin_ans() {}
};

// Export create/destroy functions for dynamic loading
extern "C" LIB_EXPORT umat_plugin_ans_api* create_api() {
    return new umat_plugin_ans();
}

extern "C" LIB_EXPORT void destroy_api(umat_plugin_ans_api* p) {
    delete p;
}
