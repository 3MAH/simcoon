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
#include <algorithm>
#include <vector>
#include <armadillo>

#include <simcoon/parameter.hpp>
#include <simcoon/Simulation/Phase/phase_characteristics.hpp>
#include <simcoon/Simulation/Phase/state_variables_M.hpp>
#include <simcoon/Continuum_mechanics/Umat/umat_smart.hpp>
#include <simcoon/Continuum_mechanics/Umat/umat_plugin_api.hpp>
#include <simcoon/Continuum_mechanics/Umat/fea_transfer.hpp>

using namespace std;
using namespace arma;
using namespace simcoon;

// LIB_EXPORT comes from umat_plugin_api.hpp
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

// simcoon calls a foreign Ansys USERMAT: the components are those of Ansys (fea_transfer.hpp),
// the matrices column-major as Fortran, and the law's statev are its own (no simcoon header)
class LIB_EXPORT umat_plugin_ans : public umat_plugin_ans_api {
public:

    std::string name() const {
        return "umans";
    }

    void umat_ansys(simcoon::phase_characteristics &rve, const arma::mat &DR, const double &Time, const double &DTime, const int &ndi, const int &nshr, bool &start, const int &solver_type, double &tnew_dt) {
        (void)solver_type;
        auto sv = std::dynamic_pointer_cast<state_variables_M>(rve.sptr_sv_local);
        const fea_components K = ansys_components(ndi, nshr);

        int matId = 1, elemId = 1, kDomIntPt = 1, kLayer = 1, kSectPt = 1;
        int ldstep = 1;
        int isubst = start ? 1 : 2;   // the USERMAT initialises its state on the first increment
        int keycut = 0;
        int nDirect = ndi, nShear = nshr, ncomp = K.ntens;
        int nStatev = static_cast<int>(sv->statev.n_elem);
        int nProp = static_cast<int>(rve.sptr_matprops->props.n_elem);
        std::vector<double> ustatev(sv->statev.begin(), sv->statev.end());
        std::vector<double> prop(rve.sptr_matprops->props.begin(), rve.sptr_matprops->props.end());

        if (start) sv->sigma.zeros();
        const vec sigma_start = sv->sigma;
        double stress[6], Strain[6], dStrain[6], epsPl[6] = {0.}, dsdde[36] = {0.};
        vec2fe(sv->sigma, K, stress);
        vec2fe(sv->Etot, K, Strain);
        vec2fe(sv->DEtot, K, dStrain);
        double rotateM[9];
        std::copy(DR.memptr(), DR.memptr() + 9, rotateM);   // column-major, as rotateM(3,3)
        double sedEl = sv->Wm(1), sedPl = sv->Wm(3), epseq = 0.;
        double Time_loc = Time, dTime_loc = DTime, Temp_loc = sv->T, dTemp_loc = sv->DT;
        double coords[3] = {0., 0., 0.};
        double defGrad_t[9] = {1., 0., 0., 0., 1., 0., 0., 0., 1.};
        double defGrad[9] = {1., 0., 0., 0., 1., 0., 0., 0., 1.};
        double tsstif[2] = {0., 0.};
        double epsZZ = 0.;
        int var1 = 0, var2 = 0, var3 = 0, var4 = 0, var5 = 0;

        usermat_(&matId, &elemId, &kDomIntPt, &kLayer, &kSectPt,
                 &ldstep, &isubst, &keycut,
                 &nDirect, &nShear, &ncomp, &nStatev, &nProp,
                 &Time_loc, &dTime_loc, &Temp_loc, &dTemp_loc,
                 stress, ustatev.data(), dsdde,
                 &sedEl, &sedPl, &epseq,
                 Strain, dStrain, epsPl,
                 prop.data(), coords, rotateM,
                 defGrad_t, defGrad,
                 tsstif, &epsZZ,
                 &var1, &var2, &var3, &var4, &var5);

        fe2vec(stress, K, sv->sigma);
        fe2mat(dsdde, K, sv->Lt);
        std::copy(ustatev.begin(), ustatev.end(), sv->statev.begin());
        // work by the trapezoidal rule; recoverable and dissipated parts from the USERMAT energies
        sv->Wm(0) += 0.5*sum((sigma_start + sv->sigma) % sv->DEtot);
        sv->Wm(1) = sedEl;
        sv->Wm(3) = sedPl;
        sv->Wm(2) = sv->Wm(0) - sedEl - sedPl;

        if (keycut == 1) {
            tnew_dt = 0.5;
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
