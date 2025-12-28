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

#include <iostream>
#include <fstream>
#include <assert.h>
#include <string.h>
#include <armadillo>
#include <dylib.hpp>

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

//declaration of the extern function to use
extern "C"{
	void umat_(double *stress, double *statev, double *ddsdde, double &sse, double &spd, double &scd, double &rpl, double *ddsddt, double *drplde, double &drpldt, const double *stran, const double *dstran, const double *time, const double &dtime, const double &temperature, const double &Dtemperature, const double &predef, const double &dpred, char *cmname, const int &ndi, const int &nshr, const int &ntens, const int &nstatev, const double *props, const int &nprops, const double &coords, const double *drot, double &pnewdt, const double &celent, const double *dfgrd0, const double *dfgrd1, const int &noel, const int &npt, const double &layer, const int &kspt, const int &kstep, const int &kinc);
	}

class LIB_EXPORT umat_plugin_aba : public umat_plugin_aba_api {
public:

    std::string name() const {
        return "umaba";
    }
    
    void umat_abaqus(simcoon::phase_characteristics &rve, const arma::mat &DR, const double &Time, const double &DTime, const int &ndi, const int &nshr, bool &start, const int &solver_type, double &tnew_dt) {
        
        ///@brief Macroscopic state variables and control increments
        double stress[6];
        double stran[6];
        double dstran[6];
        double temperature;
        double Dtemperature;
        
        ///@brief Umat variable list unused here
        double sse = 0.;
        double spd = 0.;
        double scd = 0.;
        double rpl = 0.;
        double drpldt = 0.;
        double predef = 0.;
        double dpred = 0.;
        int layer = 0;
        int kspt = 0;
        double celent = 0.;
        double dfgrd0[9];
        double dfgrd1[3];
        double drplde[6];
        double coords = 0;
        double drot[9];
            
        ///@brief Usefull UMAT variables
        int ntens = 6;
        int noel = 1;
        int npt = 1;
        int kstep = 0;
        int kinc = 0;
        double ddsdde[36];
        double ddsddt[6];
        char cmname[6];
        strcpy(cmname, rve.sptr_matprops->umat_name.c_str());
        int nprops = rve.sptr_matprops->nprops;
        double props[nprops];
        int nstatev = rve.sptr_sv_global->nstatev;
        double statev[nstatev+4];                   //+4 for a mechanical response to store Wm components
        double pnewdt = tnew_dt;
        double time[2];
        double dtime;
        
        auto rve_sv_M = std::dynamic_pointer_cast<state_variables_M>(rve.sptr_sv_local);
        
        smart2abaqus_M_full(stress, ddsdde, stran, dstran, time, dtime, temperature, Dtemperature, nprops, props, nstatev, statev, ndi, nshr, drot, rve_sv_M->sigma, rve_sv_M->Lt, rve_sv_M->Etot, rve_sv_M->DEtot, rve_sv_M->T, rve_sv_M->DT, Time, DTime, rve.sptr_matprops->props, rve_sv_M->Wm, rve_sv_M->statev, DR, start);
//        smart2abaqus_M(stress, ddsdde, statev, ndi, nshr, rve_sv_M->sigma, rve_sv_M->statev, rve_sv_M->Wm, rve_sv_M->Lt);
                
        umat_(stress, statev, ddsdde, sse, spd, scd, rpl, ddsddt, drplde, drpldt, stran, dstran, time, dtime, temperature, Dtemperature, predef, dpred, cmname, ndi, nshr, ntens, nstatev, props, nprops, coords, drot, pnewdt, celent, dfgrd0, dfgrd1, noel, npt, layer, kspt, kstep, kinc);
        
        abaqus2smart_M_light(stress, ddsdde, nstatev, statev, ndi, nshr, rve_sv_M->sigma, rve_sv_M->Lt, rve_sv_M->Wm, rve_sv_M->statev);
    }

   ~umat_plugin_aba() {}
};

// Export create/destroy functions using dylib macros
extern "C" LIB_EXPORT umat_plugin_aba_api* create_api() {
    return new umat_plugin_aba();
}

extern "C" LIB_EXPORT void destroy_api(umat_plugin_aba_api* p) {
    delete p;
}