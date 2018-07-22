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
#include <boost/config.hpp> // for BOOST_SYMBOL_EXPORT

#include <simcoon/parameter.hpp>
#include <simcoon/Simulation/Phase/phase_characteristics.hpp>
#include <simcoon/Simulation/Phase/state_variables_M.hpp>
#include <simcoon/Continuum_mechanics/Umat/umat_smart.hpp>
#include <simcoon/Continuum_mechanics/Umat/umat_plugin_api.hpp>

using namespace std;
using namespace arma;
using namespace simcoon;

class umat_plugin_aba : public umat_plugin_api {
public:
    umat_plugin_aba() {};
    
    std::string name() const {
        return "umat";
    }
    
    void umat_(double *stress, double *statev, double *ddsdde, double &sse, double &spd, double &scd, double &rpl, double *ddsddt, double *drplde, double &drpldt, const double *stran, const double *dstran, const double *time, const double &dtime, const double &temperature, const double &Dtemperature, const double &predef, const double &dpred, char *cmname, const int &ndi, const int &nshr, const int &ntens, const int &nstatev, const double *props, const int &nprops, const double &coords, const double *drot, double &pnewdt, const double &celent, const double *dfgrd0, const double *dfgrd1, const int &noel, const int &npt, const double &layer, const int &kspt, const int &kstep, const int &kinc) {
    
        UNUSED(sse);
        UNUSED(spd);
        UNUSED(scd);
        UNUSED(rpl);
        UNUSED(ddsddt);
        UNUSED(drplde);
        UNUSED(drpldt);
        UNUSED(predef);
        UNUSED(dpred);
        UNUSED(ntens);
        UNUSED(coords);
        UNUSED(celent);
        UNUSED(dfgrd0);
        UNUSED(dfgrd1);
        UNUSED(noel);
        UNUSED(npt);
        UNUSED(layer);
        UNUSED(kspt);
        UNUSED(kstep);
        UNUSED(kinc);
        
        bool start = false;
        double Time = 0.;
        double DTime = 0.;
        int solver_type = 0;
        
        int nstatev_smart = nstatev-4;
        vec props_smart = zeros(nprops);
        vec statev_smart = zeros(nstatev_smart);
        
        mat DR = zeros(3,3);
        
        string umat_name(cmname);
        umat_name = umat_name.substr(0, 5);
           
        phase_characteristics rve;
        rve.construct(0,1);
        auto rve_sv_M = std::dynamic_pointer_cast<state_variables_M>(rve.sptr_sv_global);
        rve_sv_M->resize(nstatev_smart);
        rve.sptr_matprops->resize(nprops);
        
        abaqus2smart_M(stress, ddsdde, stran, dstran, time, dtime, temperature, Dtemperature, nprops, props, nstatev, statev, ndi, nshr, drot, rve_sv_M->sigma, rve_sv_M->Lt, rve_sv_M->Etot, rve_sv_M->DEtot, rve_sv_M->T, rve_sv_M->DT, Time, DTime, props_smart, rve_sv_M->Wm, rve_sv_M->statev, DR, start);
        rve.sptr_matprops->update(0, umat_name, 1., 0., 0., 0., nprops, props_smart);
        select_umat_M(rve, DR, Time, DTime, ndi, nshr, start, solver_type, pnewdt);
        smart2abaqus_M(stress, ddsdde, statev, ndi, nshr, rve_sv_M->sigma, rve_sv_M->statev, rve_sv_M->Wm, rve_sv_M->Lt);
    }
    
    virtual ~umat_plugin_aba() {};
};

// Exporting `my_namespace::plugin` variable with alias name `external_umat`
// (Has the same effect as `BOOST_DLL_ALIAS(my_namespace::plugin, external_umat)`)
extern "C" BOOST_SYMBOL_EXPORT umat_plugin_aba external_umat;
umat_plugin_aba external_umat;
