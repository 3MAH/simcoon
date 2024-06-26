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

///@file Taba2sim.cpp
///@brief Test for the transfer from Abaqus to Simcoon Umat subroutines
///@version 1.0

#include <gtest/gtest.h>
#include <fstream>
#include <iterator>
#include <armadillo>

#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Umat/umat_smart.hpp>
#include <simcoon/Simulation/Maths/random.hpp>
#include <simcoon/Simulation/Phase/phase_characteristics.hpp>
#include <simcoon/Simulation/Phase/state_variables.hpp>
#include <simcoon/Simulation/Phase/state_variables_M.hpp>
#include <simcoon/Simulation/Phase/read.hpp>
#include <simcoon/Simulation/Solver/read.hpp>


using namespace std;
using namespace arma;
using namespace simcoon;

TEST(Taba2sim, read_write)
{

    /* initialize random seed: */
    srand(time(NULL));
    
    string path_data = "data";
    string materialfile = "material.dat";    
    string inputfile = "Nellipsoids0.dat";
    
    //double psi_rve = 0.;
    //double theta_rve = 0.;
    //double phi_rve = 0.;
    
    phase_characteristics rve;
    rve.construct(2,1);
    
    //Characterisitcs for the statev test
    double *stress = new double[6];
    double *stran = new double[6];
    double *dstran = new double[6];
    double *ddsdde = new double[36];
    double *time = new double[2];
    double dtime = 0.;
    double temperature = 273.15;
    double Dtemperature = 0.;
    unsigned int nprops = 3;
    double *props = new double[5];
    props[0] = 2;
    props[1] = 0;
    props[2] = 20;
    props[3] = 20;
    props[4] = 0;
    int ndi = 3;
    int nshr = 3;
    double *drot = new double[9];
    char *cmname = new char[5];
    cmname[0] = 'M';
    cmname[1] = 'I';
    cmname[2] = 'M';
    cmname[3] = 'T';
    cmname[4] = 'N';
    double pnewdt = 1.;
    
    bool start = false;
    double Time = 0.;
    double DTime = 0.;
    //    int solver_type = 0;
    
    dstran[0] = 0.001;
    dtime = 0.1;
    
    unsigned int nstatev_multi = 0;
    unsigned int nstatev_macro = 0;
    vec props_smart = zeros(nprops);
    for (unsigned int i=0; i<nprops; i++) {
        props_smart(i) = props[i];
    }
    int solver_type = 0;
    
    vec sigma = zeros(6);
    vec Etot = zeros(6);
    vec DEtot = zeros(6);
    mat Lt = zeros(6,6);
    mat DR = zeros(3,3);
    vec Wm = zeros(4);
    
    string umat_name(cmname);
    umat_name = umat_name.substr(0, 5);
    
    rve.sptr_matprops->update(0, umat_name, 1, 0., 0., 0., nprops, props_smart);

    read_ellipsoid(rve, path_data, inputfile);
    size_statev(rve, nstatev_multi);

    rve.sptr_matprops->umat_name = umat_name;
    
    unsigned int nstatev = nstatev_multi + nstatev_macro + 4;
    double *statev = new double[nstatev];
    for (unsigned int i=0; i<nstatev; i++) {
        statev[i] = alead(0.,100.);
    }
    
    auto rve_sv_M = std::dynamic_pointer_cast<state_variables_M>(rve.sptr_sv_global);
    rve_sv_M->resize(nstatev_macro);

    vec statev_macro = zeros(nstatev_macro);
    vec statev_multi(&statev[nstatev_macro+4], nstatev_multi, false, false);
    vec statev_multi_n = statev_multi;

    unsigned int pos = 0;
    statev_2_phases(rve, pos, statev_multi);
    
    abaqus2smart_M(stress, ddsdde, stran, dstran, time, dtime, temperature, Dtemperature, nprops, props, nstatev_macro, statev, ndi, nshr, drot, rve_sv_M->sigma, rve_sv_M->Lt, rve_sv_M->Etot, rve_sv_M->DEtot, rve_sv_M->T, rve_sv_M->DT, Time, DTime, props_smart, rve_sv_M->Wm, rve_sv_M->statev, DR, start);
    smart2abaqus_M(stress, ddsdde, statev, ndi, nshr, rve_sv_M->sigma, rve_sv_M->statev, rve_sv_M->Wm, rve_sv_M->Lt);
    
    pos = 0;
    phases_2_statev(statev_multi, pos, rve);
    
    ASSERT_TRUE(std::equal(statev_multi.begin(), statev_multi.end(), statev_multi_n.begin()));
    
    delete[] statev;
    delete[] stress;
    delete[] stran;
    delete[] dstran;
    delete[] ddsdde;
    delete[] time;
}
