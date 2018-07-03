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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "aba_2_umat"
#include <boost/test/unit_test.hpp>

#include <fstream>
#include <iterator>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Umat/umat_smart.hpp>
#include <simcoon/Simulation/Phase/phase_characteristics.hpp>
#include <simcoon/Simulation/Phase/state_variables.hpp>
#include <simcoon/Simulation/Phase/state_variables_M.hpp>
#include <simcoon/Simulation/Solver/read.hpp>

using namespace std;
using namespace arma;
using namespace simcoon;

BOOST_AUTO_TEST_CASE( read_write )
{
    
    string path_data = "data";
    string materialfile = "material.dat";    
    
    double psi_rve = 0.;
    double theta_rve = 0.;
    double phi_rve = 0.;
    
    unsigned int nstatev = 0;
    
    phase_characteristics rve;
    rve.construct(0,1);
    
    //Characterisitcs for the statecv test
    double *statev = new double[nstatev];
    double *stress = new double[6];
    double *stran = new double[6];
    double *dstran = new double[6];
    double *ddsdde = new double[36];
    double *time = new double[2];
    double dtime = 0.;
    double temperature = 273.15;
    double Dtemperature = 0.;
    unsigned int nprops = 3;
    vec props = zeros(3);
    props[0] = 70000.;
    props[1] = 0.3;
    props[2] = 1.E-5;
    int ndi = 3;
    int nshr = 3;
    double *drot = new double[9];
    char *cmname = new char[5];
    cmname[0] = 'E';
    cmname[1] = 'L';
    cmname[2] = 'I';
    cmname[3] = 'S';
    cmname[4] = 'O';
    
    bool start = false;
    double Time = 0.;
    double DTime = 0.;
    int solver_type = 0;
    
    unsigned int nstatev_smart = 0;
    vec props_smart = zeros(nprops);
    
    vec sigma = zeros(6);
    vec Etot = zeros(6);
    vec DEtot = zeros(6);
    mat Lt = zeros(6,6);
    mat DR = zeros(3,3);
    vec Wm = zeros(4);
    
    string umat_name(cmname);
    umat_name = umat_name.substr(0, 5);

    read_matprops(umat_name, nprops, props, nstatev_smart, psi_rve, theta_rve, phi_rve, path_data, materialfile);
    vec statev_smart = zeros(nstatev_smart);
    
    size_statev(rve, nstatev_smart);
    cout << "nstatev_smart = " << nstatev_smart << endl;
    
    nstatev += 4;
    cout << "nstatev = " << nstatev << endl;
    
    vec statev_smart_n = statev_smart;
//    std::copy(std::begin(statev_n), std::end(statev_n), std::begin(statev_n));
    
    auto rve_sv_M = std::dynamic_pointer_cast<state_variables_M>(rve.sptr_sv_global);
    rve_sv_M->resize(nstatev_smart);
    rve.sptr_matprops->resize(nprops);
    
    unsigned int pos = 0;
    statev_2_phases(rve, pos, statev_smart);
    
//    abaqus2smart_M(stress, ddsdde, stran, dstran, time, dtime, temperature, Dtemperature, nprops, props, nstatev, statev, ndi, nshr, drot, rve_sv_M->sigma, rve_sv_M->Lt, rve_sv_M->Etot, rve_sv_M->DEtot, rve_sv_M->T, rve_sv_M->DT, Time, DTime, props_smart, rve_sv_M->Wm, rve_sv_M->statev, DR, start);
    
//    smart2abaqus_M(stress, ddsdde, statev, ndi, nshr, rve_sv_M->sigma, rve_sv_M->statev, rve_sv_M->Wm, rve_sv_M->Lt);
    
    phases_2_statev(statev_smart, pos, rve);
    
    BOOST_CHECK_EQUAL_COLLECTIONS(statev_smart.begin(), statev_smart.end(), statev_smart_n.begin(), statev_smart_n.end());
    
    delete[] statev;
    delete[] stress;
    delete[] stran;
    delete[] dstran;
    delete[] ddsdde;
    delete[] time;
}
