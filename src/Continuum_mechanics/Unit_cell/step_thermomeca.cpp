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

///@file steps.cpp
///@brief Characteristics of a step, as an Abaqus input
///@version 1.0

#include <iostream>
#include <sstream>
#include <fstream>
#include <assert.h>
#include <math.h>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Simulation/Solver/step.hpp>
#include <simcoon/Simulation/Solver/step_thermomeca.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/steps_thermomeca.hpp>

using namespace std;
using namespace arma;

namespace simcoon{
    
//=====Private methods for ellipsoid_characteristics===================================

//=====Public methods for ellipsoid_characteristics============================================

//@brief default constructor
//-------------------------------------------------------------
aba_step_thermomeca::aba_step_thermomeca() : step_thermomeca()
//-------------------------------------------------------------
{
    name = "step_0";    //name of the step
    nlgeom = false;
    type = 0; //0 : Automatic or  1: fixed
    max_temp = 10.; //Max temperature per increment
}

/*!
 \brief Constructor with parameters
 name = mname;          Name of the step
 nlgeom = mnlgeom;      Bool for geometric nonlinearities
 BC_T = mBC_T;          target temperature
 Etot = mEtot;          current strain
 sigma = msigma;        current stress
 T = mT;                current temperature
 */    
    
//-------------------------------------------------------------
aba_step_thermomeca::aba_step_thermomeca(const std::string &mname, const bool &mnlgeom, const int &mtype, const double &mmax_temp, const int &mnumber, const double &mDn_init, const double &mDn_mini, const double &mDn_inc, const int &mmode, const Col<int> &mcBC_meca, const vec &mBC_meca, const mat &mmecas, const double &mBC_T, const int &mcBC_T, const vec &mTs) : step_thermomeca(mnumber, mDn_init, mDn_mini, mDn_inc, mmode, mcBC_meca, mBC_meca, mmecas, mBC_T, mcBC_T, mTs)
//-------------------------------------------------------------
{
    name = mname;
    nlgeom = mnlgeom;
    type = mtype;
    max_temp = mmax_temp;
}

/*!
 \brief Copy constructor
 \param bl block object to duplicate
 */

//------------------------------------------------------
aba_step_thermomeca::aba_step_thermomeca(const aba_step_thermomeca& stm) : step_thermomeca(stm)
//------------------------------------------------------
{
    name = stm.name;
    nlgeom = stm.nlgeom;
    type = stm.type;
}
    
//------------------------------------------------------
void aba_step_thermomeca::update(const step_thermomeca& stm, const string &mname, const bool &mnlgeom, const int &mtype, const int &mmax_temp)
//------------------------------------------------------
{
    name = mname;
    nlgeom = mnlgeom;
    type = mtype;
    max_temp = mmax_temp;
    
    number = stm.number;
    Dn_init = stm.Dn_init;
    Dn_mini = stm.Dn_mini;
    Dn_inc = stm.Dn_inc;
    ninc = stm.ninc;
    mode = stm.mode;
    
    BC_Time = stm.BC_Time;
    times = stm.times;
    file = stm.file;
    
    cBC_meca = stm.cBC_meca;
    BC_meca = stm.BC_meca;
    BC_T = stm.BC_T;
    cBC_T = stm.cBC_T;
    
    step::generate();
}
    

/*!
 \brief destructor
 */

aba_step_thermomeca::~aba_step_thermomeca() {}
    
/*!
 \brief Standard operator = for block
 */

//----------------------------------------------------------------------
aba_step_thermomeca& aba_step_thermomeca::operator = (const aba_step_thermomeca& stm)
//----------------------------------------------------------------------
{
//	assert(stm.ninc>0);
//	assert(stm.mode>0);

    name = stm.name;
    nlgeom = stm.nlgeom;
    type = stm.type;
    
    number = stm.number;
    Dn_init = stm.Dn_init;
    Dn_mini = stm.Dn_mini;
    Dn_inc = stm.Dn_inc;
    ninc = stm.ninc;
    mode = stm.mode;
    
    BC_Time = stm.BC_Time;
    times = stm.times;
    file = stm.file;
    
    cBC_meca = stm.cBC_meca;
    BC_meca = stm.BC_meca;
    mecas = stm.mecas;
    BC_T = stm.BC_T;
    cBC_T = stm.cBC_T;
    Ts = stm.Ts;
    
	return *this;
}
    
//-------------------------------------------------------------
void aba_step_thermomeca::write(const string &path_data, const string &inputfile, const bool &append)
//-------------------------------------------------------------
{
    std::string filename = path_data + "/" + inputfile;
    std::ofstream param_aba;
    
    if(append == true) {
        param_aba.open(filename, ios::app);
    }
    else {
        param_aba.open(filename, ios::out);
    }
    
    //Insert the step name and type, as well as the controls
    param_aba << "*Step, name=" << name << ", inc=" << ninc << "\n";
    if (type == 0) {
        param_aba << "*Coupled Temperature-displacement, creep=none, deltmx=" << max_temp << "\n" << Dn_init*Dn_inc*BC_Time << " ," << BC_Time << " ," << Dn_mini*Dn_inc*BC_Time << " ," << Dn_inc*BC_Time << "\n";
    }
    else if(type ==1) {
        param_aba << "*Coupled Temperature-displacement, creep=none" << Dn_inc*BC_Time << " ," << BC_Time << "\n";
    }
    else {
        cout << "Error, the type of the step is not recognized (0: Automatic, 1: Fixed)" << endl;
        return;
    }

    //Insert the 'Strain' and 'Stress' boundary conditions
    param_aba << "**\n";
    Col<int> CD = {11,22,33,12,13,23};

    param_aba << "**Strain - BC\n";
    param_aba << "*Boundary, op=NEW\n";
    param_aba << "CentreNode, 1, 1\n";
    param_aba << "CentreNode, 2, 2\n";
    param_aba << "CentreNode, 3, 3\n";
    for(int k = 0 ; k < 6 ; k++) {
        if(cBC_meca(k) == 0)
            param_aba << "CD" << CD(k) << " , 1, 1, " << BC_meca(k) << "\n";
    }

    param_aba << "**Stress - BC\n";
    param_aba << "*Cload, op=NEW\n";
    for(int k = 0 ; k < 6 ; k++) {
        if(cBC_meca(k) == 1)
            param_aba << "CD" << CD(k) << " , 1, 1, " << BC_meca(k) << "\n";
        else if(cBC_meca(k) == 2)
            param_aba << "CD" << CD(k) << " , 1, " << 0. << "\n";
    }
    
    param_aba << "**\n";

    param_aba << "*Output, field\n";
    param_aba << "*Element Output, directions=YES\n";
    param_aba << "S, E, SDV,\n";
    param_aba << "*Node Output, nset=AllNodes\n";
    param_aba << "U,\n";
    param_aba << "*Node Output, nset=CD_nodes\n";
    param_aba << "RF, CF, U,\n";
    param_aba << "*Node print, nset=CD_nodes, summary=no\n";
    param_aba << "RF1, CF1, U1,\n";
    param_aba << "*End Step\n";
    param_aba << "**" << endl;
    param_aba.close();
}
    
//--------------------------------------------------------------------------
ostream& operator << (ostream& s, const aba_step_thermomeca& stm)
//--------------------------------------------------------------------------
{
    s << "\tDisplay info on the step " << stm.number << "\n";
    s << "\tStep name: " << stm.name;
    s << "\tIs non-linear (NLGEOM): " << stm.nlgeom;
    s << "\tof type (0: Automatic - 1: direct): " << stm.type;
    s << "\tLoading mode: " << stm.mode << "\n";
    s << "\tThermomechanical step \n\t";
    
    Col<int> temp;
    temp = zeros<Col<int> >(6);
    temp(0) = 0;
	temp(1) = 3;
	temp(2) = 1;
	temp(3) = 4;
	temp(4) = 5;
	temp(5) = 2;
    
    if (stm.mode == 3) {
        for(int k = 0 ; k < 6 ; k++) {
            if(stm.cBC_meca(temp(k)) == 0)
                s << "E " << (((k==0)||(k==2)||(k==5)) ? "\n\t" : "\t");
            else if(stm.cBC_meca(temp(k)) == 1)
                s << "S " << (((k==0)||(k==2)||(k==5)) ? "\n\t" : "\t");
            else if(stm.cBC_meca(temp(k)) == 2)
                s << 0 << (((k==0)||(k==2)||(k==5)) ? "\n\t" : "\t");
        }
        s << "Temperature: ";
        if(stm.cBC_T == 0)
           s << "T\n";
        else if(stm.cBC_T == 1)
            s << "Q\n" << "\n";
        else if(stm.cBC_T == 2)
            s << 0 << "\n";
    }
    else {

        s << "\tTime of the step " << stm.BC_Time << " s\n\t";
        s << "\tInitial fraction: " << stm.Dn_init << "\tMinimal fraction: " << stm.Dn_mini << "\tIncrement fraction: " << stm.Dn_inc << "\n\t";
        for(int k = 0 ; k < 6 ; k++) {
            s << ((stm.cBC_meca(temp(k)) == 0) ? "\tE " : "\tS ") << stm.BC_meca(temp(k)) << (((k==0)||(k==2)||(k==5)) ? "\n\t" : "\t");
        }
        s << "Temperature at the end of step: " << stm.BC_T << "\n";
    }
	return s;
}

} //namespace simcoon
