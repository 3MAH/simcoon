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

///@file step_thermomeca.hpp
///@brief object that defines a thermomechanical step
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
#include <simcoon/Simulation/Phase/state_variables_T.hpp>

using namespace std;
using namespace arma;

namespace simcoon{
    
//=====Private methods for ellipsoid_characteristics===================================

//=====Public methods for ellipsoid_characteristics============================================

//@brief default constructor
//-------------------------------------------------------------
step_thermomeca::step_thermomeca() : step()
//-------------------------------------------------------------
{
    BC_T = 0.;
    cBC_T = 0;
}
    
step_thermomeca::step_thermomeca(const unsigned int &control_type) : step()
{
    if (control_type == 4 || control_type == 5) { //Control with F (4) or gradU (5)
        cBC_meca = zeros<Col<int> >(9);
        BC_meca = zeros(9);
    }
    else {
        cBC_meca = zeros<Col<int> >(6);
        BC_meca = zeros(6);
    }
    BC_T = 0.;
    cBC_T = 0;
    BC_R = eye(3,3);
}

/*!
 \brief Constructor with parameters
 cBC_meca = mcBC_meca;  Type of mechanical boundary conditions
 BC_meca = mBC_meca;    Target values of the mechanical BC
 cBC_meca = mcBC_meca;  Type of thermal boundary conditions
 BC_T = mBC_T;          target values of the thermal BC
 */

//-------------------------------------------------------------
step_thermomeca::step_thermomeca(const int &mnumber, const double &mDn_init, const double &mDn_mini, const double &mDn_inc, const int &mmode, const unsigned int &mcontrol_type, const Col<int> &mcBC_meca, const vec &mBC_meca, const mat &mmecas, const double &mBC_T, const int &mcBC_T, const vec &mTs) : step(mnumber, mDn_init, mDn_mini, mDn_inc, mmode, mcontrol_type)
//-------------------------------------------------------------
{
    cBC_meca = mcBC_meca;
    BC_meca = mBC_meca;
    mecas = mmecas;
    BC_T = mBC_T;
    cBC_T = mcBC_T;
    Ts = mTs;
}

/*!
 \brief Copy constructor
 \param bl block object to duplicate
 */

//------------------------------------------------------
step_thermomeca::step_thermomeca(const step_thermomeca& stm) : step(stm)
//------------------------------------------------------
{
    cBC_meca = stm.cBC_meca;
    BC_meca = stm.BC_meca;
    mecas = stm.mecas;
    BC_T = stm.BC_T;
    cBC_T = stm.cBC_T;
    Ts = stm.Ts;
}

/*!
 \brief destructor
 */

step_thermomeca::~step_thermomeca() {}

//-------------------------------------------------------------
void step_thermomeca::generate(const double &mTime, const vec &mEtot, const vec &msigma, const double &mT, const mat& mR)
//-------------------------------------------------------------
{

    assert(control_type <= 3);
    
    //This in for the case of an incremental path file, to get the number of increments
    string buffer;
    ifstream pathinc;
    if(mode == 3){
        ninc = 0;
        pathinc.open(file, ios::in);
        if(!pathinc)
        {
            cout << "Error: cannot open the file << " << file << "\n Please check if the file is correct and is you have added the extension\n";
        }
        //read the file to get the number of increments
        while (!pathinc.eof())
        {
            getline (pathinc,buffer);
            if (buffer != "") {
                ninc++;
            }
        }
        pathinc.close();
    }
    
    step::generate();
    
    Ts = zeros(ninc);
    unsigned int size_meca = BC_meca.n_elem;
    mecas = zeros(ninc, size_meca);
    Rs = zeros(ninc,9);
    vec mR_vec = vectorise(mR);    
    
    vec inc_coef = ones(ninc);
    if (mode == 2) {
        double sum_ = 0.;
        for(int k = 0 ; k < ninc ; k++){
            inc_coef(k) =  cos(sim_pi+ (k+1)*2.*sim_pi/(ninc+1))+1.;
            sum_ += inc_coef(k);
        }
        inc_coef = inc_coef*ninc/sum_;
    }
    
    if (mode < 3) {
        for (int i=0; i<ninc; i++) {
            times(i) = (BC_Time)/ninc;
            
            for(unsigned int k = 0 ; k < size_meca ; k++) {
                if (cBC_meca(k) == 1){
                    mecas(i,k) = inc_coef(i)*(BC_meca(k)-msigma(k))/ninc;
                }
                else if (cBC_meca(k) == 0){
                    mecas(i,k) = inc_coef(i)*(BC_meca(k)-mEtot(k))/ninc;
                }
            }
            for(unsigned int k = 0 ; k < 9 ; k++) {
                Rs(i,k) = inc_coef(i)*(BC_R(k)-mR_vec(k))/ninc;
            }
            
            if (cBC_T == 1) {
                Ts(i) = BC_T;                   //Note that here it is the flux that is imposed
            }
            else if(cBC_T == 0) {
                Ts(i) = inc_coef(i)*(BC_T - mT)/ninc;
            }
            
        }
    }
    else if (mode ==3){ ///Incremental loading
        
        //Look at how many cBc are present to know the size of the file (1 for time + 6 for each meca + 1 for temperature):
        unsigned int size_BC = size_meca + 2;
        for(unsigned int k = 0 ; k < size_meca ; k++) {
            if (cBC_meca(k) == 2){
                size_BC--;
            }
        }
        if (cBC_T > 2 ) {
            size_BC--;
        }
        
        vec BC_file_n = zeros(size_BC); //vector that temporarly stores the previous values
        vec BC_file = zeros(size_BC); //vector that temporarly stores the values
        
        BC_file_n(0) = mTime;
        int kT = 0;
        if (cBC_T == 0) {
            BC_file_n(kT+1) = mT;
            kT++;
        }
        else if (cBC_T == 1) {
            BC_file_n(kT+1) = 0.;    //Heat flux does not depend on any previous condition
            kT++;
        }
        
        for (unsigned int k=0; k<size_meca; k++) {
            if (cBC_meca(k) == 0) {
                BC_file_n(kT+1) = mEtot(k);
                kT++;
            }
            if (cBC_meca(k) == 1) {
                BC_file_n(kT+1) = msigma(k);
                kT++;
            }
        }
                
        //Read all the informations and fill the meca accordingly
        pathinc.open(file, ios::in);
        
        //For mode 3, no rotation is considered yet
        mR_vec = vectorise(eye(3,3));
        
        for (int i=0; i<ninc; i++) {
            
            pathinc >> buffer;
            for (unsigned int j=0; j<size_BC; j++) {
                pathinc >> BC_file(j);
            }
            
            times(i) = (BC_file(0) - BC_file_n(0));
            kT = 0;
            if (cBC_T == 0) {
                Ts(i) = BC_file(kT+1) - BC_file_n(kT+1);
                kT++;
            }
            else if (cBC_T == 1) {
                Ts(i) = BC_file(kT+1);  //Case of Heat, direct quantity
                kT++;
            }
            else if(cBC_T == 2) {
                Ts(i) = 0.;
            }
            
            for(unsigned int k = 0 ; k < size_meca ; k++) {
                if (cBC_meca(k) < 2){
                    mecas(i,k) = BC_file(kT+1) - BC_file_n(kT+1);
                    kT++;
                }
                else if (cBC_meca(k) == 2){
                    mecas(i,k) = 0.;
                }
            }
            BC_file_n = BC_file;

            for(unsigned int k = 0 ; k < 9 ; k++) {
                Rs(i,k) = mR_vec(k);
            }
        }
        //At the end, everything static becomes a stress-controlled with zeros
        for(unsigned int k = 0 ; k < size_meca ; k++) {
            if (cBC_meca(k) == 2)
                cBC_meca(k) = 1;
        }
        //And everything thermally static is an isothermal path
        if (cBC_T == 2) {
            cBC_T = 0;
        }
                
	}
	else {
		cout << "\nError: The mode of the step number " << number << " does not correspond to an existing loading mode.\n";
	}
    
}
    
/*!
 \brief Standard operator = for block
 */

//-------------------------------------------------------------
void step_thermomeca::generate_kin(const double &mTime, const mat &mF, const double &mT)
//-------------------------------------------------------------
{
    
    unsigned int size_meca = 9;
    assert(control_type > 3);
    for (unsigned int k=0; k<size_meca; k++) {
        assert(cBC_meca(k) != 1);
    }
    assert (cBC_T == 0);
    
    mat I2 = eye(3,3);
    //This in for the case of an incremental path file, to get the number of increments
    string buffer;
    ifstream pathinc;
    if(mode == 3){
        ninc = 0;
        pathinc.open(file, ios::in);
        if(!pathinc)
        {
            cout << "Error: cannot open the file " << file << "\n Please check if the file is correct and is you have added the extension\n";
        }
        //read the file to get the number of increments
        while (!pathinc.eof())
        {
            getline (pathinc,buffer);
            if (buffer != "") {
                ninc++;
            }
        }
        pathinc.close();
    }
    
    step::generate();
    Ts = zeros(ninc);
    mecas = zeros(ninc, size_meca);
    
    vec inc_coef = ones(ninc);          //If the mode is equal to 2, this is a sinuasoidal load control mode
    if (mode == 2) {
        double sum_ = 0.;
        for(int k = 0 ; k < ninc ; k++){
            inc_coef(k) =  cos(sim_pi+ (k+1)*2.*sim_pi/(ninc+1))+1.;
            sum_ += inc_coef(k);
        }
        inc_coef = inc_coef*ninc/sum_;
    }
    
    if (mode < 3) {
        for (int i=0; i<ninc; i++) {
            Ts(i) = (BC_T - mT)/ninc;
            times(i) = (BC_Time)/ninc;
            
            for(unsigned int k = 0 ; k < size_meca ; k++) {
                if (control_type == 4) {
                    mecas(i,k) = inc_coef(i)*(BC_meca(k)-mF(k/3,k%3))/ninc;
                }
                else if (control_type == 5) {
                    mecas(i,k) = inc_coef(i)*(BC_meca(k)-(mF(k/3,k%3)-I2(k/3,k%3)))/ninc;
                }
                else {
                    cout << "ERROR in function generate_kin of step_meca.cpp : control_type should take the value 4 or 5 and not " << control_type << endl;
                }
                Ts(i) = inc_coef(i)*(BC_T - mT)/ninc;
            }
        }
    }
    else if (mode ==3){ ///Incremental loading
        
        //Look at how many cBc are present to know the size of the file (1 for time + 6/9 for each meca + 1 for temperature):
        unsigned int size_BC = size_meca + 2;
        for(unsigned int k = 0 ; k < size_meca ; k++) {
            if (cBC_meca(k) == 2){
                size_BC--;
            }
        }
        if (cBC_T == 2 ) {
            size_BC--;
        }
        
        vec BC_file_n = zeros(size_BC); //vector that temporarly stores the previous values
        vec BC_file = zeros(size_BC); //vector that temporarly stores the values
        
        BC_file_n(0) = mTime;
        int kT = 0;
        if (cBC_T == 0) {
            BC_file_n(kT+1) = mT;
            kT++;
        }
        for (unsigned int k=0; k<size_meca; k++) {
            if (cBC_meca(k) == 0) {
                if (control_type == 4) {
                    BC_file_n(kT+1) = mF(k/3,k%3);
                }
                else if (control_type == 5) {
                    BC_file_n(kT+1) = mF(k/3,k%3)-I2(k/3,k%3);
                }
                else {
                    cout << "ERROR in function generate_kin of step_meca.cpp : control_type should take the value 4 or 5 and not " << control_type << endl;
                }
                kT++;
            }
        }
        
        //Read all the informations and fill the meca accordingly
        pathinc.open(file, ios::in);
        
        for (int i=0; i<ninc; i++) {
            
            pathinc >> buffer;
            for (unsigned int j=0; j<size_BC; j++) {
                pathinc >> BC_file(j);
            }
            
            times(i) = (BC_file(0) - BC_file_n(0));
            kT = 0;
            if (cBC_T == 0) {
                Ts(i) = BC_file(kT+1) - BC_file_n(kT+1);
                kT++;
            }
            else if(cBC_T == 2) {
                Ts(i) = 0.;
            }
            
            for(unsigned int k = 0 ; k < size_meca ; k++) {
                if (cBC_meca(k) < 2){
                    mecas(i,k) = BC_file(kT+1) - BC_file_n(kT+1);
                    kT++;
                }
                else if (cBC_meca(k) == 2){
                    mecas(i,k) = 0.;
                }
            }
            BC_file_n = BC_file;
        }
        //At the end, everything static becomes a deformation-controlled with zeros
        for(unsigned int k = 0 ; k < size_meca ; k++) {
            if (cBC_meca(k) == 2)
                cBC_meca(k) = 0;
        }
        //And everything thermally static is an isothermal path
        if (cBC_T == 2) {
            cBC_T = 0;
        }
        
    }
    else{
        cout << "\nError: The mode of the step number " << number << " does not correspond to an existing loading mode.\n";
    }
    
}

//----------------------------------------------------------------------
void step_thermomeca::assess_inc(const double &tnew_dt, double &tinc, const double &Dtinc, phase_characteristics &rve, double &Time, const double &DTime) {
    
    if(tnew_dt < 1.){
        rve.to_start();
    }
    else {
        tinc += Dtinc;
        Time += DTime;
        rve.set_start();
    }
}
    
//----------------------------------------------------------------------
step_thermomeca& step_thermomeca::operator = (const step_thermomeca& stm)
//----------------------------------------------------------------------
{
//	assert(stm.ninc>0);
//	assert(stm.mode>0);
    
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
    
//--------------------------------------------------------------------------
ostream& operator << (ostream& s, const step_thermomeca& stm)
//--------------------------------------------------------------------------
{
    s << "\tDisplay info on the step " << stm.number << "\n";
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
    
    unsigned int size_meca = stm.BC_meca.n_elem;
    if (stm.mode == 3) {
    
        if(stm.control_type == 4) {
            s << "Control: " << stm.control_type << " : Transformation gradient F\n";
        }
        else if (stm.control_type == 5) {
            s << "Control: " << stm.control_type << " : gradient of displacement gradU\n";
        }
        else {
            if(stm.control_type == 1) { s << "Control: " << stm.control_type << " : small strain hyp.\n";}
            if(stm.control_type == 2) { s << "Control: " << stm.control_type << " : Green-Lag / PKII\n";}
            if(stm.control_type == 3) { s << "Control: " << stm.control_type << " : True strain / stress\n";}
            for(unsigned int k = 0 ; k < size_meca ; k++) {
                if(stm.cBC_meca(temp(k)) == 0)
                    s << "E " << (((k==0)||(k==2)||(k==5)) ? "\n\t" : "\t");
                else if(stm.cBC_meca(temp(k)) == 1)
                    s << "S " << (((k==0)||(k==2)||(k==5)) ? "\n\t" : "\t");
                else if(stm.cBC_meca(temp(k)) == 2)
                    s << 0 << (((k==0)||(k==2)||(k==5)) ? "\n\t" : "\t");
            }
        }
        cout << "Temperature: ";
        if(stm.cBC_T == 0)
            s << "T\n";
        else if(stm.cBC_T == 1)
            s << "Q\n";
        else if(stm.cBC_T == 2)
            s << 0 << "\n";
        else if(stm.cBC_T == 3)
            s << "convexion with Tau = " << stm.BC_T << "\n";
    }
    else {
        
        s << "\tTime of the step " << stm.BC_Time << " s\n\t";
        s << "\tInitial fraction: " << stm.Dn_init << "\tMinimal fraction: " << stm.Dn_mini << "\tIncrement fraction: " << stm.Dn_inc << "\n\t";
        if(stm.control_type == 4) {
            s << "Control: " << stm.control_type << " : Transformation gradient F\n";
            for(unsigned int k = 0 ; k <size_meca ; k++) {
                s << stm.BC_meca(temp(k)) << (((k==0)||(k==2)||(k==5)) ? "\n\t" : "\t");
            }
        }
        else if (stm.control_type == 5) {
            s << "Control: " << stm.control_type << " : gradient of displacement gradU\n";
            for(unsigned int k = 0 ; k <size_meca ; k++) {
                s << stm.BC_meca(temp(k)) << (((k==0)||(k==2)||(k==5)) ? "\n\t" : "\t");
            }
        }
        else {
            if(stm.control_type == 1) { s << "Control: " << stm.control_type << " : small strain hyp.\n";}
            if(stm.control_type == 2) { s << "Control: " << stm.control_type << " : Green-Lag / PKII\n";}
            if(stm.control_type == 3) { s << "Control: " << stm.control_type << " : True strain / stress\n";}
            for(unsigned int k = 0 ; k <size_meca ; k++) {
                s << ((stm.cBC_meca(temp(k)) == 0) ? "\tE " : "\tS ") << stm.BC_meca(temp(k)) << (((k==0)||(k==2)||(k==5)) ? "\n\t" : "\t");
            }
        }
        cout << "Temperature: ";
        if(stm.cBC_T == 0)
            s << "Temperature at the end of step: " << stm.BC_T << "\n";
        else if(stm.cBC_T == 1)
            s << "Heat flux imposed throughout the step: " << stm.BC_T << "\n";
        else if(stm.cBC_T == 2)
            s << "Temperature is constant throughout the step" << "\n";
        else if(stm.cBC_T == 3)
            s << "convexion with Tau = " << stm.BC_T << "\n";
    }
	return s;
}

} //namespace simcoon
