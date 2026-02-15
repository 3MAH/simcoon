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

///@file step_meca.cpp
///@brief object that defines a mechanical step
///@version 1.0

#include <iostream>
#include <sstream>
#include <fstream>
#include <assert.h>
#include <math.h>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/exception.hpp>
#include <simcoon/Continuum_mechanics/Functions/kinematics.hpp>
#include <simcoon/Simulation/Solver/step.hpp>
#include <simcoon/Simulation/Solver/step_meca.hpp>
#include <simcoon/Simulation/Phase/state_variables_M.hpp>

using namespace std;
using namespace arma;

namespace simcoon{
    
//=====Private methods for ellipsoid_characteristics===================================

//=====Public methods for ellipsoid_characteristics============================================

//@brief default constructor
//-------------------------------------------------------------
step_meca::step_meca() : step()
//-------------------------------------------------------------
{
    BC_T = 0.;
    cBC_T = 0;
}

step_meca::step_meca(const unsigned int &control_type) : step()
{
    if (control_type == 5 || control_type == 6) { //Control with F (4) or gradU (5)
        cBC_meca = zeros<Col<int> >(9);
        BC_meca = zeros(9);
    }
    else {
        cBC_meca = zeros<Col<int> >(6);
        BC_meca = zeros(6);
    }
    BC_T = 0.;
    cBC_T = 0;
    BC_w = zeros(3,3);
    BC_R = eye(3,3);
}
    
    
/*!
 \brief Constructor with parameters
 cBC_meca = mcBC_meca;  Type of boundary conditions
 BC_meca = mBC_meca;    Target values of the BC
 BC_T = mBC_T;          target temperature
 Etot = mEtot;          current strain
 sigma = msigma;        current stress
 T = mT;                current temperature
 */    
    
//-------------------------------------------------------------
step_meca::step_meca(const int &mnumber, const double &mDn_init, const double &mDn_mini, const double &mDn_inc, const int &mmode, const unsigned int &mcontrol_type, const Col<int> &mcBC_meca, const vec &mBC_meca, const mat &mmecas, const mat &mBC_mecas, const double &mBC_T, const int &mcBC_T, const vec &mTs, const vec &mBC_Ts, const mat &mBC_w, const mat &mBC_R) : step(mnumber, mDn_init, mDn_mini, mDn_inc, mmode, mcontrol_type)
//-------------------------------------------------------------
{
    cBC_meca = mcBC_meca;
    BC_meca = mBC_meca;
    mecas = mmecas;
    BC_mecas = mBC_mecas;
    BC_T = mBC_T;
    cBC_T = mcBC_T;
    Ts = mTs;
    BC_Ts = mBC_Ts;
    BC_w = mBC_w;
    BC_R = mBC_R;
}

/*!
 \brief Copy constructor
 \param bl block object to duplicate
 */

//------------------------------------------------------
step_meca::step_meca(const step_meca& stm) : step(stm)
//------------------------------------------------------
{
    cBC_meca = stm.cBC_meca;
    BC_meca = stm.BC_meca;
    mecas = stm.mecas;
    BC_mecas = stm.mecas;
    BC_T = stm.BC_T;
    cBC_T = stm.cBC_T;
    Ts = stm.Ts;
    BC_Ts = stm.BC_Ts;
    BC_w = stm.BC_w;
    BC_R = stm.BC_R;
}

/*!
 \brief destructor
 */

step_meca::~step_meca() {}

//-------------------------------------------------------------
void step_meca::generate(const double &mTime, const vec &mEtot, const vec &msigma, const double &mT)
//-------------------------------------------------------------
{
    assert(control_type <= 4);
    
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
    BC_Ts = zeros(ninc);
    unsigned int size_meca = BC_meca.n_elem;
    mecas = zeros(ninc, size_meca);
    BC_mecas = zeros(ninc, size_meca);
    
    vec inc_coef = ones(ninc);          //If the mode is equal to 2, this is a sinuasoidal load control mode
    if (mode == 2) {
        double sum_ = 0.;
        for(int k = 0 ; k < ninc ; k++){
            inc_coef(k) =  cos(simcoon::pi+ (k+1)*2.*simcoon::pi/(ninc+1))+1.;
            sum_ += inc_coef(k);
        }
        inc_coef = inc_coef*ninc/sum_;
    }
    
    if (mode < 3) {
        for (int i=0; i<ninc; i++) {
            Ts(i) = (BC_T - mT)/ninc;
            BC_Ts(i) = (i+1)*((BC_T - mT)/ninc) + mT;
            times(i) = (BC_Time)/ninc;

            for(unsigned int k = 0 ; k < size_meca ; k++) {
                if (cBC_meca(k) == 1){
                    mecas(i,k) = inc_coef(i)*(BC_meca(k)-msigma(k))/ninc;
                    BC_mecas(i,k) = inc_coef(i)*(i+1)*((BC_meca(k)-msigma(k))/ninc) + msigma(k);
                }
                else if (cBC_meca(k) == 0){
                    mecas(i,k) = inc_coef(i)*(BC_meca(k)-mEtot(k))/ninc;
                    BC_mecas(i,k) = inc_coef(i)*(i+1)*((BC_meca(k)-mEtot(k))/ninc) + mEtot(k);
                }
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
        for (int i=0; i<ninc; i++) {

            pathinc >> buffer;
            for (unsigned int j=0; j<size_BC; j++) {
                pathinc >> BC_file(j);
            }
            
            times(i) = (BC_file(0) - BC_file_n(0));
            kT = 0;
            if (cBC_T == 0) {
                Ts(i) = BC_file(kT+1) - BC_file_n(kT+1);
                BC_Ts(i) = BC_file(kT+1);
                kT++;
            }
            else if(cBC_T == 2) {
                Ts(i) = 0.;
                BC_Ts(i) = 0.;
            }
            
            for(unsigned int k = 0 ; k < size_meca ; k++) {
                if (cBC_meca(k) < 2){
                    mecas(i,k) = BC_file(kT+1) - BC_file_n(kT+1);
                    BC_mecas(i,k) = BC_file(kT+1);
                    kT++;
                }
                else if (cBC_meca(k) == 2){
                    mecas(i,k) = 0.;
                    BC_mecas(i,k) = 0.;
                }
            }
            BC_file_n = BC_file;
            
        }
        //At the end, everything static becomes a stress-controlled with zeros
        for(unsigned int k = 0 ; k < size_meca ; k++) {
            if (cBC_meca(k) == 2)
                cBC_meca(k) = 1;
        }
        
	}
	else{
		cout << "\nError: The mode of the step number " << number << " does not correspond to an existing loading mode.\n";
	}
    
}

//-------------------------------------------------------------
void step_meca::generate_kin(const double &mTime, const mat &mF, const double &mT)
//-------------------------------------------------------------
{
    unsigned int size_meca = 9;
    assert(control_type > 4);
    for (unsigned int k=0; k<size_meca; k++) {
        assert(cBC_meca(k) != 1);
    }
    
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
    BC_Ts = zeros(ninc);
    mecas = zeros(ninc, size_meca);
    BC_mecas = zeros(ninc, size_meca);
    
    vec inc_coef = ones(ninc);          //If the mode is equal to 2, this is a sinuasoidal load control mode
    if (mode == 2) {
        double sum_ = 0.;
        for(int k = 0 ; k < ninc ; k++){
            inc_coef(k) =  cos(simcoon::pi+ (k+1)*2.*simcoon::pi/(ninc+1))+1.;
            sum_ += inc_coef(k);
        }
        inc_coef = inc_coef*ninc/sum_;
    }
    
    if (mode < 3) {

        // Previous deformation gradient
        arma::mat F_prev = mF;     // 3x3

        // Target deformation gradient
        arma::mat F_target(3,3);
        for (unsigned int k = 0; k < 9; ++k) {
            F_target(k/3, k%3) = BC_meca(k);
        }

        // Relative deformation over the step
        arma::mat F_prev_inv;
        bool inv_success = arma::inv(F_prev_inv, F_prev);
        if (!inv_success) {
            throw simcoon::exception_solver("Singular deformation gradient F_prev in step_meca::generate_kin.");
        }
        arma::mat F_tilde = F_target * F_prev_inv;

        // Logarithm of total transformation
        arma::cx_mat logF = arma::logmat(F_tilde);

        // Normalization of incremental weights
        double wsum = arma::accu(inc_coef);

        // Safety check (optional but recommended)
        if (wsum <= 0.0) {
            throw std::runtime_error("Sum of inc_coef must be positive");
        }

        // Initialize deformation at start of step
        arma::mat F_i = F_prev;

        for (unsigned int inc = 0; inc < ninc; ++inc) {

            // Incremental time for this sub-step
            times(inc) = BC_Time / ninc;

            // Incremental generator (normalized weights)
            arma::cx_mat delta_L = (inc_coef(inc) / wsum) * logF;

            // Exponential map
            arma::mat D_i = arma::real(arma::expmat(delta_L));

            // Multiplicative update
            F_i = D_i * F_i;

            // Store deformation history
            for (unsigned int k = 0; k < 9; ++k) {
                BC_mecas(inc,k) = F_i(k/3, k%3);

                if (inc == 0)
                    mecas(inc,k) = F_i(k/3, k%3) - F_prev(k/3, k%3);
                else
                    mecas(inc,k) = F_i(k/3, k%3) - BC_mecas(inc-1,k);
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
                if (control_type == 5) {
                    BC_file_n(kT+1) = mF(k/3,k%3);
                }
                else if (control_type == 6) {
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
                BC_Ts(i) = BC_file(kT+1);
                kT++;
            }
            else if(cBC_T == 2) {
                Ts(i) = 0.;
                BC_Ts(i) = 0.;
            }
            
            for(unsigned int k = 0 ; k < size_meca ; k++) {
                if (cBC_meca(k) < 2){
                    mecas(i,k) = BC_file(kT+1) - BC_file_n(kT+1);
                    BC_mecas(i,k) = BC_file(kT+1);
                    kT++;
                }
                else if (cBC_meca(k) == 2){
                    mecas(i,k) = 0.;
                    BC_mecas(i,k) = 0.;
                }
            }
            BC_file_n = BC_file;
        }
        //At the end, everything static becomes a deformation-controlled with zeros
        for(unsigned int k = 0 ; k < size_meca ; k++) {
            if (cBC_meca(k) == 2)
                cBC_meca(k) = 0;
        }
        
    }
    else{
        cout << "\nError: The mode of the step number " << number << " does not correspond to an existing loading mode.\n";
    }
    
}
    
    /*!
     \brief Standard operator = for block
     */
    
    //----------------------------------------------------------------------
void step_meca::assess_inc(const double &tnew_dt, double &tinc, const double &Dtinc, phase_characteristics &rve, double &Time, const double &DTime, const mat &DR, const int &corate_type) {
    
    if(tnew_dt < 1.){
        rve.to_start();
    }
    else {
        tinc += Dtinc;
        Time += DTime;
        BC_R = BC_R*DR;
        rve.set_start(corate_type);
    }
}
    
//----------------------------------------------------------------------
step_meca& step_meca::operator = (const step_meca& stm)
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
    BC_w = stm.BC_w;
    BC_R = stm.BC_R;
    
	return *this;
}
    
//--------------------------------------------------------------------------
ostream& operator << (ostream& s, const step_meca& stm)
//--------------------------------------------------------------------------
{
 
    s << "\tDisplay info on the step " << stm.number << "\n";
    s << "\tLoading mode: " << stm.mode << "\n";
    s << "\tMechanical step \n\t";
    
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
    
        if(stm.control_type == 5) {
            s << "Control: " << stm.control_type << " : Transformation gradient F\n";
        }
        else if (stm.control_type == 6) {
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
        else if(stm.cBC_T == 2)
            s << 0 << "\n";
    }
    else {

        s << "\tTime of the step " << stm.BC_Time << " s\n\t";
        s << "\tInitial fraction: " << stm.Dn_init << "\tMinimal fraction: " << stm.Dn_mini << "\tIncrement fraction: " << stm.Dn_inc << "\n\t";

        if(stm.control_type == 5) {
            s << "Control: " << stm.control_type << " : Eulerian velocity L\n";
            for(unsigned int k = 0 ; k <size_meca ; k++) {
                s << stm.BC_meca(temp(k)) << (((k==0)||(k==2)||(k==5)) ? "\n\t" : "\t");
            }
        }
        else {
            if(stm.control_type == 1) { s << "Control: " << stm.control_type << " : small strain hyp.\n";}
            if(stm.control_type == 2) { s << "Control: " << stm.control_type << " : Green-Lag / PKII\n";}
            if(stm.control_type == 3) { s << "Control: " << stm.control_type << " : True strain / stress\n";}
            if(stm.control_type == 4) { s << "Control: " << stm.control_type << " : Biot strain / stress\n";}            
            for(unsigned int k = 0 ; k <size_meca ; k++) {
                s << ((stm.cBC_meca(temp(k)) == 0) ? "\tE " : "\tS ") << stm.BC_meca(temp(k)) << (((k==0)||(k==2)||(k==5)) ? "\n\t" : "\t");
            }
        }
        s << "Temperature at the end of step: " << stm.BC_T << "\n";
    }
	return s;
}

} //namespace simcoon
