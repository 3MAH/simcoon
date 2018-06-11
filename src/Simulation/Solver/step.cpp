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

///@file step.cpp
///@brief object that defines an step
///@version 1.0

#include <iostream>
#include <sstream>
#include <fstream>
#include <assert.h>
#include <math.h>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Simulation/Solver/step.hpp>
#include <simcoon/Simulation/Solver/output.hpp>
#include <simcoon/Simulation/Phase/phase_characteristics.hpp>
#include <simcoon/Simulation/Phase/state_variables.hpp>

using namespace std;
using namespace arma;

namespace simcoon{

//=====Private methods for ellipsoid_characteristics===================================

//=====Public methods for ellipsoid_characteristics============================================

//@brief default constructor
//-------------------------------------------------------------
step::step()
//-------------------------------------------------------------
{
	number=0;
    Dn_init = 0.;
    Dn_mini = 0.;
    Dn_inc = 0.;
	ninc=0;
	mode=0;
    control_type=0;
    
    BC_Time = 0.;
    
    file = "";
}

/*!
 \brief Constructor with parameters
 \param mnumber : number of the step
 \param mninc : number of increments
 \param mmode : mode : type of loading
 */

//-------------------------------------------------------------
step::step(const int &mnumber, const double &mDn_init, const double &mDn_mini, const double &mDn_inc, const int &mmode, const unsigned int &mcontrol_type)
//-------------------------------------------------------------
{
    
    assert(mnumber>=0);
	assert(mDn_inc<=1.);
	assert(mDn_init<=mDn_inc);
	assert(mDn_mini<=mDn_init);
    assert(fmod(1.,mDn_inc)==0.);
    
	number = mnumber;
    Dn_init = mDn_init;
    Dn_mini = mDn_mini;
    Dn_inc = mDn_inc;
    ninc = std::round(1./mDn_inc);
	mode = mmode;
    control_type = mcontrol_type;
    
    times = zeros(ninc);
    BC_Time = 0.;
    
    file = "";
}

/*!
 \brief Copy constructor
 \param st step object to duplicate
 */

//------------------------------------------------------
step::step(const step& st)
//------------------------------------------------------
{    
	number = st.number;
    Dn_init = st.Dn_init;
    Dn_mini = st.Dn_mini;
    Dn_inc = st.Dn_inc;
	ninc = st.ninc;
	mode = st.mode;
    control_type = st.control_type;
    
    times = st.times;
    BC_Time = st.BC_Time;
    
    file = st.file;
}

/*!
 \brief destructor
 */

step::~step() {}

//-------------------------------------------------------------
void step::generate()
//-------------------------------------------------------------
{
	assert(number>=0);
    assert(Dn_inc<=1.);
	assert(mode>0);
    
    if(mode == 3) {
        Dn_inc = 1./ninc;
    }
    else {
        ninc = std::round(1./Dn_inc);
        assert(ninc*Dn_inc==1.);
    }
    times = zeros(ninc);
}

//----------------------------------------------------------------------
void step::compute_inc(double &tnew_dt, const int &inc, double &tinc, double &Dtinc, double &Dtinc_cur, const int &inforce_solver) {
//----------------------------------------------------------------------
    
    if((inc == 0)&&(Dtinc == 0.)){
        Dtinc_cur = Dn_init;
    }
    else {
        Dtinc_cur = Dtinc_cur*tnew_dt;
    }
    
    if (Dtinc_cur < Dn_mini) {
        if (inforce_solver == 1) {
//            cout << "Warning : The solver has been forced to continue with the minial increment at step:" << number << " inc: " << inc << " and fraction:" << tinc << "\n";
            //tnew_dt = 1.;
            Dtinc_cur = Dn_mini;
        }
        else {
//            cout << "\nThe increment size is less than the minimum specified\n";
            exit(0);
        }
        
    }
    
    if (Dtinc_cur >= 1.) {
        Dtinc_cur = 1.;
    }
    
    Dtinc = Dtinc_cur;
        
    if(tinc + Dtinc > 1.) {
        Dtinc = 1.-tinc;
    }
    
}
    
/*!
 \brief Standard operator = for block
 */

//----------------------------------------------------------------------
step& step::operator = (const step& st)
//----------------------------------------------------------------------
{
//	assert(st.ninc>0);
//	assert(st.mode>0);
    
	number = st.number;
    Dn_init = st.Dn_init;
    Dn_mini = st.Dn_mini;
    Dn_inc = st.Dn_inc;
    ninc = st.ninc;
	mode = st.mode;
    control_type = st.control_type;
    
    times = st.times;
    BC_Time = st.BC_Time;
    
    file = st.file;
    
	return *this;
}

//--------------------------------------------------------------------------
ostream& operator << (ostream& s, const step& st)
//--------------------------------------------------------------------------
{    
	s << "Display info on the step " << st.number << "\n";
	s << "Number of increments: " << st.ninc << "\twithin " << st.BC_Time << " s\n";
	s << "Loading mode: " << st.mode << "\n";
	s << "Control type: " << st.control_type << "\n";
	   
	return s;
}

} //namespace simcoon
