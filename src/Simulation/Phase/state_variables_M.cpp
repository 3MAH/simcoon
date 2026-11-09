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

///@file state_variables_M.cpp
///@brief State variables of a mechanical material, in a defined coordinate system:
///@version 1.0

#include <iostream>
#include <fstream>
#include <assert.h>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Simulation/Maths/rotation.hpp>
#include <simcoon/Simulation/Phase/state_variables.hpp>
#include <simcoon/Simulation/Phase/state_variables_M.hpp>

using namespace std;
using namespace arma;

namespace simcoon{

//=====Private methods for state_variables===================================

//=====Public methods for state_variables============================================

//-------------------------------------------------------------
state_variables_M::state_variables_M() : state_variables(), sigma_in(6), sigma_in_start(6), Wm(4), Wm_start(4), L(6,6), Lt(6,6)
//-------------------------------------------------------------
{
    sigma_in = zeros(6);
    sigma_in_start = zeros(6);
    
    Wm = zeros(4);
    Wm_start = zeros(4);
    
	L = zeros(6,6);
	Lt = zeros(6,6);
}

//-------------------------------------------------------------
state_variables_M::state_variables_M(const vec &mEtot, const vec &mDEtot, const vec &metot, const vec &mDetot, const vec &mPKII, const vec &mPKII_start, const vec &mtau, const vec &mtau_start, const vec &msigma, const vec &msigma_start, const mat &mF0, const mat &mF1, const mat &mU0, const mat &mU1, const mat &mR, const mat &mDR, const vec &msigma_in, const vec &msigma_in_start, const double &mT, const double &mDT, const int &mnstatev, const vec &mstatev, const vec &mstatev_start, const vec &mWm, const vec& mWm_start, const mat &mL, const mat &mLt) : state_variables(mEtot, mDEtot, metot, mDetot, mPKII, mPKII_start, mtau, mtau_start, msigma, msigma_start, mF0, mF1, mU0, mU1, mR, mDR, mT, mDT, mnstatev, mstatev, mstatev_start), sigma_in(6), sigma_in_start(6), Wm(4), Wm_start(4), L(6,6), Lt(6,6)
//-------------------------------------------------------------
{

    assert (msigma_in.size() == 6);
    assert (msigma_in_start.size() == 6);
    
    assert (mWm.size() == 4);
    assert (mWm_start.size() == 4);
    
	assert (mL.n_rows == 6);
	assert (mL.n_cols == 6);
	assert (mLt.n_rows == 6);
	assert (mLt.n_cols == 6);	

    sigma_in = msigma_in;
    sigma_in_start = msigma_in_start;
    
    Wm = mWm;
    Wm_start = mWm_start;
    
	L = mL;
	Lt = mLt;
}

//------------------------------------------------------
state_variables_M::state_variables_M(const state_variables_M& sv) : state_variables(sv), sigma_in(6), sigma_in_start(6), Wm(4), Wm_start(4), L(6,6), Lt(6,6)
//------------------------------------------------------
{
    sigma_in = sv.sigma_in;
    sigma_in_start = sv.sigma_in_start;
    
    Wm = sv.Wm;
    Wm_start = sv.Wm_start;
    
	L = sv.L;
	Lt = sv.Lt;
}

//-------------------------------------
state_variables_M::~state_variables_M()
//-------------------------------------
{

}

//----------------------------------------------------------------------
state_variables_M& state_variables_M::operator = (const state_variables_M& sv)
//----------------------------------------------------------------------
{
    state_variables::operator=(sv);   //common fields, U0, U1, R, DR and tangent_mode included
    sigma_in = sv.sigma_in;
    sigma_in_start = sv.sigma_in_start;
	L = sv.L;
	Lt = sv.Lt;
    Wm = sv.Wm;
    Wm_start = sv.Wm_start;

	return *this;

}

//----------------------------------------------------------------------
state_variables_M& state_variables_M::copy_fields_M(const state_variables_M& sv)
//----------------------------------------------------------------------
{
    state_variables::copy_fields(sv);   //common fields, U0, U1, R and DR included
	sigma_in = sv.sigma_in;
    sigma_in_start = sv.sigma_in_start;
	L = sv.L;
	Lt = sv.Lt;
    Wm = sv.Wm;
    Wm_start = sv.Wm_start;

	return *this;

}

//-------------------------------------------------------------
void state_variables_M::update(const vec &mEtot, const vec &mDEtot, const vec &metot, const vec &mDetot, const vec &mPKII, const vec &mPKII_start, const vec &mtau, const vec &mtau_start, const vec &msigma, const vec &msigma_start, const mat &mF0, const mat &mF1, const mat &mU0, const mat &mU1, const mat &mR, const mat &mDR, const vec &msigma_in, const vec &msigma_in_start, const double &mT, const double &mDT, const int &mnstatev, const vec &mstatev, const vec &mstatev_start, const vec &mWm, const vec &mWm_start, const mat &mL, const mat &mLt)
//-------------------------------------------------------------
{
    state_variables::update(mEtot, mDEtot, metot, mDetot, mPKII, mPKII_start, mtau, mtau_start, msigma, msigma_start, mF0, mF1, mU0, mU1, mR, mDR, mT, mDT, mnstatev, mstatev, mstatev_start);

    assert (msigma_in.size() == 6);
    assert (msigma_in_start.size() == 6);
    
    assert (mWm.size() == 4);
    assert (mWm_start.size() == 4);
    
    assert (mL.n_rows == 6);
    assert (mL.n_cols == 6);
    assert (mLt.n_rows == 6);
    assert (mLt.n_cols == 6);
    
    sigma_in = msigma_in;
    sigma_in_start = msigma_in_start;
    
    Wm = mWm;
    Wm_start = mWm_start;
    
    L = mL;
    Lt = mLt;
}
    
//-------------------------------------------------------------
void state_variables_M::to_start()
//-------------------------------------------------------------
{
    state_variables::to_start();
    sigma_in = sigma_in_start;
    Wm = Wm_start;
}

//-------------------------------------------------------------
void state_variables_M::set_start(const int &corate_type)
//-------------------------------------------------------------
{
    state_variables::set_start(corate_type);
    sigma_in_start = sigma_in;
    Wm_start = Wm;
}
        
//----------------------------------------------------------------------
state_variables_M& state_variables_M::rotate_l2g(const state_variables_M& sv, const frame_rotation &frame)
//----------------------------------------------------------------------
{

    state_variables::rotate_l2g(sv, frame);

    Wm = sv.Wm;
	Lt = sv.Lt;

    frame.rotate_stiffness(Lt);

	return *this;
}

//----------------------------------------------------------------------
state_variables_M& state_variables_M::rotate_g2l(const state_variables_M& sv, const frame_rotation &frame)
//----------------------------------------------------------------------
{

    state_variables::rotate_g2l(sv, frame);

    Wm = sv.Wm;
    L = sv.L;
    Lt = sv.Lt;

    // inverse operators: global -> local
    frame.rotate_stiffness(L, frame_rotation::inverse);
    frame.rotate_stiffness(Lt, frame_rotation::inverse);

	return *this;
}

//--------------------------------------------------------------------------
ostream& operator << (ostream& s, const state_variables_M& sv)
//--------------------------------------------------------------------------
{
	s << "Etot: \n" << sv.Etot << "\n";
	s << "DEtot: \n" << sv.DEtot << "\n";
	s << "etot: \n" << sv.etot << "\n";
	s << "Detot: \n" << sv.Detot << "\n";
	s << "PKII: \n" << sv.PKII << "\n";
	s << "PKII_start: \n" << sv.PKII_start << "\n";
	s << "tau: \n" << sv.tau << "\n";
	s << "tau_start: \n" << sv.tau_start << "\n";    
	s << "sigma: \n" << sv.sigma << "\n";
	s << "sigma_start: \n" << sv.sigma_start << "\n";
    s << "F0: \n" << sv.F0 << "\n";
    s << "F1: \n" << sv.F1 << "\n";
    s << "sigma_in: \n" << sv.sigma_in << "\n";
    s << "sigma_in_start: \n" << sv.sigma_in_start << "\n";    
    s << "T: \n" << sv.T << "\n";
    s << "DT: \n" << sv.DT << "\n";
    s << "Wm: \n" << sv.Wm << "\n";
    s << "Wm_start: \n" << sv.Wm_start << "\n";
	s << "L: \n" << sv.L << "\n";
	s << "Lt: \n" << sv.Lt << "\n";
	s << "\n";

    s << "nstatev: \n" << sv.nstatev << "\n";
    if (sv.nstatev) {
        s << "statev: \n";
        s << sv.statev.t();
        s << "\n";
        s << sv.statev_start.t();
        s << "\n";
    }
    
	return s;
}

} //namespace simcoon
