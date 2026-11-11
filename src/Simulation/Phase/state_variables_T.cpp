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

///@file state_variables_T.cpp
///@brief State variables of a phase, in a defined coordinate system:
///@version 1.0

#include <iostream>
#include <fstream>
#include <assert.h>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Simulation/Phase/state_variables.hpp>
#include <simcoon/Simulation/Phase/state_variables_T.hpp>
#include <simcoon/Simulation/Maths/rotation.hpp>

using namespace std;
using namespace arma;

namespace simcoon{

//=====Private methods for state_variables_T===================================

//=====Public methods for state_variables_T============================================

//-------------------------------------------------------------
state_variables_T::state_variables_T() : state_variables()
//-------------------------------------------------------------
{
    Q = 0.;
    r = 0.;
    r_in = 0.;
    Wm = zeros(4);
    Wm_start = zeros(4);
    Wt = zeros(3);
    Wt_start = zeros(3);
    dSdE = zeros(6,6);
    dSdEt = zeros(6,6);
    dSdT.zeros();
    drdE.zeros();
    drdT = 0.;
}

//-------------------------------------------------------------
state_variables_T::state_variables_T(const vec &mEtot, const vec &mDEtot, const vec &metot, const vec &mDetot, const vec &mPKII, const vec &mPKII_start, const vec &mtau, const vec &mtau_start, const vec &msigma, const vec &msigma_start, const mat &mF0, const mat &mF1, const mat &mU0, const mat &mU1, const mat &mR, const mat &mDR, const double &mT, const double &mDT, const int &mnstatev, const vec &mstatev, const vec &mstatev_start, const double &mQ, const double &mr, const double &mr_in, const vec &mWm, const vec &mWt, const vec &mWm_start, const vec &mWt_start, const mat &mdSdE, const mat &mdSdEt, const vec &mdSdT, const vec &mdrdE, const double &mdrdT) : state_variables(mEtot, mDEtot, metot, mDetot, mPKII, mPKII_start, mtau, mtau_start, msigma, msigma_start, mF0, mF1, mU0, mU1, mR, mDR, mT, mDT, mnstatev, mstatev, mstatev_start)
//-------------------------------------------------------------
{	

    
    assert (mWm.size() == 4);
    assert (mWt.size() == 3);
    
    assert (mWm_start.size() == 4);
    assert (mWt_start.size() == 3);
    
	assert (mdSdE.n_rows == 6);
	assert (mdSdE.n_cols == 6);

    assert (mdSdEt.n_rows == 6);
	assert (mdSdEt.n_cols == 6);
    
    assert (mdSdT.n_elem == 6);
    
    assert (mdrdE.n_elem == 6);
    
    Q = mQ;
    r = mr;
    r_in = mr_in;
    Wm = mWm;
    Wm_start = mWm_start;
    Wt = mWt;
    Wt_start = mWt_start;
    dSdE = mdSdE;
	dSdEt = mdSdEt;
    dSdT = mdSdT;
    drdE = mdrdE;
    drdT = mdrdT;
}

//------------------------------------------------------
state_variables_T::state_variables_T(const state_variables_T& sv) : state_variables(sv)
//------------------------------------------------------
{
    Q = sv.Q;
    r = sv.r;
    r_in = sv.r_in;
    Wm = sv.Wm;
    Wt = sv.Wt;
    Wm_start = sv.Wm_start;
    Wt_start = sv.Wt_start;
    dSdE = sv.dSdE;
	dSdEt = sv.dSdEt;
    dSdT = sv.dSdT;
    drdE = sv.drdE;
    drdT = sv.drdT;
}

//-------------------------------------
state_variables_T::~state_variables_T()
//-------------------------------------
{

}

//----------------------------------------------------------------------
state_variables_T& state_variables_T::operator = (const state_variables_T& sv)
//----------------------------------------------------------------------
{
    state_variables::operator=(sv);   //common fields, U0, U1, R, DR and tangent_mode included
    dSdE = sv.dSdE;
    dSdEt = sv.dSdEt;
    dSdT = sv.dSdT;
    drdE = sv.drdE;
    drdT = sv.drdT;
    Wm = sv.Wm;
    Wt = sv.Wt;
    Wm_start = sv.Wm_start;
    Wt_start = sv.Wt_start;

    Q = sv.Q;
    r = sv.r;
    r_in = sv.r_in;

	return *this;
}

//----------------------------------------------------------------------
state_variables_T& state_variables_T::copy_fields_T (const state_variables_T& sv)
//----------------------------------------------------------------------
{
    state_variables::copy_fields(sv);   //common fields, U0, U1, R and DR included
    dSdE = sv.dSdE;
    dSdEt = sv.dSdEt;
    dSdT = sv.dSdT;
    drdE = sv.drdE;
    drdT = sv.drdT;
    Wm = sv.Wm;
    Wt = sv.Wt;
    Wm_start = sv.Wm_start;
    Wt_start = sv.Wt_start;
    
    Q = sv.Q;
    r = sv.r;
    r_in = sv.r_in;

	return *this;
}

//-------------------------------------------------------------
void state_variables_T::update(const vec &mEtot, const vec &mDEtot, const vec &metot, const vec &mDetot, const vec &mPKII, const vec &mPKII_start, const vec &mtau, const vec &mtau_start, const vec &msigma, const vec &msigma_start, const mat &mF0, const mat &mF1, const mat &mU0, const mat &mU1, const mat &mR, const mat &mDR, const double &mT, const double &mDT, const int &mnstatev, const vec &mstatev, const vec &mstatev_start, const double &mQ, const double &mr, const double &mr_in, const vec &mWm, const vec &mWt, const vec &mWm_start, const vec &mWt_start, const mat &mdSdE, const mat &mdSdEt, const vec &mdSdT, const vec &mdrdE, const double &mdrdT)
//-------------------------------------------------------------
{
    state_variables::update(mEtot, mDEtot, metot, mDetot, mPKII, mPKII_start, mtau, mtau_start, msigma, msigma_start, mF0, mF1, mU0, mU1, mR, mDR, mT, mDT, mnstatev, mstatev, mstatev_start);
    
    
    assert (mWm.size() == 4);
    assert (mWt.size() == 3);
    
    assert (mWm_start.size() == 4);
    assert (mWt_start.size() == 3);
    
    assert (mdSdE.n_rows == 6);
    assert (mdSdE.n_cols == 6);
    
    assert (mdSdEt.n_rows == 6);
    assert (mdSdEt.n_cols == 6);
    
    assert (mdSdT.n_elem == 6);
    
    assert (mdrdE.n_elem == 6);
    
    Q = mQ;
    r = mr;
    r_in = mr_in;
    Wm = mWm;
    Wt = mWt;
    Wm_start = mWm_start;
    Wt_start = mWt_start;
    dSdE = mdSdE;
    dSdEt = mdSdEt;
    dSdT = mdSdT;
    drdE = mdrdE;
    drdT = mdrdT;
}
    
//-------------------------------------------------------------
void state_variables_T::to_start()
//-------------------------------------------------------------
{
    state_variables::to_start();
    Wm = Wm_start;
    Wt = Wt_start;
}

//-------------------------------------------------------------
void state_variables_T::set_start(const int &corate_type)
//-------------------------------------------------------------
{
    state_variables::set_start(corate_type);
    Wm_start = Wm;
    Wt_start = Wt;
}
    
//----------------------------------------------------------------------
state_variables_T& state_variables_T::rotate_l2g(const state_variables_T& sv, const frame_rotation &frame)
//----------------------------------------------------------------------
{

    state_variables::rotate_l2g(sv, frame);

    dSdE = sv.dSdE;
    dSdT = sv.dSdT;
    drdE = sv.drdE;
    drdT = sv.drdT;
    r = sv.r;
    Wm = sv.Wm;
    Wt = sv.Wt;

    frame.rotate_stiffness(dSdE);
    frame.rotate_stress(dSdT);   // dSdT, drdE: Voigt covectors, rotate with the stress operator (header)
    frame.rotate_stress(drdE);

	return *this;
}

//----------------------------------------------------------------------
state_variables_T& state_variables_T::rotate_g2l(const state_variables_T& sv, const frame_rotation &frame)
//----------------------------------------------------------------------
{

    state_variables::rotate_g2l(sv, frame);

    Q = sv.Q;
    Wm = sv.Wm;
    Wt = sv.Wt;

	return *this;
}

//--------------------------------------------------------------------------
ostream& operator << (ostream& s, const state_variables_T& sv)
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
    s << "T: \n" << sv.T << "\n";
    s << "DT: \n" << sv.DT << "\n";
    s << "Q: \n" << sv.Q << "\n";
    s << "r: \n" << sv.r << "\n";
    s << "r_in: \n" << sv.r_in << "\n";
    s << "Wm: \n" << sv.Wm << "\n";
    s << "Wm_start: \n" << sv.Wm_start << "\n";
    s << "Wt: \n" << sv.Wt << "\n";
    s << "Wt_start: \n" << sv.Wt_start << "\n";
	s << "dSdE: \n" << sv.dSdE << "\n";
	s << "dSdEt: \n" << sv.dSdEt << "\n";
	s << "dSdT: \n" << sv.dSdT << "\n";
	s << "drdE: \n" << sv.drdE << "\n";
	s << "drdT: \n" << sv.drdT << "\n";
    
    s << "nstatev: \n" << sv.nstatev << "\n";
    if (sv.nstatev) {
        s << "statev: \n";
        s << sv.statev.t();
        s << "\n";
        s << sv.statev_start.t();
        s << "\n";
    }    
    
	s << "\n";

	return s;
}

} //namespace simcoon
