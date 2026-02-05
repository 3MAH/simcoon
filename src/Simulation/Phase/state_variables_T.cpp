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
#include <simcoon/Simulation/Maths/rotation_class.hpp>
#include <simcoon/Continuum_mechanics/Functions/natural_basis.hpp>

using namespace std;
using namespace arma;

namespace simcoon{

//=====Private methods for state_variables_T===================================

//=====Public methods for state_variables_T============================================

/*!
  \brief default constructor
*/

//-------------------------------------------------------------
state_variables_T::state_variables_T() : state_variables(), sigma_in(6), sigma_in_start(6), Wm(4), Wt(3), Wm_start(4), Wt_start(3), dSdE(6,6), dSdEt(6,6), dSdT(1,6), drdE(1,6), drdT(1,1)
//-------------------------------------------------------------
{
    Q = 0.;
    r = 0.;
    r_in = 0.;
    sigma_in = zeros(6);
    sigma_in_start = zeros(6);
    Wt = zeros(3);
    Wt_start = zeros(3);
    dSdE = zeros(6,6);
    dSdEt = zeros(6,6);
    dSdT = zeros(1,6);
    drdE = zeros(1,6);
    drdT = zeros(1,1);
}

/*!
  \brief Constructor with parameters
  \n\n
  \f$ \textbf{Examples :} \f$ \n
*/

//-------------------------------------------------------------
state_variables_T::state_variables_T(const vec &mEtot, const vec &mDEtot, const vec &metot, const vec &mDetot, const vec &mPKII, const vec &mPKII_start, const vec &mtau, const vec &mtau_start, const vec &msigma, const vec &msigma_start, const mat &mF0, const mat &mF1, const mat &mU0, const mat &mU1, const mat &mR, const mat &mDR, const vec &msigma_in, const vec &msigma_in_start, const double &mT, const double &mDT, const int &mnstatev, const vec &mstatev, const vec &mstatev_start, const natural_basis &mnb, const double &mQ, const double &mr, const double &mr_in, const vec &mWm, const vec &mWt, const vec &mWm_start, const vec &mWt_start, const mat &mdSdE, const mat &mdSdEt, const mat &mdSdT, const mat &mdrdE, const mat &mdrdT) : state_variables(mEtot, mDEtot, metot, mDetot, mPKII, mPKII_start, mtau, mtau_start, msigma, msigma_start, mF0, mF1, mU0, mU1, mR, mDR, mT, mDT, mnstatev, mstatev, mstatev_start, mnb), sigma_in(6), sigma_in_start(6), Wm(4), Wt(3), Wm_start(4), Wt_start(3), dSdE(6,6), dSdEt(6,6), dSdT(1,6), drdE(1,6), drdT(1,1)
//-------------------------------------------------------------
{	

    assert (msigma_in.size() == 6);
    assert (msigma_in_start.size() == 6);
    
    assert (mWm.size() == 4);
    assert (mWt.size() == 3);
    
    assert (mWm_start.size() == 4);
    assert (mWt_start.size() == 3);
    
	assert (mdSdE.n_rows == 6);
	assert (mdSdE.n_cols == 6);

    assert (mdSdEt.n_rows == 6);
	assert (mdSdEt.n_cols == 6);
    
    assert (mdSdT.n_rows == 1);
    assert (mdSdT.n_cols == 6);
    
    assert (mdrdE.n_rows == 1);
    assert (mdrdE.n_cols == 6);
    
    assert (mdrdT.n_rows == 1);
    assert (mdrdT.n_cols == 1);

    Q = mQ;
    r = mr;
    r_in = mr_in;
    sigma_in = msigma_in;
    sigma_in_start = msigma_in_start;
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

/*!
  \brief Copy constructor
  \param s state_variables_T object to duplicate
*/

//------------------------------------------------------
state_variables_T::state_variables_T(const state_variables_T& sv) : state_variables(sv), sigma_in(6), sigma_in_start(6), Wm(4), Wt(3), Wm_start(4), Wt_start(3), dSdE(6,6), dSdEt(6,6), dSdT(1,6), drdE(1,6), drdT(1,1)
//------------------------------------------------------
{
    Q = sv.Q;
    r = sv.r;
    r_in = sv.r_in;
    sigma_in = sv.sigma_in;
    sigma_in_start = sv.sigma_in_start;
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

/*!
  \brief Destructor

  Deletes statev_variables, the vectors and matrix, and table statev if nstatev is not null.
*/

//-------------------------------------
state_variables_T::~state_variables_T()
//-------------------------------------
{

}

/*!
  \brief Standard operator = for state_variables
*/

//----------------------------------------------------------------------
state_variables_T& state_variables_T::operator = (const state_variables_T& sv)
//----------------------------------------------------------------------
{
	Etot = sv.Etot;
	DEtot = sv.DEtot;
	etot = sv.etot;
	Detot = sv.Detot;
	PKII = sv.PKII;
	PKII_start = sv.PKII_start;
	tau = sv.tau;
	tau_start = sv.tau_start;
    sigma = sv.sigma;
    sigma_start = sv.sigma_start;
    F0 = sv.F0;
    F1 = sv.F1;
    sigma_in = sv.sigma_in;
    sigma_in_start = sv.sigma_in_start;
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
    T = sv.T;
    DT = sv.DT;
    
    nstatev = sv.nstatev;
    statev = sv.statev;
    statev_start = sv.statev_start;
    
    nb = sv.nb;

	return *this;
}

//----------------------------------------------------------------------
state_variables_T& state_variables_T::copy_fields_T (const state_variables_T& sv)
//----------------------------------------------------------------------
{
	Etot = sv.Etot;
	DEtot = sv.DEtot;
	etot = sv.etot;
	Detot = sv.Detot;
	PKII = sv.PKII;
	PKII_start = sv.PKII_start;
	tau = sv.tau;
	tau_start = sv.tau_start;
    sigma = sv.sigma;
    sigma_start = sv.sigma_start;
    F0 = sv.F0;
    F1 = sv.F1;
    sigma_in = sv.sigma_in;
    sigma_in_start = sv.sigma_in_start;
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
    T = sv.T;
    DT = sv.DT;

	return *this;
}

//-------------------------------------------------------------
void state_variables_T::update(const vec &mEtot, const vec &mDEtot, const vec &metot, const vec &mDetot, const vec &mPKII, const vec &mPKII_start, const vec &mtau, const vec &mtau_start, const vec &msigma, const vec &msigma_start, const mat &mF0, const mat &mF1, const mat &mU0, const mat &mU1, const mat &mR, const mat &mDR, const vec &msigma_in, const vec &msigma_in_start, const double &mT, const double &mDT, const int &mnstatev, const vec &mstatev, const vec &mstatev_start, const natural_basis &mnb, const double &mQ, const double &mr, const double &mr_in, const vec &mWm, const vec &mWt, const vec &mWm_start, const vec &mWt_start, const mat &mdSdE, const mat &mdSdEt, const mat &mdSdT, const mat &mdrdE, const mat &mdrdT)
//-------------------------------------------------------------
{
    state_variables::update(mEtot, mDEtot, metot, mDetot, mPKII, mPKII_start, mtau, mtau_start, msigma, msigma_start, mF0, mF1, mU0, mU1, mR, mDR, mT, mDT, mnstatev, mstatev, mstatev_start, mnb);
    
    assert (msigma_in.size() == 6);
    assert (msigma_in_start.size() == 6);
    
    assert (mWm.size() == 4);
    assert (mWt.size() == 3);
    
    assert (mWm_start.size() == 4);
    assert (mWt_start.size() == 3);
    
    assert (mdSdE.n_rows == 6);
    assert (mdSdE.n_cols == 6);
    
    assert (mdSdEt.n_rows == 6);
    assert (mdSdEt.n_cols == 6);
    
    assert (mdSdT.n_rows == 1);
    assert (mdSdT.n_cols == 6);
    
    assert (mdrdE.n_rows == 1);
    assert (mdrdE.n_cols == 6);
    
    assert (mdrdT.n_rows == 1);
    assert (mdrdT.n_cols == 1);
    
    Q = mQ;
    r = mr;
    r_in = mr_in;
    sigma_in = msigma_in;
    sigma_in_start = msigma_in_start;
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
    sigma_in = sigma_in_start;
    Wm = Wm_start;
    Wt = Wt_start;
}

//-------------------------------------------------------------
void state_variables_T::set_start(const int &corate_type)
//-------------------------------------------------------------
{
    state_variables::set_start(corate_type);
    sigma_in_start = sigma_in;    
    Wm_start = Wm;
    Wt_start = Wt;
}
    
//----------------------------------------------------------------------
state_variables_T& state_variables_T::rotate_l2g(const state_variables_T& sv, const double &psi, const double &theta, const double &phi)
//----------------------------------------------------------------------
{

    state_variables::rotate_l2g(sv, psi, theta, phi);

    dSdE = sv.dSdE;
    dSdEt = sv.dSdEt;
    dSdT = sv.dSdT;
    drdE = sv.drdE;
    drdT = sv.drdT;
    Q = sv.Q;
    r = sv.r;
    r_in = sv.r_in;
    sigma_in = sv.sigma_in;
    sigma_in_start = sv.sigma_in_start;
    Wm = sv.Wm;
    Wt = sv.Wt;
    Wm_start = sv.Wm_start;
    Wt_start = sv.Wt_start;

    Rotation rot = Rotation::from_euler(psi, theta, phi, "zxz");
    if (!rot.is_identity()) {
        sigma_in = rot.apply_stress(sigma_in);
        sigma_in_start = rot.apply_stress(sigma_in_start);
        dSdE = rot.apply_stiffness(dSdE);
        dSdEt = rot.apply_stiffness(dSdEt);
        vec dSdT_v = rot.apply_stress(vec(dSdT.as_col()));
        dSdT = dSdT_v.as_row();
        vec drdE_v = rot.apply_strain(vec(drdE.as_col()));
        drdE = drdE_v.as_row();
    }

	return *this;
}

//----------------------------------------------------------------------
state_variables_T& state_variables_T::rotate_g2l(const state_variables_T& sv, const double &psi, const double &theta, const double &phi)
//----------------------------------------------------------------------
{

    state_variables::rotate_g2l(sv, psi, theta, phi);

    sigma_in = sv.sigma_in;
    sigma_in_start = sv.sigma_in_start;

    dSdE = sv.dSdE;
    dSdEt = sv.dSdEt;
    dSdT = sv.dSdT;
    drdE = sv.drdE;
    drdT = sv.drdT;
    Q = sv.Q;
    r = sv.r;
    Wm = sv.Wm;
    Wt = sv.Wt;
    Wm_start = sv.Wm_start;
    Wt_start = sv.Wt_start;

    Rotation rot = Rotation::from_euler(psi, theta, phi, "zxz").inv();
    if (!rot.is_identity()) {
        sigma_in = rot.apply_stress(sigma_in);
        sigma_in_start = rot.apply_stress(sigma_in_start);
        dSdE = rot.apply_stiffness(dSdE);
        dSdEt = rot.apply_stiffness(dSdEt);
        vec dSdT_v = rot.apply_stress(vec(dSdT.as_col()));
        dSdT = dSdT_v.as_row();
        vec drdE_v = rot.apply_strain(vec(drdE.as_col()));
        drdE = drdE_v.as_row();
    }

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
    s << "sigma_in: \n" << sv.sigma_in << "\n";
    s << "sigma_in_start: \n" << sv.sigma_in_start << "\n";
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
