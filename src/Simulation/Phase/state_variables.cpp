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

///@file state_variables.cpp
///@brief State variables of a phase, in a defined coordinate system:
///@version 1.0

#include <iostream>
#include <fstream>
#include <assert.h>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Simulation/Maths/rotation.hpp>
#include <simcoon/Simulation/Phase/state_variables.hpp>
#include <simcoon/Continuum_mechanics/Functions/stress.hpp>
#include <simcoon/Continuum_mechanics/Functions/transfer.hpp>
#include <simcoon/Continuum_mechanics/Functions/natural_basis.hpp>

using namespace std;
using namespace arma;

namespace simcoon{

//=====Private methods for state_variables===================================

//=====Public methods for state_variables============================================

/*!
  \brief default constructor
*/

//-------------------------------------------------------------
state_variables::state_variables() : Etot(6), DEtot(6), etot(6), Detot(6), PKII(6), PKII_start(6), tau(6), tau_start(6), sigma(6), sigma_start(6), F0(3,3), F1(3,3), R(3,3), DR(3,3)
//-------------------------------------------------------------
{
	Etot = zeros(6);
	DEtot = zeros(6);
    etot = zeros(6);
    Detot = zeros(6);
	PKII = zeros(6);
	PKII_start = zeros(6);
    tau = zeros(6);
    tau_start = zeros(6);
	sigma = zeros(6);
	sigma_start = zeros(6);
    F0 = zeros(3,3);
    F1 = zeros(3,3);
    R = zeros(3,3);
    DR = zeros(3,3);
    T = 0.;
    DT = 0.;
	nstatev=0;
}

/*!
  \brief Constructor with parameters
  \n\n
  \f$ \textbf{Examples :} \f$ \n
*/

//-------------------------------------------------------------
state_variables::state_variables(const int &m, const bool &init, const double &value) : Etot(6), DEtot(6), etot(6), Detot(6), PKII(6), PKII_start(6), tau(6), tau_start(6), sigma(6), sigma_start(6), F0(3,3), F1(3,3), R(3,3), DR(3,3)
//-------------------------------------------------------------
{
    
    Etot = zeros(6);
    DEtot = zeros(6);
    etot = zeros(6);
    Detot = zeros(6);
	PKII = zeros(6);
	PKII_start = zeros(6);
    tau = zeros(6);
    tau_start = zeros(6);
    sigma = zeros(6);
    sigma_start = zeros(6);
    F0 = zeros(3,3);
    F1 = zeros(3,3);
    R = zeros(3,3);
    DR = zeros(3,3);
    T = 0.;
    DT = 0.;
    
    assert(m>=0);    
    nstatev = m;
    if (m>0) {
        if (init) {
            statev = value*ones(m);
            statev_start = value*ones(m);
        }
        else {
            statev = zeros(m);
            statev_start = zeros(m);
        }
    }
}
    
//-------------------------------------------------------------
state_variables::state_variables(const vec &mEtot, const vec &mDEtot, const vec &metot, const vec &mDetot, const vec &mPKII, const vec &mPKII_start, const vec &mtau, const vec &mtau_start, const vec &msigma, const vec &msigma_start, const mat &mF0, const mat &mF1, const mat &mR, const mat &mDR, const double &mT, const double &mDT, const int &mnstatev, const vec &mstatev, const vec &mstatev_start, const natural_basis &mnb) : Etot(6), DEtot(6), etot(6), Detot(6), PKII(6), PKII_start(6), tau(6), tau_start(6), sigma(6), sigma_start(6), F0(3,3), F1(3,3), R(3,3), DR(3,3)
//-------------------------------------------------------------
{	
	assert (mEtot.size() == 6);
	assert (mDEtot.size() == 6);
	assert (metot.size() == 6);
	assert (mDetot.size() == 6);
	assert (mPKII.size() == 6);
	assert (mtau.size() == 6);
	assert (msigma.size() == 6);
    assert (mF0.n_rows == 3);
    assert (mF0.n_cols == 3);
    assert (mF1.n_rows == 3);
    assert (mF1.n_cols == 3);
    assert (mR.n_rows == 3);
    assert (mR.n_cols == 3);
    assert (mDR.n_rows == 3);
    assert (mDR.n_cols == 3);
	
	Etot = mEtot;
	DEtot = mDEtot;
	etot = metot;
	Detot = mDetot;
	PKII = mPKII;
	PKII_start = mPKII_start;
	tau = mtau;
	tau_start = mtau_start;
	sigma = msigma;
	sigma_start = msigma_start;
    F0 = mF0;
    F1 = mF1;
    R = mR;
    DR = mDR;
    T = mT;
    DT = mDT;
    
    nstatev = mnstatev;
    statev = mstatev;
    statev_start = mstatev_start;
    
    nb = mnb;
}

/*!
  \brief Copy constructor
  \param s state_variables object to duplicate
*/

//------------------------------------------------------
state_variables::state_variables(const state_variables& sv) : Etot(6), DEtot(6), etot(6), Detot(6), PKII(6), PKII_start(6), tau(6), tau_start(6), sigma(6), sigma_start(6), F0(3,3), F1(3,3), R(3,3), DR(3,3)
//------------------------------------------------------
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
    R = sv.R;
    DR = sv.DR;
    T = sv.T;
    DT = sv.DT;
    
    nstatev = sv.nstatev;
    statev = sv.statev;
    statev_start = sv.statev_start;
    
    nb = sv.nb;
}

/*!
  \brief Destructor

  Deletes statev_variables, the vectors and matrix, and table statev if nstatev is not null.
*/

//-------------------------------------
state_variables::~state_variables()
//-------------------------------------
{

}

/*!
  \brief Standard operator = for state_variables
*/

//----------------------------------------------------------------------
state_variables& state_variables::operator = (const state_variables& sv)
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
    R = sv.R;
    DR = sv.DR;
    T = sv.T;
    DT = sv.DT;

    nstatev = sv.nstatev;
    statev = sv.statev;
    statev_start = sv.statev_start;
    
    nb = sv.nb;
    
	return *this;
}

state_variables& state_variables::copy_fields(const state_variables& sv)
//-------------------------------------------------------------
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
    R = sv.R;
    DR = sv.DR;
    T = sv.T;
    DT = sv.DT;
    
    return *this;
}

    
//-------------------------------------------------------------
void state_variables::resize()
//-------------------------------------------------------------
{
    assert(nstatev > 0);
    statev = zeros(nstatev);
    statev_start = zeros(nstatev);
}

//-------------------------------------------------------------
void state_variables::resize(const int &m, const bool &init, const double &value)
//-------------------------------------------------------------
{
    assert(m>=0);
    nstatev = m;
    if (m>0) {
        if (init) {
            statev = value*ones(m);
            statev_start = value*ones(m);
        }
        else {
            statev = zeros(m);
            statev_start = zeros(m);
        }
    }
}
    
    
//-------------------------------------------------------------
void state_variables::update(const vec &mEtot, const vec &mDEtot, const vec &metot, const vec &mDetot, const vec &mPKII, const vec &mPKII_start, const vec &mtau, const vec &mtau_start, const vec &msigma, const vec &msigma_start, const mat &mF0, const mat &mF1, const mat &mR, const mat &mDR, const double &mT, const double &mDT, const int &mnstatev, const vec &mstatev, const vec &mstatev_start, const natural_basis &mnb)
//-------------------------------------------------------------
{
    assert (mEtot.size() == 6);
    assert (mDEtot.size() == 6);
    assert (metot.size() == 6);
    assert (mDetot.size() == 6);
    assert (mPKII.size() == 6);
    assert (msigma.size() == 6);
    assert (mtau.size() == 6);
    assert (mF0.n_rows == 3);
    assert (mF0.n_cols == 3);
    assert (mF1.n_rows == 3);
    assert (mF1.n_cols == 3);
    assert (mR.n_rows == 3);
    assert (mR.n_cols == 3);
    assert (mDR.n_rows == 3);
    assert (mDR.n_cols == 3);
    
    Etot = mEtot;
    DEtot = mDEtot;
    etot = metot;
    Detot = mDetot;
    PKII = mPKII;
    PKII_start = mPKII_start;
    tau = mtau;
    tau_start = mtau_start;
    sigma = msigma;
    sigma_start = msigma_start;
    F0 = mF0;
    F1 = mF1;
    R = mR;
    DR = mDR;
    T = mT;
    DT = mDT;
    
    nstatev = mnstatev;
    statev = mstatev;
    statev_start = mstatev_start;
    
    nb = mnb;
}
    
//-------------------------------------------------------------
void state_variables::to_start()
//-------------------------------------------------------------
{
    PKII = PKII_start;
    tau = tau_start;
    sigma = sigma_start;
    statev = statev_start;
//    F1 = F0;
}
    
//-------------------------------------------------------------
void state_variables::set_start(const int &corate_type)
//-------------------------------------------------------------
{
    PKII_start = PKII;
    tau_start = rotate_stress(tau,DR);
    sigma_start = rotate_stress(sigma,DR);
    statev_start = statev;
    Etot += DEtot;
    etot = rotate_strain(etot,DR) + Detot;
    T += DT;
    F0 = F1;
//    R = R*DR;
    R = DR*R;
    nb.from_F(F1);
}
    
//----------------------------------------------------------------------
state_variables& state_variables::rotate_l2g(const state_variables& sv, const double &psi, const double &theta, const double &phi)
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
    R = sv.R;
    DR = sv.DR;
    T = sv.T;
    DT = sv.DT;
    
    nstatev = sv.nstatev;
    statev = sv.statev;
    statev_start = sv.statev_start;
    
  	if(fabs(phi) > sim_iota) {
		Etot = rotate_strain(Etot, -phi, axis_phi);
		DEtot = rotate_strain(DEtot, -phi, axis_phi);
		etot = rotate_strain(etot, -phi, axis_phi);
		Detot = rotate_strain(Detot, -phi, axis_phi);
		PKII = rotate_stress(PKII, -phi, axis_phi);
		PKII_start = rotate_stress(PKII_start, -phi, axis_phi);
		tau = rotate_stress(tau, -phi, axis_phi);
		tau_start = rotate_stress(tau_start, -phi, axis_phi);
		sigma = rotate_stress(sigma, -phi, axis_phi);
		sigma_start = rotate_stress(sigma_start, -phi, axis_phi);
        F0 = rotate_mat(F0, -phi, axis_phi);
        F1 = rotate_mat(F1, -phi, axis_phi);
        R = rotate_mat(R, -phi, axis_phi);
        DR = rotate_mat(DR, -phi, axis_phi);
	}
  	if(fabs(theta) > sim_iota) {
		Etot = rotate_strain(Etot, -theta, axis_theta);
		DEtot = rotate_strain(DEtot, -theta, axis_theta);
		etot = rotate_strain(etot, -theta, axis_theta);
		Detot = rotate_strain(Detot, -theta, axis_theta);
		PKII = rotate_stress(PKII, -theta, axis_theta);
		PKII_start = rotate_stress(PKII_start, -theta, axis_theta);
		tau = rotate_stress(tau, -theta, axis_theta);
		tau_start = rotate_stress(tau_start, -theta, axis_theta);
		sigma = rotate_stress(sigma, -theta, axis_theta);
		sigma_start = rotate_stress(sigma_start, -theta, axis_theta);
        F0 = rotate_mat(F0, -theta, axis_theta);
        F1 = rotate_mat(F1, -theta, axis_theta);
        R = rotate_mat(R, -theta, axis_theta);
        DR = rotate_mat(DR, -theta, axis_theta);
	}
	if(fabs(psi) > sim_iota) {
		Etot = rotate_strain(Etot, -psi, axis_psi);
		DEtot = rotate_strain(DEtot, -psi, axis_psi);
		etot = rotate_strain(etot, -psi, axis_psi);
		Detot = rotate_strain(Detot, -psi, axis_psi);
		PKII = rotate_stress(PKII, -psi, axis_psi);
		PKII_start = rotate_stress(PKII_start, -psi, axis_psi);
		tau = rotate_stress(tau, -psi, axis_psi);
		tau_start = rotate_stress(tau_start, -psi, axis_psi);
		sigma = rotate_stress(sigma, -psi, axis_psi);
		sigma_start = rotate_stress(sigma_start, -psi, axis_psi);
        F0 = rotate_mat(F0, -psi, axis_psi);
        F1 = rotate_mat(F1, -psi, axis_psi);
        R = rotate_mat(R, -psi, axis_psi);
        DR = rotate_mat(DR, -psi, axis_psi);
	}
    
	return *this;
}
    
//----------------------------------------------------------------------
state_variables& state_variables::rotate_g2l(const state_variables& sv, const double &psi, const double &theta, const double &phi)
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
    R = sv.R;
    DR = sv.DR;
    T = sv.T;
    DT = sv.DT;
    
    nstatev = sv.nstatev;
    statev = sv.statev;
    statev_start = sv.statev_start;    
    
  	if(fabs(psi) > sim_iota) {
		Etot = rotate_strain(Etot, psi, axis_psi);
		DEtot = rotate_strain(DEtot, psi, axis_psi);
		etot = rotate_strain(etot, psi, axis_psi);
		Detot = rotate_strain(Detot, psi, axis_psi);
		PKII = rotate_stress(PKII, psi, axis_psi);
		PKII_start = rotate_stress(PKII_start, psi, axis_psi);
		tau = rotate_stress(tau, psi, axis_psi);
		tau_start = rotate_stress(tau_start, psi, axis_psi);
		sigma = rotate_stress(sigma, psi, axis_psi);
		sigma_start = rotate_stress(sigma_start, psi, axis_psi);
        F0 = rotate_mat(F0, psi, axis_psi);
        F1 = rotate_mat(F1, psi, axis_psi);
        R = rotate_mat(R, psi, axis_psi);
        DR = rotate_mat(DR, psi, axis_psi);
	}
	if(fabs(theta) > sim_iota) {
		Etot = rotate_strain(Etot, theta, axis_theta);
		DEtot = rotate_strain(DEtot, theta, axis_theta);
		etot = rotate_strain(etot, theta, axis_theta);
		Detot = rotate_strain(Detot, theta, axis_theta);
		PKII = rotate_stress(PKII, theta, axis_theta);
		PKII_start = rotate_stress(PKII_start, theta, axis_theta);
		tau = rotate_stress(tau, theta, axis_theta);
		tau_start = rotate_stress(tau_start, theta, axis_theta);
		sigma = rotate_stress(sigma, theta, axis_theta);
		sigma_start = rotate_stress(sigma_start, theta, axis_theta);
        F0 = rotate_mat(F0, theta, axis_theta);
        F1 = rotate_mat(F1, theta, axis_theta);
        R = rotate_mat(R, theta, axis_theta);
        DR = rotate_mat(DR, theta, axis_theta);
	}
	if(fabs(phi) > sim_iota) {
		Etot = rotate_strain(Etot, phi, axis_phi);
		DEtot = rotate_strain(DEtot, phi, axis_phi);
		etot = rotate_strain(etot, phi, axis_phi);
		Detot = rotate_strain(Detot, phi, axis_phi);
		PKII = rotate_stress(PKII, phi, axis_phi);
		PKII_start = rotate_stress(PKII_start, phi, axis_phi);
		tau = rotate_stress(tau, phi, axis_phi);
		tau_start = rotate_stress(tau_start, phi, axis_phi);
		sigma = rotate_stress(sigma, phi, axis_phi);
		sigma_start = rotate_stress(sigma_start, phi, axis_phi);
        F0 = rotate_mat(F0, phi, axis_phi);
        F1 = rotate_mat(F1, phi, axis_phi);
        R = rotate_mat(R, phi, axis_phi);
        DR = rotate_mat(DR, phi, axis_phi);
    }
    
	return *this;
}

//--------------------------------------------------------------------------
ostream& operator << (ostream& s, const state_variables& sv)
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
    s << "R: \n" << sv.R << "\n";
    s << "DR: \n" << sv.DR << "\n";
    s << "T: \n" << sv.T << "\n";
    s << "DT: \n" << sv.DT << "\n";
    
    s << "nstatev: \n" << sv.nstatev << "\n";
    if (sv.nstatev) {
        s << "statev: \n";
        s << sv.statev.t();
        s << "\n";
        s << "statev_start: \n";
        s << sv.statev_start.t();
        s << "\n";
    }

    s << "natural_basis: \n" << sv.nb << "\n";
    
	return s;
}

} //namespace simcoon
