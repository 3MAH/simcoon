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

///@file ellipsoid_multi.cpp
///@brief Micromechanical characteristics of a phase
///@version 1.0

#include <iostream>
#include <fstream>
#include <string>
#include <assert.h>
#include <armadillo>
#include <simcoon/Continuum_Mechanics/Functions/constitutive.hpp>
#include <simcoon/Simulation/Maths/rotation.hpp>
#include <simcoon/Simulation/Geometry/ellipsoid.hpp>
#include <simcoon/Continuum_Mechanics/Homogenization/ellipsoid_multi.hpp>
#include <simcoon/Continuum_Mechanics/Homogenization/eshelby.hpp>

using namespace std;
using namespace arma;

namespace simcoon{

//Definition of the static variables
int ellipsoid_multi::mp;
int ellipsoid_multi::np;
vec ellipsoid_multi::x;
vec ellipsoid_multi::wx;
vec ellipsoid_multi::y;
vec ellipsoid_multi::wy;
    
    
//=====Private methods for ellipsoid_multi===================================

//=====Public methods for ellipsoid_multi====================================

/*!
  \brief default constructor
*/
    
//-------------------------------------------------------------
ellipsoid_multi::ellipsoid_multi() : phase_multi(), S_loc(6,6), P_loc(6,6), T_loc(6,6), T(6,6), T_in_loc(6,6), T_in(6,6)
//-------------------------------------------------------------
{
    //This calls only the constructor of the two matrix A & B
    
}

/*!
  \brief Constructor with parameters
*/

//-------------------------------------------------------------
ellipsoid_multi::ellipsoid_multi(const mat &mA, const mat &mA_start, const mat &mB, const mat &mB_start, const vec &mA_in, const mat &mS_loc, const mat &mP_loc, const mat &mT_loc, const mat &mT, const mat &mT_in_loc, const mat &mT_in) : phase_multi(mA, mA_start, mB, mB_start, mA_in), S_loc(6,6), P_loc(6,6), T_loc(6,6), T(6,6), T_in_loc(6,6), T_in(6,6)
//-------------------------------------------------------------
{
    S_loc = mS_loc;
    P_loc = mP_loc;
    T_loc = mT_loc;
    T = mT;
    T_in_loc = mT_in_loc;
    T_in = mT_in;
}

/*!
  \brief Copy constructor
  \param s ellipsoid_multi object to duplicate
*/
    
//------------------------------------------------------
ellipsoid_multi::ellipsoid_multi(const ellipsoid_multi& pc) : phase_multi(pc), S_loc(6,6), P_loc(6,6), T_loc(6,6), T(6,6), T_in_loc(6,6), T_in(6,6)
//------------------------------------------------------
{
    S_loc = pc.S_loc;
    P_loc = pc.P_loc;
    T_loc = pc.T_loc;
    T = pc.T;
    T_in_loc = pc.T_in_loc;
    T_in = pc.T_in;
}

/*!
  \brief Destructor

  Deletes phase_multi (the arma::mat).
*/

//-------------------------------------
ellipsoid_multi::~ellipsoid_multi() {}
//-------------------------------------

/*!
  \brief Standard operator = for phase_multi
*/

//-------------------------------------
void ellipsoid_multi::fillS_loc(const mat& Lt_m, const ellipsoid &ell)
//-------------------------------------
{
    mat Ltm_local_geom = rotate_g2l_L(Lt_m, ell.psi_geom, ell.theta_geom, ell.phi_geom);
    S_loc = Eshelby(Ltm_local_geom, ell.a1, ell.a2, ell.a3, x, wx, y, wy, mp, np);
}
    
//-------------------------------------
void ellipsoid_multi::fillP_loc(const mat& Lt_m, const ellipsoid &ell)
//-------------------------------------
{
    mat Ltm_local_geom = rotate_g2l_L(Lt_m, ell.psi_geom, ell.theta_geom, ell.phi_geom);
    P_loc = T_II(Ltm_local_geom, ell.a1, ell.a2, ell.a3, x, wx, y, wy, mp, np);
}
    

//-------------------------------------
void ellipsoid_multi::fillT(const mat& Lt_m, const mat& Lt, const ellipsoid &ell)
//This method correspond to the classical Eshelby method
//-------------------------------------
{
    mat Lt_m_local_geom = rotate_g2l_L(Lt_m, ell.psi_geom, ell.theta_geom, ell.phi_geom);
    S_loc = Eshelby(Lt_m_local_geom, ell.a1, ell.a2, ell.a3, x, wx, y, wy, mp, np);
    mat Lt_local_geom = rotate_g2l_L(Lt, ell.psi_geom, ell.theta_geom, ell.phi_geom);
    
    T_loc = inv(eye(6,6) + S_loc*inv(Lt_m_local_geom)*(Lt_local_geom - Lt_m_local_geom));
    
    T = rotate_l2g_A(T_loc, ell.psi_geom, ell.theta_geom, ell.phi_geom);
}

//-------------------------------------
void ellipsoid_multi::fillT_iso(const mat& Lt_m, const mat& Lt, const ellipsoid &ell)
//This method corresponf to the isotropization method
//-------------------------------------
{
    mat Lt_m_iso = Isotropize(Lt_m);
    S_loc = Eshelby(Lt_m_iso, ell.a1, ell.a2, ell.a3, x, wx, y, wy, mp, np);
    mat Lt_local_geom = rotate_g2l_L(Lt, ell.psi_geom, ell.theta_geom, ell.phi_geom);
    
    T_loc = inv(eye(6,6) + S_loc*inv(Lt_m_iso)*(Lt_local_geom - Lt_m_iso));
    
    T = rotate_l2g_A(T_loc, ell.psi_geom, ell.theta_geom, ell.phi_geom);
}

//-------------------------------------
void ellipsoid_multi::fillT_mec_in(const mat& L_m, const mat& L, const ellipsoid &ell)
//This method corresponds to the Eshelby inhomogeneous inhomogeneity - inelasticity problem. It calculates
//the interaction tensors T for the elastic and the inelastic part.
//-------------------------------------
{
    mat L_m_local_geom = rotate_g2l_L(L_m, ell.psi_geom, ell.theta_geom, ell.phi_geom);
    S_loc = Eshelby(L_m_local_geom, ell.a1, ell.a2, ell.a3, x, wx, y, wy, mp, np);
    mat L_local_geom = rotate_g2l_L(L, ell.psi_geom, ell.theta_geom, ell.phi_geom);
    
    T_loc = inv(eye(6,6) + S_loc*inv(L_m_local_geom)*(L_local_geom - L_m_local_geom));
    T = rotate_l2g_A(T_loc, ell.psi_geom, ell.theta_geom, ell.phi_geom);
    
    T_in_loc = (eye(6,6)-T_loc)*inv(L_m_local_geom - L_local_geom);
    T_in = rotate_l2g_M(T_in_loc, ell.psi_geom, ell.theta_geom, ell.phi_geom);
}
    
//----------------------------------------------------------------------
ellipsoid_multi& ellipsoid_multi::operator = (const ellipsoid_multi& pc)
//----------------------------------------------------------------------
{
    A = pc.A;
    B = pc.B;

    A_start = pc.A_start;
    B_start = pc.B_start;
    
    S_loc = pc.S_loc;
    P_loc = pc.P_loc;
    T_loc = pc.T_loc;
    T = pc.T;
    T_in_loc = pc.T_in_loc;
    T_in = pc.T_in;
    
	return *this;
}
    
//--------------------------------------------------------------------------
ostream& operator << (ostream& s, const ellipsoid_multi& pc)
//--------------------------------------------------------------------------
{
	s << "Display phase multi:\n";
	s << "Display strain concentration tensor:\n";
    s << pc.A;
    s << "Display stress concentration tensor:\n";
    s << pc.B;

    s << "Display Eshelby tensor (local coordinates):\n";
    s << pc.S_loc;
    s << "Display Polarization tensor (local coordinates):\n";
    s << pc.P_loc;
    s << "Display Interaction concentration tensor (local coordinates):\n";
    s << pc.T_loc;
    s << "Display Interaction concentration tensor (global coordinates):\n";
    s << pc.T;
    s << "Display Inelastic Interaction concentration tensor (local coordinates):\n";
    s << pc.T_in_loc;
    s << "Display Inelastic Interaction concentration tensor (global coordinates):\n";
    s << pc.T_in;
    
    s << "\n\n";

	return s;
}

} //namespace simcoon
