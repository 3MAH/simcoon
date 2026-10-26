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
#include <cstring>
#include <vector>
#include <armadillo>
#include <simcoon/exception.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Simulation/Maths/rotation.hpp>
#include <simcoon/Simulation/Geometry/ellipsoid.hpp>
#include <simcoon/Continuum_mechanics/Homogenization/ellipsoid_multi.hpp>
#include <simcoon/Continuum_mechanics/Homogenization/eshelby.hpp>

using namespace std;
using namespace arma;

namespace simcoon{

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
const mat* eshelby_memo::find(const mat &Lt_local, const double &a1, const double &a2, const double &a3) const
//-------------------------------------
{
    for (const entry &e : entries) {
        // bitwise key: same medium in the same frame, same semi-axes
        if (e.a1 == a1 && e.a2 == a2 && e.a3 == a3
            && std::memcmp(e.Lt.memptr(), Lt_local.memptr(), 36 * sizeof(double)) == 0)
            return &e.S;
    }
    return nullptr;
}

//-------------------------------------
const mat& eshelby_memo::store(const mat &Lt_local, const double &a1, const double &a2, const double &a3, const mat &S)
//-------------------------------------
{
    entries.push_back(entry{arma::mat::fixed<6,6>(Lt_local), a1, a2, a3, S});
    return entries.back().S;
}

//-------------------------------------
void eshelby_memo::clear()
//-------------------------------------
{
    entries.clear();
}

namespace {

// Eshelby tensors integrated by this thread, kept across calls: the pass memo and the
// inclusion's own memo live only as long as a pass or a phase object, which every L_eff and
// every solve rebuilds. Exact key (bitwise medium, semi-axes, quadrature tables), so a hit is
// the tensor a fresh integration gives; a ring of 32 entries.
struct eshelby_cache_entry {
    arma::mat::fixed<6,6> Lt;
    double a1, a2, a3;
    const eshelby_quadrature *quad;
    arma::mat::fixed<6,6> S;
};
const std::size_t eshelby_cache_size = 32;
thread_local std::vector<eshelby_cache_entry> eshelby_cache;
thread_local std::size_t eshelby_cache_next = 0;

const arma::mat::fixed<6,6> *eshelby_cache_find(const mat &Lt_local, const ellipsoid &ell, const eshelby_quadrature &quad)
{
    for (const eshelby_cache_entry &e : eshelby_cache) {
        if (e.quad == &quad && e.a1 == ell.a1 && e.a2 == ell.a2 && e.a3 == ell.a3
            && std::memcmp(e.Lt.memptr(), Lt_local.memptr(), 36 * sizeof(double)) == 0)
            return &e.S;
    }
    return nullptr;
}

void eshelby_cache_store(const mat &Lt_local, const ellipsoid &ell, const eshelby_quadrature &quad, const mat &S)
{
    eshelby_cache_entry e{arma::mat::fixed<6,6>(Lt_local), ell.a1, ell.a2, ell.a3, &quad, arma::mat::fixed<6,6>(S)};
    if (eshelby_cache.size() < eshelby_cache_size) {
        eshelby_cache.push_back(e);
    } else {
        eshelby_cache[eshelby_cache_next] = e;
        eshelby_cache_next = (eshelby_cache_next + 1) % eshelby_cache_size;
    }
}

// True when L is isotropic up to round-off: every entry within 1e-14 (relative to the largest)
// of the isotropic tensor of lambda = L(0,1), mu = L(3,3). A rotation would give it back but
// for round-off.
bool isotropic_up_to_roundoff(const mat &L)
{
    const double lambda = L(0,1);
    const double mu = L(3,3);
    const double tol = 1.E-14*abs(L).max();
    for (int j=0; j<6; j++) {
        for (int i=0; i<6; i++) {
            double L_iso_ij = 0.;
            if (i < 3 && j < 3)
                L_iso_ij = (i == j) ? lambda + 2.*mu : lambda;
            else if (i == j)
                L_iso_ij = mu;
            if (fabs(L(i,j) - L_iso_ij) > tol)
                return false;
        }
    }
    return true;
}

}  // namespace

// Eshelby tensor of ell in the medium Lt_local (inclusion frame): the inclusion's own last
// tensor when the medium has not changed, else the pass memo when given, else the closed form
// when the settings allow it and one applies, else the tensors this thread integrated lately,
// else integrated.
mat ellipsoid_multi::eshelby_memoized(eshelby_memo *memo, const mat &Lt_local, const ellipsoid &ell, const eshelby_quadrature &quad)
{
    if (const mat *own = S_loc_memo.find(Lt_local, ell.a1, ell.a2, ell.a3))
        return *own;
    mat S;
    const mat *hit = (memo != nullptr) ? memo->find(Lt_local, ell.a1, ell.a2, ell.a3) : nullptr;
    if (hit != nullptr) {
        S = *hit;
    } else {
        if (quad.closed_form && Eshelby_closed_form(Lt_local, ell.a1, ell.a2, ell.a3, S)) {
            //isotropic medium, sphere or spheroid: no integration
        } else if (const arma::mat::fixed<6,6> *cached = eshelby_cache_find(Lt_local, ell, quad)) {
            S = *cached;
        } else {
            S = Eshelby(Lt_local, ell.a1, ell.a2, ell.a3, quad.x, quad.wx, quad.y, quad.wy, quad.mp, quad.np);
            eshelby_cache_store(Lt_local, ell, quad, S);
        }
        if (memo != nullptr) memo->store(Lt_local, ell.a1, ell.a2, ell.a3, S);
    }
    S_loc_memo.clear();                                  // one entry: the current medium
    S_loc_memo.store(Lt_local, ell.a1, ell.a2, ell.a3, S);
    return S;
}

//-------------------------------------
void ellipsoid_multi::fillS_loc(const mat& Lt_m, const ellipsoid &ell, const eshelby_quadrature &quad)
//-------------------------------------
{
    Rotation rot = Rotation::from_euler(ell.psi_geom, ell.theta_geom, ell.phi_geom, "zxz");
    mat Ltm_local_geom = rot.apply_stiffness(Lt_m, false);
    if (!(quad.closed_form && Eshelby_closed_form(Ltm_local_geom, ell.a1, ell.a2, ell.a3, S_loc)))
        S_loc = Eshelby(Ltm_local_geom, ell.a1, ell.a2, ell.a3, quad.x, quad.wx, quad.y, quad.wy, quad.mp, quad.np);
}

//-------------------------------------
void ellipsoid_multi::fillP_loc(const mat& Lt_m, const ellipsoid &ell, const eshelby_quadrature &quad)
//-------------------------------------
{
    Rotation rot = Rotation::from_euler(ell.psi_geom, ell.theta_geom, ell.phi_geom, "zxz");
    mat Ltm_local_geom = rot.apply_stiffness(Lt_m, false);
    P_loc = T_II(Ltm_local_geom, ell.a1, ell.a2, ell.a3, quad.x, quad.wx, quad.y, quad.wy, quad.mp, quad.np);
}


//-------------------------------------
void ellipsoid_multi::fillT(const mat& Lt_m, const mat& Lt, const ellipsoid &ell, const eshelby_quadrature &quad, eshelby_memo *memo)
//This method correspond to the classical Eshelby method
//-------------------------------------
{
    const double geometry[6] = {ell.a1, ell.a2, ell.a3, ell.psi_geom, ell.theta_geom, ell.phi_geom};
    if (T_memo.valid && T_memo.quad == &quad
        && std::memcmp(T_memo.geometry, geometry, sizeof(geometry)) == 0
        && std::memcmp(T_memo.Lt_m.memptr(), Lt_m.memptr(), 36 * sizeof(double)) == 0
        && std::memcmp(T_memo.Lt.memptr(), Lt.memptr(), 36 * sizeof(double)) == 0) {
        S_loc = T_memo.S_loc;
        T_loc = T_memo.T_loc;
        T = T_memo.T;
        return;
    }

    //A geometry frame that is the global one needs no rotation (from_euler(0,0,0) is the identity)
    const bool global_frame = (ell.psi_geom == 0. && ell.theta_geom == 0. && ell.phi_geom == 0.);
    Rotation rot = global_frame ? Rotation::identity() : Rotation::from_euler(ell.psi_geom, ell.theta_geom, ell.phi_geom, "zxz");
    //An isotropic medium is not rotated either: the rotation would only add round-off, and the
    //same medium in every frame shares one Eshelby tensor between orientations (memo, cache)
    mat Lt_m_local_geom = (global_frame || isotropic_up_to_roundoff(Lt_m)) ? Lt_m : rot.apply_stiffness(Lt_m, false);
    S_loc = eshelby_memoized(memo, Lt_m_local_geom, ell, quad);
    mat Lt_local_geom = global_frame ? Lt : rot.apply_stiffness(Lt, false);

    try {
      T_loc = inv(eye(6,6) + S_loc*inv(Lt_m_local_geom)*(Lt_local_geom - Lt_m_local_geom));
    } catch (const std::runtime_error &e) {
      cerr << "Error in inv: " << e.what() << endl;
      throw simcoon::exception_inv("Error in inv function inside ellipsoid_multi::fillT.");
    }
    T = global_frame ? T_loc : rot.apply_strain_concentration(T_loc, true);

    T_memo.Lt_m = Lt_m;
    T_memo.Lt = Lt;
    std::memcpy(T_memo.geometry, geometry, sizeof(geometry));
    T_memo.quad = &quad;
    T_memo.S_loc = S_loc;
    T_memo.T_loc = T_loc;
    T_memo.T = T;
    T_memo.valid = true;
}

//-------------------------------------
void ellipsoid_multi::fillT_iso(const mat& Lt_m, const mat& Lt, const ellipsoid &ell, const eshelby_quadrature &quad, eshelby_memo *memo)
//This method corresponf to the isotropization method
//-------------------------------------
{
    Rotation rot = Rotation::from_euler(ell.psi_geom, ell.theta_geom, ell.phi_geom, "zxz");
    mat Lt_m_iso = Isotropize(Lt_m);
    S_loc = eshelby_memoized(memo, Lt_m_iso, ell, quad);
    mat Lt_local_geom = rot.apply_stiffness(Lt, false);

    try {
      T_loc = inv(eye(6,6) + S_loc*inv(Lt_m_iso)*(Lt_local_geom - Lt_m_iso));
    } catch (const std::runtime_error &e) {
      cerr << "Error in inv: " << e.what() << endl;
      throw simcoon::exception_inv("Error in inv function inside ellipsoid_multi::fillT_is.");
    }

    T = rot.apply_strain_concentration(T_loc, true);
}

//-------------------------------------
void ellipsoid_multi::fillT_mec_in(const mat& L_m, const mat& L, const ellipsoid &ell, const eshelby_quadrature &quad, eshelby_memo *memo)
//This method corresponds to the Eshelby inhomogeneous inhomogeneity - inelasticity problem. It calculates
//the interaction tensors T for the elastic and the inelastic part.
//-------------------------------------
{
    Rotation rot = Rotation::from_euler(ell.psi_geom, ell.theta_geom, ell.phi_geom, "zxz");
    mat L_m_local_geom = rot.apply_stiffness(L_m, false);
    S_loc = eshelby_memoized(memo, L_m_local_geom, ell, quad);
    mat L_local_geom = rot.apply_stiffness(L, false);

    try {
      T_loc = inv(eye(6,6) + S_loc*inv(L_m_local_geom)*(L_local_geom - L_m_local_geom));
      T_in_loc = (eye(6,6)-T_loc)*inv(L_m_local_geom - L_local_geom);
    } catch (const std::runtime_error &e) {
      cerr << "Error in inv: " << e.what() << endl;
      throw simcoon::exception_inv("Error in inv function inside ellipsoid_multi::fillT_mec_in.");
    }

    T = rot.apply_strain_concentration(T_loc, true);
    T_in = rot.apply_compliance(T_in_loc, true);
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
