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

#include <iostream>
#include <fstream>
#include <assert.h>
#include <armadillo>
#include <simcoon/exception.hpp>
#include <simcoon/Continuum_mechanics/Functions/natural_basis.hpp>

using namespace std;
using namespace arma;

namespace simcoon{

// default constructor
//-------------------------------------------------------------
natural_basis::natural_basis()
//-------------------------------------------------------------
{
    g_i.resize(3);  // Covariant Vectors
    g0i.resize(3);  // Contravariant Vectors

    // Default to the undeformed (Cartesian) natural basis: g_i = e_i, metric = identity.
    // This makes U() = R() = I, so any consumer that does not call from_F (small strain,
    // log_R) sees no convection -- a safe, neutral default rather than a singular zero metric.
    for (unsigned int i=0; i<3; i++) {
        g_i[i] = zeros(3);
        g_i[i](i) = 1.;
        g0i[i] = g_i[i];
    }

    g_ij = eye(3,3); // Covariant components of the metric tensor
    g0ij = eye(3,3); // Contravariant components of the metric tensor
}

//Constructor with parameters
//-------------------------------------------------------------
natural_basis::natural_basis(const std::vector<arma::vec> &mg_i)
//-------------------------------------------------------------
{
    g_i.resize(3);  // Covariant Vectors
    g0i.resize(3);  // Contravariant Vectors
    
    for (unsigned int i=0; i<3; i++) {
        g_i[i] = zeros(3);
        g0i[i] = zeros(3);
    }
    
    g_ij = zeros(3,3); // Covariant components of the metric tensor
    g0ij = zeros(3,3); // Contravariant components of the metric tensor
    
    for (unsigned int i=0; i<3; i++) {
        g_i[i] = mg_i[i];
    }

    //
    for (unsigned int i=0; i<3; i++) {
        for (unsigned int j=0; j<3; j++) {
            g_ij(i,j) = sum(g_i[i]%g_i[j]);
        }
    }
    
    try {
        g0ij = inv(g_ij);
    } catch (const std::runtime_error &e) {
        cerr << "Error in inv: " << e.what() << endl;
        throw simcoon::exception_inv("Error in inv function inside natural_basis constructor.");
    }   
    
    for (unsigned int i=0; i < 3; i++) {
        g0i[i] = zeros(3);
        for (unsigned int j=0; j<3; j++) {
            g0i[i] += g0ij(i,j)*g_i[j];   // g^i = g^{ij} g_j  (sum over j)
        }
    }
}

// copy constructor
//------------------------------------------------------
natural_basis::natural_basis(const natural_basis &nb)
//------------------------------------------------------
{

    g_i = nb.g_i;  // Covariant Vectors
    g0i = nb.g0i;  // Contravariant Vectors
        
    g_ij = nb.g_ij; // Covariant components of the metric tensor
    g0ij = nb.g0ij;
}

// destructor
//-------------------------------------
natural_basis::~natural_basis()
//-------------------------------------
{

}

// affectation = operator
//----------------------------------------------------------------------
natural_basis& natural_basis::operator = (const natural_basis& nb)
//----------------------------------------------------------------------
{
    g_i = nb.g_i;  // Covariant Vectors
    g0i = nb.g0i;  // Contravariant Vectors
        
    g_ij = nb.g_ij; // Covariant components of the metric tensor
    g0ij = nb.g0ij;
    
	return *this;
}
    
// update operator
//-------------------------------------------------------------
void natural_basis::update(const std::vector<arma::vec> &mg_i)
//-------------------------------------------------------------
{
    for (unsigned int i=0; i < 3; i++) {
        g_i[i] = mg_i[i];
    }

    //
    for (unsigned int i=0; i<3; i++) {
        for (unsigned int j=0; j<3; j++) {
            g_ij(i,j) = sum(g_i[i]%g_i[j]);
        }
    }
    
    try {
        g0ij = inv(g_ij);
    } catch (const std::runtime_error &e) {
        cerr << "Error in inv: " << e.what() << endl;
        throw simcoon::exception_inv("Error in inv function inside natural_basis update.");
    } 
    
    for (unsigned int i=0; i<3; i++) {
        g0i[i] = zeros(3);
        for (unsigned int j=0; j<3; j++) {
            g0i[i] += g0ij(i,j)*g_i[j];   // g^i = g^{ij} g_j  (sum over j)
        }
    }
}

// update operator from F
//-------------------------------------------------------------
void natural_basis::from_F(const arma::mat &F)
//-------------------------------------------------------------
{

    for (unsigned int i=0; i<3; i++) {
        for (unsigned int j=0; j<3; j++) {
            g_i[i](j) = F(j,i); //each transform/image of the basis vector are the columns of F
        }
    }

    //
    for (unsigned int i=0; i<3; i++) {
        for (unsigned int j=0; j<3; j++) {
            g_ij(i,j) = sum(g_i[i]%g_i[j]);
        }
    }
    
    try {
        g0ij = inv(g_ij);
    } catch (const std::runtime_error &e) {
        cerr << "Error in inv: " << e.what() << endl;
        throw simcoon::exception_inv("Error in inv function inside natural_basis from_F.");
    } 
    
    for (unsigned int i=0; i<3; i++) {
        g0i[i] = zeros(3);
        for (unsigned int j=0; j<3; j++) {
            g0i[i] += g0ij(i,j)*g_i[j];   // g^i = g^{ij} g_j  (sum over j)
        }
    }
}

// transformation gradient F = [g_1 g_2 g_3] (covariant vectors as columns)
//-------------------------------------------------------------
arma::mat natural_basis::F() const
//-------------------------------------------------------------
{
    mat Fmat(3,3);
    for (unsigned int i=0; i<3; i++) {
        Fmat.col(i) = g_i[i];
    }
    return Fmat;
}

// inverse transformation gradient F^{-1} (contravariant vectors g^i as rows)
//-------------------------------------------------------------
arma::mat natural_basis::F_inv() const
//-------------------------------------------------------------
{
    mat Finv(3,3);
    for (unsigned int i=0; i<3; i++) {
        Finv.row(i) = g0i[i].t();   // g^i = F^{-T} e_i is the i-th row of F^{-1}
    }
    return Finv;
}

// right stretch U = (F^T F)^{1/2} = (g_ij)^{1/2}
//-------------------------------------------------------------
arma::mat natural_basis::U() const
//-------------------------------------------------------------
{
    return sqrtmat_sympd(g_ij);
}

// polar rotation R = F U^{-1}
//-------------------------------------------------------------
arma::mat natural_basis::R() const
//-------------------------------------------------------------
{
    mat Uinv;
    if (!inv(Uinv, U())) {
        throw simcoon::exception_inv("Error inverting the right stretch U in natural_basis::R.");
    }
    return F()*Uinv;
}

// left (spatial / Eulerian) stretch V = (F F^T)^{1/2} = b^{1/2} = R U R^T
//-------------------------------------------------------------
arma::mat natural_basis::V() const
//-------------------------------------------------------------
{
    mat Fmat = F();
    return sqrtmat_sympd(Fmat*Fmat.t());
}

// contravariant components of a spatial (stress-like) tensor on the natural basis:
// sigma^{ij} = g^i . sigma . g^j = F^{-1} sigma F^{-T}
//-------------------------------------------------------------
arma::mat natural_basis::contravariant(const arma::mat &sigma) const
//-------------------------------------------------------------
{
    mat Finv = F_inv();
    return Finv*sigma*Finv.t();
}

// covariant components of a spatial (strain-like) tensor on the natural basis:
// eps_{ij} = g_i . eps . g_j = F^T eps F
//-------------------------------------------------------------
arma::mat natural_basis::covariant(const arma::mat &eps) const
//-------------------------------------------------------------
{
    mat Fmat = F();
    return Fmat.t()*eps*Fmat;
}

// stream operator
//--------------------------------------------------------------------------
ostream& operator << (ostream& s, const natural_basis& nb)
//--------------------------------------------------------------------------
{
    
    for (unsigned int i=0; i<3; i++) {
        s << "g_i[" << i << "] = " << nb.g_i[i].t()  << "\n";
    }
    for (unsigned int i=0; i<3; i++) {
        s << "g0i[" << i << "] = " << nb.g0i[i].t()  << "\n";
    }
    
	s << "g_ij: \n" << nb.g_ij << "\n";
	s << "g0ij: \n" << nb.g0ij << "\n";

	return s;
}

} //namespace simcoon
