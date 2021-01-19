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

///@file natural_basis.cpp
///@brief State variables of a phase, in a defined coordinate system:
///@version 1.0

#include <iostream>
#include <fstream>
#include <assert.h>
#include <armadillo>
#include <simcoon/Continuum_mechanics/Functions/natural_basis.hpp>

using namespace std;
using namespace arma;

namespace simcoon{

//Fonctions that utilise natural basis objects
//=====Private methods for state_variables===================================

//=====Public methods for state_variables============================================

/*!
  \brief default constructor
*/

//-------------------------------------------------------------
natural_basis::natural_basis()
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
}

/*!
  \brief Constructor with parameters
  \n\n
  \f$ \textbf{Examples :} \f$ \n
*/

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
    
    g0ij = inv(g_ij);
    
    for (unsigned int i=0; i < 3; i++) {
        for (unsigned int j=0; j<3; j++) {
            g0i[i] = g0ij(i,j)*g_i[j];
        }
    }
}

/*!
  \brief Copy constructor
  \param s state_variables object to duplicate
*/

//------------------------------------------------------
natural_basis::natural_basis(const natural_basis &nb)
//------------------------------------------------------
{

    g_i = nb.g_i;  // Covariant Vectors
    g0i = nb.g0i;  // Contravariant Vectors
        
    g_ij = nb.g_ij; // Covariant components of the metric tensor
    g0ij = nb.g0ij;
}

/*!
  \brief Destructor

  Nothing to delete, classical objects
*/

//-------------------------------------
natural_basis::~natural_basis()
//-------------------------------------
{

}

/*!
  \brief Standard operator = for natural_basis
*/

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
    
    g0ij = inv(g_ij);
    
    for (unsigned int i=0; i<3; i++) {
        for (unsigned int j=0; j<3; j++) {
            g0i[i] = g0ij(i,j)*g_i[j];
        }
    }
}

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
    
    g0ij = inv(g_ij);
    
    for (unsigned int i=0; i<3; i++) {
        for (unsigned int j=0; j<3; j++) {
            g0i[i] = g0ij(i,j)*g_i[j];
        }
    }
}

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
