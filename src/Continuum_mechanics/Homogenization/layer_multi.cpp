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

///@file layer_multi.cpp
///@brief Micromechanical characteristics of a layer
///@version 1.0

#include <iostream>
#include <fstream>
#include <string>
#include <assert.h>
#include <armadillo>
#include <simcoon/Continuum_Mechanics/Homogenization/layer_multi.hpp>

using namespace std;
using namespace arma;

namespace simcoon{

//=====Private methods for layer_multi===================================

//=====Public methods for layer_multi====================================

/*!
  \brief default constructor
*/

//-------------------------------------------------------------
layer_multi::layer_multi() : phase_multi(), Dnn(3,3), Dnt(3,3), dXn(3,3), dXt(3,3), sigma_hat(3), dzdx1(3)
//-------------------------------------------------------------
{

}

/*!
  \brief Constructor with parameters
*/

//-------------------------------------------------------------
layer_multi::layer_multi(const mat &mA, const mat &mA_start, const mat &mB, const mat &mB_start, const vec &mA_in, const mat &mDnn, const mat &mDnt, const mat &mdXn, const mat &mdXt, const vec &msigma_hat, const vec &mdzdx1) : phase_multi(mA, mA_start, mB, mB_start, mA_in), Dnn(3,3), Dnt(3,3), dXn(3,3), dXt(3,3), sigma_hat(3), dzdx1(3)
//-------------------------------------------------------------
{
    Dnn = mDnn;
    Dnt = mDnt;
    dXn = mdXn;
    dXt = mdXt;
    sigma_hat = msigma_hat;
    dzdx1 = mdzdx1;
}

/*!
  \brief Copy constructor
  \param s phase_multi object to duplicate
*/
    
//------------------------------------------------------
layer_multi::layer_multi(const layer_multi& pc) : phase_multi(pc)
//------------------------------------------------------
{
    Dnn = pc.Dnn;
    Dnt = pc.Dnt;
    dXn = pc.dXn;
    dXt = pc.dXt;
    sigma_hat = pc.sigma_hat;
    dzdx1 = pc.dzdx1;
}

/*!
  \brief Destructor

  Deletes phase_multi (the arma::mat).
*/

//-------------------------------------
layer_multi::~layer_multi() {}
//-------------------------------------

/*!
  \brief Standard operator = for phase_multi
*/

//----------------------------------------------------------------------
layer_multi& layer_multi::operator = (const layer_multi& pc)
//----------------------------------------------------------------------
{
    A = pc.A;
    B = pc.B;
    
    A_start = pc.A_start;
    B_start = pc.B_start;
    
    Dnn = pc.Dnn;
    Dnt = pc.Dnt;
    dXn = pc.dXn;
    dXt = pc.dXt;
    sigma_hat = pc.sigma_hat;
    dzdx1 = pc.dzdx1;
    
    return *this;
}
    
//--------------------------------------------------------------------------
ostream& operator << (ostream& s, const layer_multi& pc)
//--------------------------------------------------------------------------
{
	s << "Display layer multi:\n";
	s << "Display strain concentration tensor:\n";
    s << pc.A;
    s << "Display stress concentration tensor:\n";
    s << pc.B;

    s << "Display normal, tangent part of the tangent modulus:\n";
    s << pc.Dnn << "\n" << pc.Dnt << "\n";
    s << "Display dXn, dXt\n";
    s << pc.dXn << "\n" << pc.dXt << "\n";
    
    s << "\n\n";

	return s;
}

} //namespace simcoon
