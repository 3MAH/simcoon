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

///@file ellipsoid_multi.hpp
///@brief characteristics of an ellipsoid for homogenization purposes
///@version 1.0

#pragma once

#include <iostream>
#include <string>
#include <armadillo>
#include <simcoon/Simulation/Geometry/ellipsoid.hpp>
#include <simcoon/Continuum_mechanics/Homogenization/phase_multi.hpp>

namespace simcoon{

/**
 * @file ellipsoid_multi.hpp
 * @brief Ellipsoidal inclusion homogenization functions.
 */

/** @addtogroup homogenization
 *  @{
 */


//======================================
class ellipsoid_multi : public phase_multi
//======================================
{
private:
    
protected:
    
	public :
    
    arma::mat S_loc;
    arma::mat P_loc;
    arma::mat T_loc;
    arma::mat T;

    arma::mat T_in_loc;
    arma::mat T_in;
    
    static int mp;
    static int np;
    static arma::vec x;
    static arma::vec wx;
    static arma::vec y;
    static arma::vec wy;
    
    ellipsoid_multi(); //default constructor
    ellipsoid_multi(const arma::mat&, const arma::mat&, const arma::mat&, const arma::mat&, const arma::vec&, const arma::mat&, const arma::mat&, const arma::mat&, const arma::mat&, const arma::mat&, const arma::mat&); //Constructor with parameters
    ellipsoid_multi(const ellipsoid_multi&);	//Copy constructor
    ~ellipsoid_multi();
    
    virtual void fillS_loc(const arma::mat&, const ellipsoid &); //need the L_global of the matrix
    virtual void fillP_loc(const arma::mat&, const ellipsoid &); //need the L_global of the matrix
    virtual void fillT(const arma::mat&, const arma::mat&, const ellipsoid &); //need the L_global of the matrix
    virtual void fillT_iso(const arma::mat&, const arma::mat&, const ellipsoid &); //need the L_global
    virtual void fillT_mec_in(const arma::mat&, const arma::mat&, const ellipsoid &); //need the L_global
    
//    virtual void l2g_T();
//    virtual void g2l_T();
    
	virtual ellipsoid_multi& operator = (const ellipsoid_multi&);
    
    friend std::ostream& operator << (std::ostream&, const ellipsoid_multi&);
    
};


/** @} */ // end of homogenization group

} //namespace simcoon
