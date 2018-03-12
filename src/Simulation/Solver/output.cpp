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

///@file output.cpp
///@brief object that defines the output
///@version 1.0

#include <iostream>
#include <sstream>
#include <fstream>
#include <assert.h>
#include <math.h>
#include <armadillo>
#include <simcoon/Simulation/Solver/output.hpp>

using namespace std;
using namespace arma;

namespace simcoon{

//=====Private methods for ellipsoid_characteristics===================================

//=====Public methods for ellipsoid_characteristics============================================

//@brief default constructor
//-------------------------------------------------------------
solver_output::solver_output()
//-------------------------------------------------------------
{
    o_nb_meca = 0;
    o_nb_T = 0;
    o_nw_statev = 0;
}

/*!
 \brief Constructor with parameters
 \param nblock : number of blocks
 */

//-------------------------------------------------------------
solver_output::solver_output(const int &nblock)
//-------------------------------------------------------------
{
    o_nb_meca = 0;
    o_nb_T = 0;
    o_nw_statev = 0;
    
    o_type.zeros(nblock);
    o_nfreq.zeros(nblock);
    o_tfreq.zeros(nblock);
}

/*!
 \brief Copy constructor
 \param so solver_output object to duplicate
 */

//------------------------------------------------------
solver_output::solver_output(const solver_output& so)
//------------------------------------------------------
{
    
    o_nb_meca = so.o_nb_meca;
    o_meca = so.o_meca;
    o_nb_T = so.o_nb_T;
    o_nw_statev = so.o_nw_statev;
	o_wanted_statev = so.o_wanted_statev;
	o_range_statev = so.o_range_statev;
    o_type = so.o_type;
    o_nfreq = so.o_nfreq;
    o_tfreq = so.o_tfreq;
}

/*!
 \brief destructor
 */

solver_output::~solver_output() {}

/*!
 \brief Standard operator = for solver_output objects
 */

//----------------------------------------------------------------------
solver_output& solver_output::operator = (const solver_output& so)
//----------------------------------------------------------------------
{

    o_nb_meca = so.o_nb_meca;
    o_meca = so.o_meca;
    o_nb_T = so.o_nb_T;
    o_nw_statev = so.o_nw_statev;
	o_wanted_statev = so.o_wanted_statev;
	o_range_statev = so.o_range_statev;
    o_type = so.o_type;
    o_nfreq = so.o_nfreq;
    o_tfreq = so.o_tfreq;
    
	return *this;
}

//--------------------------------------------------------------------------
ostream& operator << (ostream& s, const solver_output& so)
//--------------------------------------------------------------------------
{    
	s << "Display info on the output:\n ";
	s << "Number of meca: " << so.o_nb_meca << "\t heat" << so.o_nb_T << " \n";
	s << "meca\n" << so.o_meca << "\n";
    
    s << "statev:\n";
    s << "number wanted statev = " << so.o_nw_statev << "\n";
    s << "standalone statev\n" << so.o_wanted_statev << "\n";
    s << "range statev\n" << so.o_range_statev << "\n";
    
    s << "block\t type\t every\n";
    for (unsigned int i=0; i<so.o_type.n_elem; i++) {
        s << i+1 << "\t" << so.o_type(i) << "\t";
        if (so.o_type(i) == 1)
            s << so.o_nfreq(i);
        else if (so.o_type(i) ==2)
            s << so.o_tfreq(i);
        else
            s << "N/A";
        s << "\n\n";
        
    }
    
	return s;
}

} //namespace simcoon
