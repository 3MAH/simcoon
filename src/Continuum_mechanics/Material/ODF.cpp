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

///@file phase_characteristics.cpp
///@brief Characteristics of a phase, which hereditates from:
///- material_characteristics
///@version 1.0

#include <iostream>
#include <fstream>
#include <string>
#include <assert.h>
#include <armadillo>
#include <memory>
#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Material/peak.hpp>
#include <simcoon/Continuum_mechanics/Material/ODF.hpp>


using namespace std;
using namespace arma;

namespace simcoon{
    
    //=====Private methods for ODF===================================
    
    //=====Public methods for ODF============================================
    
/*!
 \brief default constructor
 */

//-------------------------------------------------------------
ODF::ODF() : limits(2)
//-------------------------------------------------------------
{

    Nphases = 0;
    Angle = 0;
    radian = false;
    n_densities = 0;
    norm = 0.;
}

//-------------------------------------------------------------
ODF::ODF(const int &nAngle, const bool &nradian, const double &angle_min, const double &angle_max) : limits(2)
//-------------------------------------------------------------
{
    
    Nphases = 0;
    Angle = nAngle; //0 : Psi, 1: Theta, 2: Phi
    radian = nradian;
    n_densities = 0;
    norm = 0.;
    limits(0) = angle_min;
    limits(1) = angle_max;
    
    if (radian == false) {
        limits(0) *= (pi/180.);
        limits(1) *= (pi/180.);
    }
}
    
/*!
 \brief Constructor with parameters
 \f$ \textbf{Examples :} \f$ \n
 */
        
//-------------------------------------------------------------
ODF::ODF(const int &mNphases, const int &mAngle, const int &mn_densities, const std::vector<peak> &mpeaks, const double &mnorm, const vec &mdensities, const vec &mlimits, const bool &mradian) : limits(2)
//-------------------------------------------------------------
{
    Nphases = mNphases;
    Angle = mAngle;
    radian = mradian;
    n_densities = mn_densities;
    
    peaks = mpeaks;
    norm = mnorm;
    
    densities = mdensities;
    limits = mlimits;
}
    
    /*!
     \brief Copy constructor
     \param s phase_characteristics object to duplicate
     */
    
//------------------------------------------------------
ODF::ODF(const ODF& pc)
//------------------------------------------------------
{
    Nphases = pc.Nphases;
    Angle = pc.Angle;
    radian = pc.radian;
    n_densities = pc.n_densities;
    
    peaks = pc.peaks;
    norm = pc.norm;
    
    densities = pc.densities;
    limits = pc.limits;
}
    
/*!
 \brief Destructor
 
 Deletes ODF, 
 */
    
//-------------------------------------
ODF::~ODF() {}
//-------------------------------------
    
/*!
 \brief Standard operator = for ODF
 */

//----------------------------------------------------------------------
void ODF::construct(const int &npeaks)
//----------------------------------------------------------------------
{
    for(int i=0; i<npeaks; i++) {
        peaks.push_back(peak());
    }

}
    
//----------------------------------------------------------------------
double ODF::density(const double &alpha)
//----------------------------------------------------------------------
{
    double density = 0.;
    
    for(auto p : peaks) {
        density += p.get_density_ODF(alpha);
    }
    return density;
}
    
    
//----------------------------------------------------------------------
ODF& ODF::operator = (const ODF& pc)
//----------------------------------------------------------------------
{
    Nphases = pc.Nphases;
    Angle = pc.Angle;
    radian = pc.radian;
    n_densities = pc.n_densities;
    
    peaks = pc.peaks;
    norm = pc.norm;
    
    densities = pc.densities;
    limits = pc.limits;
    
    return *this;
}

//--------------------------------------------------------------------------
ostream& operator << (ostream& s, const ODF& pc)
//--------------------------------------------------------------------------
{
    s << "Display ODF characteristics:\n";
    s << "Number of phases Nphases:\t" << pc.Nphases << "\n";
    s << "Angle:\t" << pc.Angle << "\n";
    s << "radian:\t" << pc.radian << "\n";

    s << "number of peaks:" << pc.peaks.size() << "\n";    
    s << "peaks:\n";
    
    for(auto p : pc.peaks) {
        s << p << "\n";
    }
//    s << "n_densities:\t" << pc.n_densities << "\n";
    s << "norm:\t" << pc.norm << "\n";
    
    return s;
}
    
} //namespace simcoon
