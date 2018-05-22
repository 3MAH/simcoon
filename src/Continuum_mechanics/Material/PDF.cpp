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
#include <simcoon/Continuum_mechanics/Material/PDF.hpp>


using namespace std;
using namespace arma;

namespace simcoon{
    
    //=====Private methods for ODF===================================
    
    //=====Public methods for ODF============================================
    
/*!
 \brief default constructor
 */

//-------------------------------------------------------------
PDF::PDF() : limits(2)
//-------------------------------------------------------------
{

    Nphases = 0;
    Parameter = 0;
    n_densities = 0;
    norm = 0.;
}

//-------------------------------------------------------------
PDF::PDF(const int &nParameter, const double &param_min, const double &param_max) : limits(2)
//-------------------------------------------------------------
{
    
    Nphases = 0;
    Parameter = nParameter; //Index of the material parameter to change
    n_densities = 0;
    norm = 0.;
    limits(0) = param_min;
    limits(1) = param_max;
    
}
    
/*!
 \brief Constructor with parameters
 \f$ \textbf{Examples :} \f$ \n
 */
        
//-------------------------------------------------------------
PDF::PDF(const int &mNphases, const int &mParameter, const int &mn_densities, const std::vector<peak> &mpeaks, const double &mnorm, const vec &mdensities, const vec &mlimits) : limits(2)
//-------------------------------------------------------------
{
    Nphases = mNphases;
    Parameter = mParameter;
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
PDF::PDF(const PDF& pc)
//------------------------------------------------------
{
    Nphases = pc.Nphases;
    Parameter = pc.Parameter;
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
PDF::~PDF() {}
//-------------------------------------
    
/*!
 \brief Standard operator = for ODF
 */

//----------------------------------------------------------------------
void PDF::construct(const int &npeaks)
//----------------------------------------------------------------------
{
    for(int i=0; i<npeaks; i++) {
        peaks.push_back(peak());
    }

}
    
//----------------------------------------------------------------------
double PDF::density(const double &alpha)
//----------------------------------------------------------------------
{
    double density = 0.;
    
    for(auto p : peaks) {
        density += p.get_density_PDF(alpha);
    }
    return density;
}
    
    
//----------------------------------------------------------------------
PDF& PDF::operator = (const PDF& pc)
//----------------------------------------------------------------------
{
    Nphases = pc.Nphases;
    Parameter = pc.Parameter;
    n_densities = pc.n_densities;
    
    peaks = pc.peaks;
    norm = pc.norm;
    
    densities = pc.densities;
    limits = pc.limits;
    
    return *this;
}

//--------------------------------------------------------------------------
ostream& operator << (ostream& s, const PDF& pc)
//--------------------------------------------------------------------------
{
    s << "Display PDF characteristics:\n";
    s << "Number of phases Nphases:\t" << pc.Nphases << "\n";
    s << "Index of parameter:\t" << pc.Parameter << "\n";

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
