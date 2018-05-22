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

///@file peak.cpp
///@brief peak charcteristics from a density function:
///@version 0.9

#include <iostream>
#include <fstream>
#include <string>
#include <assert.h>
#include <armadillo>
#include <memory>
#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Material/peak.hpp>
#include <simcoon/Simulation/Maths/stats.hpp>

using namespace std;
using namespace arma;

namespace simcoon{
    
//=====Private methods for peak===================================

//=====Public methods for peak============================================

/*!
 \brief default constructor
 */
    
//-------------------------------------------------------------
peak::peak()
//-------------------------------------------------------------
{
    number = 0;
    method = 0;
    mean = 0.;
    s_dev = 0.;
    width = 0.;
    ampl = 0.;
}
    
/*!
 \brief Constructor with parameters
 \f$ \textbf{Examples :} \f$ \n
 */
    
//-------------------------------------------------------------
peak::peak(const int &mnumber, const int &mmethod, const double &mmean, const double &ms_dev, const double &mwidth, const double &mampl, const vec &mparams, const vec &mdensities)
//-------------------------------------------------------------
{
    number = mnumber;
    method = mmethod;
    mean = mmean;
    s_dev = ms_dev;
    width = mwidth;
    ampl = mampl;
    params = mparams;
    densities = mdensities;
}
    
/*!
 \brief Copy constructor
 \param s phase_characteristics object to duplicate
 */

//------------------------------------------------------
peak::peak(const peak &pc)
//------------------------------------------------------
{
    number = pc.number;
    method = pc.method;
    mean = pc.mean;
    s_dev = pc.s_dev;
    width = pc.width;
    ampl = pc.ampl;
    params = pc.params;
    densities = pc.densities;
}

/*!
 \brief Destructor
 
 Deletes peak
 */
    
//-------------------------------------
peak::~peak() {}
//-------------------------------------
    
/*!
 \brief Standard operator = for peak
 */
    
//-------------------------------------------------------------
double peak::get_density_ODF(const double &theta)
//-------------------------------------------------------------
{
    
    switch (method) {
        case 1: {
            return ODF_sd(theta, mean, params);
            break;
        }
        case 2: {
            return ODF_hard(theta, mean, s_dev, ampl) + ODF_hard(theta - pi, mean, s_dev, ampl) + ODF_hard(theta + pi, mean, s_dev, ampl);
            break;
        }
        case 3: {
            return Gaussian(theta, mean, s_dev, ampl) + Gaussian(theta - pi, mean, s_dev, ampl) + Gaussian(theta + pi, mean, s_dev, ampl);
            break;
        }
        case 4: {
            return Lorentzian(theta, mean, width, ampl) + Lorentzian(theta - pi, mean, width, ampl) + Lorentzian(theta + pi, mean, width, ampl);
            break;
        }
        case 5: {
            return PseudoVoigt(theta, mean, s_dev, width, ampl, params) + PseudoVoigt(theta - pi, mean, s_dev, width, ampl, params) + PseudoVoigt(theta + pi, mean, s_dev, width, ampl, params);
            break;
        }
        case 6: {
            assert(width > 0.);
            double inv_width = 1./width;
            return Pearson7(theta, mean, inv_width, params) + Pearson7(theta - pi, mean, inv_width, params) + Pearson7(theta + pi, mean, inv_width, params);
            break;
        }
        case 7: {
            return 1.;
            break;
        }
        default : {
            cout << "Error: The peak type specified is not recognized" << endl;
            return 0;
            break;            
        }
    }
}

//-------------------------------------------------------------
double peak::get_density_PDF(const double &theta)
//-------------------------------------------------------------
{
    
    switch (method) {
        case 1: {
            return ODF_sd(theta, mean, params);
            break;
        }
        case 2: {
            return ODF_hard(theta, mean, s_dev, ampl);
            break;
        }
        case 3: {
            return Gaussian(theta, mean, s_dev, ampl);
            break;
        }
        case 4: {
            return Lorentzian(theta, mean, width, ampl);
            break;
        }
        case 5: {
            return PseudoVoigt(theta, mean, s_dev, width, ampl, params);
            break;
        }
        case 6: {
            assert(width > 0.);
            double inv_width = 1./width;
            return Pearson7(theta, mean, inv_width, params);
            break;
        }
        case 7: {
            return 1.;
            break;
        }
        default : {
            cout << "Error: The peak type specified is not recognized" << endl;
            return 0;
            break;            
        }
    }
}
    
//----------------------------------------------------------------------
peak& peak::operator = (const peak& pc)
//----------------------------------------------------------------------
{
    number = pc.number;
    method = pc.method;
    mean = pc.mean;
    s_dev = pc.s_dev;
    width = pc.width;
    ampl = pc.ampl;
    params = pc.params;
    densities = pc.densities;
    
    return *this;
}
    
    
//--------------------------------------------------------------------------
ostream& operator << (ostream& s, const peak& pc)
//--------------------------------------------------------------------------
{
    s << "Display peak characteristics:\n";
    s << "number: \t" << pc.number << "\n";
    s << "method: \t" << pc.method << "\n";
    s << "mean = \t" << pc.mean << "\n";
    s << "standard deviation = \t" << pc.s_dev << "\n";
    s << "width = \t" << pc.width << "\n";
    s << "params = \t" << pc.params.t() << "\n";
    
    s << "\n\n";
    
    return s;
}
    
} //namespace simcoon
