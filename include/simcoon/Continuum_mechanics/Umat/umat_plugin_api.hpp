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

#pragma once
#include <string>
#include <armadillo>
#include <simcoon/Simulation/Phase/phase_characteristics.hpp>

class BOOST_SYMBOL_VISIBLE umat_plugin_aba_api {
public:

//    umat_plugin_aba_api();
    virtual std::string name() const = 0;
    
    virtual void umat_abaqus(simcoon::phase_characteristics &, const arma::mat &, const double &, const double &, const int &, const int &, bool &, const int &, double &) = 0;
    
    virtual ~umat_plugin_aba_api() {}
};

class BOOST_SYMBOL_VISIBLE umat_plugin_ext_api {
public:

//    umat_plugin_ext_api();
    virtual std::string name() const = 0;
    
    virtual void umat_external_M(const arma::vec &, const arma::vec &, arma::vec &, arma::mat &, arma::mat &, arma::vec &, const arma::mat &, const int &, const arma::vec &, const int &, arma::vec &, const double &, const double &,const double &,const double &, double &, double &, double &, double &, const int &, const int &, const bool &, const int &, double &) = 0;
    
    virtual ~umat_plugin_ext_api() {}
};
