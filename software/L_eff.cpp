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

///@file L_eff.cpp
///@brief solver: Determination of the effective elastic properties	of a composite
///@version 1.9

#include <iostream>
#include <fstream>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Umat/umat_L_elastic.hpp>
#include <simcoon/Continuum_mechanics/Functions/natural_basis.hpp>
#include <simcoon/Simulation/Solver/read.hpp>
#include <simcoon/Simulation/Phase/state_variables_M.hpp>
#include <simcoon/Simulation/Phase/phase_characteristics.hpp>

using namespace std;
using namespace arma;
using namespace simcoon;

ofstream output("L.txt");

int main() {
    
	///Material properties reading, use "material.dat" to specify parameters values
    string umat_name;
    string path_data = "data";
    string materialfile = "material.dat";
	
    unsigned int nprops = 0;
    unsigned int nstatev = 0;
    vec props;
    
    double psi_rve = 0.;
    double theta_rve = 0.;
    double phi_rve = 0.;
    
    double T_init = 273.15;

    read_matprops(umat_name, nprops, props, nstatev, psi_rve, theta_rve, phi_rve, path_data, materialfile);
    phase_characteristics rve;
    
    rve.construct(0,1);
    natural_basis nb;
    rve.sptr_matprops->update(0, umat_name, 1, psi_rve, theta_rve, phi_rve, props.n_elem, props);
    rve.sptr_sv_global->update(zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), eye(3,3), eye(3,3), eye(3,3), eye(3,3), eye(3,3), eye(3,3),T_init, 0., nstatev, zeros(nstatev), zeros(nstatev), nb);
    
    auto sv_M = std::dynamic_pointer_cast<state_variables_M>(rve.sptr_sv_global);

    //Second we call a recursive method that find all the elastic moduli iof the phases
    get_L_elastic(rve);
    output << sv_M->Lt << "\n";
    output.close();
    
	return 0;
}
