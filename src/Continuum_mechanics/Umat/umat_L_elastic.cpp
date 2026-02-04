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

///@file umat_L_elastic.cpp
///@brief elastic properties of composite materials
///@version 1.0

#include <iostream>
#include <map>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Umat/umat_L_elastic.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Simulation/Phase/phase_characteristics.hpp>
#include <simcoon/Simulation/Phase/state_variables_M.hpp>
#include <simcoon/Simulation/Phase/read.hpp>
#include <simcoon/Continuum_mechanics/Homogenization/ellipsoid_multi.hpp>
#include <simcoon/Continuum_mechanics/Homogenization/eshelby.hpp>
#include <simcoon/Continuum_mechanics/Micromechanics/schemes.hpp>

using namespace std;
using namespace arma;

namespace simcoon{
    
void get_L_elastic(phase_characteristics &rve)
{
    
    string path_data = "data";
    string inputfile; //file # that stores the microstructure properties
    
    std::map<string, int> list_umat;
    list_umat = {{"ELISO",1},{"ELIST",2},{"ELORT",3},{"MIHEN",100},{"MIMTN",101},{"MISCN",103},{"MIPLN",104}};
    
    int method = list_umat[rve.sptr_matprops->umat_name];
    
    //first we read the behavior of the phases & we construct the tensors if necessary
    switch (method) {
            
        case 100: case 101: case 103: {
            //Definition of the static vectors x,wx,y,wy
            ellipsoid_multi::mp = rve.sptr_matprops->props(2);
            ellipsoid_multi::np = rve.sptr_matprops->props(3);
            ellipsoid_multi::x.set_size(ellipsoid_multi::mp);
            ellipsoid_multi::wx.set_size(ellipsoid_multi::mp);
            ellipsoid_multi::y.set_size(ellipsoid_multi::np);
            ellipsoid_multi::wy.set_size(ellipsoid_multi::np);
            points(ellipsoid_multi::x, ellipsoid_multi::wx, ellipsoid_multi::y, ellipsoid_multi::wy,ellipsoid_multi::mp, ellipsoid_multi::np);
            
            inputfile = "Nellipsoids" + to_string(int(rve.sptr_matprops->props(1))) + ".dat";
            read_ellipsoid(rve, path_data, inputfile);
            break;
        }
        case 104: {
            inputfile = "Nlayers" + to_string(int(rve.sptr_matprops->props(1))) + ".dat";
            read_layer(rve, path_data, inputfile);
            break;
        }
    }
    
    rve.global2local();
    auto umat_M = std::dynamic_pointer_cast<state_variables_M>(rve.sptr_sv_local);
    shared_ptr<state_variables_M> umat_sub_phases_M; //shared_ptr on state variables
    
    switch (method) {
            
        case 1: {
            double E = rve.sptr_matprops->props(0);
            double nu = rve.sptr_matprops->props(1);
            
            umat_M->Lt = L_iso(E, nu, "Enu");
            break;
        }
        case 2: {
            double axis = rve.sptr_matprops->props(0);
            double EL = rve.sptr_matprops->props(1);
            double ET = rve.sptr_matprops->props(2);
            double nuTL = rve.sptr_matprops->props(3);
            double nuTT = rve.sptr_matprops->props(4);
            double GLT = rve.sptr_matprops->props(5);
            
            umat_M->Lt = L_isotrans(EL, ET, nuTL, nuTT, GLT, axis);
            break;
        }
        case 3: {
            double Ex = rve.sptr_matprops->props(0);
            double Ey = rve.sptr_matprops->props(1);
            double Ez = rve.sptr_matprops->props(2);
            double nuxy = rve.sptr_matprops->props(3);
            double nuxz = rve.sptr_matprops->props(4);
            double nuyz = rve.sptr_matprops->props(5);
            double Gxy = rve.sptr_matprops->props(6);
            double Gxz = rve.sptr_matprops->props(7);
            double Gyz = rve.sptr_matprops->props(8);
            
            umat_M->Lt = L_ortho(Ex,Ey,Ez,nuxy,nuxz,nuyz,Gxy,Gxz,Gyz, "EnuG");
            break;
        }
        case 100: {
            for (unsigned int i=0; i<rve.sub_phases.size(); i++) {
                get_L_elastic(rve.sub_phases[i]);
            }
            Lt_Homogeneous_E(rve);
            break;
        }
        case 101: {
            for (unsigned int i=0; i<rve.sub_phases.size(); i++) {
                get_L_elastic(rve.sub_phases[i]);
            }
            int n_matrix = rve.sptr_matprops->props(4);
            Lt_Mori_Tanaka(rve, n_matrix);
            break;	
        }
        case 103: {
            for (unsigned int i=0; i<rve.sub_phases.size(); i++) {
                get_L_elastic(rve.sub_phases[i]);
            }
            int n_matrix = rve.sptr_matprops->props(4);
            Lt_Self_Consistent(rve, n_matrix, true, 1);
            
            mat Lt_n = zeros(6,6);
            int nbiter=0;
            double error = 1.;
            
            while ((error > simcoon::precision_micro)&&(nbiter <= simcoon::maxiter_micro)) {
                Lt_n = umat_M->Lt;
                for (unsigned int i=0; i<rve.sub_phases.size(); i++) {
                    get_L_elastic(rve.sub_phases[i]);
                }
                Lt_Self_Consistent(rve, n_matrix, false, 1);
                umat_M->Lt = zeros(6,6);
                for (auto r : rve.sub_phases) {
                    umat_sub_phases_M = std::dynamic_pointer_cast<state_variables_M>(r.sptr_sv_global);
                    umat_M->Lt += r.sptr_shape->concentration*(umat_sub_phases_M->Lt*r.sptr_multi->A);
                }
                error = norm(umat_M->Lt - Lt_n,2.);
                nbiter++;
             }
            break;
        }
        case 104: {
            for (unsigned int i=0; i<rve.sub_phases.size(); i++) {
                get_L_elastic(rve.sub_phases[i]);
            }
            Lt_Periodic_Layer(rve);
            break;
        }
        default: {
            cout << "Error: The choice of Cnstitutive model is not purely linear elastic or could not be found in the umat library :" << rve.sptr_matprops->umat_name << "\n";
            return;
        }
    }
    
    switch (method) {
            
	    case 100: case 101: case 103: case 104: {

            // Compute the effective tangent modulus, and the effective stress
            umat_M->Lt = zeros(6,6);
            for (auto r : rve.sub_phases) {
                umat_sub_phases_M = std::dynamic_pointer_cast<state_variables_M>(r.sptr_sv_global);
                umat_M->Lt += r.sptr_shape->concentration*(umat_sub_phases_M->Lt*r.sptr_multi->A);
            }
            break;
        }
        default: break;
    }
    rve.local2global();
    
}
	
} //namespace simcoon
