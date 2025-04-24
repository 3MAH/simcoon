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

///@file schemes.cpp
///@brief micromechanical schemes for non-linear N-phases heterogeneous materials:
///@brief Mori-Tanaka scheme
///@version 1.0

#include <iostream>
#include <fstream>
#include <sstream>
#include <assert.h>
#include <armadillo>
#include <memory>
#include <simcoon/Simulation/Maths/rotation.hpp>
#include <simcoon/Simulation/Phase/phase_characteristics.hpp>
#include <simcoon/Simulation/Phase/state_variables.hpp>
#include <simcoon/Simulation/Phase/state_variables_M.hpp>
#include <simcoon/Simulation/Geometry/geometry.hpp>
#include <simcoon/Simulation/Geometry/layer.hpp>
#include <simcoon/Simulation/Geometry/ellipsoid.hpp>
// #include <simcoon/Simulation/Geometry/cylinder.hpp>
#include <simcoon/Continuum_mechanics/Homogenization/phase_multi.hpp>
#include <simcoon/Continuum_mechanics/Homogenization/layer_multi.hpp>
#include <simcoon/Continuum_mechanics/Homogenization/ellipsoid_multi.hpp>
// #include <smartplus/Libraries/Homogenization/cylinder_multi.hpp>
#include <simcoon/Continuum_mechanics/Homogenization/eshelby.hpp>

using namespace std;
using namespace arma;

namespace simcoon{
  
void Lt_Homogeneous_E(phase_characteristics &phase) {
    
    //Compute the strain concentration tensor A
    for(auto r : phase.sub_phases) {
        r.sptr_multi->A = eye(6,6);
    }
}

void DE_Homogeneous_E(phase_characteristics &phase) {
    
    std::shared_ptr<state_variables_M> sv_r;
    std::shared_ptr<state_variables_M> sv_eff = std::dynamic_pointer_cast<state_variables_M>(phase.sptr_sv_local);
    
    //Compute the strain concentration tensor A
    for(auto r : phase.sub_phases) {
        r.sptr_multi->A = eye(6,6);
        sv_r->DEtot = r.sptr_multi->A*sv_eff->DEtot; //Recall that the global coordinates of subphases is the local coordinates of the generic phase
    }
}

void Lt_Mori_Tanaka(phase_characteristics &phase, const int &n_matrix) {
    
    mat sumT = zeros(6,6);
    mat inv_sumT = zeros(6,6);

    std::shared_ptr<ellipsoid_multi> elli_multi;
    std::shared_ptr<ellipsoid> elli;
    //ptr on the matrix properties
    std::shared_ptr<state_variables_M> sv_0 = std::dynamic_pointer_cast<state_variables_M>(phase.sub_phases[n_matrix].sptr_sv_global);
    std::shared_ptr<state_variables_M> sv_eff = std::dynamic_pointer_cast<state_variables_M>(phase.sptr_sv_local);
    std::shared_ptr<state_variables_M> sv_r;
    
    //Compute the Eshelby tensor and the interaction tensor for each phase
    for(auto r : phase.sub_phases) {
        elli_multi = std::dynamic_pointer_cast<ellipsoid_multi>(r.sptr_multi);
        elli = std::dynamic_pointer_cast<ellipsoid>(r.sptr_shape);
        sv_r = std::dynamic_pointer_cast<state_variables_M>(r.sptr_sv_global);
        
        //Note The tangent modulus are turned in the coordinate system of the ellipspoid in the fillT function
        //Note The tangent modulus are turned in the coordinate system of the ellipspoid in the fillT function
        if (r.sptr_matprops->number == n_matrix)
            elli_multi->T = eye(6,6);
        else
            elli_multi->fillT(sv_0->Lt, sv_r->Lt, *elli);

        //Compute the normalization interaction tensir sumT
		sumT += elli->concentration*elli_multi->T;
    }

    try {
        inv_sumT = inv(sumT);
      } catch (const std::runtime_error &e) {
        cerr << "Error in inv: " << e.what() << endl;
        throw simcoon::exception_inv("Error in inv function inside Lt_Mori_Tanaka.");
      }    

    //Compute the strain concentration tensor A
    for(auto r : phase.sub_phases) {
        elli_multi = std::dynamic_pointer_cast<ellipsoid_multi>(r.sptr_multi);
        elli_multi->A = elli_multi->T*inv_sumT;
    }
}

void Lt_Mori_Tanaka_iso(phase_characteristics &phase, const int &n_matrix) {
    
    mat sumT = zeros(6,6);
    mat inv_sumT = zeros(6,6);
    
    std::shared_ptr<ellipsoid_multi> elli_multi;
    std::shared_ptr<ellipsoid> elli;
    //ptr on the matrix properties
    std::shared_ptr<state_variables_M> sv_0 = std::dynamic_pointer_cast<state_variables_M>(phase.sub_phases[n_matrix].sptr_sv_global);
    std::shared_ptr<state_variables_M> sv_eff = std::dynamic_pointer_cast<state_variables_M>(phase.sptr_sv_local);
    std::shared_ptr<state_variables_M> sv_r;
    
    //Compute the Eshelby tensor and the interaction tensor for each phase
    for(auto r : phase.sub_phases) {
        elli_multi = std::dynamic_pointer_cast<ellipsoid_multi>(r.sptr_multi);
        elli = std::dynamic_pointer_cast<ellipsoid>(r.sptr_shape);
        sv_r = std::dynamic_pointer_cast<state_variables_M>(r.sptr_sv_global);
        
        //Note The tangent modulus are turned in the coordinate system of the ellipspoid in the fillT function
        //Note The tangent modulus are turned in the coordinate system of the ellipspoid in the fillT function
        if (r.sptr_matprops->number == n_matrix)
            elli_multi->T = eye(6,6);
        else
            elli_multi->fillT_iso(sv_0->Lt, sv_r->Lt, *elli);
        
        //Compute the normalization interaction tensir sumT
        sumT += elli->concentration*elli_multi->T;
    }

    try {
        inv_sumT = inv(sumT);
      } catch (const std::runtime_error &e) {
        cerr << "Error in inv: " << e.what() << endl;
        throw simcoon::exception_inv("Error in inv function inside Lt_Mori_Tanaka_iso.");
      }        
    
    //Compute the strain concentration tensor A
    for(auto r : phase.sub_phases) {
        elli_multi = std::dynamic_pointer_cast<ellipsoid_multi>(r.sptr_multi);
        elli_multi->A = elli_multi->T*inv_sumT;
    }
}
    
void DE_Mori_Tanaka(phase_characteristics &phase, const int &n_matrix) {
    
    mat sumT = zeros(6,6);
    mat inv_sumT = zeros(6,6);
    
    std::shared_ptr<ellipsoid_multi> elli_multi;
    std::shared_ptr<ellipsoid> elli;
    //ptr on the matrix properties
    std::shared_ptr<state_variables_M> sv_0 = std::dynamic_pointer_cast<state_variables_M>(phase.sub_phases[n_matrix].sptr_sv_global);
    std::shared_ptr<state_variables_M> sv_eff = std::dynamic_pointer_cast<state_variables_M>(phase.sptr_sv_local);
    std::shared_ptr<state_variables_M> sv_r;
    
    //Compute the Eshelby tensor and the interaction tensor for each phase
    for(auto r : phase.sub_phases) {
        elli_multi = std::dynamic_pointer_cast<ellipsoid_multi>(r.sptr_multi);
        elli = std::dynamic_pointer_cast<ellipsoid>(r.sptr_shape);
        sv_r = std::dynamic_pointer_cast<state_variables_M>(r.sptr_sv_global);
        
        //Note The tangent modulus are turned in the coordinate system of the ellipspoid in the fillT function
        if (r.sptr_matprops->number == n_matrix)
            elli_multi->T = eye(6,6);
        else
            elli_multi->fillT(sv_0->Lt, sv_r->Lt, *elli);

        //Compute the normalization interaction tensir sumT
        sumT += elli->concentration*elli_multi->T;
    }
    
    try {
        inv_sumT = inv(sumT);
      } catch (const std::runtime_error &e) {
        cerr << "Error in inv: " << e.what() << endl;
        throw simcoon::exception_inv("Error in inv function inside DE_Mori_Tanaka.");
      }   
    
    //Compute the strain concentration tensor A
    for(auto r : phase.sub_phases) {
        elli_multi = std::dynamic_pointer_cast<ellipsoid_multi>(r.sptr_multi);
        sv_r = std::dynamic_pointer_cast<state_variables_M>(r.sptr_sv_global);
        
        elli_multi->A = elli_multi->T*inv_sumT;
        sv_r->DEtot = elli_multi->A*sv_eff->DEtot; //Recall that the global coordinates of subphases is the local coordinates of the generic phase
    }
}

void DE_Mori_Tanaka_iso(phase_characteristics &phase, const int &n_matrix) {
    
    mat sumT = zeros(6,6);
    mat inv_sumT = zeros(6,6);
    
    std::shared_ptr<ellipsoid_multi> elli_multi;
    std::shared_ptr<ellipsoid> elli;
    //ptr on the matrix properties
    std::shared_ptr<state_variables_M> sv_0 = std::dynamic_pointer_cast<state_variables_M>(phase.sub_phases[n_matrix].sptr_sv_global);
    std::shared_ptr<state_variables_M> sv_eff = std::dynamic_pointer_cast<state_variables_M>(phase.sptr_sv_local);
    std::shared_ptr<state_variables_M> sv_r;
    
    //Compute the Eshelby tensor and the interaction tensor for each phase
    for(auto r : phase.sub_phases) {
        elli_multi = std::dynamic_pointer_cast<ellipsoid_multi>(r.sptr_multi);
        elli = std::dynamic_pointer_cast<ellipsoid>(r.sptr_shape);
        sv_r = std::dynamic_pointer_cast<state_variables_M>(r.sptr_sv_global);
        
        //Note The tangent modulus are turned in the coordinate system of the ellipspoid in the fillT function
        if (r.sptr_matprops->number == n_matrix)
            elli_multi->T = eye(6,6);
        else
            elli_multi->fillT_iso(sv_0->Lt, sv_r->Lt, *elli);
        
        //Compute the normalization interaction tensir sumT
        sumT += elli->concentration*elli_multi->T;
    }
    
    try {
        inv_sumT = inv(sumT);
      } catch (const std::runtime_error &e) {
        cerr << "Error in inv: " << e.what() << endl;
        throw simcoon::exception_inv("Error in inv function inside DE_Mori_Tanaka_iso.");
      }   
    
    //Compute the strain concentration tensor A
    for(auto r : phase.sub_phases) {
        elli_multi = std::dynamic_pointer_cast<ellipsoid_multi>(r.sptr_multi);
        sv_r = std::dynamic_pointer_cast<state_variables_M>(r.sptr_sv_global);
        
        elli_multi->A = elli_multi->T*inv_sumT;
        sv_r->DEtot = elli_multi->A*sv_eff->DEtot; //Recall that the global coordinates of subphases is the local coordinates of the generic phase
    }
}
    
void Lt_Self_Consistent(phase_characteristics &phase, const int &n_matrix, const bool &start, const int &option_start) {
    
    std::shared_ptr<ellipsoid_multi> elli_multi;
    std::shared_ptr<ellipsoid> elli;
    //ptr on the matrix properties
    std::shared_ptr<state_variables_M> sv_r;
    std::shared_ptr<state_variables_M> sv_eff = std::dynamic_pointer_cast<state_variables_M>(phase.sptr_sv_local);
    
    //In the self_consistent scheme we need to have the effective tangent modulus first, based on some guessed initial concentration tensor.
    if(start) {
        //ensure that the Lt of phase is actually 0
        if((option_start == 0)&&(n_matrix < 0))
            Lt_Homogeneous_E(phase);
        else if(option_start == 1)
            Lt_Mori_Tanaka(phase, n_matrix);
        else {
            cout << "error , option is not valid for the start option of Self-Consistent scheme (0 : h_E, 1 : MT)";
        }
        
        //Compute the effective tensor from the previous strain localization tensors
        mat Lt_eff = zeros(6,6);
        for(unsigned int i=0; i<phase.sub_phases.size(); i++) {
            elli_multi = std::dynamic_pointer_cast<ellipsoid_multi>(phase.sub_phases[i].sptr_multi);
            elli = std::dynamic_pointer_cast<ellipsoid>(phase.sub_phases[i].sptr_shape);
            sv_r = std::dynamic_pointer_cast<state_variables_M>(phase.sub_phases[i].sptr_sv_global);
            Lt_eff += elli->concentration*elli_multi->A*sv_r->Lt;
        }
        sv_eff->Lt = Lt_eff;
    }
    
    mat sumA = zeros(6,6);
    //Compute the Eshelby tensor and the interaction tensor for each phase
    for(unsigned int i=0; i<phase.sub_phases.size(); i++) {
        elli_multi = std::dynamic_pointer_cast<ellipsoid_multi>(phase.sub_phases[i].sptr_multi);
        elli = std::dynamic_pointer_cast<ellipsoid>(phase.sub_phases[i].sptr_shape);
        sv_r = std::dynamic_pointer_cast<state_variables_M>(phase.sub_phases[i].sptr_sv_global);
        
        //Note The tangent modulus are turned in the coordinate system of the ellipspoid in the fillT function
        if (phase.sub_phases[i].sptr_matprops->number == n_matrix)
            elli_multi->T = eye(6,6);
        else {
            elli_multi->fillT(sv_eff->Lt, sv_r->Lt, *elli);
            sumA += elli->concentration*elli_multi->T;
        }

    }
    
    for(unsigned int i=0; i<phase.sub_phases.size(); i++) {
        //Note The tangent modulus are turned in the coordinate system of the ellipspoid in the fillT function
        elli_multi = std::dynamic_pointer_cast<ellipsoid_multi>(phase.sub_phases[i].sptr_multi);
        elli = std::dynamic_pointer_cast<ellipsoid>(phase.sub_phases[i].sptr_shape);
        sv_r = std::dynamic_pointer_cast<state_variables_M>(phase.sub_phases[i].sptr_sv_global);
        if (phase.sub_phases[i].sptr_matprops->number == n_matrix)
            elli_multi->A = (eye(6,6) - sumA)*(1./elli->concentration);
        else {
            elli_multi->A = elli_multi->T;
        }
    }    
}
    
void DE_Self_Consistent(phase_characteristics &phase, const int &n_matrix, const bool &start, const int &option_start) {
    
    std::shared_ptr<ellipsoid_multi> elli_multi;
    std::shared_ptr<ellipsoid> elli;
    //ptr on the matrix properties
    std::shared_ptr<state_variables_M> sv_r;
    std::shared_ptr<state_variables_M> sv_eff = std::dynamic_pointer_cast<state_variables_M>(phase.sptr_sv_local);

    //In the self_consistent scheme we need to have the effective tangent modulus first, based on some guessed initial concentration tensor.
    if(start) {
        //ensure that the Lt of phase is actually 0
        if((option_start == 0)&&(n_matrix < 0))
            Lt_Homogeneous_E(phase);
        else if(option_start == 1)
            Lt_Mori_Tanaka(phase, n_matrix);
        else {
            cout << "error , option is not valid for the start option of Self-Consistent scheme (0 : h_E, 1 : MT)";
        }
        
        //Compute the effective tensor from the previous strain localization tensors
        mat Lt_eff = zeros(6,6);
        for(auto r : phase.sub_phases) {
            sv_r = std::dynamic_pointer_cast<state_variables_M>(r.sptr_sv_global);
            Lt_eff += r.sptr_shape->concentration*r.sptr_multi->A*sv_r->Lt;
        }
        sv_eff->Lt = Lt_eff;
    }
    
    mat sumA = zeros(6,6);
    //Compute the Eshelby tensor and the interaction tensor for each phase
    for(unsigned int i=0; i<phase.sub_phases.size(); i++) {
        elli_multi = std::dynamic_pointer_cast<ellipsoid_multi>(phase.sub_phases[i].sptr_multi);
        elli = std::dynamic_pointer_cast<ellipsoid>(phase.sub_phases[i].sptr_shape);
        sv_r = std::dynamic_pointer_cast<state_variables_M>(phase.sub_phases[i].sptr_sv_global);
        
        //Note The tangent modulus are turned in the coordinate system of the ellipspoid in the fillT function
        if (phase.sub_phases[i].sptr_matprops->number == n_matrix)
            elli_multi->T = eye(6,6);
        else {
            elli_multi->fillT(sv_eff->Lt, sv_r->Lt, *elli);
            sumA += elli->concentration*elli_multi->T;
        }
        
    }

    for(unsigned int i=0; i<phase.sub_phases.size(); i++) {
        //Note The tangent modulus are turned in the coordinate system of the ellipspoid in the fillT function
        elli_multi = std::dynamic_pointer_cast<ellipsoid_multi>(phase.sub_phases[i].sptr_multi);
        elli = std::dynamic_pointer_cast<ellipsoid>(phase.sub_phases[i].sptr_shape);
        sv_r = std::dynamic_pointer_cast<state_variables_M>(phase.sub_phases[i].sptr_sv_global);
        if (phase.sub_phases[i].sptr_matprops->number == n_matrix)
            elli_multi->A = (eye(6,6) - sumA)*(1./elli->concentration);
        else {
            elli_multi->A = elli_multi->T;
        }

        sv_r->DEtot = elli_multi->A*sv_eff->DEtot; //Recall that the global coordinates of subphases is the local coordinates of the generic phase
    }
}
    
void dE_Periodic_Layer(phase_characteristics &phase, const int &nbiter) {
    
    std::shared_ptr<layer_multi> lay_multi;
    std::shared_ptr<layer> lay;
    //ptr on the matrix properties
    std::shared_ptr<state_variables_M> sv_r;
    std::shared_ptr<state_variables_M> sv_eff = std::dynamic_pointer_cast<state_variables_M>(phase.sptr_sv_local);
    
    mat Lt_loc = zeros(6,6);
    
    if (nbiter == 0) {
        for (auto r : phase.sub_phases) {
            r.sptr_sv_global->DEtot = sv_eff->DEtot;
        }
    }
    
    //Compute the increment of strain
    mat sumDnn = zeros(3,3);
	vec sumcDsig = zeros(3);
    for(auto r : phase.sub_phases) {
        
        sv_r = std::dynamic_pointer_cast<state_variables_M>(r.sptr_sv_global);
        lay_multi = std::dynamic_pointer_cast<layer_multi>(r.sptr_multi);
        lay = std::dynamic_pointer_cast<layer>(r.sptr_shape);
        Lt_loc = rotate_g2l_L(sv_r->Lt, lay->psi_geom, lay->theta_geom, lay->phi_geom);
        
        lay_multi->Dnn(0,0) = Lt_loc(0,0);
        lay_multi->Dnn(0,1) = Lt_loc(0,3);
        lay_multi->Dnn(0,2) = Lt_loc(0,4);
        lay_multi->Dnn(1,0) = Lt_loc(3,0);
        lay_multi->Dnn(1,1) = Lt_loc(3,3);
        lay_multi->Dnn(1,2) = Lt_loc(3,4);
        lay_multi->Dnn(2,0) = Lt_loc(4,0);
        lay_multi->Dnn(2,1) = Lt_loc(4,3);
        lay_multi->Dnn(2,2) = Lt_loc(4,4);
        
        mat sigma_local = rotate_g2l_stress(sv_r->sigma, lay->psi_geom, lay->theta_geom, lay->phi_geom);
        lay_multi->sigma_hat(0) = sigma_local(0);
        lay_multi->sigma_hat(1) = sigma_local(3);
        lay_multi->sigma_hat(2) = sigma_local(4);
    }

    try {
        for(auto r : phase.sub_phases) {
            lay_multi = std::dynamic_pointer_cast<layer_multi>(r.sptr_multi);
            lay = std::dynamic_pointer_cast<layer>(r.sptr_shape);
            sumDnn += lay->concentration*inv(lay_multi->Dnn);
            sumcDsig += lay->concentration*inv(lay_multi->Dnn)*lay_multi->sigma_hat;
        }
    } catch (const std::runtime_error &e) {
        cerr << "Error in inv: " << e.what() << endl;
        throw simcoon::exception_inv("Error in inv function inside dE_Periodic_Layer.");
    }       

    vec m;
    try {
        m = inv(sumDnn)*sumcDsig;
      } catch (const std::runtime_error &e) {
        cerr << "Error in inv: " << e.what() << endl;
        throw simcoon::exception_inv("Error in inv function inside dE_Periodic_Layer.");
      }       
    
      try {
        for(auto r : phase.sub_phases) {
            lay_multi = std::dynamic_pointer_cast<layer_multi>(r.sptr_multi);
            lay = std::dynamic_pointer_cast<layer>(r.sptr_shape);
            sumDnn += lay->concentration*inv(lay_multi->Dnn);
            sumcDsig += lay->concentration*inv(lay_multi->Dnn)*lay_multi->sigma_hat;
        }
      } catch (const std::runtime_error &e) {
        cerr << "Error in inv: " << e.what() << endl;
        throw simcoon::exception_inv("Error in inv function inside dE_Periodic_Layer.");
      }   
      
    for(auto r : phase.sub_phases) {
        
        sv_r = std::dynamic_pointer_cast<state_variables_M>(r.sptr_sv_global);
        lay_multi = std::dynamic_pointer_cast<layer_multi>(r.sptr_multi);

        try {
            lay_multi->dzdx1 = inv(lay_multi->Dnn)*(m-lay_multi->sigma_hat);
        } catch (const std::runtime_error &e) {
            cerr << "Error in inv: " << e.what() << endl;
            throw simcoon::exception_inv("Error in inv function inside dE_Periodic_Layer.");
        }        
        lay = std::dynamic_pointer_cast<layer>(r.sptr_shape);
        
        mat dEtot_local = zeros(6);
        dEtot_local(0) = lay_multi->dzdx1(0);
        dEtot_local(3) = lay_multi->dzdx1(1);
        dEtot_local(4) = lay_multi->dzdx1(2);
        mat dEtot_global = rotate_l2g_strain(dEtot_local, lay->psi_geom, lay->theta_geom, lay->phi_geom);
        
        sv_r->DEtot += dEtot_global;
    }
}
    
void Lt_Periodic_Layer(phase_characteristics &phase) {
    
    std::shared_ptr<layer_multi> lay_multi;
    std::shared_ptr<layer> lay;
    //ptr on the matrix properties
    std::shared_ptr<state_variables_M> sv_r;
    
    mat Lt_loc = zeros(6,6);
    mat A_loc = zeros(6,6);
    
    //Compute the strain concentration tensor A
    for(auto r : phase.sub_phases) {
        
        sv_r = std::dynamic_pointer_cast<state_variables_M>(r.sptr_sv_global);
        lay_multi = std::dynamic_pointer_cast<layer_multi>(r.sptr_multi);
        lay = std::dynamic_pointer_cast<layer>(r.sptr_shape);
        Lt_loc = rotate_g2l_L(sv_r->Lt, lay->psi_geom, lay->theta_geom, lay->phi_geom);
        
        lay_multi->Dnn(0,0) = Lt_loc(0,0);
        lay_multi->Dnn(0,1) = Lt_loc(0,3);
        lay_multi->Dnn(0,2) = Lt_loc(0,4);
        lay_multi->Dnn(1,0) = Lt_loc(3,0);
        lay_multi->Dnn(1,1) = Lt_loc(3,3);
        lay_multi->Dnn(1,2) = Lt_loc(3,4);
        lay_multi->Dnn(2,0) = Lt_loc(4,0);
        lay_multi->Dnn(2,1) = Lt_loc(4,3);
        lay_multi->Dnn(2,2) = Lt_loc(4,4);
        
        lay_multi->Dnt(0,0) = Lt_loc(0,1);
        lay_multi->Dnt(0,1) = Lt_loc(0,2);
        lay_multi->Dnt(0,2) = Lt_loc(0,5);
        lay_multi->Dnt(1,0) = Lt_loc(3,1);
        lay_multi->Dnt(1,1) = Lt_loc(3,2);
        lay_multi->Dnt(1,2) = Lt_loc(3,5);
        lay_multi->Dnt(2,0) = Lt_loc(4,1);
        lay_multi->Dnt(2,1) = Lt_loc(4,2);
        lay_multi->Dnt(2,2) = Lt_loc(4,5);
    }
    
    mat sumDnn = zeros(3,3);
    mat sumDnt = zeros(3,3);

    try {
        for(auto r : phase.sub_phases) {
            lay_multi = std::dynamic_pointer_cast<layer_multi>(r.sptr_multi);
            lay = std::dynamic_pointer_cast<layer>(r.sptr_shape);
            sumDnn += lay->concentration*inv(lay_multi->Dnn);
            sumDnt += lay->concentration*inv(lay_multi->Dnn)*lay_multi->Dnt;
        }
    } catch (const std::runtime_error &e) {
        cerr << "Error in inv: " << e.what() << endl;
        throw simcoon::exception_inv("Error in inv function inside Lt_Periodic_Layer.");
    }    

    mat m_n;
    try {
        m_n = inv(sumDnn);
    } catch (const std::runtime_error &e) {
        cerr << "Error in inv: " << e.what() << endl;
        throw simcoon::exception_inv("Error in inv function inside Lt_Periodic_Layer.");
    }       

    mat m_t = m_n*sumDnt;

    for(auto r : phase.sub_phases) {
        lay_multi = std::dynamic_pointer_cast<layer_multi>(r.sptr_multi);
        lay = std::dynamic_pointer_cast<layer>(r.sptr_shape);   

        try {
            lay_multi->dXn = inv(lay_multi->Dnn)*(m_n-lay_multi->Dnn);
            lay_multi->dXt = inv(lay_multi->Dnn)*(m_t-lay_multi->Dnt);
        } catch (const std::runtime_error &e) {
            cerr << "Error in inv: " << e.what() << endl;
            throw simcoon::exception_inv("Error in inv function inside Lt_Periodic_Layer.");
        }           
        
        A_loc = eye(6,6);
        
        A_loc(0,0) += lay_multi->dXn(0,0);
        A_loc(0,1) += lay_multi->dXt(0,0);
        A_loc(0,2) += lay_multi->dXt(0,1);
        A_loc(0,3) += lay_multi->dXn(0,1);
        A_loc(0,4) += lay_multi->dXn(0,2);
        A_loc(0,5) += lay_multi->dXt(0,2);
        
        A_loc(3,0) += lay_multi->dXn(1,0);
        A_loc(3,1) += lay_multi->dXt(1,0);
        A_loc(3,2) += lay_multi->dXt(1,1);
        A_loc(3,3) += lay_multi->dXn(1,1);
        A_loc(3,4) += lay_multi->dXn(1,2);
        A_loc(3,5) += lay_multi->dXt(1,2);
        
        A_loc(4,0) += lay_multi->dXn(2,0);
        A_loc(4,1) += lay_multi->dXt(2,0);
        A_loc(4,2) += lay_multi->dXt(2,1);
        A_loc(4,3) += lay_multi->dXn(2,1);
        A_loc(4,4) += lay_multi->dXn(2,2);
        A_loc(4,5) += lay_multi->dXt(2,2);

        lay_multi->A = rotate_l2g_A(A_loc, lay->psi_geom, lay->theta_geom, lay->phi_geom);
    }
    
}
    
    
} //namespace simcoon
