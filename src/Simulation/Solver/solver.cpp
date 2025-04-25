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

///@file constitutive.hpp
///@brief solver: solve the mechanical thermomechanical equilibrium			//
//	for a homogeneous loading path, allowing repeatable steps
///@version 1.9

#include <iostream>
#include <fstream>
#include <string>
#include <assert.h>
#include <math.h>
#include <memory>
#include <filesystem>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Simulation/Phase/material_characteristics.hpp>
#include <simcoon/Simulation/Phase/state_variables.hpp>
#include <simcoon/Simulation/Phase/state_variables_M.hpp>
#include <simcoon/Simulation/Phase/state_variables_T.hpp>
#include <simcoon/Simulation/Maths/rotation.hpp>
#include <simcoon/Continuum_mechanics/Functions/transfer.hpp>
#include <simcoon/Continuum_mechanics/Functions/kinematics.hpp>
#include <simcoon/Continuum_mechanics/Functions/stress.hpp>
#include <simcoon/Continuum_mechanics/Functions/objective_rates.hpp>
#include <simcoon/Continuum_mechanics/Functions/natural_basis.hpp>
#include <simcoon/Continuum_mechanics/Umat/umat_smart.hpp>
#include <simcoon/Simulation/Solver/read.hpp>
#include <simcoon/Simulation/Solver/block.hpp>
#include <simcoon/Simulation/Solver/step.hpp>
#include <simcoon/Simulation/Solver/step_meca.hpp>
#include <simcoon/Simulation/Solver/step_thermomeca.hpp>

using namespace std;
using namespace arma;

namespace simcoon{

void solver(const string &umat_name, const vec &props, const unsigned int &nstatev, const double &psi_rve, const double &theta_rve, const double &phi_rve, const int &solver_type, const int &corate_type, const double &div_tnew_dt_solver, const double &mul_tnew_dt_solver, const int &miniter_solver, const int &maxiter_solver, const int &inforce_solver, const double &precision_solver, const double &lambda_solver, const std::string &path_data, const std::string &path_results, const std::string &pathfile, const std::string &outputfile) {

    //Check if the required directories exist:
    if(!filesystem::is_directory(path_data)) {
        cout << "error: the folder for the data, " << path_data << ", is not present" << endl;
        return;
    }
    if(!filesystem::is_directory(path_results)) {
        cout << "The folder for the results, " << path_results << ", is not present and has been created" << endl;
        filesystem::create_directory(path_results);
    }
    
    std::string ext_filename = outputfile.substr(outputfile.length()-4,outputfile.length());
    std::string filename = outputfile.substr(0,outputfile.length()-4); //to remove the extension
    
    std::string outputfile_global = filename + "_global" + ext_filename;
    std::string outputfile_local = filename + "_local" + ext_filename;
    
    std::string output_info_file = "output.dat";
    
	///Usefull UMAT variables
	int ndi = 3;
	int nshr = 3;    
    std::vector<block> blocks;  //loading blocks
    phase_characteristics rve;  // Representative volume element
    
    unsigned int size_meca = 0; //6 for small perturbation, 9 for finite deformation
	bool start = true;
	double Time = 0.;
	double DTime = 0.;
    double T_init = 0.;
    double tnew_dt = 1.;
    
    mat C = zeros(6,6); //Stiffness dS/dE
    mat c = zeros(6,6); //stifness dtau/deps
    mat DR = eye(3,3);
    mat R = eye(3,3);
    
//    mat dSdE = zeros(6,6);
//    mat dSdT = zeros(1,6);
    mat dQdE = zeros(6,1);
    mat dQdT = zeros(1,1);
    
    //read the material properties
    //Read the loading path
    read_path(blocks, T_init, path_data, pathfile);
    
    ///Material properties reading, use "material.dat" to specify parameters values
    rve.sptr_matprops->update(0, umat_name, 1, psi_rve, theta_rve, phi_rve, props.n_elem, props);
    
    //Output
    int o_ncount = 0;
    double o_tcount = 0.;
    
    solver_output so(blocks.size());
    read_output(so, blocks.size(), nstatev, path_data, output_info_file);
    
    //Check output and step files
    check_path_output(blocks, so);

    double error = 0.;
    vec residual;
    vec Delta;
    int nK = 0; // The size of the problem to solve
    mat K;
    mat invK;
    int compteur = 0.;
    
    int inc = 0.;
    double tinc=0.;
    double Dtinc=0.;
    double Dtinc_cur=0.;
    double q_conv = 0.;        //q_conv parameter for 0D convexion, Q_conv = qconv (T-T_init), with q_conv = rho*c_p\tau, tau being a time constant for convexion thermal mechanical conditions
    
    /// Block loop
    for(unsigned int i = 0 ; i < blocks.size() ; i++){

        switch(blocks[i].type) {
            case 1: { //Mechanical
                
                /// resize the problem to solve
                residual = zeros(6);
                Delta = zeros(6);
                K = zeros(6,6);
                invK = zeros(6,6);

                if(blocks[i].control_type <= 4) {
                    size_meca = 6;
                }
                else {
                    size_meca = 9;
                }
//                if((blocks[i].control_type == 1)||(blocks[i].control_type == 2))
//                    size_meca = 6;
//                else if(blocks[i].control_type == 3)
//                    size_meca = 9;
                
                shared_ptr<state_variables_M> sv_M;
                
                if(start) {
                    rve.construct(0,blocks[i].type);
                    natural_basis nb;
                    rve.sptr_sv_global->update(zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), eye(3,3), eye(3,3), eye(3,3), eye(3,3), eye(3,3), eye(3,3), T_init, 0., nstatev, zeros(nstatev), zeros(nstatev), nb);
                    sv_M = std::dynamic_pointer_cast<state_variables_M>(rve.sptr_sv_global);
                }
                else {
                    //sv_M is reassigned properly
                    sv_M = std::dynamic_pointer_cast<state_variables_M>(rve.sptr_sv_global);
                }
                sv_M->L = zeros(6,6);
                sv_M->Lt = zeros(6,6);
                
                //At start, the rotation increment is null
                DTime = 0.;
                sv_M->DEtot = zeros(6);
                sv_M->DT = 0.;
                
                //Run the umat for the first time in the block. So that we get the proper tangent properties
                run_umat_M(rve, DR, Time, DTime, ndi, nshr, start, solver_type, blocks[i].control_type, tnew_dt);
                
                shared_ptr<step_meca> sptr_meca;
                if(solver_type == 1) {
                    //RNL
                    sptr_meca = std::dynamic_pointer_cast<step_meca>(blocks[0].steps[0]);
                    assert(blocks[i].control_type == 1);
                    sptr_meca->generate(Time, sv_M->Etot, sv_M->sigma, sv_M->T);
                    
                    Lt_2_K(sv_M->Lt, K, sptr_meca->cBC_meca, lambda_solver);
                    
                    //jacobian inversion
                    invK = inv(K);
                }
                else if ((solver_type < 0)||(solver_type > 2)) {
                    cout << "Error, the solver type is not properly defined";
                    return;
                }
                
                if(start) {
                    //Use the number of phases saved to define the files
                    rve.define_output(path_results, outputfile_global, "global");
                    rve.define_output(path_results, outputfile_local, "local");
                    //Write the initial results
//                    rve.output(so, -1, -1, -1, -1, Time, "global");
//                    rve.output(so, -1, -1, -1, -1, Time, "local");
                }
                //Set the start values of sigma_start=sigma and statev_start=statev for all phases
                rve.set_start(corate_type); //DEtot = 0 and DT = 0 and DR = 0 so we can use it safely here
                start = false;
                
                /// Cycle loop
                for(unsigned int n = 0; n < blocks[i].ncycle; n++){
                    
                    /// Step loop
                    for(unsigned int j = 0; j < blocks[i].nstep; j++){
                    
                        sptr_meca = std::dynamic_pointer_cast<step_meca>(blocks[i].steps[j]);
                        if (blocks[i].control_type == 1) {
                            sptr_meca->generate(Time, sv_M->Etot, sv_M->sigma, sv_M->T);
                        }
                        else if (blocks[i].control_type == 2) {
                            sptr_meca->generate(Time, sv_M->Etot, sv_M->PKII, sv_M->T);
                        }
                        else if (blocks[i].control_type == 3) {
                            sptr_meca->generate(Time, sv_M->etot, sv_M->sigma, sv_M->T);
//                            sptr_meca->generate(Time, sv_M->etot, sv_M->tau, sv_M->T);
                        }
                        else if (blocks[i].control_type == 4) {
                            vec Biot_vec = t2v_stress(sv_M->Biot_stress());
                            sptr_meca->generate(Time, t2v_strain(sv_M->U0), Biot_vec, sv_M->T); 
                        }                        
                        else if((blocks[i].control_type == 5)||(blocks[i].control_type == 6)) {
                            sptr_meca->generate_kin(Time, sv_M->F0, sv_M->T);
                        }
                        else {
                            cout << "error in Simulation/Solver/solver.cpp: control_type should be a int value in a range of 1 to 5" << endl;
                            exit(0);
                        }
                    
                        nK = sum(sptr_meca->cBC_meca);
                        
                        inc = 0;
                        while(inc < sptr_meca->ninc) {
                            
                            if(error > precision_solver) {
                                for(int k = 0 ; k < 6 ; k++)
                                {
                                    if(sptr_meca->cBC_meca(k)) {
                                        sptr_meca->mecas(inc,k) -= residual(k);
                                    }
                                }
                            }
                            
                            while (tinc<1.) {
                                
                                sptr_meca->compute_inc(tnew_dt, inc, tinc, Dtinc, Dtinc_cur, inforce_solver);
                                                             
                                if(nK == 0){
                                    
                                    if (blocks[i].control_type == 1) {
                                        sv_M->DEtot = Dtinc*sptr_meca->mecas.row(inc).t();
                                        sv_M->DT = Dtinc*sptr_meca->Ts(inc);
                                        sv_M->DR = eye(3,3);
                                        DTime = Dtinc*sptr_meca->times(inc);
                                    }
                                    else if (blocks[i].control_type == 2) {
                                        sv_M->DEtot = Dtinc*sptr_meca->mecas.row(inc).t();
                                        sv_M->DT = Dtinc*sptr_meca->Ts(inc);
                                        //Application of the Hughes-Winget (1980) algorithm
                                        DTime = Dtinc*sptr_meca->times(inc);
                                        DR = inv(eye(3,3)-0.5*DTime*sptr_meca->BC_w)*(eye(3,3) + 0.5*sptr_meca->BC_w*DTime);
                                        
                                        sv_M->F0 = ER_to_F(v2t_strain(sv_M->Etot), sptr_meca->BC_R);
                                        sv_M->F1 = ER_to_F(v2t_strain(sv_M->Etot + sv_M->DEtot), sptr_meca->BC_R*DR);
                                        
                                        mat D = zeros(3,3);
                                        mat Omega = zeros(3,3);
                                        if(corate_type == 0) {
                                            Jaumann(sv_M->DR, D, Omega, DTime, sv_M->F0, sv_M->F1);
                                        }
                                        if(corate_type == 1) {
                                            Green_Naghdi(sv_M->DR, D, Omega, DTime, sv_M->F0, sv_M->F1);
                                        }
                                        if(corate_type == 2) {
                                            logarithmic(sv_M->DR, D, Omega, DTime, sv_M->F0, sv_M->F1);
                                        }

                                        sv_M->Detot = t2v_strain(Delta_log_strain(D, Omega, DTime));
                                        //mat e_tot_log = t2v_strain(0.5*logmat_sympd(L_Cauchy_Green(sv_M->F1)));
                                        //mat E_dot2 = (1./DTime)*v2t_strain(sv_M->DEtot);
                                    }
                                    else if (blocks[i].control_type == 3) {
                                        sv_M->Detot = Dtinc*sptr_meca->mecas.row(inc).t();
                                        sv_M->DT = Dtinc*sptr_meca->Ts(inc);
                                        //Application of the Hughes-Winget (1980) algorithm
                                        DTime = Dtinc*sptr_meca->times(inc);

                                        DR = inv(eye(3,3)-0.5*DTime*sptr_meca->BC_w)*(eye(3,3) + 0.5*sptr_meca->BC_w*DTime);
                                        
                                        sv_M->F0 = eR_to_F(v2t_strain(sv_M->etot), sptr_meca->BC_R);
                                        sv_M->F1 = eR_to_F(v2t_strain(sv_M->etot + sv_M->Detot), sptr_meca->BC_R*DR);

                                        mat D = zeros(3,3);
                                        mat Omega = zeros(3,3);
                                        if(corate_type == 0) {
                                            Jaumann(sv_M->DR, D, Omega, DTime, sv_M->F0, sv_M->F1);
                                        }
                                        if(corate_type == 1) {
                                            Green_Naghdi(sv_M->DR, D, Omega, DTime, sv_M->F0, sv_M->F1);
                                        }
                                        if(corate_type == 2) {
                                            logarithmic(sv_M->DR, D, Omega, DTime, sv_M->F0, sv_M->F1);
                                        }

                                        sv_M->DEtot = t2v_strain(Green_Lagrange(sv_M->F1)) - sv_M->Etot;
                                        
                                        if (DTime > simcoon::iota)
                                            D = sv_M->Detot/DTime;
                                        else
                                            D = zeros(3,3);
                                    }
                                    else if (blocks[i].control_type == 4) {

                                        sv_M->U1 = sv_M->U0 + v2t_strain(Dtinc*sptr_meca->mecas.row(inc).t());
                                        sv_M->DT = Dtinc*sptr_meca->Ts(inc);
                                        //Application of the Hughes-Winget (1980) algorithm

                                        DTime = Dtinc*sptr_meca->times(inc);
                                        DR = inv(eye(3,3)-0.5*DTime*sptr_meca->BC_w)*(eye(3,3) + 0.5*sptr_meca->BC_w*DTime);
                                        
                                        sv_M->F0 = sptr_meca->BC_R*sv_M->U0;
                                        sv_M->F1 = (sptr_meca->BC_R*DR)*(sv_M->U1);
                                        sv_M->DEtot = t2v_strain(Green_Lagrange(sv_M->F1)) - sv_M->Etot;
                                                                                
                                        mat D = zeros(3,3);
                                        mat Omega = zeros(3,3);
                                        if(corate_type == 0) {
                                            Jaumann(sv_M->DR, D, Omega, DTime, sv_M->F0, sv_M->F1);
                                        }
                                        if(corate_type == 1) {
                                            Green_Naghdi(sv_M->DR, D, Omega, DTime, sv_M->F0, sv_M->F1);
                                        }
                                        if(corate_type == 2) {
                                            logarithmic(sv_M->DR, D, Omega, DTime, sv_M->F0, sv_M->F1);
                                        }

                                        sv_M->Detot = t2v_strain(Delta_log_strain(D, Omega, DTime));
                                    }                                    
                                    else {
                                        sv_M->F1 = sv_M->F0 + Dtinc*v2t(sptr_meca->mecas.row(inc).t());
                                        sv_M->DT = Dtinc*sptr_meca->Ts(inc);
                                        DTime = Dtinc*sptr_meca->times(inc);
                                        
                                        mat D = zeros(3,3);
                                        mat Omega = zeros(3,3);
                                        mat Omega2 = zeros(3,3);
                                        mat Omega3 = zeros(3,3);
                                        if(corate_type == 0) {
                                            Jaumann(sv_M->DR, D, Omega, DTime, sv_M->F0, sv_M->F1);
                                            sv_M->Detot = t2v_strain(Delta_log_strain(D, Omega, DTime));
                                        }
                                        if(corate_type == 1) {
                                            Green_Naghdi(sv_M->DR, D, Omega, DTime, sv_M->F0, sv_M->F1);
                                            sv_M->Detot = t2v_strain(Delta_log_strain(D, Omega, DTime));
                                        }
                                        if(corate_type == 2) {
                                            logarithmic(sv_M->DR, D, Omega, DTime, sv_M->F0, sv_M->F1);
                                            sv_M->Detot = t2v_strain(Delta_log_strain(D, Omega, DTime));
                                        }
                                        mat N_1 = zeros(3,3);
                                        mat N_2 = zeros(3,3);
                                        if(corate_type == 3) {
                                            logarithmic_R(sv_M->DR, N_1, N_2, D, Omega, DTime, sv_M->F0, sv_M->F1);
                                            mat I = eye(3,3);
                                            mat DR_N = (inv(I-0.5*DTime*(N_1-N_2)))*(I+0.5*DTime*(N_1-N_2));
                                            
                                            sv_M->Detot = t2v_strain(Delta_log_strain(D, Omega, DTime));
                                            sv_M->etot = rotate_strain(sv_M->etot, DR_N);
                                            sv_M->sigma_start = rotate_stress(sv_M->sigma_start, DR_N);                                            
                                            sv_M->Detot = rotate_strain(sv_M->Detot, DR_N);
                                        }
                                        if(corate_type == 4) {
                                            Truesdell(sv_M->DR, D, Omega, DTime, sv_M->F0, sv_M->F1);
                                            sv_M->Detot = t2v_strain(Delta_log_strain(D, Omega, DTime));
    //                                            log_modified2(sv_M->DR, N_1, N_2, D, Omega, DTime, sv_M->F0, sv_M->F1);
                                        }
                                        if(corate_type == 5) {
                                            logarithmic_F(sv_M->DR, N_1, N_2, D, Omega, DTime, sv_M->F0, sv_M->F1);
                                            mat I = eye(3,3);
                                            mat DR_N = (inv(I-0.5*DTime*(N_1-D)))*(I+0.5*DTime*(N_1-D));
                                            
//                                            cout << "DR_N = \n"  << DR_N << endl;

                                            mat Detot_nat = Delta_log_strain(D, Omega, DTime);
                                            sv_M->etot = t2v_strain(DR_N*v2t_strain(sv_M->etot)*inv(DR_N));
                                            sv_M->sigma_start = t2v_stress(DR_N*v2t_stress(sv_M->sigma_start)*inv(DR_N));                                            
                                            sv_M->Detot = t2v_strain(DR_N*Detot_nat*inv(DR_N));
                                            
/*                                            sv_M->Detot = t2v_strain(Delta_log_strain(D, Omega, DTime));
                                            sv_M->etot = t2v_strain()rotate_strain(sv_M->etot, DR_N);
                                            sv_M->Detot = rotate_strain(sv_M->Detot, DR_N);
*/                                            
                                        }

                                        sv_M->DEtot = t2v_strain(Green_Lagrange(sv_M->F1)) - sv_M->Etot;

                                    }
                                    rve.to_start();
                                    run_umat_M(rve, sv_M->DR, Time, DTime, ndi, nshr, start, solver_type, blocks[i].control_type, tnew_dt);
                                }
                                else{
                                    /// ********************** SOLVING THE MIXED PROBLEM NRSTRUCT ***********************************
                                    ///Saving stress and stress set point at the beginning of the loop
                                
                                    error = 1.;
                                    
                                    if (blocks[i].control_type == 1) {
                                    
                                        sv_M->DEtot = zeros(6);
                                        for(int k = 0 ; k < 6 ; k++)
                                        {
                                            if (sptr_meca->cBC_meca(k)) {
                                                residual(k) = sv_M->sigma(k) - sv_M->sigma_start(k) - Dtinc*sptr_meca->mecas(inc,k);
                                            }
                                            else {
                                                residual(k) = lambda_solver*(sv_M->DEtot(k) - Dtinc*sptr_meca->mecas(inc,k));
                                            }
                                        }
                                    }
                                    else if (blocks[i].control_type == 2) {
                                        sv_M->DEtot = zeros(6);
                                        for(int k = 0 ; k < 6 ; k++)
                                        {
                                            if (sptr_meca->cBC_meca(k)) {
                                                residual(k) = sv_M->PKII(k) - sv_M->PKII_start(k) - Dtinc*sptr_meca->mecas(inc,k);
                                            }
                                            else {
                                                residual(k) = lambda_solver*(sv_M->DEtot(k) - Dtinc*sptr_meca->mecas(inc,k));
                                            }
                                        }
                                    }
                                    else if (blocks[i].control_type == 3) {
                                        sv_M->Detot = zeros(6);
                                        for(int k = 0 ; k < 6 ; k++)
                                        {
                                            if (sptr_meca->cBC_meca(k)) {
//                                                residual(k) = sv_M->tau(k) - sv_M->tau_start(k) - Dtinc*sptr_meca->mecas(inc,k);
                                                residual(k) = sv_M->sigma(k) - sv_M->sigma_start(k) - Dtinc*sptr_meca->mecas(inc,k);
                                            }
                                            else {
                                                residual(k) = lambda_solver*(sv_M->Detot(k) - Dtinc*sptr_meca->mecas(inc,k));
                                            }
                                        }
                                    }
                                    else if (blocks[i].control_type == 4) {
                                        vec Biot_stress = t2v_stress(sv_M->Biot_stress());
                                        vec Biot_stress_start = t2v_stress(sv_M->Biot_stress_start());
                                        vec DU_vec = t2v_strain(sv_M->U1 - sv_M->U0);                                        
                                        for(int k = 0 ; k < 6 ; k++)
                                        {
                                            if (sptr_meca->cBC_meca(k)) {
                                                residual(k) = Biot_stress(k) - Biot_stress_start(k) - Dtinc*sptr_meca->mecas(inc,k);
                                            }
                                            else {
                                                residual(k) = lambda_solver*(DU_vec(k) - Dtinc*sptr_meca->mecas(inc,k));
                                            }
                                        }
                                    }                                    
                                    else {
                                        cout << "error , Those control types are inteded for use in strain-controlled loading only" << endl;
                                        exit(0);
                                    }
                                    while((error > precision_solver)&&(compteur < maxiter_solver)) {

                                        if(solver_type != 1){
                                            // classic
                                            ///Prediction of the strain increment using the tangent modulus given from the umat_ function
                                            //we use the ddsdde (Lt) from the previous increment
                                            if (blocks[i].control_type == 1) {
                                                Lt_2_K(sv_M->Lt, K, sptr_meca->cBC_meca, lambda_solver);
                                            }
                                            else if (blocks[i].control_type == 2) {
                                                C = Dsigma_LieDD_2_DSDE(sv_M->Lt, sv_M->F1);
                                                Lt_2_K(C, K, sptr_meca->cBC_meca, lambda_solver);                                             
                                            }
                                            else if (blocks[i].control_type == 3) {

                                                if(corate_type == 0) {
                                                    C = Dsigma_LieDD_Dsigma_JaumannDD(sv_M->Lt, v2t_stress(sv_M->sigma));
                                                    Lt_2_K(C, K, sptr_meca->cBC_meca, lambda_solver);
                                                }
                                                if(corate_type == 1) {
                                                    C = Dsigma_LieDD_Dsigma_GreenNaghdiDD(sv_M->Lt, sv_M->F1, v2t_stress(sv_M->sigma));
                                                    Lt_2_K(C, K, sptr_meca->cBC_meca, lambda_solver);
                                                }
                                                if(corate_type == 2) {
                                                    C = Dsigma_LieDD_Dsigma_logarithmicDD(sv_M->Lt, sv_M->F1, v2t_stress(sv_M->sigma));
                                                    Lt_2_K(C, K, sptr_meca->cBC_meca, lambda_solver);
                                                }
                                            }
                                            else if (blocks[i].control_type == 4) {
                                                mat DSDE = Dsigma_LieDD_2_DSDE(sv_M->Lt, sv_M->F1);    
                                                mat R = zeros(3,3);
                                                mat U = zeros(3,3);                                                
                                                RU_decomposition(R,U,sv_M->F1);
                                                C = DSDE_DBiotStressDU(DSDE, U, v2t_stress(sv_M->PKII));
                                                Lt_2_K(C, K, sptr_meca->cBC_meca, lambda_solver);
                                            }                                            
                                            
                                            ///jacobian inversion
                                            invK = inv(K);
                                            
                                            /// Prediction of the component of the strain tensor
                                            Delta = -invK * residual;
                                        }
                                        else if(solver_type == 1) {
                                            //RNL
                                            vec sigma_in_red = zeros(6);
                                            for(int k = 0 ; k < 6 ; k++)
                                            {
                                                if (sptr_meca->cBC_meca(k)) {
                                                    sigma_in_red(k) = sv_M->sigma_in(k) - sv_M->sigma_in_start(k);
                                                }
                                                else {
                                                    sigma_in_red(k) = 0.;
                                                }
                                            }
                                            Delta = -invK * residual;
                                        }
                                        
                                        if (blocks[i].control_type == 1) {
                                            sv_M->DR = eye(3,3);
                                            sv_M->DEtot += Delta;
                                            sv_M->DT = Dtinc*sptr_meca->Ts(inc);
                                            DTime = Dtinc*sptr_meca->times(inc);
                                        }
                                        else if (blocks[i].control_type == 2) {
                                        
                                            sv_M->DEtot += Delta;
                                            sv_M->DT = Dtinc*sptr_meca->Ts(inc);
                                            //Application of the Hughes-Winget (1980) algorithm
                                            DTime = Dtinc*sptr_meca->times(inc);
                                            DR = inv(eye(3,3)-0.5*DTime*sptr_meca->BC_w)*(eye(3,3) + 0.5*sptr_meca->BC_w*DTime);
                                            
                                            sv_M->F0 = ER_to_F(v2t_strain(sv_M->Etot), sptr_meca->BC_R);
                                            sv_M->F1 = ER_to_F(v2t_strain(sv_M->Etot + sv_M->DEtot), sptr_meca->BC_R*DR);
                                    
                                            mat D = zeros(3,3);
                                            mat Omega = zeros(3,3);
                                            if(corate_type == 0) {
                                                Jaumann(sv_M->DR, D, Omega, DTime, sv_M->F0, sv_M->F1);
                                            }
                                            if(corate_type == 1) {
                                                Green_Naghdi(sv_M->DR, D, Omega, DTime, sv_M->F0, sv_M->F1);
                                            }
                                            if(corate_type == 2) {
                                                logarithmic(sv_M->DR, D, Omega, DTime, sv_M->F0, sv_M->F1);
                                            }
                                            if(corate_type == 4) {
                                                Truesdell(sv_M->DR, D, Omega, DTime, sv_M->F0, sv_M->F1);
                                            }                                            
                                            sv_M->Detot = t2v_strain(Delta_log_strain(D, Omega, DTime));
                                        }
                                        else if (blocks[i].control_type == 3) {

                                            sv_M->Detot += Delta;
                                            sv_M->DT = Dtinc*sptr_meca->Ts(inc);
                                            //Application of the Hughes-Winget (1980) algorithm
                                            DTime = Dtinc*sptr_meca->times(inc);                                            
                                            DR = inv(eye(3,3)-0.5*DTime*sptr_meca->BC_w)*(eye(3,3) + 0.5*sptr_meca->BC_w*DTime);
                                            
                                            sv_M->F0 = eR_to_F(v2t_strain(sv_M->etot), sptr_meca->BC_R);
                                            sv_M->F1 = eR_to_F(v2t_strain(sv_M->etot + sv_M->Detot), sptr_meca->BC_R*DR);

                                            sv_M->DEtot = t2v_strain(Green_Lagrange(sv_M->F1)) - sv_M->Etot;

                                            mat D = zeros(3,3);
                                            mat Omega = zeros(3,3);
                                            if(corate_type == 0) {
                                                Jaumann(sv_M->DR, D, Omega, DTime, sv_M->F0, sv_M->F1);
                                            }
                                            if(corate_type == 1) {
                                                Green_Naghdi(sv_M->DR, D, Omega, DTime, sv_M->F0, sv_M->F1);
                                            }
                                            if(corate_type == 2) {
                                                logarithmic(sv_M->DR, D, Omega, DTime, sv_M->F0, sv_M->F1);
                                            }
                                            if (DTime > simcoon::iota)
                                                D = sv_M->Detot/DTime;
                                            else
                                                D = zeros(3,3);
                                        }
                                        else if (blocks[i].control_type == 4) {
                                    
                                            sv_M->U1 += v2t_strain(Delta);
                                            sv_M->DT = Dtinc*sptr_meca->Ts(inc);
                                            //Application of the Hughes-Winget (1980) algorithm
                                            DTime = Dtinc*sptr_meca->times(inc);
                                            DR = inv(eye(3,3)-0.5*DTime*sptr_meca->BC_w)*(eye(3,3) + 0.5*sptr_meca->BC_w*DTime);
                                            sv_M->F0 = sptr_meca->BC_R*sv_M->U0;
                                            sv_M->F1 = (DR*sptr_meca->BC_R)*sv_M->U1;
                                            sv_M->DEtot = t2v_strain(Green_Lagrange(sv_M->F1)) - sv_M->Etot;;
                                            mat D = zeros(3,3);
                                            mat Omega = zeros(3,3);
                                            if(corate_type == 0) {
                                                Jaumann(sv_M->DR, D, Omega, DTime, sv_M->F0, sv_M->F1);
                                            }
                                            if(corate_type == 1) {
                                                Green_Naghdi(sv_M->DR, D, Omega, DTime, sv_M->F0, sv_M->F1);
                                            }
                                            if(corate_type == 2) {
                                                logarithmic(sv_M->DR, D, Omega, DTime, sv_M->F0, sv_M->F1);
                                            }
                                            if(corate_type == 4) {
                                                Truesdell(sv_M->DR, D, Omega, DTime, sv_M->F0, sv_M->F1);
                                            }    
                                            sv_M->Detot = t2v_strain(Delta_log_strain(D, Omega, DTime));
                                        }      
                                        rve.to_start();
                                        run_umat_M(rve, sv_M->DR, Time, DTime, ndi, nshr, start, solver_type, blocks[i].control_type, tnew_dt);

                                        if (blocks[i].control_type == 1) {
                                        
                                            //sv_M->DEtot = zeros(6);
                                            for(int k = 0 ; k < 6 ; k++)
                                            {
                                                if (sptr_meca->cBC_meca(k)) {
                                                    residual(k) = sv_M->sigma(k) - sv_M->sigma_start(k) - Dtinc*sptr_meca->mecas(inc,k);
                                                }
                                                else {
                                                    residual(k) = lambda_solver*(sv_M->DEtot(k) - Dtinc*sptr_meca->mecas(inc,k));
                                                }
                                            }
                                        }
                                        else if (blocks[i].control_type == 2) {
                                            for(int k = 0 ; k < 6 ; k++)
                                            {
                                                if (sptr_meca->cBC_meca(k)) {
                                                    residual(k) = sv_M->PKII(k) - sv_M->PKII_start(k) - Dtinc*sptr_meca->mecas(inc,k);
                                                }
                                                else {
                                                    residual(k) = lambda_solver*(sv_M->DEtot(k) - Dtinc*sptr_meca->mecas(inc,k));
                                                }
                                            }
                                        }
                                        else if (blocks[i].control_type == 3) {
                                            //sv_M->DEtot = zeros(6);
                                            for(int k = 0 ; k < 6 ; k++)
                                            {
                                                if (sptr_meca->cBC_meca(k)) {
//                                                    residual(k) = sv_M->tau(k) - sv_M->tau_start(k) - Dtinc*sptr_meca->mecas(inc,k);
                                                    residual(k) = sv_M->sigma(k) - sv_M->sigma_start(k) - Dtinc*sptr_meca->mecas(inc,k);
                                                }
                                                else {
                                                    residual(k) = lambda_solver*(sv_M->Detot(k) - Dtinc*sptr_meca->mecas(inc,k));
                                                }
                                            }
                                        }
                                        else if (blocks[i].control_type == 4) {

                                            vec Biot_stress = t2v_stress(sv_M->Biot_stress());
                                            vec Biot_stress_start = t2v_stress(sv_M->Biot_stress_start());                                        
                                            vec DU_vec = t2v_strain(sv_M->U1 - sv_M->U0);
                                            for(int k = 0 ; k < 6 ; k++)
                                            {
                                                if (sptr_meca->cBC_meca(k)) {
                                                    residual(k) = Biot_stress(k) - Biot_stress_start(k) - Dtinc*sptr_meca->mecas(inc,k);
                                                }
                                                else {
                                                    residual(k) = lambda_solver*(DU_vec(k) - Dtinc*sptr_meca->mecas(inc,k));
                                                }
                                            }
                                        }                                        
                                        compteur++;
                                        error = norm(residual, 2.);
                                        
                                        if(tnew_dt < 1.) {
                                            if((fabs(Dtinc_cur - sptr_meca->Dn_mini) > simcoon::iota)||(inforce_solver == 0)) {
                                                compteur = maxiter_solver;
                                            }
                                        }
                                        
                                    }
                                    
                                }
/*                                if((fabs(Dtinc_cur - sptr_meca->Dn_mini) < simcoon::iota)&&(tnew_dt < 1.)) {
//                                    cout << "The subroutine has required a step reduction lower than the minimal indicated at" << sptr_meca->number << " inc: " << inc << " and fraction:" << tinc << "\n";
                                    //The solver has been inforced!
                                    return;
                                }
                                
                                if((error > 1000.*precision_solver)&&(Dtinc_cur == sptr_meca->Dn_mini)) {
//                                    cout << "The error has exceeded 100 times the precision, the simulation has stopped at " << sptr_meca->number << " inc: " << inc << " and fraction:" << tinc << "\n";
                                    //The solver has been inforced!
                                    return;
                                }*/
                                
                                if(error > precision_solver) {
                                    if(Dtinc_cur == sptr_meca->Dn_mini) {
                                        if(inforce_solver == 1) {
                                            
                                            cout << "The solver has been inforced to proceed (Solver issue) at step:" << sptr_meca->number << " inc: " << inc << " and fraction:" << tinc << ", with the error: " << error << "\n";
//                                            cout << "The next increment has integrated the error to avoid propagation\n";
                                            //The solver has been inforced!
                                            tnew_dt = 1.;
                                            
                                            if (inc+1<sptr_meca->ninc) {
                                                for(int k = 0 ; k < 6 ; k++)
                                                {
                                                    if(sptr_meca->cBC_meca(k)) {
                                                        sptr_meca->mecas(inc+1,k) -= residual(k);
                                                    }
                                                }
                                            }
                                        }
                                        else if (inforce_solver == 2) {
                                            tnew_dt = 1.;
                                            
                                            if (inc+1<sptr_meca->ninc) {
                                                for(int k = 0 ; k < 6 ; k++)
                                                {
                                                    if(sptr_meca->cBC_meca(k)) {
                                                        sptr_meca->mecas(inc+1,k) -= residual(k);
                                                    }
                                                }
                                            }
                                        }
                                        else return;
                                        
                                    }
                                    else {
                                        tnew_dt = div_tnew_dt_solver;
                                    }
                                }
                                
                                if((compteur < miniter_solver)&&(tnew_dt >= 1.)) {
                                    tnew_dt = mul_tnew_dt_solver;
                                }
                                compteur = 0;
                                
                                sptr_meca->assess_inc(tnew_dt, tinc, Dtinc, rve ,Time, DTime, DR, corate_type);
                                //start variables ready for the next increment
                                
                            }
                            
                            //At the end of each increment, check if results should be written
                            if (so.o_type(i) == 1) {
                                o_ncount++;
                            }
                            if (so.o_type(i) == 2) {
                                o_tcount+=DTime;
                            }
                            
                            //Write the results
                            if (((so.o_type(i) == 1)&&(o_ncount == so.o_nfreq(i)))||(((so.o_type(i) == 2)&&(fabs(o_tcount - so.o_tfreq(i)) < 1.E-12)))) {

                                rve.output(so, i, n, j, inc, Time, "global");
                                rve.output(so, i, n, j, inc, Time, "local");
                                
                                if (so.o_type(i) == 1) {
                                    o_ncount = 0;
                                }
                                if (so.o_type(i) == 2) {
                                    o_tcount = 0.;
                                }
                            }
                            
                            tinc = 0.;
                            inc++;
                         }
                                                
                    }
                        
                }
                break;
            }
            case 2: { //Thermomechanical
                
                /// resize the problem to solve
                residual = zeros(7);
                Delta = zeros(7);
                K = zeros(7,7);
                invK = zeros(7,7);
                
                shared_ptr<state_variables_T> sv_T;
                
                if(start) {
                    rve.construct(0,blocks[i].type);
                    natural_basis nb;
                    rve.sptr_sv_global->update(zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), eye(3,3), eye(3,3), eye(3,3), eye(3,3), eye(3,3), eye(3,3), T_init, 0., nstatev, zeros(nstatev), zeros(nstatev), nb);
                    sv_T = std::dynamic_pointer_cast<state_variables_T>(rve.sptr_sv_global);
                }
                else {
                    //Dynamic cast from some other (possible state_variable_M/T)
                    /*sv_M = std::dynamic_pointer_cast<state_variables_M>(rve.sptr_sv_global);
                     rve.construct(0,blocks[i].type);
                     rve.sptr_sv_global->update(sv_M->Etot, sv_M->DEtot, sv_M->sigma, sv_M->sigma_start, sv_M->T, sv_M->DT, sv_M->sse, sv_M->spd, nstatev, sv_M->statev, sv_M->statev_start);*/
                    //sv_M is reassigned properly
                    sv_T = std::dynamic_pointer_cast<state_variables_T>(rve.sptr_sv_global);
                }
                
                sv_T->dSdE = zeros(6,6);
                sv_T->dSdT = zeros(6,1);
                dQdE = zeros(1,6);
                dQdT = zeros(1,1);
                
                DR = eye(3,3);
                DTime = 0.;
                sv_T->DEtot = zeros(6);
                sv_T->DT = 0.;
                
                //Run the umat for the first time in the block. So that we get the proper tangent properties
                run_umat_T(rve, DR, Time, DTime, ndi, nshr, start, solver_type, blocks[i].control_type, tnew_dt);
                
                sv_T->Q = -1.*sv_T->r;    //Since DTime=0;
                dQdT = lambda_solver;  //To avoid any singularity in the system                
                
                shared_ptr<step_thermomeca> sptr_thermomeca;
                if(solver_type == 1) {
                    //RNL
                    sptr_thermomeca = std::dynamic_pointer_cast<step_thermomeca>(blocks[0].steps[0]);
                    sptr_thermomeca->generate(Time, sv_T->Etot, sv_T->sigma, sv_T->T);
                    
                    Lth_2_K(sv_T->dSdE, sv_T->dSdT, dQdE, dQdT, K, sptr_thermomeca->cBC_meca, sptr_thermomeca->cBC_T, lambda_solver);
                    
                    //jacobian inversion
                    invK = inv(K);
                }
                else if ((solver_type < 0)||(solver_type > 2)) {
                    cout << "Error, the solver type is not properly defined";
                    return;
                }
                
                if(start) {
                    //Use the number of phases saved to define the files
                    rve.define_output(path_results, outputfile_global, "global");
                    rve.define_output(path_results, outputfile_local, "local");
                    //Write the initial results
//                    rve.output(so, -1, -1, -1, -1, Time, "global");
//                    rve.output(so, -1, -1, -1, -1, Time, "local");
                }
                //Set the start values of sigma_start=sigma and statev_start=statev for all phases
                rve.set_start(corate_type); //DEtot = 0 and DT = 0 so we can use it safely here
                start = false;
                
                /// Cycle loop
                for(unsigned int n = 0; n < blocks[i].ncycle; n++){
                    
                    /// Step loop
                    for(unsigned int j = 0; j < blocks[i].nstep; j++){
                        
                        
                        shared_ptr<step_thermomeca> sptr_thermomeca = std::dynamic_pointer_cast<step_thermomeca>(blocks[i].steps[j]);
                        sptr_thermomeca->generate(Time, sv_T->Etot, sv_T->sigma, sv_T->T);
                        
                        nK = sum(sptr_thermomeca->cBC_meca);
                        
                        inc = 0;
                        if(sptr_thermomeca->cBC_T == 3)
                            q_conv = sptr_thermomeca->BC_T;
                        
                        while(inc < sptr_thermomeca->ninc) {
                            
                            
                            if(error > precision_solver) {
                                for(int k = 0 ; k < 6 ; k++)
                                {
                                    if (sptr_thermomeca->cBC_meca(k)) {
                                        sptr_thermomeca->mecas(inc,k) -= residual(k);
                                    }
                                }
                                if (sptr_thermomeca->cBC_T) {
                                    sptr_thermomeca->Ts(inc) -= residual(6);
                                }
                            }
                            
                            while (tinc<1.) {
                                
                                sptr_thermomeca->compute_inc(tnew_dt, inc, tinc, Dtinc, Dtinc_cur, inforce_solver);
                                
                                if(nK + sptr_thermomeca->cBC_T == 0){
                                    
                                    sv_T->DEtot = Dtinc*sptr_thermomeca->mecas.row(inc).t();
                                    sv_T->DT = Dtinc*sptr_thermomeca->Ts(inc);
                                    DTime = Dtinc*sptr_thermomeca->times(inc);
                                    
                                    run_umat_T(rve, DR, Time, DTime, ndi, nshr, start, solver_type, blocks[i].control_type, tnew_dt);
                                    sv_T->Q = -1.*sv_T->r;
                                    
                                }
                                else{
                                    /// ********************** SOLVING THE MIXED PROBLEM NRSTRUCT ***********************************
                                    ///Saving stress and stress set point at the beginning of the loop
                                    
                                    error = 1.;
                                    
                                    sv_T->DEtot = zeros(6);
                                    sv_T->DT = 0.;
                                    
                                    //Construction of the initial residual
                                    for(int k = 0 ; k < 6 ; k++)
                                    {
                                        if (sptr_thermomeca->cBC_meca(k)) {
                                            residual(k) = sv_T->sigma(k) - sv_T->sigma_start(k) - Dtinc*sptr_thermomeca->mecas(inc,k);
                                        }
                                        else {
                                            residual(k) = lambda_solver*(sv_T->DEtot(k) - Dtinc*sptr_thermomeca->mecas(inc,k));
                                        }
                                    }
                                    if (sptr_thermomeca->cBC_T == 1) {
                                        residual(6) = sv_T->Q - sptr_thermomeca->Ts(inc);
                                    }
                                    else if(sptr_thermomeca->cBC_T == 0) {
                                        residual(6) = lambda_solver*(sv_T->DT - Dtinc*sptr_thermomeca->Ts(inc));
                                    }
                                    else if(sptr_thermomeca->cBC_T == 3) { //Special case of 0D convexion that depends on temperature assumption
                                        residual(6) = sv_T->Q + q_conv*(sv_T->T-T_init);
                                    }
                                    else {
                                        cout << "error : The Thermal BC is not recognized\n";
                                        return;
                                    }
                                    
                                    while((error > precision_solver)&&(compteur < maxiter_solver)) {
                                        
                                        if(solver_type != 1){
                                            // classic
                                            ///Prediction of the strain increment using the tangent modulus given from the umat_ function
                                            //we use the ddsdde (Lt) from the previous increment
                                            Lth_2_K(sv_T->dSdE, sv_T->dSdT, dQdE, dQdT, K, sptr_thermomeca->cBC_meca, sptr_thermomeca->cBC_T, lambda_solver);
                                            
                                            ///jacobian inversion
                                            invK = inv(K);
                                            
                                            /// Prediction of the component of the strain tensor
                                            Delta = -invK * residual;
                                        }
                                        else if(solver_type == 1) {
                                            //RNL
                                            vec sigma_in_red = zeros(7);
                                            for(int k = 0 ; k < 6 ; k++)
                                            {
                                                if (sptr_thermomeca->cBC_meca(k)) {
                                                    sigma_in_red(k) = sv_T->sigma_in(k) - sv_T->sigma_in_start(k);
                                                }
                                                else {
                                                    sigma_in_red(k) = 0.;
                                                }
                                            }
                                            sigma_in_red(6) = -1.*sv_T->r_in;
                                            Delta = -invK * residual;
                                        }
                                        
                                        for(int k = 0 ; k < 6 ; k++)
                                        {
                                            sv_T->DEtot(k) += Delta(k);
                                        }
                                        sv_T->DT += Delta(6);
                                        DTime = Dtinc*sptr_thermomeca->times(inc);
                                        
                                        rve.to_start();
                                        run_umat_T(rve, DR, Time, DTime, ndi, nshr, start, solver_type, blocks[i].control_type, tnew_dt);
                                        
                                        if (DTime < 1.E-12) {
                                            sv_T->Q = -1.*sv_T->r;    //Since DTime=0;
                                            
                                            dQdE = -1.*sv_T->drdE.t();
                                            dQdT = lambda_solver;  //To avoid any singularity in the system
                     
                                        }
                                        else{
                                            //Attention here, we solve Phi = Q - Q_conv = Q - q_conv*(sv_T->T-T_init)) = Q - (rho*c_p*(1./tau)*(sv_T->T-T_init)) = 0
                                            //So the derivative / T has the extra q_conv = rho*c_p*(1./tau)
                                            //It is actually icluded here in dQdT
                                            sv_T->Q = -1.*sv_T->r;
                                            
                                            dQdE = -sv_T->drdE.t();
                                            
                                            if (sptr_thermomeca->cBC_T < 3) {
                                                dQdT = -1.*sv_T->drdT;
                                            }
                                            else if(sptr_thermomeca->cBC_T == 3) {
                                                dQdT = -1.*sv_T->drdT + q_conv;
                                            }
                                            
                                        }
                                                                                
                                        for(int k = 0 ; k < 6 ; k++)
                                        {
                                            if (sptr_thermomeca->cBC_meca(k)) {
                                                residual(k) = sv_T->sigma(k) - sv_T->sigma_start(k) - Dtinc*sptr_thermomeca->mecas(inc,k);
                                            }
                                            else {
                                                residual(k) = lambda_solver*(sv_T->DEtot(k) - Dtinc*sptr_thermomeca->mecas(inc,k));
                                            }
                                        }
                                        if (sptr_thermomeca->cBC_T == 1) {
                                            residual(6) = sv_T->Q - sptr_thermomeca->Ts(inc);
                                        }
                                        else if(sptr_thermomeca->cBC_T == 0) {
                                            residual(6) = lambda_solver*(sv_T->DT - Dtinc*sptr_thermomeca->Ts(inc));
                                        }
                                        else if(sptr_thermomeca->cBC_T == 3) { //Special case of 0D convexion that depends on temperature assumption
                                            residual(6) = sv_T->Q + q_conv*(sv_T->T-T_init);
                                        }
                                        else {
                                            cout << "error : The Thermal BC is not recognized\n";
                                            return;
                                        }
                                        
                                        compteur++;
                                        error = norm(residual, 2.);
                                        
                                        if(tnew_dt < 1.) {
                                            if((fabs(Dtinc_cur - sptr_thermomeca->Dn_mini) > simcoon::iota)||(inforce_solver == 0)) {
                                                compteur = maxiter_solver;
                                            }
                                        }
                                        
                                    }
                                    
                                }
                                
/*                                if((fabs(Dtinc_cur - sptr_thermomeca->Dn_mini) < simcoon::iota)&&(tnew_dt < 1.)) {
                                    cout << "The subroutine has required a step reduction lower than the minimal indicated at" << sptr_thermomeca->number << " inc: " << inc << " and fraction:" << tinc << "\n";
                                    //The solver has been inforced!
                                    return;
                                }
                                
                                if((error > 1000.*precision_solver)&&(Dtinc_cur == sptr_thermomeca->Dn_mini)) {
                                    cout << "The error has exceeded 1000 times the precision, the simulation has stopped at " << sptr_thermomeca->number << " inc: " << inc << " and fraction:" << tinc << "\n";
                                    //The solver has been inforced!
                                    return;
                                }
                                */
                                
                                if(error > precision_solver) {
                                    if(Dtinc_cur == sptr_thermomeca->Dn_mini) {
                                        if(inforce_solver == 1) {
                                        
                                            cout << "The solver has been inforced to proceed (Solver issue) at step:" << sptr_thermomeca->number << " inc: " << inc << " and fraction:" << tinc << ", with the error: " << error << "\n";
                                            cout << "The next increment has integrated the error to avoid propagation\n";
                                            //The solver has been inforced!
                                            tnew_dt = 1.;
                                        
                                            if (inc+1<sptr_thermomeca->ninc) {
                                                for(int k = 0 ; k < 6 ; k++)
                                                {
                                                    if(sptr_thermomeca->cBC_meca(k)) {
                                                        sptr_thermomeca->mecas(inc+1,k) -= residual(k);
                                                    }
                                                    if (sptr_thermomeca->cBC_T) {
                                                        sptr_thermomeca->Ts(inc+1) -= residual(6);
                                                    }
                                                    
                                                }
                                            }
                                        }
                                        else return;
                                        
                                    }
                                    else {
                                        tnew_dt = div_tnew_dt_solver;
                                    }
                                }
                                
                                if((compteur < miniter_solver)&&(tnew_dt >= 1.)) {
                                    tnew_dt = mul_tnew_dt_solver;
                                }
                                compteur = 0;
                                
                                sptr_thermomeca->assess_inc(tnew_dt, tinc, Dtinc, rve ,Time, DTime, DR, corate_type);
                                //start variables ready for the next increment
                                
                            }
                            
                            //At the end of each increment, check if results should be written
                            if (so.o_type(i) == 1) {
                                o_ncount++;
                            }
                            if (so.o_type(i) == 2) {
                                o_tcount+=DTime;
                            }
                            
                            //Write the results
                            if (((so.o_type(i) == 1)&&(o_ncount == so.o_nfreq(i)))||(((so.o_type(i) == 2)&&(fabs(o_tcount - so.o_tfreq(i)) < 1.E-12)))) {
                    
                                rve.output(so, i, n, j, inc, Time, "global");
                                rve.output(so, i, n, j, inc, Time, "local");
                                if (so.o_type(i) == 1) {
                                    o_ncount = 0;
                                }
                                if (so.o_type(i) == 2) {
                                    o_tcount = 0.;
                                }
                            }
                            
                            tinc = 0.;
                            inc++;
                        }
                        
                    }
                    
                }
                break;
            }
            default: {
                cout << "the block type is not defined!\n";
                break;
            }
        }
        //end of blocks loops
    }
    
}
    
} //namespace simcoon
