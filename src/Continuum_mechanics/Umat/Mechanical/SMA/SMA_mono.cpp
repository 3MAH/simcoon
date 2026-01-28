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

///@file sma_mono.cpp
///@brief User subroutine for monocrystalline SMA with multiple elastic symmetry options
///@author Chemisky

#include <iostream>
#include <fstream>
#include <string>
#include <armadillo>

#include <simcoon/parameter.hpp>
#include <simcoon/Simulation/Maths/lagrange.hpp>
#include <simcoon/Simulation/Maths/rotation.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Continuum_mechanics/Material/variant.hpp>

using namespace std;
using namespace arma;

namespace simcoon {

///@brief Unified SMA monocrystal UMAT supporting multiple elastic symmetries
///@brief The elastic symmetry is determined by umat_name:
///@brief   SMAMO (isotropic): props(0)=E, props(1)=nu, then common params from props(2)
///@brief   SMAMC (cubic): props(0)=C11, props(1)=C12, props(2)=C44, then common params from props(3)
///@brief
///@brief Common parameters (offset depends on symmetry):
///@brief   alpha_iso : Coefficient of thermal expansion
///@brief   b         : Slope parameter (corresponds to g)
///@brief   g         : Shear strain magnitude
///@brief   Ms        : Martensite start temperature
///@brief   Af        : Austenite finish temperature
///@brief   nvariants : Number of martensite variants
///@brief   c_lambda0, p0_lambda0, n_lambda0, alpha_lambda0 : Lagrange parameters (variant)
///@brief   c_lambda1, p0_lambda1, n_lambda1, alpha_lambda1 : Lagrange parameters (total)

/*void Amortissement(double &xi, const double &dxi, const double &xi_start)
{
		
	double limit1 = 1. - limit;
	double limit2 = limit;

		
    if (xi >= limit1)
    {
		xi = (xi - dxi + limit1)*0.5;
	}
    
    if (xi <= limit2)
    {
        xi = (xi - dxi + limit2)*0.5;
    }
}*/

void umat_sma_mono(const string &umat_name, const vec &Etot, const vec &DEtot, vec &sigma, mat &Lt, mat &L, const mat &DR, const int &nprops, const vec &props, const int &nstatev, vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, double &tnew_dt) {

    UNUSED(nprops);
    UNUSED(nstatev);
    UNUSED(Time);
    UNUSED(DTime);
    UNUSED(nshr);
    UNUSED(tnew_dt);

	// ######################  Props #################################
	//From the props to the material properties
    // Property offset depends on elastic symmetry type
    int offset = 0;
    double alpha_iso = 0.;
    double b = 0.;
    double g = 0.;
    double Ms = 0.;
    double Af = 0.;
    int nvariants = 0;
    double c_lambda0 = 0.;
    double p0_lambda0 = 0.;
    double n_lambda0 = 0.;
    double alpha_lambda0 = 0.;
    double c_lambda1 = 0.;
    double p0_lambda1 = 0.;
    double n_lambda1 = 0.;
    double alpha_lambda1 = 0.;

    // Build elastic stiffness tensor based on umat_name
    if (umat_name == "SMAMO") {
        // Isotropic elasticity: E, nu
        double E = props(0);
        double nu = props(1);
        L = L_iso(E, nu, "Enu");
        offset = 2;
    }
    else if (umat_name == "SMAMC") {
        // Cubic elasticity: C11, C12, C44
        double C11 = props(0);
        double C12 = props(1);
        double C44 = props(2);
        L = L_cubic(C11, C12, C44, "Cii");
        offset = 3;
    }
    else {
        cout << "Error: Unknown umat_name in umat_sma_mono: " << umat_name << "\n";
        exit(0);
    }

    // Extract common parameters using offset
    alpha_iso = props(offset);
    b = props(offset + 1);
    g = props(offset + 2);
    Ms = props(offset + 3);
    Af = props(offset + 4);
    nvariants = int(props(offset + 5));
    c_lambda0 = props(offset + 6);
    p0_lambda0 = props(offset + 7);
    n_lambda0 = props(offset + 8);
    alpha_lambda0 = props(offset + 9);
    c_lambda1 = props(offset + 10);
    p0_lambda1 = props(offset + 11);
    n_lambda1 = props(offset + 12);
    alpha_lambda1 = props(offset + 13);

    double roS0 = -1.*b*g;
    double romu0 = 0.5*roS0*(Ms + Af);
    double Y0 = 0.5*roS0*(Ms - Af);

    //definition of the CTE tensor
    vec alpha = alpha_iso*Ith();
    std::string data_path= std::getenv("SIMCOON_DATA_PATH") ? std::getenv("SIMCOON_DATA_PATH") : "data" ;
	mat Hnm = zeros(nvariants, nvariants);
	Hnm.load(data_path+"/Hnm.inp", raw_ascii);
    
	// ######################  Statev #################################
	
	///@brief Temperature initialization
	double T_init = statev(0);
    
	std::vector<variant> var(nvariants);
	
	///@brief Properties of the variants, use "test.dat" to specify the parameters (for now)
	ifstream paramvariant;
	paramvariant.open(data_path+"/variant.inp", ios::in);
	if(paramvariant) {
		string chaine1;
		for(int i=0; i<nvariants; i++) {
			paramvariant >> chaine1 >> var[i].n(0) >> var[i].n(1) >> var[i].n(2) >> var[i].m(0) >> var[i].m(1) >> var[i].m(2);
			var[i].build(g);
		}		
	}
	else {
		cout << "Error: cannot open .dat file \n";
	}
	paramvariant.close();
    
	vec xin(nvariants);
	vec xin_start(nvariants);
    
    for(int i=0; i<nvariants; i++) {
		xin(i) = statev(i+1);
	}
	
	xin_start = xin;
    
	vec ET = zeros(6);
	ET(0) = statev(nvariants+1);
	ET(1) = statev(nvariants+2);
	ET(2) = statev(nvariants+3);
	ET(3) = statev(nvariants+4);
	ET(4) = statev(nvariants+5);
	ET(5) = statev(nvariants+6);
	
	vec ET_start = ET;
    
	double xi = statev(nvariants+7);
    double xi_start = xi;
    double dxi = 0.;
    
	//Rotation of internal variables (tensors)
    ET = rotate_strain(ET, DR);

    // Note: Elastic stiffness tensor L is already constructed based on umat_name

    ///@brief Initialization
    if(start)
    {
        T_init = T;
        vec vide = zeros(6);
        sigma = vide;
        ET = vide;
        xi = 0.;
        for(int i=0; i<nvariants; i++) {
            xin(i) = sim_limit;
            xi += xin(i);
        }
        
        Wm = 0.;
        Wm_r = 0.;
        Wm_ir = 0.;
        Wm_d = 0.;
    }
	
    vec sigma_start = sigma;
    
	vec PhiF(nvariants);	
	vec PhiF_start(nvariants);			
	vec PhiR(nvariants);
	vec PhiR_start(nvariants);
   	
   	std::vector<int> transfo_actif(24);
    
    ///Elastic prediction - Accounting for the thermal prediction
    vec Eel = Etot + DEtot - alpha*(T+DT-T_init) - ET;
    sigma = el_pred(L, Eel, ndi);
    
    //Need to define the thermodynamic parameters:
	vec lambda0 = zeros(nvariants);
  
	double lambda1 = 0.;
    
	mat dlambda0dxin = zeros(nvariants, nvariants);
  
	double dlambda1dxin = 0.;  
    
    vec Hnmxim = Hnm*xin;
    vec Hnmxim_start = Hnm*xin_start;
    
    lambda1 = lagrange_pow_1(xi, c_lambda1, p0_lambda1, n_lambda1, alpha_lambda1);
    
    //Set the Lagrange multiplliers due to Physical limitations for each Phi
    for(int i=0; i<nvariants; i++) {
		lambda0(i) = -1.*lagrange_pow_0(xin(i), c_lambda0, p0_lambda0, n_lambda0, alpha_lambda0);
        
		//Set the thermo forces  
		PhiF(i) = sum(sigma%var[i].ETn) + roS0*(T+DT) - romu0 - Hnmxim(i) - lambda0(i) - lambda1 - Y0;
		PhiR(i) = sum(sigma%var[i].ETn) + roS0*(T+DT) - romu0 - Hnmxim(i) - lambda0(i) - lambda1 + Y0;
		PhiF_start(i) = sum(sigma_start%var[i].ETn) + roS0*(T) - romu0 - Hnmxim_start(i) - lambda0(i) - lambda1 - Y0;
		PhiR_start(i) = sum(sigma_start%var[i].ETn) + roS0*(T) - romu0 - Hnmxim_start(i) - lambda0(i) - lambda1 + Y0;  
        
		//Define which variant system is active 
		transfo_actif[i] = 0;
	
/*
		if((PhiF(i) > 0)&&(PhiF(i) > PhiF_start(i)))	{
			transfo_actif[i] = 1;
		}
		else if((PhiR(i) < 0)&&(PhiR(i) < PhiR_start(i))) {
			transfo_actif[i] = -1;
		}*/

		if(PhiF(i) > 0)	{
			transfo_actif[i] = 1;
		}
		else if(PhiR(i) < 0) {
			transfo_actif[i] = -1;
		}
		
    }
	
	//Collect the number of potential transforming systems
	//Note : It is possible that some systems have a reverse transformation while others have a forward one ..
	int nactive = 0;
	for(int i=0; i<nvariants; i++) {
		if(abs(transfo_actif[i]) > 0) {
            nactive++;
        }
	}
	
	vec dxin = zeros(nactive);
	vec Phi = zeros(nactive);
	mat dPhidxi = zeros(nactive, nactive);
    mat dPhidsigma = zeros(nactive, 6);
    mat lambda = zeros(6, nactive);
    
	std::vector<int> active;
	std::vector<int> tactive;

	double sumPhi = 0.;
    for(int i=0; i<nvariants; i++) {
        
		if(transfo_actif[i] > 0) {
			sumPhi +=  fabs(PhiF(i));
			active.push_back(i);
			tactive.push_back(1);
		}
		else if (transfo_actif[i] < 0) {
			sumPhi +=  fabs(PhiR(i));
			active.push_back(i);
			tactive.push_back(-1);
        }
	}
    
    int compteur = 0;
    
	if(nactive > 0)
	{
        
		for(compteur = 0; ((compteur < maxiter_umat) && (sumPhi/Y0 > precision_umat)); compteur++)
		{
            dPhidxi = zeros(nactive, nactive);

			lambda1 = lagrange_pow_1(xi, c_lambda1, p0_lambda1, n_lambda1, alpha_lambda1);        
			dlambda1dxin = dlagrange_pow_1(xi, c_lambda1, p0_lambda1, n_lambda1, alpha_lambda1);    

			for(int i=0; i<nactive; i++) {

				lambda0(active[i]) = -1.*lagrange_pow_0(xin(active[i]), c_lambda0, p0_lambda0, n_lambda0, alpha_lambda0);
                dlambda0dxin(active[i],active[i]) = -1.*dlagrange_pow_0(xin(active[i]), c_lambda0, p0_lambda0, n_lambda0, alpha_lambda0);          
                
				if(tactive[i] == 1) {

                    Hnmxim = Hnm*xin;
					Phi(i) = sum(sigma%var[active[i]].ETn) + roS0*(T+DT) - romu0 - Hnmxim(active[i]) - lambda0(active[i]) - lambda1 - Y0;

                    for (int j=0; j<nactive; j++) {
                        dPhidxi(i,j) += -1.*dlambda0dxin(active[i],active[j]) - dlambda1dxin;
                    }
				}
				else if(tactive[i] == -1) {
					
                    Hnmxim = Hnm*xin;
					Phi(i) = sum(sigma%var[active[i]].ETn) + roS0*(T+DT) - romu0 - Hnmxim(active[i]) - lambda0(active[i]) - lambda1 + Y0;
                  
                    for (int j=0; j<nactive; j++) {
                        dPhidxi(i,j) += -1.*dlambda0dxin(active[i],active[j]) - dlambda1dxin;
                    }

				}
                
                for(int j=0; j<6; j++) {
                    dPhidsigma(i,j) = var[active[i]].ETn(j);
                    lambda(j,i) = var[active[i]].ETn(j);
                }
                  
                for(int j=0; j<nactive; j++) {
                    dPhidxi(i,j) += -1.*Hnm(active[i], active[j]);
                }
                
			}
            
            // Assemble the derivatives (CCP algorithm)
            dPhidxi += dPhidsigma*(-1.*L*lambda);
            
			if(fabs(det(dPhidxi)) > 0.)
                dxin = -1.*inv(dPhidxi)*Phi;
            else {
                //pnewdt = 0.1;
//                Phi = zeros(nactive);
                dxin = zeros(nactive);
            }
            
            dxi = 0.;
            for(int i=0; i<nactive; i++) {

				dxi += dxin(i);

                xin(active[i]) += dxin(i);
//				Amortissement(xin(active[i]), dxin(i), xin_start(active[i]));
                
            }       
            
            ET = zeros(6);
            sumPhi = 0.;
            xi = 0.;
            
            for(int i=0; i<nactive; i++) {
                sumPhi +=  fabs(Phi(i));
            }

            for(int i=0; i<nvariants; i++) {
                xi += xin(i);
            }
			
//			Amortissement(xi, dxi, xi_start);            
            for(int i=0; i<nvariants; i++) {
                ET += var[i].ETn*xin(i);
            }

            //the stress is now computed using the relationship sigma = L(E-Ep)
            Eel = Etot + DEtot - alpha*(T + DT - T_init) - ET;
            sigma = el_pred(L, Eel, ndi);
        }
                       
		vec lambda_eff = zeros(6);
		mat Lt_eff = zeros(6,6);
		double xi_path = 0.;		
		double H_eff = 0.;
		
		vec B1(6);
		vec B2(6);
		double At2;		
		
		if(fabs(xi - xi_start) > 1.E-8) {
			
			xi_path = xi - xi_start;
			
	        for(int i=0; i<nactive; i++) {
				lambda_eff +=  1./(xi_path)*(xin(active[i])-xin_start(active[i]))*var[i].ETn;
				for(int j=0; j<nactive; j++) {
					H_eff +=  1./pow(xi_path,2.)*(fabs(xin(active[i])-xin_start(active[i])))*(Hnm(active[i], active[j]) + dlambda0dxin(active[i],active[j]) + dlambda1dxin)*fabs((xin(active[j])-xin_start(active[j])));
				}
	        }

			B1 = L*lambda_eff;
			B2 = L*lambda_eff;
			At2 = -1.*H_eff - sum(lambda_eff%B1);
		
			Lt_eff = L + (B1*trans(B2))/At2;        
		}
		else
			Lt_eff = L;
		
        //Computation of the tangent modulus !
//		if((fabs(det(dPhidxi)) > 0.)&&(fabs(xi - xi_start) > 1.E-5)) {
		if(fabs(det(dPhidxi)) > 0.) {
            Lt = L + ((L*lambda)*(inv(dPhidxi))*(dPhidsigma*L));
        }
        else
            Lt = L;
       
		vec eigval = eig_sym(Lt);
		
/*		///Note : To use with Self-consistent micromechanics only!
		while(eigval(0) < 0.) {
				Lt = 0.99*Lt + 0.01*Lt_eff;
				eigval = eig_sym(Lt);			
		}*/
        
    }
	else {
        Eel = Etot + DEtot - alpha*(T+DT-T_init) - ET;
        sigma = el_pred(L, Eel, ndi);
        
        ///Computation of the tangent modulus
        Lt = L;
    }
    
	vec Dsigma = sigma - sigma_start;
	vec eigval = eig_sym(Lt);
    
/*	if (compteur < 20) {
        tnew_dt = 2.;
    }
    
	if (sum(Dsigma%DEtot) < 0.) {
        //        This constraint is a form of the Drucker Postulate : Stability of hardening consitions
		tnew_dt = 0.2;
        Lt = L;
    }
    
	if (compteur == maxitNewton) {
        tnew_dt = 0.2;
        Lt = L;
    }
    
	if (Mises_stress(Dsigma) > 20.) {
		tnew_dt = 0.2;
        Lt = L;
    }
    
	if (isnan(Mises_stress(sigma))) {
		tnew_dt = 0.2;
        Lt = L;
    }*/
    
    statev(0) = T_init;
    
    for(int i=0; i<nvariants; i++) {
		statev(i+1) = xin(i);
	}
    
	statev(nvariants+1) = ET(0);
	statev(nvariants+2) = ET(1);
	statev(nvariants+3) = ET(2);
	statev(nvariants+4) = ET(3);
	statev(nvariants+5) = ET(4);
	statev(nvariants+6) = ET(5);
	
    statev(nvariants+7) = xi;
}
    
} //namespace simcoon
