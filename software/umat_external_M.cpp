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

///@file umat_singleM.cpp
///@brief umat template to run smart subroutines using Abaqus
///@brief Implemented in 1D-2D-3D
///@author Chemisky & Despringre
///@version 1.0
///@date 12/04/2013

#include <iostream>
#include <fstream>
#include <assert.h>
#include <string.h>
#include <armadillo>

#include <simcoon/parameter.hpp>
#include <simcoon/Simulation/Phase/phase_characteristics.hpp>
#include <simcoon/Simulation/Phase/state_variables_M.hpp>
#include <simcoon/Continuum_mechanics/Umat/umat_smart.hpp>

#include <simcoon/Simulation/Maths/lagrange.hpp>
#include <simcoon/Simulation/Maths/rotation.hpp>
#include <simcoon/Simulation/Maths/num_solve.hpp>
#include <simcoon/Continuum_Mechanics/Functions/contimech.hpp>
#include <simcoon/Continuum_Mechanics/Functions/constitutive.hpp>
#include <simcoon/Continuum_Mechanics/Functions/recovery_props.hpp>
#include <simcoon/Continuum_Mechanics/Functions/criteria.hpp>

///@param stress array containing the components of the stress tensor (dimension ntens)
///@param statev array containing the evolution variables (dimension nstatev)
///@param ddsdde array containing the mechanical tangent operator (dimension ntens*ntens)
///@param sse unused
///@param spd unused
///@param scd unused
///@param rpl unused
///@param ddsddt array containing the thermal tangent operator
///@param drple unused
///@param drpldt unused
///@param stran array containing total strain component (dimension ntens) at the beginning of increment
///@param dstran array containing the component of total strain increment (dimension ntens)
///@param time two compoenent array : first component is the value of step time at the beginning of the current increment and second component is the value of total time at the beginning of the current increment
///@param dtime time increment
///@param temperature temperature avlue at the beginning of increment
///@param Dtemperature temperature increment
///@param predef unused
///@param dpred unused
///@param cmname user-defined material name
///@param ndi number of direct stress components
///@param nshr number of shear stress components
///@param ntens number stress and strain components
///@param nstatev number of evolution variables
///@param props array containing material properties
///@param nprops number of material properties
///@param coords coordinates of the considered point
///@param drot rotation increment matrix (dimension 3*3)
///@param pnewdt ratio of suggested new time increment
///@param celent characteristic element length
///@param dfgrd0 array containing the deformation gradient at the beginning of increment (dimension 3*3)
///@param dfgrd1 array containing the deformation gradient at the end of increment (dimension 3*3)
///@param noel element number
///@param npt integration point number
///@param layer layer number - not used
///@param kspt section point number within the current layer - not used
///@param kstep step number
///@param kinc increment number

using namespace std;
using namespace arma;
using namespace simcoon;

void umat_sma_unified_TD(const arma::vec &, const arma::vec &, arma::vec &, arma::mat &, const arma::mat &, const int &, const arma::vec &, const int &, arma::vec &, const double &, const double &,const double &,const double &, double &, double &, double &, double &, const int &, const int &, const bool &, double &);

extern "C" void umat_(double *stress, double *statev, double *ddsdde, double &sse, double &spd, double &scd, double &rpl, double *ddsddt, double *drplde, double &drpldt, const double *stran, const double *dstran, const double *time, const double &dtime, const double &temperature, const double &Dtemperature, const double &predef, const double &dpred, char *cmname, const int &ndi, const int &nshr, const int &ntens, const int &nstatev, const double *props, const int &nprops, const double &coords, const double *drot, double &pnewdt, const double &celent, const double *dfgrd0, const double *dfgrd1, const int &noel, const int &npt, const double &layer, const int &kspt, const int &kstep, const int &kinc)
{
    UNUSED(sse);
    UNUSED(spd);
	UNUSED(scd);
	UNUSED(rpl);
	UNUSED(ddsddt);
	UNUSED(drplde);
	UNUSED(drpldt);
	UNUSED(predef);
	UNUSED(dpred);
	UNUSED(ntens);
	UNUSED(coords);
	UNUSED(celent);
	UNUSED(dfgrd0);
	UNUSED(dfgrd1);
	UNUSED(noel);
	UNUSED(npt);
	UNUSED(layer);
	UNUSED(kspt);
	UNUSED(kstep);
	UNUSED(kinc);
	
	bool start = false;
	double Time = 0.;
	double DTime = 0.;
    int solver_type = 0;
    
	int nstatev_smart = nstatev-4;
    vec props_smart = zeros(nprops);
    vec statev_smart = zeros(nstatev_smart);
    
	vec sigma = zeros(6);
	vec Etot = zeros(6);
	vec DEtot = zeros(6);
	mat Lt = zeros(6,6);
	mat DR = zeros(3,3);
    vec Wm = zeros(4);
    double T = 0.;
    double DT = 0.;
    
	abaqus2smart_M(stress, ddsdde, stran, dstran, time, dtime, temperature, Dtemperature, nprops, props, nstatev, statev, ndi, nshr, drot, sigma, Lt, Etot, DEtot, T, DT, Time, DTime, props_smart, Wm, statev_smart, DR, start);
    umat_sma_unified_TD(Etot, DEtot, sigma, Lt, DR, nprops, props_smart, nstatev_smart, statev_smart, T, DT, Time, DTime, Wm(0), Wm(1), Wm(2), Wm(3), ndi, nshr, start, pnewdt);
	smart2abaqus_M(stress, ddsdde, statev, ndi, nshr, sigma, statev_smart, Wm, Lt);
}

void umat_sma_unified_TD(const vec &Etot, const vec &DEtot, vec &sigma, mat &Lt, const mat &DR, const int &nprops, const vec &props, const int &nstatev, vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, double &tnew_dt) {
    
    /*    UNUSED(nprops);
     UNUSED(nstatev);
     UNUSED(Time);
     UNUSED(DTime);
     UNUSED(nshr);
     UNUSED(tnew_dt);*/
    
    ///@brief The mechanical transformation UMAT for SMAs has the following statev and props
    
    ///@brief props[0] : flagT: 0 transformation temperatures linearly extrapolated; 1 : smooth temperatures
    ///@brief props[1] : EA: Young's modulus of Austenite
    ///@brief props[2] : EM: Young's modulus of Martensite
    ///@brief props[3] : nuA : Poisson's ratio of Austenite
    ///@brief props[4] : nuM : Poisson's ratio of Martensite
    ///@brief props[5] : alphaA_iso : CTE of Austenite
    ///@brief props[6] : alphaM_iso : CTE of Martensite
    ///@brief props[7] : Hmin : Minimal transformation strain magnitude
    ///@brief props[8] : Hmax : Maximal transformation strain magnitude
    ///@brief props[9] : k1 : Exponential evolution of transformation strain magnitude
    ///@brief props[10] : sigmacrit : Critical stress for change of transformation strain magnitude
    ///@brief props[11]: C_A : Slope of martesnite -> austenite parameter
    ///@brief props[12]: C_M : Slope of austenite -> martensite parameter
    ///@brief props[13]: Ms0 : Martensite start at zero stress
    ///@brief props[14]: Mf0 : Martensite finish at zero stress
    ///@brief props[15]: As0 : Austenite start at zero stress
    ///@brief props[16]: Af0 : Austenite finish at zero stress
    ///@brief props[17]: n1 : Martensite start smooth exponent
    ///@brief props[18]: n2 : Martensite finish smooth exponent
    ///@brief props[19]: n3 : Austenite start smooth exponent
    ///@brief props[20]: n4 : Austenite finish smooth exponent
    
    ///@brief props[21]: c_lambda : penalty function exponent start point
    ///@brief props[22]: p0_lambda : penalty function exponent limit penalty value
    ///@brief props[23]: n_lambda : penalty function power law exponent
    ///@brief props[24]: alpha_lambda : penalty function power law parameter
    
    ///@brief The elastic-plastic UMAT with isotropic hardening requires 14 statev:
    ///@brief statev[0] : T_init : Initial temperature
    ///@brief statev[1] : xi : MVF (Martensitic volume fraction)
    ///@brief statev[2] : Transformation strain 11: ET(0,0)
    ///@brief statev[3] : Transformation strain 22: ET(1,1)
    ///@brief statev[4] : Transformation strain 33: ET(2,2)
    ///@brief statev[5] : Transformation strain 12: ET(0,1) (*2)
    ///@brief statev[6] : Transformation strain 13: ET(0,2) (*2)
    ///@brief statev[7] : Transformation strain 23: ET(1,2) (*2)
    
    ///@brief statev[8] : xiF : forward MVF
    ///@brief statev[9] : xiR : reverse MVF
    ///@brief statev[10] : rhoDs0 difference in entropy for the phases (M - A)
    ///@brief statev[11] : rhoDs0 difference in internal energy for the phases (M - A)
    ///@brief statev[12] : parameter for the stress dependance of transformation limits
    ///@brief statev[13] : a1 : forward hardening parameter
    ///@brief statev[14] : a2 : reverse hardening parameter
    ///@brief statev[15] : a3 : Equilibrium hardening parameter
    ///@brief statev[16] : Y0t : Initial transformation critical value
    
    // Taking the values from the material created
    int flagT = props(0);
    double E_A = props(1);
    double E_M = props(2);
    double nu_A = props(3);
    double nu_M= props(4);
    double alphaA_iso = props(5);
    double alphaM_iso = props(6);
    //parameters for Hcur
    double Hmin = props(7);
    double Hmax = props(8);
    double k1 = props(9);
    double sigmacrit = props(10);
    //parameters for phase transformation
    double C_A = props(11);
    double C_M = props(12);
    double Ms0 = props(13);
    double Mf0 = props(14);
    double As0 = props(15);
    double Af0 = props(16);
    //Additional parameters
    double n1 = props(17);
    double n2 = props(18);
    double n3 = props(19);
    double n4 = props(20);
    double sigmacaliber = props(21);
    double prager_b = props(22);
    double prager_n = props(23);
    double sigmastar = 0.;
    
    //Set the Lagrange multipliers coefficients
    double c_lambda = props(24);
    double p0_lambda = props(25);
    double n_lambda = props(26);
    double alpha_lambda = props(27);
    
    double wtp = props(28);
    double C0tp = props(29);
    double C1tp = props(30);
    double Cp = props(31);
    double gammap = props(32);
    double C2tp = props(33);
    double sigmaY_ratch = props(34);
    double alphap = props(35);
    double p0_coa = props(36);
    double Rtp = props(37);
    
    double D_crit = props(38);
    double D_coa = props(39);
    double gammad = props(40);
    double Cd = props(41);
    double Nf_0 = props(42);
    
    ///@brief Temperature initialization
    double T_init = statev(0);
    ///@brief Martensite volume fraction initialization
    double xi = statev(1);
    ///@brief mean strain tensor creation
    vec ET = zeros(6);
    ET(0) = statev(2);
    ET(1) = statev(3);
    ET(2) = statev(4);
    ET(3) = statev(5);
    ET(4) = statev(6);
    ET(5) = statev(7);
    
    vec ETT = zeros(6);
    ETT(0) = statev(8);
    ETT(1) = statev(9);
    ETT(2) = statev(10);
    ETT(3) = statev(11);
    ETT(4) = statev(12);
    ETT(5) = statev(13);
    
    vec ETP(6);
    ETP(0) = statev(14);
    ETP(1) = statev(15);
    ETP(2) = statev(16);
    ETP(3) = statev(17);
    ETP(4) = statev(18);
    ETP(5) = statev(19);
    
    vec ETD(6);
    ETD(0) = statev(20);
    ETD(1) = statev(21);
    ETD(2) = statev(22);
    ETD(3) = statev(23);
    ETD(4) = statev(24);
    ETD(5) = statev(25);
    
    ///@brief ETMax allow the definition of the lambdaTTR
    vec DETF = zeros(6);
    vec DETR = zeros(6);
    vec ETMean = zeros(6);
    
    double xiF = statev(26);
    double xiR = statev(27);
    
    double p = statev(28);
    double pF = statev(29);
    double pR = statev(30);
    
    double d = statev(31);
    double dF = statev(32);
    double dR = statev(33);
    
    double xi_IR = statev(34);
    
    double rhoDs0 = statev(35);
    double rhoDE0 = statev(36);
    double Dtr = statev(37);
    double a1 = statev(38);
    double a2 = statev(39);
    double a3 = statev(40);
    double Y0t = statev(41);
    
    // ######################  Elastic compliance and stiffness #################################
    //defines K and mu explicitely
    //Find the elastic stiffness tensor that is dependent on fraction volume of martensite
    double K_A = E_A/(3.*(1.-2*nu_A));
    double mu_A = E_A/(2.*(1.+nu_A));
    
    double K_M = E_M/(3.*(1.-2*nu_M));
    double mu_M = E_M/(2.*(1.+nu_M));
    
    double K_eff = (K_A*K_M)/(xi*K_A + (1. - xi)*K_M);
    double mu_eff = (mu_A*mu_M)/(xi*mu_A + (1. - xi)*mu_M);
    
    //defines M_A and M_M
    mat M_A = M_iso(K_A, mu_A, "Kmu");
    mat M_M = M_iso(K_M, mu_M, "Kmu");
    mat M = M_iso(K_eff, mu_eff, "Kmu");
    mat L = L_iso(K_eff, mu_eff, "Kmu");
    mat DM = M_M - M_A;
    mat L_tilde = (1.-d)*L;
    
    //definition of the CTE tensor
    vec alpha = (alphaM_iso*xi + alphaA_iso*(1.-xi))*Ith();
    vec Dalpha = (alphaM_iso - alphaA_iso)*Ith();
    
    ///@brief Initialization
    if(start) {
        
        T_init = T;
        vec vide = zeros(6);
        sigma = zeros(6);
        ET = zeros(6);
        ETT = zeros(6);
        ETP = zeros(6);
        ETD = zeros(6);
        xiF = limit;
        xiR = 0.;
        xi = xiF;
        
        pF = limit;
        pR = 0.;
        p = limit;
        
        dF = limit;
        dR = 0.;
        d = dF;
        
        xi_IR = 0.;
        
        Wm = 0.;
        Wm_r = 0.;
        Wm_ir = 0.;
        Wm_d = 0.;
        
        //Get Hcurstar and its derivative to determine rhoDs0
        if (sigmacaliber > sigmacrit)
            sigmastar = sigmacaliber - sigmacrit;
        else
            sigmastar = 0.;
        
        double Hcurstar = Hmin + (Hmax - Hmin)*(1. - exp(-k1*sigmastar));
        assert(Hcurstar > 1E-12);
        
        double dHcurstar = (Hmax - Hmin)*k1*exp(-k1*sigmastar);
        rhoDs0 = -2.*C_M*C_A*(Hcurstar + sigmacaliber*(dHcurstar + (1/E_M-1/E_A)))/(C_M+C_A);
        
        //Determination of the Smooth temperatures
        double MsSmooth = 0.;
        double MfSmooth = 0.;
        double AsSmooth = 0.;
        double AfSmooth = 0.;
        
        if(flagT == 0) {
            MsSmooth = 0.5*Ms0*(1.+ (n1+1)*pow(2.,-n1) + (n2-1)*pow(2.,-n2))/(n1*pow(2.,-n1) + n2*pow(2.,-n2)) + 0.5*Mf0*(-1.+ (n1-1)*pow(2.,-n1) + (n2+1)*pow(2.,-n2))/(n1*pow(2.,-n1) + n2*pow(2.,-n2));
            MfSmooth = 0.5*Ms0*(-1.+ (n1+1)*pow(2.,-n1) + (n2-1)*pow(2.,-n2))/(n1*pow(2.,-n1) + n2*pow(2.,-n2)) + 0.5*Mf0*(1.+ (n1-1)*pow(2.,-n1) + (n2+1)*pow(2.,-n2))/(n1*pow(2.,-n1) + n2*pow(2.,-n2));
            AsSmooth = 0.5*As0*(1.+ (n3-1)*pow(2.,-n3) + (n4+1)*pow(2.,-n4))/(n3*pow(2.,-n3) + n4*pow(2.,-n4)) + 0.5*Af0*(-1.+ (n3+1)*pow(2.,-n3) + (n4-1)*pow(2.,-n4))/(n3*pow(2.,-n3) + n4*pow(2.,-n4));
            AfSmooth = 0.5*As0*(-1.+ (n3-1)*pow(2.,-n3) + (n4+1)*pow(2.,-n4))/(n3*pow(2.,-n3) + n4*pow(2.,-n4)) + 0.5*Af0*(1.+ (n3+1)*pow(2.,-n3) + (n4-1)*pow(2.,-n4))/(n3*pow(2.,-n3) + n4*pow(2.,-n4));
        }
        else {
            MsSmooth = Ms0;
            MfSmooth = Mf0;
            AsSmooth = As0;
            AfSmooth = Af0;
        }
        
        rhoDE0 = 0.5*rhoDs0*(MsSmooth + AfSmooth);
        Dtr = (C_M-C_A)*(Hcurstar + sigmacaliber*(dHcurstar + (1/E_M-1/E_A)))/((C_A+C_M)*(Hcurstar + sigmacaliber*dHcurstar));
        a1 = rhoDs0*(MfSmooth - MsSmooth);
        a2 = rhoDs0*(AsSmooth - AfSmooth);
        a3 = -0.25*a1*(1 + 1./(n1 + 1.) - 1./(n2 + 1.))+0.25*a2*(1. + 1./(n3 + 1.) - 1./(n4 + 1.));
        Y0t = 0.5*rhoDs0*(MsSmooth - AfSmooth) - a3;
    }
    
    //Rotation of internal variables (tensors)
    rotate_strain(ET, DR);
    rotate_strain(ETT, DR);
    rotate_strain(ETP, DR);
    rotate_strain(ETD, DR);
    
    //Variables values at the start of the increment
    vec	sigma_start = sigma;
    vec ET_start = ET;
    
    double pF_start = pF;
    double pR_start = pR;
    double dF_start = dF;
    double dR_start = dR;
    
    // Find Hcur explicit
    if (Mises_stress(sigma) > sigmacrit)
        sigmastar = Mises_stress(sigma) - sigmacrit;
    else
        sigmastar = 0.;
    
    double Hcur = Hmin + (Hmax - Hmin)*(1. - exp(-1.*k1*sigmastar));
    
    //definition of Lambdas associated to transformation
    vec lambdaTTF = Hcur*dPrager_stress(sigma, prager_b, prager_n);
    
    if (Mises_strain(ETT) > 1E-6)
        ETMean = dev(ETT) / (xi);
    else if (Mises_stress(sigma) < 1.E-6)
        ETMean = lambdaTTF;
    else
        ETMean = 0.*Ith();
    
    vec lambdaTTR = -1.*ETMean;
    
    vec lambdaTPF = zeros(6);
    vec lambdaTPR = zeros(6);
    
    vec lambdaTDF = zeros(6);
    vec lambdaTDR = zeros(6);
    
    //Definition of the modified Y function
    double YtF = Y0t + Dtr*Hcur*Mises_stress(sigma);
    double YtR = Y0t + Dtr*sum(sigma%ETMean); //*D Changed by Dtr.
    
    double HfF = 0.;
    double HfR = 0.;
    
    //Hardening function definition (Smooth hardening functions)
    if ((xi > 0.)&&((1. - xi) > 0.)) {
        HfF = 0.5*a1*(1. + pow(xi,n1) - pow(1. - xi,n2)) + a3;
        HfR = 0.5*a2*(1. + pow(xi,n3) - pow((1. - xi),n4)) - a3;
    }
    else if ((xi <= 0.)&&((1. - xi) > 0.)) {
        HfF = 0.5*a1*(1. - pow(1. - xi,n2)) + a3;
        HfR = 0.5*a2*(1. - pow((1. - xi),n4)) - a3;
    }
    else if ((xi > 0.)&&((1. - xi) <= 0.)) {
        HfF = 0.5*a1*(1. + pow(xi,n1)) + a3;
        HfR = 0.5*a2*(1. + pow(xi,n3)) - a3;
    }
    else {
        HfF =  a3;
        HfR = - a3;
    }
    
    double HtpF = Rtp;
    double HtpR = Rtp;
    
    //Set the Lagrange multiplliers due to Physical limitations
    double lambda0 = -1.*lagrange_pow_0(xi, c_lambda, p0_lambda, n_lambda, alpha_lambda);
    double lambda1 = lagrange_pow_1(xi, c_lambda, p0_lambda, n_lambda, alpha_lambda);
    
    //Define the value of DM_sig and M_sig;
    vec DM_sig = (DM*sigma_start);
    vec M_sig = (M*sigma_start);
    //Define the value of Dalpha_T
    vec Dalpha_T = Dalpha*(T+DT);
    vec alpha_T = alpha*(T+DT);
    
    //Set the thermo forces
    double A_xiF = rhoDs0*(T+DT) - rhoDE0 + 0.5*sum(sigma%DM_sig) + sum(sigma%Dalpha)*(T+DT-T_init);
    double A_xiF_start = rhoDs0*(T) - rhoDE0 + 0.5*sum(sigma_start%DM_sig) + sum(sigma_start%Dalpha)*(T-T_init);
    double A_xiR = -1.*rhoDs0*(T+DT) + rhoDE0 - 0.5*sum(sigma%DM_sig) - sum(sigma%Dalpha)*(T+DT-T_init);
    double A_xiR_start = -1.*rhoDs0*(T) + rhoDE0 - 0.5*sum(sigma_start%DM_sig) - sum(sigma_start%Dalpha)*(T-T_init);
    double A_gtF = -1.;
    double A_gtR = 1.;
    
    double A_pF = 0.;
    double A_pF_start = 0.;
    double A_pR = 0.;
    double A_pR_start = 0.;
    double A_gpF = -1;
    double A_gpR = -1;
    
    double A_dF = 1./(2.*pow(1.-d,2.))*sum(sigma%M_sig);
    double A_dF_start = A_dF;
    double A_dR = 1./(2.*pow(1.-d,2.))*sum(sigma%M_sig);
    double A_dR_start = A_dR;
    
    //Transformation criteria
    double PhihatF = Hcur*Prager_stress(sigma, prager_b, prager_n);
    double PhihatR = sum(sigma%ETMean);
    
    //Set the plastic hardening functions
    double ftpF = 0.;
    double ftpR = 0.;
    
    //Set the damage accumulation functions
    double ftdF = 0.;
    double ftdR = 0.;
    
    //Set the damage accumulation functions
    if (fabs(PhihatF) > iota) {
        ftdF = (D_crit/2.)*(1./(pow((PhihatF/Cd),-1.*gammad)-Nf_0));
        double fact_ratch = Macaulay_p((Mises_stress(sigma_start)-sigmaY_ratch)/sigmaY_ratch);
        if ((d > iota)&&(fact_ratch > 0.)) {
            ftpF = wtp*C0tp*(pow(PhihatF/Cp,gammap)*(C1tp*p+exp(-1.*p/C2tp)) + pow(fact_ratch,alphap)*lagrange_pow_1(d/D_crit, 1.-D_coa/D_crit, p0_coa, 0., 2.));
            
        }
        else {
            ftpF = wtp*C0tp*(pow(PhihatF/Cp,gammap)*(C1tp*p+exp(-1.*p/C2tp)));
        }
    }
    else {
        ftpF = 0.;
        ftdF = 0.;
    }
    
    if (fabs(PhihatR) > iota) {
        ftdR = (D_crit/2.)*(1./(pow((PhihatR/Cd),-1.*gammad)-Nf_0));
        double fact_ratch = Macaulay_p((Mises_stress(sigma_start)-sigmaY_ratch)/sigmaY_ratch);
        if ((d > iota)&&(fact_ratch > 0.)) {
            ftpR = (1.-wtp)*C0tp*(pow(Macaulay_p(PhihatR)/Cp,gammap)*(C1tp*p+exp(-1.*p/C2tp))+ pow(fact_ratch,alphap)*lagrange_pow_1(d/D_crit, 1.-D_coa/D_crit, p0_coa, 0., 2.));
        }
        else {
            ftpR = (1.-wtp)*C0tp*(pow(Macaulay_p(PhihatR)/Cp,gammap)*(C1tp*p+exp(-1.*p/C2tp)));
        }
    }
    else {
        ftpR = 0.;
        ftdR = 0.;
    }
    
    //Variables required for the loop
    vec s_j = zeros(2);
    s_j(0) = xiF;
    s_j(1) = xiR;
    vec Ds_j = zeros(2);
    vec ds_j = zeros(2);
    
    ///Elastic prediction - Accounting for the thermal prediction
    vec Eel = Etot + DEtot - alpha*(T+DT-T_init) - ET;
    //sigma = el_pred(L, Eel, ndi);
    sigma = el_pred(L_tilde, Eel, ndi);
    
    //Define the functions for the system to solve
    vec Phi = zeros(2);
    mat B = zeros(2,2);
    vec Y_crit = zeros(2);
    
    //Define the function for the system to solve
    double dHfF = 0.;
    double dHfR = 0.;
    double dHtpF = 0.;
    double dHtpR = 0.;
    vec dHcurdsigma = zeros(6);
    //Relative to forward transformation
    vec dPhihatFdsigma = zeros(6);
    double dPhihatFdxiF = 0.;
    double dPhihatFdxiR = 0.;
    
    vec dA_xiFdsigma = zeros(6);
    double dA_xiFdxiF = 0.;
    double dA_xiFdxiR = 0.;
    
    vec dlambda1dsigma = zeros(6);
    double dlambda1dxiF = 0.;
    double dlambda1dxiR = 0.;
    
    vec dYtFdsigma = zeros(6);
    double dYtFdxiF = 0.;
    double dYtFdxiR = 0.;
    
    vec dPhiFdsigma = zeros(6);
    double dPhiFdxiF = 0.;
    double dPhiFdxiR = 0.;
    double dPhiFdpF = 0.;
    double dPhiFdpR = 0.;
    double dPhiFddF = 0.;
    double dPhiFddR = 0.;
    
    //Relative to reverse transformation
    vec dPhihatRdsigma = zeros(6);
    double dPhihatRdxiF = 0.;
    double dPhihatRdxiR = 0.;
    vec dPhihatRdETF = zeros(6);
    vec dPhihatRdETR = zeros(6);
    
    vec dA_xiRdsigma = zeros(6);
    double dA_xiRdxiF = 0.;
    double dA_xiRdxiR = 0.;
    
    vec dlambda0dsigma = zeros(6);
    double dlambda0dxiF = 0.;
    double dlambda0dxiR = 0.;
    
    vec dYtRdsigma = zeros(6);
    double dYtRdxiF = 0.;
    double dYtRdxiR = 0.;
    vec dYtRdETF = zeros(6);
    vec dYtRdETR = zeros(6);
    
    vec dPhiRdsigma = zeros(6);
    double dPhiRdxiF = 0.;
    double dPhiRdxiR = 0.;
    vec dPhiRdETTF = zeros(6);
    vec dPhiRdETTR = zeros(6);
    double dPhiRdpF = 0.;
    double dPhiRdpR = 0.;
    double dPhiRddF = 0.;
    double dPhiRddR = 0.;
    
    vec dftpFdsigma = zeros(6);
    vec dftpRdsigma = zeros(6);
    
    //Set the derivative of damage accumulation functions
    vec dftdFdDsig = zeros(6);
    vec dftdRdDsig = zeros(6);
    
    //Compute the explicit flow direction
    std::vector<vec> kappa_j(2);
    mat K = zeros(2,2);
    
    //Loop parameters
    int compteur = 0;
    double error = 1.;
    
    //Loop
    for (compteur = 0; ((compteur < maxiter_umat) && (error > precision_umat)); compteur++) {
        
        K_eff = (K_A*K_M) / (xi*K_A + (1. - xi)*K_M);
        mu_eff = (mu_A*mu_M) / (xi*mu_A + (1. - xi)*mu_M);
        L = L_iso(K_eff, mu_eff, "Kmu");
        M = M_iso(K_eff, mu_eff, "Kmu");
        
        DM_sig = DM*sigma; //With respect to xi
        M_sig = M*sigma;    // With respect to d
        Dalpha_T = Dalpha*(T+DT);
        
        //Hardening function definition (Smooth hardening functions)
        if ((xi > 0.) && ((1. - xi) > 0.)) {
            HfF = 0.5*a1*(1. + pow(xi, n1) - pow(1. - xi, n2)) + a3;
        }
        else if ((xi <= 0.) && ((1. - xi) > 0.)) {
            HfF = 0.5*a1*(1. - pow(1. - xi, n2)) + a3;
        }
        else if ((xi > 0.) && ((1. - xi) <= 0.)) {
            HfF = 0.5*a1*(1. + pow(xi, n1)) + a3;
        }
        else  {
            HfF = 0.5*a1 + a3;
        }
        
        if ((xi > 0.) && ((1. - xi) > 0.)) {
            HfR = 0.5*a2*(1. + pow(xi, n3) - pow((1. - xi), n4)) - a3;
        }
        else if ((xi <= 0.) && ((1. - xi) > 0.)) {
            HfR = 0.5*a2*(1. - pow((1. - xi), n4)) - a3;
        }
        else if ((xi > 0.) && ((1. - xi) <= 0.)) {
            HfR = 0.5*a2*(1. + pow(xi, n3)) - a3;
        }
        else
            HfR = 0.5*a2* -a3;
        
        HtpF = Rtp;
        HtpR = Rtp;
        
        // Find Hcur explicit
        if (Mises_stress(sigma) > sigmacrit)
            sigmastar = Mises_stress(sigma) - sigmacrit;
        else
            sigmastar = 0.;
        
        Hcur = Hmin + (Hmax - Hmin)*(1. - exp(-1.*k1*sigmastar));
        
        //Forward transformation thermodynamic force
        A_gtF = -1.;
        A_gtR = 1.;
        A_gpF = -1.;
        A_gpR = -1.;
        
        PhihatF = Hcur*Prager_stress(sigma, prager_b, prager_n);
        A_xiF = rhoDs0*(T + DT) - rhoDE0 + 0.5*sum(sigma%DM_sig) + sum(sigma%Dalpha)*(T + DT - T_init);
        A_pF = 0.;
        A_dF = 1./(2.*pow(1.-d,2.))*sum(sigma%M_sig);
        lambda1 = lagrange_pow_1(xi, c_lambda, p0_lambda, n_lambda, alpha_lambda);
        YtF = Y0t + Dtr*Hcur*Mises_stress(sigma);
        Phi(0) = PhihatF + A_xiF + A_gtF*HfF + A_gpF*HtpF*p - lambda1 - YtF; //A_gtF = A_gt
        
        //Reverse transformation thermodynamic force
        PhihatR = sum(sigma%ETMean);
        A_xiR = -1.*rhoDs0*(T + DT) + rhoDE0 - 0.5*sum(sigma%DM_sig) - sum(sigma%Dalpha)*(T + DT - T_init);
        A_pR = 0.;
        A_dR = 1./(2.*pow(1.-d,2.))*sum(sigma%M_sig);
        lambda0 = -1.*lagrange_pow_0(xi, c_lambda, p0_lambda, n_lambda, alpha_lambda);
        YtR = Y0t + Dtr*sum(sigma%ETMean);
        
        Phi(1) = -1.*PhihatR + A_xiR + A_gtR*HfR - A_gpR*HtpR*p + lambda0 - YtR;  // PhiR < 0. ; A_gtR = -1.*A_gt
        
        mat dMddF = (1./pow(1.-d,2.))*M;
        mat dMddR = (1./pow(1.-d,2.))*M;
        vec dalphaddF = (1./pow(1.-d,2.))*alpha;
        vec dalphaddR = (1./pow(1.-d,2.))*alpha;
        
        //Set the evolution tensors
        lambdaTTF = Hcur * dPrager_stress(sigma, prager_b, prager_n);
        lambdaTTR = -1. * ETMean;
        
        lambdaTPF = dPrager_stress(sigma, prager_b, prager_n);
        lambdaTPR = dPrager_stress(sigma, prager_b, prager_n);
        
        lambdaTDF = zeros(6);
        lambdaTDR = zeros(6);
        
        kappa_j[0] = L_tilde*(lambdaTTF + lambdaTPF*ftpF + lambdaTDF*ftdF + DM_sig + (dMddF*sigma)*ftdF + Dalpha_T + dalphaddF*(T+DT-T_init)*ftdF);
        kappa_j[1] = L_tilde*(lambdaTTR + lambdaTPR*ftpR + lambdaTDR*ftdR - DM_sig + (dMddR*sigma)*ftdR - Dalpha_T + dalphaddR*(T+DT-T_init)*ftdR); //derivative
        
        if ((xi > 0.) && ((1. - xi) > 0.)) {
            dHfF = 0.5*a1*(n1*pow(xi, n1 - 1.) + n2*pow(1. - xi, n2 - 1.));
        }
        else if ((xi <= 0.) && ((1. - xi) > 0.)) {
            dHfF = 0.5*a1*(n2*pow(1. - xi, n2 - 1.));
        }
        else if ((xi > 0.) && ((1. - xi) <= 0.)) {
            dHfF = 0.5*a1*(n1*pow(xi, n1 - 1.));
        }
        else
        {
            dHfF = 0.;
        }
        
        if ((xi > 0.) && ((1. - xi) > 0.)) {
            dHfR = 0.5*a2*(n3*pow(xi, n3 - 1.) + n4*pow(1. - xi, n4 - 1.));
        }
        else if ((xi <= 0.) && ((1. - xi) > 0.)) {
            dHfR = 0.5*a2*(n4*pow(1. - xi, n4 - 1.));
        }
        else if ((xi > 0.) && ((1. - xi) <= 0.)) {
            dHfR = 0.5*a2*(n3*pow(xi, n3 - 1.));
        }
        else
        {
            dHfR = 0.;
        }
        
        dHtpF = 0.;
        dHtpR = 0.;
        
        dHcurdsigma = k1*(Hmax - Hmin)*exp(-1.*k1*sigmastar)*eta_stress(sigma);
        
        //Related to forward transformation
        dPhihatFdsigma = dHcurdsigma * Prager_stress(sigma, prager_b, prager_n) + Hcur * dPrager_stress(sigma, prager_b, prager_n);
        dPhihatFdxiF = 0.;
        dPhihatFdxiR = 0.;
        
        dA_xiFdsigma = DM_sig + Dalpha*(T+DT);
        dA_xiFdxiF = 0.;
        dA_xiFdxiR = 0.;
        
        dlambda1dsigma = zeros(6);
        dlambda1dxiF = dlagrange_pow_1(xi, c_lambda, p0_lambda, n_lambda, alpha_lambda);
        dlambda1dxiR = -1.*dlagrange_pow_1(xi, c_lambda, p0_lambda, n_lambda, alpha_lambda);
        
        dYtFdsigma = Dtr*(dHcurdsigma * Mises_stress(sigma) + Hcur * eta_stress(sigma));
        dYtFdxiF = 0.;
        dYtFdxiR = 0.;
        
        dPhiFdsigma = dPhihatFdsigma + dA_xiFdsigma - dlambda1dsigma - dYtFdsigma;
        dPhiFdxiF = dPhihatFdxiF + dA_xiFdxiF - dHfF - dlambda1dxiF - dYtFdxiF;
        dPhiFdxiR = dPhihatFdxiR + dA_xiFdxiR + dHfF - dlambda1dxiR - dYtFdxiR;
        
        dPhiFdpF = A_gpF*HtpF;
        dPhiFdpR = A_gpF*HtpF;
        dPhiFddF = 0.;
        dPhiFddR = 0.;
        
        //Relative to reverse transformation
        dPhihatRdsigma = ETMean;
        dPhihatRdxiF = (-1./xi)*sum(sigma%ETMean); //Need to add the derivative with respect to
        dPhihatRdxiR = (1./xi)*sum(sigma%ETMean);
        dPhihatRdETF = sigma/xi;
        dPhihatRdETR = sigma/xi;
        
        dA_xiRdsigma = -1.*DM_sig -1.*Dalpha*(T+DT-T_init);
        dA_xiRdxiF = 0.;
        dA_xiRdxiR = 0.;
        
        dlambda0dsigma = zeros(6);
        dlambda0dxiF = -1.*dlagrange_pow_0(xi, c_lambda, p0_lambda, n_lambda, alpha_lambda);
        dlambda0dxiR = dlagrange_pow_0(xi, c_lambda, p0_lambda, n_lambda, alpha_lambda);
        
        dYtRdsigma = 1.*Dtr*ETMean;
        dYtRdxiF = (-Dtr/xi)*sum(sigma%ETMean);
        dYtRdxiR = (Dtr/xi)*sum(sigma%ETMean);
        dYtRdETF = Dtr*sigma/xi;
        dYtRdETR = Dtr*sigma/xi;
        
        dPhiRdsigma = -1.*dPhihatRdsigma + dA_xiRdsigma + dlambda0dsigma - dYtRdsigma;
        
        dPhiRdxiF = -1.*dPhihatRdxiF + dA_xiRdxiF + dHfR + dlambda0dxiF - dYtRdxiF;
        dPhiRdxiR = -1.*dPhihatRdxiR + dA_xiRdxiR - dHfR + dlambda0dxiR - dYtRdxiR;
        
        dPhiRdpF = A_gpR*HtpR;
        dPhiRdpR = A_gpR*HtpR;
        
        dPhiRddF = 0.;
        dPhiRddR = 0.;
        
        dPhiRdETTF = -1.*dPhihatRdETF - dYtRdETF;
        dPhiRdETTR = -1.*dPhihatRdETF - dYtRdETF;
        
        K(0,0) = dPhiFdxiF + dPhiFdpF*ftpF + dPhiFddF*ftdF;
        K(0,1) = dPhiFdxiR + dPhiFdpR*ftpR + dPhiFddR*ftdR;
        K(1,0) = dPhiRdxiF + dPhiRdpF*ftpF + dPhiRddF*ftdF + sum(dPhiRdETTF%lambdaTTF);
        K(1,1) = dPhiRdxiR + dPhiRdpR*ftpR + dPhiRddR*ftdR + sum(dPhiRdETTR%lambdaTTR);
        
        B(0,0) = -1.*sum(dPhiFdsigma%kappa_j[0]) + K(0,0);
        B(0,1) = -1.*sum(dPhiFdsigma%kappa_j[1]) + K(0,1);
        B(1,0) = -1.*sum(dPhiRdsigma%kappa_j[0]) + K(1,0);
        B(1,1) = -1.*sum(dPhiRdsigma%kappa_j[1]) + K(1,1);
        
        Y_crit(0) = YtF;
        Y_crit(1) = YtR;
        
        Fischer_Burmeister_m(Phi, Y_crit, B, Ds_j, ds_j, error);
        
        s_j(0) += ds_j(0);
        s_j(1) += ds_j(1);
        
        xiF = s_j(0);
        xiR = s_j(1);
        
        xi = xi + ds_j(0) - 1.*ds_j(1);
        pF = pF + ftpF*ds_j(0);
        pR = pR + ftpR*ds_j(1);
        p = pF + pR;
        
        dF = dF + ftdF*ds_j(0);
        dR = dR + ftdR*ds_j(1);
        d = dF + dR;
        
        if (d>D_crit) {
            d = D_crit;
        }
        
        L_tilde = (1.-d)*L;
        
        ETT = ETT + ds_j(0)*lambdaTTF + ds_j(1)*lambdaTTR;
        ETP = ETP + ftpF*ds_j(0)*lambdaTPF + ftpR*ds_j(1)*lambdaTPR;
        ETD = ETD + ftdF*ds_j(0)*lambdaTDF + ftdR*ds_j(1)*lambdaTDR;
        
        ET = ETT + ETP + ETD;
        
        DETF += ds_j(0)*(lambdaTTF + lambdaTPF*ftpF + lambdaTDF*ftdF);
        DETR += ds_j(1)*(lambdaTTR + lambdaTPR*ftpR + lambdaTDR*ftdR);
        
        if((Mises_strain(ETT) > precision_umat)&&(xi > precision_umat))
        {
            ETMean = dev(ETT) / (xi);
        }
        else {
            ETMean = lambdaTTF;
        }
        
        //the stress is now computed using the relationship sigma = L(E-Ep)
        Eel = Etot + DEtot - alpha*(T + DT - T_init) - ET;
        //sigma = el_pred(L, Eel, ndi);
        sigma = el_pred(L_tilde, Eel, ndi);
    }
    
    //Computation of the increments of variables
    vec Dsigma = sigma - sigma_start;
    vec DET = ET - ET_start;
    double DxiF = Ds_j[0];
    double DxiR = Ds_j[1];
    
    double DpF = pF - pF_start;
    double DpR = pR - pR_start;
    double DdF = dF - dF_start;
    double DdR = dR - dR_start;
    
    //Computation of the tangent modulus
    mat Bhat = zeros(2, 2);
    Bhat(0,0) = sum(dPhiFdsigma%kappa_j[0]) - K(0,0);
    Bhat(0,1) = sum(dPhiFdsigma%kappa_j[1]) - K(0,1);
    Bhat(1,0) = sum(dPhiRdsigma%kappa_j[0]) - K(1,0);
    Bhat(1,1) = sum(dPhiRdsigma%kappa_j[1]) - K(1,1);
    
    vec op = zeros(2);
    mat delta = eye(2,2);
    
    for (int i=0; i<2; i++) {
        if(Ds_j[i] > iota)
            op(i) = 1.;
    }
    
    mat Bbar = zeros(2,2);
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            Bbar(i, j) = op(i)*op(j)*Bhat(i, j) + delta(i,j)*(1-op(i)*op(j));
        }
    }
    
    mat invBbar = zeros(2, 2);
    mat invBhat = zeros(2, 2);
    invBbar = inv(Bbar);
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            invBhat(i, j) = op(i)*op(j)*invBbar(i, j);
        }
    }
    
    std::vector<vec> P_epsilon(2);
    P_epsilon[0] = invBhat(0, 0)*(L_tilde*dPhiFdsigma) + invBhat(1, 0)*(L_tilde*dPhiRdsigma);
    P_epsilon[1] = invBhat(0, 1)*(L_tilde*dPhiFdsigma) + invBhat(1, 1)*(L_tilde*dPhiRdsigma);
    
    Lt = L_tilde - (kappa_j[0]*P_epsilon[0].t() + kappa_j[1]*P_epsilon[1].t());
    
    //Preliminaries for the computation of mechanical work
    
    double Dgamma_loc = 0.5*sum((sigma_start+sigma)%(DETF-DETR)) + 0.5*(A_xiF_start + A_xiF)*DxiF + 0.5*(A_xiR_start + A_xiR)*DxiR + 0.5*(A_pF_start + A_pF)*DpF + 0.5*(A_pR_start + A_pR)*DpR + 0.5*(A_dF_start + A_dF)*DdF + 0.5*(A_dR_start + A_dR)*DdR;
    
    //Computation of the mechanical and thermal work quantities
    Wm += 0.5*sum((sigma_start+sigma)%DEtot);
    Wm_r += 0.5*sum((sigma_start+sigma)%(DEtot-DETF+DETR))- 0.5*(A_xiF_start + A_xiF)*DxiF - 0.5*(A_xiR_start + A_xiR)*DxiR;
    Wm_ir += -0.5*(A_pF_start + A_pF)*DpF - 0.5*(A_pR_start + A_pR)*DpR - 0.5*(A_dF_start + A_dF)*DdF - 0.5*(A_dR_start + A_dR)*DdR;
    Wm_d += Dgamma_loc;
    
    ///@brief statev evolving variables
    statev(0) = T_init;
    statev(1) = xi;
    
    statev(2) = ET(0);
    statev(3) = ET(1);
    statev(4) = ET(2);
    statev(5) = ET(3);
    statev(6) = ET(4);
    statev(7) = ET(5);
    
    statev(8) = ETT(0);
    statev(9) = ETT(1);
    statev(10) = ETT(2);
    statev(11) = ETT(3);
    statev(12) = ETT(4);
    statev(13) = ETT(5);
    
    statev(14) = ETP(0);
    statev(15) = ETP(1);
    statev(16) = ETP(2);
    statev(17) = ETP(3);
    statev(18) = ETP(4);
    statev(19) = ETP(5);
    
    statev(20) = ETD(0);
    statev(21) = ETD(1);
    statev(22) = ETD(2);
    statev(23) = ETD(3);
    statev(24) = ETD(4);
    statev(25) = ETD(5);
    
    statev(26) = xiF;
    statev(27) = xiR;
    
    statev(28)= p;
    statev(29) = pF;
    statev(30) = pR;
    
    statev(31) = d;
    statev(32) = dF;
    statev(33) = dR;
    
    statev(34) = xi_IR;
    
    if(start) {
        statev(35) = rhoDs0;
        statev(36) = rhoDE0;
        statev(37) = Dtr;
        
        statev(38) = a1;
        statev(39) = a2;
        statev(40) = a3;
        
        statev(41) = Y0t;
    }
    
    statev(42) = PhihatF;
    statev(43) = PhihatR;
    statev(44) = Macaulay_p((Mises_stress(sigma)-sigmaY_ratch)/sigmaY_ratch);
    statev(45) = pow(Mises_stress(sigma)/sigmaY_ratch,alphap);
    statev(46) = lagrange_pow_1(d/D_crit, 1.-D_coa/D_crit, p0_coa, 0., 2.);
}

