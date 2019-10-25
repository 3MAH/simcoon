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

///@file umat_smart.hpp
///@brief Selection of constitutive laws and transfer to between Abaqus and simcoon formats
///@brief Implemented in 1D-2D-3D
///@version 1.0

#pragma once
#include <armadillo>
#include <simcoon/Simulation/Phase/phase_characteristics.hpp>

extern "C"{
	void umat_(double *, double *, double *, double &, double &, double &, double &, double *, double *, double &, const double *, const double *, const double *, const double &, const double &, const double &, const double &, const double &, char *, const int &, const int &, const int &, const int &, const double *, const int &, const double &, const double *, double &, const double &, const double *, const double *, const int &, const int &, const double &, const int &, const int &, const int &);
}

namespace simcoon{

///@param stress array containing the components of the stress tensor (dimension ntens)
///@param stran array containing total strain component (dimension ntens) at the beginning of increment
///@param dstran array containing the component of total strain increment (dimension ntens)
///@param time two compoenent array : first component is the value of step time at the beginning of the current increment and second component is the value of total time at the beginning of the current increment
///@param dtime time increment
///@param temperature temperature avlue at the beginning of increment
///@param Dtemperature temperature increment
///@param ndi number of direct stress components
///@param nshr number of shear stress components
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

void size_statev(phase_characteristics &, unsigned int &);
    
void statev_2_phases(phase_characteristics &, unsigned int &, const arma::vec &);

void phases_2_statev(arma::vec &, unsigned int &, const phase_characteristics &);
    
void abaqus2smart_M_light(const double *, const double *, const int &, double *, const int &, const int &, arma::vec &, arma::mat &, arma::vec &, arma::vec &);
    
void abaqus2smart_M(const double *, const double *, const double *, const double *, const double *, const double &, const double &, const double &, const int &,const double *, const int &, double *, const int &, const int &, const double *, arma::vec &, arma::mat &, arma::vec &, arma::vec &, double &, double &, double &, double &, arma::vec &, arma::vec &, arma::vec &, arma::mat &, bool &);

void abaqus2smart_M_full(const double *, const double *, const double *, const double *, const double *, const double &, const double &, const double &, const int &,const double *, const int &, double *, const int &, const int &, const double *, arma::vec &, arma::mat &, arma::vec &, arma::vec &, double &, double &, double &, double &, arma::vec &, arma::vec &, arma::vec &, arma::mat &, bool &);

void abaqus2smart_T(const double *, const double *, const double *, const double *, const double &, const double *, const double *, const double *, const double &, const double &, const double &, const int &,const double *, const int &, double *, const int &, const int &, const double *, arma::vec &, arma::mat &, arma::mat &, arma::mat &, arma::mat &, arma::vec &, arma::vec &, double &, double &, double &, double &, arma::vec &, arma::vec &, arma::vec &, arma::vec &, arma::mat &, bool &);

void run_umat_M_aba(phase_characteristics &, const arma::mat &, const double &, const double &, const int &, const int &, bool &, const int &, double &, const std::string &);
    
void select_umat_T(phase_characteristics &, const arma::mat &, const double &, const double &, const int &, const int &, bool &, const int &, double &);

void select_umat_M_finite(phase_characteristics &, const arma::mat &, const double &, const double &, const int &, const int &, bool &, const int &, double &);
    
void select_umat_M(phase_characteristics &, const arma::mat &, const double &, const double &, const int &, const int &, bool &, const int &, double &);
    
void run_umat_T(phase_characteristics &, const arma::mat &, const double &, const double &, const int &, const int &, bool &, const int &, const unsigned int &, double &);

void run_umat_M(phase_characteristics &, const arma::mat &, const double &, const double &, const int &, const int &, bool &, const int &, const unsigned int &, double &);
    
void smart2abaqus_M(double *, double *, double *, const int &, const int &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::mat &);
   
void smart2abaqus_M_full(double *, double *, double *, double *, double *, double &, double &, double &, int &, double *, int &, double *, const int &, const int &, double *, const arma::vec &, const arma::mat &, const arma::vec &, const arma::vec &, const double &, const double &, const double &, const double &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::mat &, bool &);
    
void smart2abaqus_T(double *, double *, double *, double *, double &, double &, double *, const int &, const int &, const arma::vec &, const arma::vec &, const double &, const arma::vec &, const arma::vec &, const arma::mat &, const arma::mat &, const arma::mat &, const arma::mat &);
            
} //namespace simcoon
