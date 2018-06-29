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

///@file materials.cpp
///@brief Characteristics of a material, as an Abaqus input
///@version 1.0

#include <iostream>
#include <fstream>
#include <string>
#include <assert.h>
#include <armadillo>
#include <simcoon/Continuum_mechanics/Unit_cell/materials.hpp>

using namespace std;
using namespace arma;

namespace simcoon{

//=====Private methods for material_characteristics===================================

//=====Public methods for material_characteristics============================================

/*!
  \brief default constructor
*/

//-------------------------------------------------------------
aba_material::aba_material() : material_characteristics()
//-------------------------------------------------------------
{
	id=0;
    nstatev = 0;
    density = 0.;
    conductivity = 0.;
}

/*!
  \brief Constructor
  \param nprops : size of the vector props
  \param nstatev : size of the vector statev
  \param init boolean that indicates if the constructor has to initialize the material characeristics vector (default value is true) \n
  \n\n
  \f$ \textbf{Examples :} \f$ \n
*/

//-------------------------------------------------------------
aba_material::aba_material(const int &n, const bool &init, const double &value)
//-------------------------------------------------------------
{

	assert(n>0);

    id=0;
    nstatev = 0;
    
	number = -1;
    save = 0;
    
	psi_mat=0.;
	theta_mat=0.;
	phi_mat=0.;
	
    density = 0.;
    conductivity = 0.;
    
    nprops = n;
    if (init) {
        props = value*ones(n);
    }
    else{
        props = zeros(n);
    }
}

/*!
  \brief Constructor with parameters
  \param nstatev : size of the table statev 
  \param init boolean that indicates if the constructor has to initialize the state variables vectors and statev table (default value is true) \n
  \n\n
  \f$ \textbf{Examples :} \f$ \n
*/

//-------------------------------------------------------------
aba_material::aba_material(const int &mnumber, const int &mid, const string &mumat_name, const int &msave, const double &mpsi_mat, const double &mtheta_mat, const double &mphi_mat, const double &mconductivity, const double &mdensity, const int &mnprops, const int &mnstatev, const vec &mprops)
//-------------------------------------------------------------
{	
	assert(mnprops);
	
	number = mnumber;
    id = mid;
	umat_name = mumat_name;
    save = msave;
    
    psi_mat = mpsi_mat;
	theta_mat = mtheta_mat;
	phi_mat = mphi_mat;
    
    conductivity = mconductivity;
    density = mdensity;
    
	nprops = mnprops;
	nstatev = mnstatev;
	props = mprops;
}

/*!
  \brief Copy constructor
  \param mc aba_material object to duplicate
*/

//------------------------------------------------------
aba_material::aba_material(const aba_material& mc) : material_characteristics(mc)
//------------------------------------------------------
{
    id = mc.id;
	nstatev = mc.nstatev;
    
    conductivity = mc.conductivity;
    density = mc.density;
}

/*!
  \brief Destructor

  Deletes aba_material, the vectors and matrix.
*/

//-------------------------------------
aba_material::~aba_material() {}
//-------------------------------------

/*!
  \brief Standard operator = for phase_characteristics
*/
    
//-------------------------------------------------------------
void aba_material::update(const int &mnumber, const int &mid, const string &mumat_name, const int &msave, const double &mpsi_mat, const double &mtheta_mat, const double &mphi_mat, const double &mconductivity, const double &mdensity, const int &mnprops, const int &mnstatev, const vec &mprops)
//-------------------------------------------------------------
{
    
    material_characteristics::update(mnumber, mumat_name, msave, mpsi_mat, mtheta_mat, mphi_mat, mnprops, mprops);
    id = mid;
	nstatev = mnstatev;
    conductivity = mconductivity;
    density = mdensity;
}
    
//-------------------------------------------------------------
    void aba_material::update(const material_characteristics &mc, const double &mconductivity, const double &mdensity, const state_variables &sv, const int &mid)
//-------------------------------------------------------------
{

    number = mc.number;
    id = mid;
    umat_name = mc.umat_name;
    save = mc.save;
    
    psi_mat = mc.psi_mat;
    theta_mat = mc.theta_mat;
    phi_mat = mc.phi_mat;
    
    conductivity = mconductivity;
    density = mdensity;
    
    nprops = mc.nprops;
    nstatev = sv.nstatev;
    props = mc.props;
}
    
//----------------------------------------------------------------------
aba_material& aba_material::operator = (const aba_material& mc)
//----------------------------------------------------------------------
{
	assert(mc.nprops);

	number = mc.number;
	id = mc.id;
	umat_name = mc.umat_name;
    save = mc.save;
    
    psi_mat = mc.psi_mat;
	theta_mat = mc.theta_mat;
	phi_mat = mc.phi_mat;
    
    conductivity = mc.conductivity;
    density = mc.density;
		
	nprops = mc.nprops;
    nstatev = mc.nstatev;
	props = mc.props;
    
	return *this;
}
    
//-------------------------------------------------------------
void aba_material::write(const unsigned int &loading_type, const string &path_data, const string &inputfile, const bool &append)
//-------------------------------------------------------------
{
    std::string filename = path_data + "/" + inputfile;
    std::ofstream param_aba;
    
    if(append == true) {
        param_aba.open(filename, ios::app);
    }
    else {
        param_aba.open(filename, ios::out);
    }
    
    param_aba << "*Material, name=" << umat_name << "-" << id << "\n";
    if (loading_type==1) {
        param_aba << "*Depvar\n     ";
        param_aba << nstatev+4 << "\n";
    }
    else if (loading_type==2){
        
        param_aba << "*Depvar\n     ";
        param_aba << nstatev+7 << "\n";
    }
    else {
        cout << "error in aba_material::write (/Continuum_mechanics/Unit_cell/materials.cpp) : loading_type should take a value that is either 1 (Mechanical loading) or 2 (Thermomechanical loading)" << endl;
    }
    param_aba << "*User Material, constants=" << nprops << "\n";
    param_aba << "**" << endl;
    
    int z=0;
    for (int i=0; i<nprops; i++) {
        param_aba << props(i);
        z++;
        if ((z==8)||(i+1==nprops)) {
            param_aba << "\n";
            z=0;
        }
        else {
            param_aba << ",";
        }
    }
    param_aba << "**" << endl;
    param_aba.close();
}

//--------------------------------------------------------------------------
ostream& operator << (ostream& s, const aba_material& mc)
//--------------------------------------------------------------------------
{
	assert(mc.nprops);

	s << "Display state variables\n";
	s << "Number of the phase: " << mc.number << "\n";
	s << "Name of the material = " << mc.umat_name << "\n";
	s << "local material orientation: psi = " << mc.psi_mat << "\t theta = " << mc.theta_mat << "\t phi = " << mc.phi_mat << "\n";
	s << "nprops: \n" << mc.nprops << "\n";
	s << "nstatev: \n" << mc.nstatev << "\n";
	s << "props: \n";
    s << mc.props.t();
	s << "\n";

	s << "\n\n";

	return s;
}

} //namespace simcoon
