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

///@file section_characteristics.cpp
///@brief Characteristics of a phase, which hereditates from:
///- material_characteristics
///@version 1.0

#include <iostream>
#include <fstream>
#include <string>
#include <assert.h>
#include <armadillo>
#include <memory>
#include <simcoon/Simulation/Phase/phase_characteristics.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/section_characteristics.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/materials.hpp>

using namespace std;
using namespace arma;

namespace simcoon{

//=====Private methods for section_characteristics===================================

//=====Public methods for section_characteristics============================================

/*!
  \brief default constructor
*/

//-------------------------------------------------------------
section_characteristics::section_characteristics()
//-------------------------------------------------------------
{

}

/*!
  \brief Constructor with parameters
  \f$ \textbf{Examples :} \f$ \n
*/

//-------------------------------------------------------------
section_characteristics::section_characteristics(const std::string &melset_name, const int &mid, const aba_material &mabamat, const std::string &msub_sections_file)
//-------------------------------------------------------------
{
    elset_name = melset_name;
    id = mid;
    abamat = mabamat;
    sub_sections_file = msub_sections_file;
}
    
/*!
  \brief Copy constructor
  \param s section_characteristics object to duplicate
*/
    
//------------------------------------------------------
section_characteristics::section_characteristics(const section_characteristics& sc)
//------------------------------------------------------
{
    elset_name = sc.elset_name;
    id = sc.id;
    abamat = sc.abamat;
    
    sub_sections = sc.sub_sections;
    sub_sections_file = sc.sub_sections_file;
}

/*!
  \brief Destructor

  Deletes section_characteristics, the shared ptr related object will be detroyed automatically if they is no pointer to it.
*/

//-------------------------------------
section_characteristics::~section_characteristics() {}
//-------------------------------------

/*!
  \brief Standard operator = for section_characteristics
*/
  
    
//-------------------------------------------------------------
void section_characteristics::update(const std::string &melset_name, const int &mid, const phase_characteristics &pc)
//-------------------------------------------------------------
{
    elset_name = melset_name;
    id = mid;
    //Note : there is no density nor conductivity here since here section_characteristics inheritates from phases characteristics
    abamat.update(*pc.sptr_matprops, 0., 0., *pc.sptr_sv_global, mid);
    
    int nphases = pc.sub_phases.size();
    sub_sections.resize(nphases);
    
    for (int i=0; i<nphases; i++) {
        string melset_name_sub = melset_name + to_string(i);
        int mid = 100000 + pc.sptr_matprops->props(1)*1000 + pc.sub_phases[i].sptr_matprops->number;
        sub_sections[i].update(melset_name_sub, mid, pc.sub_phases[i]);
    }
    sub_sections_file = pc.sub_phases_file;
}

//----------------------------------------------------------------------
section_characteristics& section_characteristics::operator = (const section_characteristics& sc)
//----------------------------------------------------------------------
{
    elset_name = sc.elset_name;
    id = sc.id;
    abamat = sc.abamat;
    
    sub_sections = sc.sub_sections;
    sub_sections_file = sc.sub_sections_file;
    
    return *this;
}
    
//--------------------------------------------------------------------------
ostream& operator << (ostream& s, const section_characteristics& sc)
//--------------------------------------------------------------------------
{
	s << "Display section characteristics:\n";
    s << "display Element set name: " << sc.elset_name << "\n";
    s << "display section id: " << sc.id << "\n";
	s << "Display material characteristics:\n";
    s << sc.abamat;
    
    for(auto r : sc.sub_sections) {
        s << r;
    }
    
    s << "\n\n";

	return s;
}

} //namespace simcoon
