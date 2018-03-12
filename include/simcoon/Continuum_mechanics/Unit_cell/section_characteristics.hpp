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

///@file phase_characteristics.hpp
///@brief Characteristics of a phase, which hereditates from:
///-material characteristics
///@version 1.0

#pragma once

#include <iostream>
#include <string>
#include <memory>
#include <armadillo>
#include <simcoon/Simulation/Phase/phase_characteristics.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/materials.hpp>

namespace simcoon{

//======================================
class section_characteristics
//======================================
{
	private:

	protected:

	public :

        std::string elset_name;
        int id;
    
        aba_material abamat;
        std::vector<section_characteristics> sub_sections;
        std::string sub_sections_file;
    
		section_characteristics(); 	//default constructor
    
        section_characteristics(const std::string &, const int &, const aba_material &, const std::string &);
    
		section_characteristics(const section_characteristics&);	//Copy constructor
        virtual ~section_characteristics();
    
        virtual void update(const std::string &, const int &, const phase_characteristics &);

		virtual section_characteristics& operator = (const section_characteristics&);
        friend std::ostream& operator << (std::ostream&, const section_characteristics&);
};

} //namespace simcoon
