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

///@file generation.cpp
///@brief generation for genetic algorithm (among others)
///@author Chemisky & Despringre
///@version 1.0

#include <iostream>
#include <assert.h>
#include <math.h>
#include <armadillo>
#include <simcoon/Simulation/Identification/individual.hpp>
#include <simcoon/Simulation/Identification/generation.hpp>

using namespace std;
using namespace arma;

namespace simcoon{

//=====Private methods for phases_characteristics===================================

//=====Public methods for phases_characteristics============================================

///@brief default constructor
//----------------------------------------------------------------------
generation::generation()
//----------------------------------------------------------------------
{

}

///@brief Constructor
///@param n : number of individuals
///@param init boolean that indicates if the constructor has to initialize  (default value is true)
//----------------------------------------------------------------------
generation::generation(const int &n, const int &m, int &idnumber, const double &mlambda)
//----------------------------------------------------------------------
{
	assert(n>0);
	
    for (int i=0; i<n; i++) {
        pop.push_back(individual(m, idnumber, mlambda));
        idnumber++;
    }
}

///@brief Copy constructor
///@param gp generation object to duplicate
//----------------------------------------------------------------------
generation::generation(const generation& gp)
//----------------------------------------------------------------------
{
    pop = gp.pop;
}

///@brief destructor
//----------------------------------------------------------------------
generation::~generation() {}
//----------------------------------------------------------------------

///@brief Construct : A method to construct the generation after its initial construction
//----------------------------------------------------------------------
void generation::construct(const int &n, const int &m, int &idnumber, const double &mlambda)
//----------------------------------------------------------------------
{
	assert(n>0);

    for (int i=0; i<n; i++) {
        pop.push_back(individual(m, idnumber, mlambda));
        idnumber++;
    }
}

///@brief classify : A method to classify the individuals in a generation
//----------------------------------------------------------------------
void generation::classify()
//----------------------------------------------------------------------
{
    
	assert(pop.size()>0);

	double mini = -1.;
	int posmini = 0;
	individual temp;
    unsigned int number_not_NaN = pop.size();
	
    for(unsigned int i=0; i < number_not_NaN; i++) {
        int check_nan = 0;
        if (std::isnan(pop[i].cout)) {
            check_nan++;
        }
        
        if(check_nan > 0) {
            temp = pop[i];
            pop.erase(pop.begin() + i);
            pop.push_back(temp);
            number_not_NaN--;
            i--;
        }
    }
    
    for(unsigned int i=0; i < pop.size(); i++) {
        pop[i].rank=i+1;
    }
    
	for(unsigned int i=0; i < number_not_NaN-1; i++) {
		mini=pop[i].cout;
		posmini=i;
		
		for(unsigned int j=i; j < number_not_NaN; j++) {
            if(pop[j].cout < mini) {
				mini=pop[j].cout;
				posmini=j;
			}
		}
		temp=pop[posmini];
		pop[posmini]=pop[i];
		pop[i]=temp;
	} 	

	for(unsigned int i=0; i < number_not_NaN; i++) {
		pop[i].rank=i+1;
	}

}

///@brief newid : A method to assign new id's to a generation
//----------------------------------------------------------------------
void generation::newid(int &idnumber)
//----------------------------------------------------------------------
{
    
    assert(pop.size()>0);

    for (unsigned int i=0; i<pop.size(); i++) {
		pop[i].id = idnumber;
        idnumber++;
    }
}

///@brief newid : A method to destruct and free memory of a generation
//----------------------------------------------------------------------
void generation::destruct() {}
//----------------------------------------------------------------------

///@brief Standard operator = for generation
//----------------------------------------------------------------------
generation& generation::operator = (const generation& gp)
//----------------------------------------------------------------------
{
    pop = gp.pop;    

	return *this;
}

///@brief Standard operator << for generation
//--------------------------------------------------------------------------
ostream& operator << (ostream& s, const generation& gp)
//--------------------------------------------------------------------------
{
	s << "Display info on the generation\n";
	s << "Number of individuals: " << gp.pop.size() << "\n";
	
	s << "Characteristics of each individual: \n";
    for(auto ind : gp.pop) {
		s << ind;
	}
	s << "\n\n";

	return s;
}

} //namespace simcoon