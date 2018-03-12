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

///@file random.cpp
///@brief random number generators
///@version 1.0

#include <iostream>
#include <math.h>
#include <assert.h>
#include <armadillo>

using namespace std;
using namespace arma;

namespace simcoon{

//This function returns a random in number between 0 and a
int alea(const int &n)
{ 
	assert (0 < n && n <= RAND_MAX); 
	int partSize = n == RAND_MAX ? 1 : 1 + (RAND_MAX-n)/(n+1); 
	int maxUsefull = partSize * n + (partSize-1); 
	int draw; 
 
	do { 
		draw = rand(); 
	} while (draw > maxUsefull);
		 
	return draw/partSize; 
}

//This function returns a random in number between a and b
int aleab(const int &a, const int &b)
{ 
	int n=b-a;
	return alea(n) + a; 
}

//This function returns a random double number between a and b
double alead(const double &a, const double &b){
  
    return ((double)rand()/((double)RAND_MAX)) * (b-a) + a;
}

} //namespace simcoon
