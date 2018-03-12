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

///@file output.cpp
///@brief Function to output phases data

#include <armadillo>

using namespace std;
using namespace arma;

namespace simcoon{

void output_modulus(const mat &L, const int &i0, vec &statev) {
	int iterator =0 ;
	for (unsigned int i = 0; i < L.n_rows; i++) {
		for (unsigned int j = 0; j < L.n_cols; j++) {
			statev(i0 + iterator) = L(i,j);
            iterator += 1;
		}
	}
}

void output_moduli(const vector<mat> &list, const int &i0, vec &statev) {
	int iterator =0 ;
	mat L;
	for (unsigned int k = 0; k < list.size(); k++) {
		L = list[k];
		for (unsigned int i = 0; i < L.n_rows; i++) {
			for (unsigned int j = 0; j < L.n_cols; j++) {
				statev(i0 + iterator) = L(i,j);
				iterator += 1;
			}
		}
	}
	
}

} //namespace simcoon
