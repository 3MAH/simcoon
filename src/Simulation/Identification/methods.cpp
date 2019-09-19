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

///@file methods.cpp
///@brief methods for genetic algorithm (among others)
///@author Chemisky & Despringre
///@version 1.0

#include <iostream>
#include <fstream>
#include <assert.h>
#include <math.h>
#include <armadillo>
#include <simcoon/Simulation/Maths/random.hpp>
#include <simcoon/Simulation/Maths/lagrange.hpp>
#include <simcoon/Simulation/Identification/parameters.hpp>
#include <simcoon/Simulation/Identification/methods.hpp>
#include <simcoon/Simulation/Identification/generation.hpp>
#include <simcoon/Simulation/Identification/optimize.hpp>

using namespace std;
using namespace arma;

namespace simcoon{
    
//Genetic method
void genetic(generation &gen_g, generation &gensons, int &idnumber, const double &probaMut, const double &pertu, const vector<parameters> &params){
    
    int n_param = params.size();
    int maxpop = gensons.size();
    
    //Generate two genitoers
	individual dad(n_param, 0, 0.);
	individual mom(n_param, 0, 0.);

    int chromosome = 0;
    //Very small pertubation
    
    gensons.newid(idnumber);
    for(int i=0; i<maxpop; i++) {
        /// Random determination of "father" and "mother"
        dad = gen_g.pop[alea(maxpop-1)];
        mom = gen_g.pop[alea(maxpop-1)];
        while(dad.id==mom.id)
            mom = gen_g.pop[alea(maxpop-1)];
            
        for(int j=0; j<n_param; j++) {
            chromosome = alea(1);
            if(chromosome==0) {
                gensons.pop[i].p(j)=dad.p(j)*alead(1.-pertu,1.+pertu);
            }
            else {
                gensons.pop[i].p(j)=mom.p(j)*alead(1.-pertu,1.+pertu);
            }
                
            if (gensons.pop[i].p(j) > params[j].max_value)
                gensons.pop[i].p(j) = params[j].max_value;
            if (gensons.pop[i].p(j) < params[j].min_value)
                gensons.pop[i].p(j) = params[j].min_value;
            
            ///Apply a mutation
            if (alea(99)<probaMut)
                gensons.pop[i].p(j) = alead(params[j].min_value, params[j].max_value);
        }
    }
    
}

void find_best(generation &gen_cur, generation &gboys_cur, const generation &gen_old, const generation &gboys_old, const generation &gensons, const int &maxpop, const int &n_param, int& id0) {
    
    generation genall((maxpop > 1) ?  2*maxpop : maxpop, n_param, id0);
    
    for(int i=0; i<gboys_old.size(); i++) {
        genall.pop[i] = gboys_old.pop[i];
    }
    for(int i=gboys_old.size(); i<maxpop; i++) {
        genall.pop[i]=gen_old.pop[i];
    }
    
    if (maxpop > 1) {
        for(int i=0; i<gensons.size(); i++) {
            genall.pop[i+maxpop]=gensons.pop[i];
        }
        genall.classify();
    }

    gen_cur.construct(maxpop, n_param, id0, 0.);
    
    if(gboys_old.size()) {
        gboys_cur.construct(gboys_old.size(), n_param, id0, 0.);
    }
    
    for(int i=0; i<maxpop; i++) {
        gen_cur.pop[i] = genall.pop[i];
    }
    
    for(int i=0; i<gboys_cur.size(); i++) {
        gboys_cur.pop[i] = genall.pop[i];
    }    
}
    
void write_results(ofstream &result, const string &outputfile, const generation &gen_cur, const int &g, const int &maxpop, const int &n_param) {
    
    result.open(outputfile,  ios::out | ios::app);
    for(int i=0; i<maxpop; i++) {
        
        result << g << "\t" << gen_cur.pop[i].id << "\t" << gen_cur.pop[i].cout << "\t";
        for(int j=0; j<n_param;j++) {
            result << gen_cur.pop[i].p(j);
            if(j==n_param-1)
                result << "\n";
            else
                result << "\t";
        }
    }
    result.close();
}

} //namespace simcoon
