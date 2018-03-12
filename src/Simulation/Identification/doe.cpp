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

///@file doe.cpp
///@brief Design of Experiments library

#include <iostream>
#include <fstream>
#include <assert.h>
#include <math.h>
#include <armadillo>
#include <simcoon/Simulation/Maths/random.hpp>
#include <simcoon/Simulation/Identification/doe.hpp>
#include <simcoon/Simulation/Identification/parameters.hpp>
#include <simcoon/Simulation/Identification/generation.hpp>
#include <simcoon/Simulation/Identification/read.hpp>

using namespace std;
using namespace arma;

namespace simcoon{

mat doe_uniform(const int &spop, const int &n_param, const vector<parameters> &params) {

    double pfactor = 0.;
    int pinc = 1;
    int z = 1;
    int pcol = 0;

    int n_samples = (int)pow(spop,n_param);
    mat doe = zeros(n_samples, n_param);

    ///Determination of parameters_equally_spaced
    for(int j=0; j<n_param; j++) {
        
        pcol = pow(spop,j);
        pinc=1;
        z=1;
        
        for(int i=0; i<n_samples; i++) {
            
            pfactor = (double) pinc/(spop+1);
            
            doe(i,j) = params[j].min_value + pfactor*(params[j].max_value-params[j].min_value);
            
            if (z==pcol) {
                if (pinc==spop) {
                    pinc = 0;
                }	  
                pinc++;
                z=0;
            }	
            z++;
        }
    }
    return doe;
}

mat doe_uniform_limit(const int &spop, const int &n_param, const vector<parameters> &params) {

    //This doe can only work if spop >=2
    assert(spop >= 2);
    
    double pfactor = 0.;
    int pinc = 0;
    int z = 1;
    int pcol = 0;

    int n_samples = (int)pow(spop,n_param);
    mat doe = zeros(n_samples, n_param);
    
     ///Determination of parameters_equally_spaced
     for(int j=0; j<n_param; j++) {
     
         pcol = pow(spop,j);
         pinc=0;
         z=1;
     
         for(int i=0; i<n_samples; i++) {
     
             pfactor = (double) pinc/(spop-1);
     
             doe(i,j) = params[j].min_value + pfactor*(params[j].max_value-params[j].min_value);
     
             if (z==pcol) {
					if (pinc==spop-1) {
                        pinc = -1;
					}
					pinc++;
					z=0;
             }
             z++;
         }
    }
    return doe;
}
    
mat doe_random(const int &n_samples, const int &n_param, const vector<parameters> &params) {
    
    mat doe = zeros(n_samples, n_param);
    
    for(int j=0; j<n_param; j++) {
        for(int i=0; i<n_samples; i++) {
            doe(i,j) = alead(params[j].min_value, params[j].max_value);
        }
    }
    
    return doe;
}
    
void gen_initialize(generation &geninit, int &spop, int &apop, int &idnumber, const int &aleaspace, const int &n_param, const vector<parameters> &params, const double &lambda) {
    
    if(aleaspace==0) {
        
        ///Populate the space with equidistant values. First generation
        int geninit_nindividuals = (int)pow(spop,n_param);
        geninit.construct(geninit_nindividuals, n_param, idnumber, lambda);
        
        ///Determination of parameters_equally_spaced
        mat samples = doe_uniform(spop, n_param, params);
        
        for(int j=0; j<n_param; j++) {
            for(int i=0; i<geninit.size(); i++) {
                geninit.pop[i].p(j) = samples(i,j);
            }
        }
        
    }
    if(aleaspace==1) {
        
        ///Populate the space with equidistant values. First generation
        int geninit_nindividuals = (int)pow(spop,n_param);
        geninit.construct(geninit_nindividuals, n_param, idnumber, lambda);
        
        ///Determination of parameters_equally_spaced
        mat samples = doe_uniform_limit(spop, n_param, params);
        
        for(int j=0; j<n_param; j++) {
            for(int i=0; i<geninit.size(); i++) {
                geninit.pop[i].p(j) = samples(i,j);
            }
        }
    }
    else if(aleaspace==2) {
        
        int geninit_nindividuals=apop;
        geninit.construct(geninit_nindividuals, n_param, idnumber, lambda);
                
        ///Determination of random values
        ///Determination of parameters_equally_spaced
        mat samples = doe_random(geninit_nindividuals, n_param, params);
        
        for(int j=0; j<n_param; j++) {
            for(int i=0; i<geninit.size(); i++) {
                geninit.pop[i].p(j) = samples(i,j);
            }
        }
    }
    else if(aleaspace==3) {
        
        mat samples;
        read_gen(apop, samples, n_param);
        int geninit_nindividuals=apop;
        geninit.construct(geninit_nindividuals, n_param, idnumber, lambda);
        
        for(int j=0; j<n_param; j++) {
            for(int i=0; i<geninit.size(); i++) {
                geninit.pop[i].p(j) = samples(i,j);
            }
        }
    }    
}
    
} //namespace simcoon
