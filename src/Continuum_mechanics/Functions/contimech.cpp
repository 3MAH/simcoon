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

///@file contimech.cpp
///@brief Functions that computes Mises stress/strains, directions, etc
///@version 1.0

#include <iostream>
#include <fstream>
#include <assert.h>
#include <math.h>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Continuum_mechanics/Functions/transfer.hpp>

using namespace std;
using namespace arma;

namespace simcoon{

//This function returns the deviatoric part of m
mat dev(const mat &m) {
    assert((m.n_cols==3)&&(m.n_rows==3));

    return m - sph(m);
}

//This function returns the spherical part of m
mat sph(const mat &m) {
    assert((m.n_cols==3)&&(m.n_rows==3));
    
    return (1./3.)*(trace(m))*eye(3,3);
}

//This function returns the trace of the tensor v
double tr(const vec &v) {
	assert(v.size()==6);
	double tr = v(0) + v(1) + v(2);
	
	return tr;
}

//This function returns the deviatoric part of v
vec dev(const vec &v) {
	assert(v.size()==6);

	vec vdev(6);
	double sph = (1./3.)*(v(0) + v(1) + v(2));
	for (int i=0; i<3; i++) {
		vdev(i) = v(i) - sph; 
	}
	for (int i=3; i<6; i++) {
		vdev(i) = v(i); 
	}	

	return vdev;
}

//This function determines the Mises equivalent of a stress tensor, according to the Voigt convention for stress 
double Mises_stress(const vec &v) {
	assert(v.size()==6);

	vec vdev = dev(v);
	vec vdev2 = vdev;
	for (int i=3; i<6; i++)
		vdev2(i) = 2.*vdev2(i);	

	return sqrt(3./2.*sum(vdev%vdev2));
}

//This function determines the strain flow (direction) from a stress tensor (Mises convention), according to the Voigt convention for strains
vec eta_stress(const vec &v) {
	assert(v.size()==6);
	
	vec vdev = dev(v);
	vec vdev2 = vdev;	
	for (int i=3; i<6; i++)
		vdev2(i) = 2.*vdev2(i);	
		
	double n = sqrt(3./2.*sum(vdev%vdev2));

	if (n > 0.) {		
		return (3./2.)*vdev2*(1./n);
	}
	else {
		return zeros(6);
	}

}
    
//This function determines the strain flow (direction) from a stress tensor, according to the Voigt convention for strains
vec eta_norm_stress(const vec &v) {
    assert(v.size()==6);
    
    vec v2 = v;
    for (int i=3; i<6; i++)
    v2(i) = 2.*v2(i);
    
    double n = sqrt(sum(v%v2));
    
    if (n > 0.) {
        return v2*(1./n);
    }
    else {
        return zeros(6);
    }
}

//This function determines the strain flow (direction) from a stress tensor, according to the Voigt convention for strains
vec eta_norm_strain(const vec &v) {
    assert(v.size()==6);
    
    vec v2 = v;
    for (int i=3; i<6; i++)
        v2(i) = 0.5*v2(i);
    
    double n = sqrt(sum(v%v2));
    
    if (n > 0.) {
        return v*(1./n);
    }
    else {
        return zeros(6);
    }
}
    
double norm_stress(const vec &v) {
    
    assert(v.size()==6);
    
    vec v2 = v;
    for (int i=3; i<6; i++)
    v2(i) = 2.*v2(i);
    
    return sqrt(sum(v%v2));
}

double norm_strain(const vec &v) {
    
    assert(v.size()==6);
    
    vec v2 = v;
    for (int i=3; i<6; i++)
        v2(i) = 0.5*v2(i);
    
    return sqrt(sum(v%v2));
}

//This function determines the Mises equivalent of a strain tensor, according to the Voigt convention for strains 
double Mises_strain(const vec &v) {
	assert(v.size()==6);

	vec vdev = dev(v);
	vec vdev2 = vdev;
	for (int i=3; i<6; i++)
		vdev2(i) = 0.5*vdev2(i);	

	return sqrt(2./3.*sum(vdev%vdev2));
}

//This function determines the strain flow (direction) from a strain tensor, according to the Voigt convention for strains
vec eta_strain(const vec &v) {
	assert(v.size()==6);
	
	vec vdev = dev(v);
	vec vdev2 = vdev;	
	for (int i=3; i<6; i++)
		vdev2(i) = 0.5*vdev2(i);	
		
	double n = sqrt(2./3.*sum(vdev%vdev2));

	if (n > 0.) {
		return (2./3.)*vdev*(1./n);
	}
	else {
		return zeros(6);
	}

}

//Returns the second invariant of the deviatoric part of a second order stress tensor written as a Voigt vector
double J2_stress(const vec &v) {
	assert(v.size()==6);

	vec vdev = dev(v);
	vec vdev2 = vdev;
	for (int i=3; i<6; i++)
		vdev2(i) = 2.*vdev2(i);	

	return 0.5*sum(vdev%vdev2);
}


//Returns the second invariant of the deviatoric part of a second order strain tensor written as a Voigt vector
double J2_strain(const vec &v) {
	assert(v.size()==6);

	vec vdev = dev(v);
	vec vdev2 = vdev;
	for (int i=3; i<6; i++)
		vdev2(i) = 0.5*vdev2(i);	

	return 0.5*sum(vdev%vdev2);
}

//Returns the third invariant of the deviatoric part of a second order stress tensor written as a Voigt vector
double J3_stress(const vec &v) {
    vec vdev = dev(v);
    mat mat1 = v2t_stress(vdev);
    
    return (1./3.)*accu(mat1%(mat1*mat1));
}

//Returns the third invariant of the deviatoric part of a second order stress tensor written as a Voigt vector
double J3_strain(const vec &v) {
    vec vdev = dev(v);
    mat mat1 = v2t_strain(vdev);
    
    return (1./3.)*accu(mat1%(mat1*mat1));
}
    
//This function returns the value if it's positive, zero if it's negative (Macaulay brackets <>+)
double Macaulay_p(const double &d) {

    if (d >= 0)
        return d;
    else
        return 0.;
}
    
//This function returns the value if it's negative, zero if it's positive (Macaulay brackets <>-)
double Macaulay_n(const double &d) {

    if (d <= 0)
        return d;
    else
        return 0.;
}
    
//This function returns the value if it's negative, zero if it's positive (Macaulay brackets <>-)
double sign(const double &d) {
    
    if ((d < sim_iota)&&(fabs(d) > sim_iota))
        return -1.;
    else if(d > sim_iota)
        return 1.;
    else
        return 0.;
}    

//Returns the normalized vector normal to an ellipsoid with semi-principal axes of length a1, a2, a3. The direction of the normalized vector is set by angles u
vec normal_ellipsoid(const double &u, const double &v, const double &a1, const double &a2, const double &a3) {

	vec normal = zeros(3);
	double x0 = a1 * cos(u)*sin(v);
	double y0 = a2 * sin(u)*sin(v);
	double z0 = a3 * cos(v);
	
	normal(0) = x0 / (a1*a1);
	normal(1) = y0 / (a2*a2);
	normal(2) = z0 / (a3*a3);
	
	double 	norme = norm(normal,2);
	for (int i=0; i<3; i++){
		normal(i)= normal(i)/norme;
	}
	return normal;
}

//-----------------------------------------------------------------------
double curvature_ellipsoid(const double &u, const double &v, const double &a1, const double &a2, const double &a3)
{
	double curvature = (a1*a1*a2*a2*a3*a3)/(pow((a1*a1*a2*a2*cos(v)*cos(v) + a3*a3*sin(v)*sin(v) * (a2*a2*cos(u)*cos(u) + a1*a1*sin(u)*sin(u))),2));
	
	return curvature;
}

//Returns the normal and tangent components of the stress vector in the normal direction n to an ellipsoid with axes a1, a2, a3. The direction of the normalized vector is set by angles u
vec sigma_int(const vec &sigma_in, const double &a1, const double &a2, const double &a3, const double &u, const double &v)
//-----------------------------------------------------------------------
{   
	mat s_in = v2t_stress(sigma_in);
	vec s_in_i = zeros(3);
	vec s_inter = zeros(2);
	
	vec normal = normal_ellipsoid(u,v,a1,a2,a3);
	s_in_i = s_in * normal;
	s_inter(0)= sum(s_in_i % normal);
	s_inter(1) = sqrt(fabs(norm(s_in_i,2)*norm(s_in_i,2)-s_inter(0)*s_inter(0)));
			
    return s_inter;
}

///This computes the Hill interfacial operator according to a normal a (see papers of Siredey and Entemeyer phD dissertation)
mat p_ikjl(const vec &a) {

	mat A = (a)*trans(a);
	mat F = zeros(6,6);

	int ij=0;
	int kl=0;
	int ik=0;
	int jl=0;
	int il=0;
	int jk=0;
    
    umat Id(3,3);
    Id(0,0) = 0;
    Id(0,1) = 3;
    Id(0,2) = 4;
    Id(1,0) = 3;
    Id(1,1) = 1;
    Id(1,2) = 5;
    Id(2,0) = 4;
    Id(2,1) = 5;
    Id(2,2) = 2;
	
	for (int i=0; i<3; i++) {
		for (int j=i; j<3; j++) {
			ij = Id(i,j);
			for (int k=0; k<3; k++) {
				for (int l=k; l<3; l++) {
					kl = Id(k,l);
					ik = Id(i,k);
					jl = Id(j,l);
					il = Id(i,l);
					jk = Id(j,k);
					
					F(ij,kl) += 0.5*(A(ik,jl)+A(jk,il));
				}
			}
		}
	}
	
	return F;
}

} //namespace simcoon
