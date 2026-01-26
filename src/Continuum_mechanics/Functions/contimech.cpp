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
#include <simcoon/FTensor.hpp>
#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Continuum_mechanics/Functions/transfer.hpp>

using namespace std;
using namespace arma;
using namespace FTensor;

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
    
    if ((d < simcoon::iota)&&(fabs(d) > simcoon::iota))
        return -1.;
    else if(d > simcoon::iota)
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
vec sigma_int(const vec &sigma_in, const double &u, const double &v, const double &a1, const double &a2, const double &a3)
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

mat auto_sym_dyadic(const mat &A) {
//T
    vec A_v = t2v_sym(A);

	return A_v * A_v.t();
}

mat sym_dyadic(const mat &A, const mat &B) {

    vec A_v = t2v_sym(A);
    vec B_v = t2v_sym(B);
	
	return A_v * B_v.t();
}

mat auto_dyadic(const mat &A) {

    Tensor2<double,3,3> A_ = mat_FTensor2(A);
    Tensor4<double,3,3,3,3> C_;
    
    Index<'i', 3> i;
    Index<'j', 3> j;
    Index<'k', 3> k;
    Index<'l', 3> l;
        
    C_(i,j,k,l) = A_(i,j)*A_(k,l);
    return FTensor4_mat(C_);
    
}

mat dyadic_4vectors_sym(const vec &n_a, const vec &n_b, const std::string  &conv) {

	mat C = zeros(6,6);

	int ij=0;
	int kl=0;
    
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

	if (conv == "aabb") {

		for (int i=0; i<3; i++) {
			for (int j=i; j<3; j++) {
				ij = Id(i,j);
				for (int k=0; k<3; k++) {
					for (int l=k; l<3; l++) {
						kl = Id(k,l);
						C(ij,kl) = n_a(i)*n_a(j)*n_b(k)*n_b(l);
					}
				}
			}
		}
		return C;
	}
	else if (conv == "abab") {

		for (int i=0; i<3; i++) {
			for (int j=i; j<3; j++) {
				ij = Id(i,j);
				for (int k=0; k<3; k++) {
					for (int l=k; l<3; l++) {
						kl = Id(k,l);
						C(ij,kl) = 0.5*(n_a(i)*n_b(j)*n_a(k)*n_b(l) + n_a(i)*n_b(j)*n_b(k)*n_a(l));
					}
				}
			}
		}
		return C;		
		
	}
	else {
	    throw std::invalid_argument("conv string must be either aabb or abab");
	}
}

mat dyadic(const mat &A, const mat &B) {
            
    Tensor2<double,3,3> A_ = mat_FTensor2(A);
    Tensor2<double,3,3> B_ = mat_FTensor2(B);
    Tensor4<double,3,3,3,3> C_;
    
    Index<'i', 3> i;
    Index<'j', 3> j;
    Index<'k', 3> k;
    Index<'l', 3> l;
        
    C_(i,j,k,l) = A_(i,j)*B_(k,l);
    return FTensor4_mat(C_);    
}

///This computes the symmetric 4th-order dyadic product A o A = 0.5*(A(i,k)*A(j,l) + A(i,l)*A(j,k));
mat auto_sym_dyadic_operator(const mat &A) {

	mat C = zeros(6,6);

	int ij=0;
	int kl=0;
    
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
					C(ij,kl) = 0.5*(A(i,k)*A(j,l) + A(i,l)*A(j,k));
				}
			}
		}
	}
	
	return C;
}

///This computes the symmetric 4th-order dyadic product A o B = 0.5*(A(i,k)*B(j,l) + A(i,l)*B(j,k));
mat sym_dyadic_operator(const mat &A, const mat &B) {

	mat C = zeros(6,6);

	int ij=0;
	int kl=0;
    
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
					C(ij,kl) = 0.5*(A(i,k)*B(j,l) + A(i,l)*B(j,k));
				}
			}
		}
	}
	return C;
}

///This computes the operator BBBB such that B_i x D x B_j = BBBB x D, considering B_i and B_j are projection tensors of the vectors b_i and b_j
mat B_klmn(const vec &b_i, const vec &b_j) {

	mat Bij = b_i*b_j.t();
	mat BBBB = zeros(6,6);

	int ij=0;
	int kl=0;
    
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
					BBBB(ij,kl) = 0.5*(Bij(i,j)*Bij(k,l) + Bij(i,j)*Bij(l,k));
				}
			}
		}
	}
	
	return BBBB;
}


/*mat eulerian_determinant(const mat &A) {
	mat Id = eye(3,3);
	vec A_Id = A*Ith();
	double Id_A_Id = sum(Ith()%A_Id);
	return A - (1./3.)*(sym_dyadic(Id,A_Id)+sym_dyadic(A_Id,Id))+(1./9.)*Id_A_Id*auto_sym_dyadic(Id);

}*/

} //namespace simcoon
