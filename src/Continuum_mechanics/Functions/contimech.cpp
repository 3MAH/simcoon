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

#include <iostream>
#include <fstream>
#include <assert.h>
#include <math.h>
#include <armadillo>
#include <Fastor/Fastor.h>
#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Continuum_mechanics/Functions/transfer.hpp>
#include <simcoon/Continuum_mechanics/Functions/fastor_bridge.hpp>

using namespace std;
using namespace arma;

namespace simcoon{

namespace {

// dev(vec) of a Voigt 6-vector and its partner with shears scaled by s_shear (2 stress, 0.5 strain)
inline void dev_voigt(const vec &v, const double &s_shear, double *d, double *d2) {
    const double sph = (1./3.)*(v(0) + v(1) + v(2));
    for (int i=0; i<3; i++) {
        d[i] = v(i) - sph;
        d2[i] = d[i];
    }
    for (int i=3; i<6; i++) {
        d[i] = v(i);
        d2[i] = s_shear*v(i);
    }
}

// v and its partner with shears scaled by s_shear (2 stress, 0.5 strain)
inline void shear_scaled(const vec &v, const double &s_shear, double *v1, double *v2) {
    for (int i=0; i<6; i++) {
        v1[i] = v(i);
        v2[i] = (i < 3) ? v1[i] : s_shear*v1[i];
    }
}

// sum(a % b) on 6 entries, in the accumulation order of Armadillo's accu (two interleaved partial sums)
inline double dot_voigt(const double *a, const double *b) {
    double val1 = 0.;
    double val2 = 0.;
    for (int i=0; i<6; i+=2) {
        val1 += (a[i] * b[i]);
        val2 += (a[i+1] * b[i+1]);
    }
    return (val1 + val2);
}

} // namespace

mat dev(const mat &m) {
    if ((m.n_rows != 3) || (m.n_cols != 3)) {
        return m - sph(m);
    }
    // m - sph(m) entry by entry; off-diagonals subtract 0*tr/3 to keep signed zeros and non-finite traces
    double val1 = 0.;
    double val2 = 0.;
    val1 += m(0,0);
    val2 += m(1,1);
    val1 += m(2,2);
    const double k = (1./3.)*(val1 + val2);
    const double zero_k = 0.*k;
    mat d(3, 3, fill::none);
    for (int j=0; j<3; j++) {
        for (int i=0; i<3; i++) {
            d(i,j) = m(i,j) - ((i == j) ? k : zero_k);
        }
    }
    return d;
}

mat sph(const mat &m) {
    assert((m.n_cols==3)&&(m.n_rows==3));
    
    return (1./3.)*(trace(m))*eye(3,3);
}

double tr(const vec &v) {
	assert(v.size()==6);
	double tr = v(0) + v(1) + v(2);
	
	return tr;
}

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

double Mises_stress(const vec &v) {
	assert(v.size()==6);

	double vdev[6], vdev2[6];
	dev_voigt(v, 2., vdev, vdev2);

	return sqrt(3./2.*dot_voigt(vdev, vdev2));
}

vec eta_stress(const vec &v) {
	assert(v.size()==6);
	
	double vdev[6], vdev2[6];
	dev_voigt(v, 2., vdev, vdev2);
		
	double n = sqrt(3./2.*dot_voigt(vdev, vdev2));

	if (n > 0.) {
		// (3/2) vdev2 (1/n), evaluated as the expression template did
		const double inv_n = 1./n;
		vec eta(6, fill::none);
		for (int i=0; i<6; i++)
			eta(i) = ((3./2.)*vdev2[i])*inv_n;
		return eta;
	}
	else {
		return zeros(6);
	}

}
    
mat deta_stress(const vec &v) {
	assert(v.size()==6);

	double vdev[6], vdev2[6];
	dev_voigt(v, 2., vdev, vdev2);
	double sigeq = sqrt(3./2.*dot_voigt(vdev, vdev2));
	if (sigeq < simcoon::iota) {
		return zeros(6,6);
	}
	double eta[6];
	const double inv_n = 1./sigeq;
	for (int i=0; i<6; i++)
		eta[i] = ((3./2.)*vdev2[i])*inv_n;
	// P_Mises in simcoon's engineering-shear Voigt basis: normal block = (3/2) dev,
	// shear diagonal = 3 (the factor 2 is already embedded, as in eta_stress / Mises_stress).
	// the product stays its own statement: rounded before the subtraction, as in the matrix temporary
	static const double Pm[6][6] = { {  1., -0.5, -0.5, 0., 0., 0.},
	                                 {-0.5,   1., -0.5, 0., 0., 0.},
	                                 {-0.5, -0.5,   1., 0., 0., 0.},
	                                 {  0.,   0.,   0., 3., 0., 0.},
	                                 {  0.,   0.,   0., 0., 3., 0.},
	                                 {  0.,   0.,   0., 0., 0., 3.} };
	mat deta(6, 6, fill::none);
	for (int j=0; j<6; j++) {
		for (int i=0; i<6; i++) {
			const double eta_ij = eta[i]*eta[j];
			deta(i,j) = (Pm[i][j] - eta_ij)/sigeq;
		}
	}
	return deta;
}

vec eta_norm_stress(const vec &v) {
    assert(v.size()==6);
    
    double v1[6], v2[6];
    shear_scaled(v, 2., v1, v2);
    
    double n = sqrt(dot_voigt(v1, v2));
    
    if (n > 0.) {
        const double inv_n = 1./n;
        vec eta(6, fill::none);
        for (int i=0; i<6; i++)
            eta(i) = v2[i]*inv_n;
        return eta;
    }
    else {
        return zeros(6);
    }
}

vec eta_norm_strain(const vec &v) {
    assert(v.size()==6);
    
    double v1[6], v2[6];
    shear_scaled(v, 0.5, v1, v2);
    
    double n = sqrt(dot_voigt(v1, v2));
    
    if (n > 0.) {
        const double inv_n = 1./n;
        vec eta(6, fill::none);
        for (int i=0; i<6; i++)
            eta(i) = v1[i]*inv_n;
        return eta;
    }
    else {
        return zeros(6);
    }
}
    
double norm_stress(const vec &v) {
    
    assert(v.size()==6);
    
    double v1[6], v2[6];
    shear_scaled(v, 2., v1, v2);
    
    return sqrt(dot_voigt(v1, v2));
}

double norm_strain(const vec &v) {
    
    assert(v.size()==6);
    
    double v1[6], v2[6];
    shear_scaled(v, 0.5, v1, v2);
    
    return sqrt(dot_voigt(v1, v2));
}

double Mises_strain(const vec &v) {
	assert(v.size()==6);

	double vdev[6], vdev2[6];
	dev_voigt(v, 0.5, vdev, vdev2);

	return sqrt(2./3.*dot_voigt(vdev, vdev2));
}

vec eta_strain(const vec &v) {
	assert(v.size()==6);
	
	double vdev[6], vdev2[6];
	dev_voigt(v, 0.5, vdev, vdev2);
		
	double n = sqrt(2./3.*dot_voigt(vdev, vdev2));

	if (n > 0.) {
		// (2/3) vdev (1/n), evaluated as the expression template did
		const double inv_n = 1./n;
		vec eta(6, fill::none);
		for (int i=0; i<6; i++)
			eta(i) = ((2./3.)*vdev[i])*inv_n;
		return eta;
	}
	else {
		return zeros(6);
	}

}

double J2_stress(const vec &v) {
	assert(v.size()==6);

	double vdev[6], vdev2[6];
	dev_voigt(v, 2., vdev, vdev2);

	return 0.5*dot_voigt(vdev, vdev2);
}


double J2_strain(const vec &v) {
	assert(v.size()==6);

	double vdev[6], vdev2[6];
	dev_voigt(v, 0.5, vdev, vdev2);

	return 0.5*dot_voigt(vdev, vdev2);
}

double J3_stress(const vec &v) {
    vec vdev = dev(v);
    mat mat1 = v2t_stress(vdev);
    
    return (1./3.)*accu(mat1%(mat1*mat1));
}

double J3_strain(const vec &v) {
    vec vdev = dev(v);
    mat mat1 = v2t_strain(vdev);
    
    return (1./3.)*accu(mat1%(mat1*mat1));
}
    
double Macaulay_p(const double &d) {

    if (d >= 0)
        return d;
    else
        return 0.;
}
    
double Macaulay_n(const double &d) {

    if (d <= 0)
        return d;
    else
        return 0.;
}
    
double sign(const double &d) {
    
    if ((d < simcoon::iota)&&(fabs(d) > simcoon::iota))
        return -1.;
    else if(d > simcoon::iota)
        return 1.;
    else
        return 0.;
}    

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

double curvature_ellipsoid(const double &u, const double &v, const double &a1, const double &a2, const double &a3)
{
	double curvature = (a1*a1*a2*a2*a3*a3)/(pow((a1*a1*a2*a2*cos(v)*cos(v) + a3*a3*sin(v)*sin(v) * (a2*a2*cos(u)*cos(u) + a1*a1*sin(u)*sin(u))),2));
	
	return curvature;
}

vec sigma_int(const vec &sigma_in, const double &u, const double &v, const double &a1, const double &a2, const double &a3)
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

// rank-one product by hand: no BLAS dispatch for 36 multiplications
static mat outer_voigt(const vec &a, const vec &b) {
    mat C(6, 6, fill::none);
    for (int j=0; j<6; j++) {
        for (int i=0; i<6; i++) {
            C(i,j) = a(i)*b(j);
        }
    }
    return C;
}

mat auto_sym_dyadic(const mat &A) {
    vec A_v = t2v_sym(A);

	return outer_voigt(A_v, A_v);
}

mat sym_dyadic(const mat &A, const mat &B) {

    vec A_v = t2v_sym(A);
    vec B_v = t2v_sym(B);

	return outer_voigt(A_v, B_v);
}

mat auto_dyadic(const mat &A) {
    return dyadic(A, A);
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
    // C_ijkl = A_ij * B_kl — full outer product, then extract Voigt 6x6
    auto A_ = arma_to_fastor2(mat::fixed<3,3>(A));  // auto-detects symmetry
    auto B_ = arma_to_fastor2(mat::fixed<3,3>(B));
    enum {i,j,k,l};
    auto C = Fastor::einsum<Fastor::Index<i,j>, Fastor::Index<k,l>>(A_, B_);
    return fastor4_to_voigt(C);
}

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

void add_linearop_eigsym(mat::fixed<6,6> &acc, const double &c, const vec::fixed<3> &b_i, const vec::fixed<3> &b_j) {
	static const int Id[3][3] = {{0, 3, 4}, {3, 1, 5}, {4, 5, 2}};
	double Bij[3][3];
	for (int i=0; i<3; i++)
		for (int j=0; j<3; j++)
			Bij[i][j] = b_i(i)*b_j(j);
	for (int i=0; i<3; i++)
		for (int j=i; j<3; j++)
			for (int k=0; k<3; k++)
				for (int l=k; l<3; l++)
					acc(Id[i][j], Id[k][l]) += c*(0.5*(Bij[i][j]*Bij[k][l] + Bij[i][j]*Bij[l][k]));
}

mat linearop_eigsym(const vec &b_i, const vec &b_j) {
	mat::fixed<6,6> BBBB(fill::zeros);
	add_linearop_eigsym(BBBB, 1., b_i, b_j);
	return BBBB;
}



} //namespace simcoon
