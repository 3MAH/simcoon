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
#include <simcoon/Simulation/Maths/small_dense.hpp>
#include <simcoon/Simulation/Maths/small_voigt.hpp>
#include <simcoon/Continuum_mechanics/Functions/tensor.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Continuum_mechanics/Functions/transfer.hpp>
#include <simcoon/Continuum_mechanics/Functions/fastor_bridge.hpp>

using namespace std;
using namespace arma;

namespace simcoon{

namespace {

// flow c x / n with the shear entries scaled by shear (2 for a strain-like output), zero at n = 0
inline vec::fixed<6> flow(const double *x, const double n, const double c, const double shear) {
    if (n <= 0.) return vec::fixed<6>(fill::zeros);
    const double k = c/n;
    vec::fixed<6> eta;
    for (int i=0; i<3; i++) eta[i] = k*x[i];
    for (int i=3; i<6; i++) eta[i] = shear*k*x[i];
    return eta;
}

// P_Mises = 3/2 Idev2 in engineering Voigt (rows strain-like, columns stress-like), as constants: deta_stress runs per Newton iteration
constexpr double P_Mises[6][6] = {{  1., -0.5, -0.5, 0., 0., 0.},
                                  {-0.5,   1., -0.5, 0., 0., 0.},
                                  {-0.5, -0.5,   1., 0., 0., 0.},
                                  {  0.,   0.,   0., 3., 0., 0.},
                                  {  0.,   0.,   0., 0., 3., 0.},
                                  {  0.,   0.,   0., 0., 0., 3.}};

inline double norm_voigt(const vec &v, Tensor2Type t) {
    return sqrt(contract_voigt(v.memptr(), v.memptr(), voigt_shear_weight(t)));
}

inline double J2_voigt(const vec &v, Tensor2Type t) {
    double d[6];
    dev_voigt(v.memptr(), d);
    return 0.5*contract_voigt(d, d, voigt_shear_weight(t));
}

// J3 = det(dev) for a traceless tensor
inline double J3_voigt(const vec &v, Tensor2Type t) {
    double d[6];
    dev_voigt(v.memptr(), d);
    return det_sym_voigt(d, voigt_shear_in(t));
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

vec::fixed<6> dev(const vec &v) {
	assert(v.size()==6);
	vec::fixed<6> d;
	dev_voigt(v.memptr(), d.memptr());
	return d;
}

double Mises_stress(const vec &v) {
	assert(v.size()==6);
	double d[6];
	return mises_voigt(v.memptr(), Tensor2Type::stress, d);
}

vec::fixed<6> eta_stress(const vec &v) {
	assert(v.size()==6);
	double d[6];
	const double n = mises_voigt(v.memptr(), Tensor2Type::stress, d);
	return flow(d, n, 3./2., 2.);   // strain-like output
}

mat::fixed<6,6> deta_stress(const vec &v) {
	assert(v.size()==6);
	double d[6];
	const double sigeq = mises_voigt(v.memptr(), Tensor2Type::stress, d);
	mat::fixed<6,6> deta(fill::zeros);
	if (sigeq < simcoon::iota) return deta;
	const vec::fixed<6> eta = flow(d, sigeq, 3./2., 2.);
	const double inv_sigeq = 1./sigeq;
	for (int j=0; j<6; j++)
		for (int i=0; i<6; i++)
			deta.at(i,j) = (P_Mises[i][j] - eta[i]*eta[j])*inv_sigeq;
	return deta;
}

vec::fixed<6> eta_norm_stress(const vec &v) {
    assert(v.size()==6);
    return flow(v.memptr(), norm_voigt(v, Tensor2Type::stress), 1., 2.);   // strain-like output
}

vec::fixed<6> eta_norm_strain(const vec &v) {
    assert(v.size()==6);
    return flow(v.memptr(), norm_voigt(v, Tensor2Type::strain), 1., 1.);   // kept in the strain convention of the input
}

double norm_stress(const vec &v) {
    assert(v.size()==6);
    return norm_voigt(v, Tensor2Type::stress);
}

double norm_strain(const vec &v) {
    assert(v.size()==6);
    return norm_voigt(v, Tensor2Type::strain);
}

double Mises_strain(const vec &v) {
	assert(v.size()==6);
	double d[6];
	return mises_voigt(v.memptr(), Tensor2Type::strain, d);
}

vec::fixed<6> eta_strain(const vec &v) {
	assert(v.size()==6);
	double d[6];
	const double n = mises_voigt(v.memptr(), Tensor2Type::strain, d);
	return flow(d, n, 2./3., 1.);   // kept in the strain convention of the input
}

double J2_stress(const vec &v) {
	assert(v.size()==6);
	return J2_voigt(v, Tensor2Type::stress);
}

double J2_strain(const vec &v) {
	assert(v.size()==6);
	return J2_voigt(v, Tensor2Type::strain);
}

double J3_stress(const vec &v) {
    assert(v.size()==6);
    return J3_voigt(v, Tensor2Type::stress);
}

double J3_strain(const vec &v) {
    assert(v.size()==6);
    return J3_voigt(v, Tensor2Type::strain);
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

mat auto_sym_dyadic(const mat &A) {
    vec A_v = t2v_sym(A);

	return outer6(A_v.memptr(), A_v.memptr());
}

mat sym_dyadic(const mat &A, const mat &B) {

    vec A_v = t2v_sym(A);
    vec B_v = t2v_sym(B);

	return outer6(A_v.memptr(), B_v.memptr());
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

mat::fixed<6,6> auto_sym_dyadic_operator(const mat::fixed<3,3> &A) {
    mat::fixed<6,6> C(fill::zeros);
    for (int J=0; J<6; J++) {
        const int k = voigt_i[J], l = voigt_j[J];
        for (int I=0; I<6; I++) {
            const int i = voigt_i[I], j = voigt_j[I];
            C.at(I,J) = 0.5*(A.at(i,k)*A.at(j,l) + A.at(i,l)*A.at(j,k));
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
