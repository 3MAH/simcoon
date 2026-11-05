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

///@file eshelby.cpp
///@brief Definition of the Eshelby tensor for ellipsoidal inclusions with
// Parts of this methods are copyrighted by Gavazzi & Lagoudas 1992 - Fair use only
///@version 1.0

#include <math.h>
#include <map>
#include <mutex>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Homogenization/eshelby.hpp>

using namespace std;
using namespace arma;

namespace simcoon{
    
//Eshelby tensor for a sphere
mat Eshelby_sphere(const double &nu) {

	mat S = zeros(6,6);
	
	S(0,0) = (7.-5.*nu)/(15.*(1.-nu));
	S(0,1) = ((5.*nu-1.)/(15.*(1-nu)));
	S(0,2) = S(0,1);
	S(1,0) = S(0,1);
	S(1,1) = S(0,0);
	S(1,2) = S(0,1);
	S(2,0) = S(0,1);
	S(2,1) = S(0,1);
	S(2,2) = S(0,0);
	S(3,3) = 2.* ((4.-5.*nu)/(15.*(1.-nu)));
	S(4,4) = S(3,3);
	S(5,5) = S(3,3);
    
    return S;
}

//	Eshelby tensor determination. The cylinder is oriented in such a way that the axis direction is the 1 direction. a2=a3 here
mat Eshelby_cylinder(const double &nu) {
    
    mat S = zeros(6,6);
    
	 S(0,0) = 0.;
	 S(0,1) = 0.;
	 S(0,2) = 0.;
	 S(1,0) = nu/(2.*(1.-nu));
	 S(1,1) = (5.-4.*nu)/(8.*(1.-nu));
	 S(1,2) = (4.*nu-1.)/(8.*(1.-nu));
	 S(2,0) = nu/(2.*(1.-nu));
	 S(2,1) = (4.*nu-1.)/(8.*(1.-nu));
	 S(2,2) = (5.-4.*nu)/(8.*(1.-nu));
	 S(3,3) = 2.*(1./4.);
	 S(4,4) = 2.*(1./4.);
	 S(5,5) = 2.*(3.-4.*nu)/(8.*(1.-nu));
	
	return S;
}

//	Eshelby tensor determination. The prolate shape is oriented in such a way that the axis direction is the 1 direction. a1>a2=a3 here
mat Eshelby_prolate(const double &nu, const double &ar) {
    
    mat S = zeros(6,6);
    double g = ar/pow(ar*ar-1,1.5)*(ar*pow((ar*ar-1),0.5)-acosh(ar));
    
	 S(0,0) = 1./(2*(1-nu))*(1-2*nu+(3*ar*ar-1)/(ar*ar-1)-(1-2*nu+3*ar*ar/(ar*ar-1))*g);
	 S(0,1) = -1./(2*(1-nu))*(1-2*nu+1/(ar*ar-1))+g/(2*(1-nu))*(1-2*nu+3/(2*(ar*ar-1)));
	 S(0,2) = S(0,1);
	 S(1,0) = -ar*ar/(2*(1-nu)*(ar*ar-1))+g/(4*(1-nu))*(3*ar*ar/(ar*ar-1)-(1-2*nu));
	 S(1,1) = 3*ar*ar/(8*(1-nu)*(ar*ar-1))+1./(4*(1-nu))*(1-2*nu-9/(4*(ar*ar-1)))*g;
	 S(1,2) = 1./(4.*(1-nu))*(ar*ar/(2*(ar*ar-1))-(1-2*nu+3/(4*(ar*ar-1)))*g);
	 S(2,0) = S(1,0);
	 S(2,1) = S(1,2);
	 S(2,2) = S(1,1);
	 S(3,3) = 2./(4.*(1-nu))*(1-2*nu-(ar*ar+1)/(ar*ar-1)-g/2.*(1-2*nu-3*(ar*ar+1)/(ar*ar-1)));
	 S(4,4) = S(3,3);
	 S(5,5) = 2./(4.*(1-nu))*(ar*ar/(2*(ar*ar-1))+g*(1-2*nu-3/(4*(ar*ar-1))));
	
	return S;
}

//	Eshelby tensor determination. The oblate shape is oriented in such a way that the axis direction is the 1 direction. a1<a2=a3 here
mat Eshelby_oblate(const double &nu, const double &ar) {
    
	mat S = zeros(6,6);
    double g = ar/pow(1-ar*ar,1.5)*(-ar*pow((1-ar*ar),0.5)+acos(ar));
    
	 S(0,0) = 1./(2*(1-nu))*(1-2*nu+(3*ar*ar-1)/(ar*ar-1)-(1-2*nu+3*ar*ar/(ar*ar-1))*g);
	 S(0,1) = -1./(2*(1-nu))*(1-2*nu+1/(ar*ar-1))+g/(2*(1-nu))*(1-2*nu+3/(2*(ar*ar-1)));
	 S(0,2) = S(0,1);
	 S(1,0) = -ar*ar/(2*(1-nu)*(ar*ar-1))+g/(4*(1-nu))*(3*ar*ar/(ar*ar-1)-(1-2*nu));
	 S(1,1) = 3*ar*ar/(8*(1-nu)*(ar*ar-1))+1./(4*(1-nu))*(1-2*nu-9/(4*(ar*ar-1)))*g;
	 S(1,2) = 1./(4.*(1-nu))*(ar*ar/(2*(ar*ar-1))-(1-2*nu+3/(4*(ar*ar-1)))*g);
	 S(2,0) = S(1,0);
	 S(2,1) = S(1,2);
	 S(2,2) = S(1,1);
	 S(3,3) = 2./(4.*(1-nu))*(1-2*nu-(ar*ar+1)/(ar*ar-1)-g/2.*(1-2*nu-3*(ar*ar+1)/(ar*ar-1)));
	 S(4,4) = S(3,3);
	 S(5,5) = 2./(4.*(1-nu))*(ar*ar/(2*(ar*ar-1))+g*(1-2*nu-3/(4*(ar*ar-1))));
	
	return S;
}

//	Eshelby tensor determination for a penny-shaped crack (limit of oblate as ar->0). The crack normal is the 1 direction.
mat Eshelby_penny(const double &nu) {
    
	mat S = zeros(6,6);
    
	// Limit values as ar -> 0 (penny-shaped crack / flat disc)
	// S_1111 -> 1
	// S_1122 = S_1133 -> nu/(1-nu)
	// S_2211 = S_3311 -> 0
	// S_2222 = S_3333 -> 0
	// S_2233 = S_3322 -> 0
	// S_2323 -> 0
	// S_1212 = S_1313 -> 1/2
	 S(0,0) = 1.;
	 S(0,1) = nu/(1.-nu);
	 S(0,2) = nu/(1.-nu);
	 S(1,0) = 0.;
	 S(1,1) = 0.;
	 S(1,2) = 0.;
	 S(2,0) = 0.;
	 S(2,1) = 0.;
	 S(2,2) = 0.;
	 S(3,3) = 0.;
	 S(4,4) = 2.*(1./2.);  // S_1313 in Voigt: factor 2 for shear
	 S(5,5) = 2.*(1./2.);  // S_1212 in Voigt: factor 2 for shear
	
	return S;
}

//This methods is using the Voigt notations for the tensors.
void calG(const double &pt, const double &a1, const double &a2, const double &a3, const double &x3, const Mat<int> &Id, const mat &Lt, mat &G)
{
    
	vec X(3);
	double x1 = sqrt(1.-pow(x3,2.));
	X(0) = x1*cos(pt)*(1./a1);
	X(1) = x1*sin(pt)*(1./a2);
	X(2) = x3*(1./a3);
	int m,n,ij,kl;
	
	mat rk = zeros(3,3);
	mat rn = zeros(3,3);
	
	for (int i=0; i<3; i++) {
		for (int j=0; j<3; j++) {
			m = Id(i,j);
			for (int k=i; k<3; k++) {
				for (int l=0; l<3; l++) {
					n = Id(k,l);
					rk(i,k) += Lt(m,n)*X(j)*X(l);                 		//to check!
				}
			}
		}	
	}
	
	rn(0,0) = rk(1,1)*rk(2,2) - rk(1,2)*rk(1,2);
	rn(0,1) = -1.*(rk(0,1)*rk(2,2) - rk(1,2)*rk(0,2));
	rn(0,2) = rk(0,1)*rk(1,2) - rk(1,1)*rk(0,2);
	rn(1,1) = rk(0,0)*rk(2,2) - rk(0,2)*rk(0,2);
	rn(1,2) = -1.*(rk(0,0)*rk(1,2) - rk(0,1)*rk(0,2));
	rn(2,2) = rk(0,0)*rk(1,1) - rk(0,1)*rk(0,1);                        //to check!
	
	double D = rn(0,0)*rk(0,0) + rn(0,1)*rk(0,1) + rn(0,2)*rk(0,2);     //to check!
    
	for (int i=0; i<3; i++) {
		for (int j=i; j<3; j++) {
			rn(i,j) = rn(i,j)*(1./D);
		}	
	} 
		 
	for (int i=0; i<3; i++) {
		for (int j=i; j<3; j++) {
			ij = Id(i,j);
			for (int k=0; k<3; k++) {
				for (int l=k; l<3; l++) {
					kl = Id(k,l);
					G(ij,kl) = rn(i,j)*X(l)*X(k);
				}
			}
		}	
	}

}

//Voigt index of the pair (i,j): the table that Eshelby and T_II hand to Gauss as Id
static const int Voigt_index[3][3] = {{0,3,4},{3,1,5},{4,5,2}};

// sum of calG in the per-entry operation order of the point loop (bitwise); sin/cos hoisted
void Gauss(Mat<int> &Id, const mat &Lt, mat &G, const double &a1, const double &a2, const double &a3, const vec &x, const vec &wx, const vec &y, const vec &wy, const int &mp, const int &np)
{
    (void)Id;   //always the Voigt table (Voigt_index)
    double L[6][6];
    for (int i=0; i<6; i++)
        for (int j=0; j<6; j++)
            L[i][j] = Lt(i,j);
    double Gl[6][6];
    for (int i=0; i<6; i++)
        for (int j=0; j<6; j++)
            Gl[i][j] = G(i,j);

    std::vector<double> cos_y(np), sin_y(np);
    for (int i=0; i<np; i++) {
        cos_y[i] = cos(y(i));
        sin_y[i] = sin(y(i));
    }
    const double inv_a1 = 1./a1;
    const double inv_a2 = 1./a2;
    const double inv_a3 = 1./a3;

    //per point of the ring: X and the upper triangle of inv(K), K_ik = L_ijkl X_j X_l
    std::vector<double> ring(9*np);
    double *X[3] = {&ring[0], &ring[np], &ring[2*np]};
    double *R00 = &ring[3*np], *R01 = &ring[4*np], *R02 = &ring[5*np];
    double *R11 = &ring[6*np], *R12 = &ring[7*np], *R22 = &ring[8*np];
    const double *R[3][3] = {{R00, R01, R02}, {nullptr, R11, R12}, {nullptr, nullptr, R22}};

    for (int l=0; l<mp; l++) {
        const double x3 = x(l);
        const double x1 = sqrt(1.-pow(x3,2.));
        const double X2 = x3*inv_a3;
        for (int p=0; p<np; p++) {
            const double Xp[3] = {x1*cos_y[p]*inv_a1, x1*sin_y[p]*inv_a2, X2};
            double rk[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
            for (int i=0; i<3; i++) {
                for (int j=0; j<3; j++) {
                    const int m = Voigt_index[i][j];
                    for (int k=i; k<3; k++) {
                        for (int ll=0; ll<3; ll++) {
                            rk[i][k] += L[m][Voigt_index[k][ll]]*Xp[j]*Xp[ll];
                        }
                    }
                }
            }
            const double rn00 = rk[1][1]*rk[2][2] - rk[1][2]*rk[1][2];
            const double rn01 = -1.*(rk[0][1]*rk[2][2] - rk[1][2]*rk[0][2]);
            const double rn02 = rk[0][1]*rk[1][2] - rk[1][1]*rk[0][2];
            const double rn11 = rk[0][0]*rk[2][2] - rk[0][2]*rk[0][2];
            const double rn12 = -1.*(rk[0][0]*rk[1][2] - rk[0][1]*rk[0][2]);
            const double rn22 = rk[0][0]*rk[1][1] - rk[0][1]*rk[0][1];
            const double D = rn00*rk[0][0] + rn01*rk[0][1] + rn02*rk[0][2];
            const double inv_D = 1./D;
            X[0][p] = Xp[0];
            X[1][p] = Xp[1];
            X[2][p] = Xp[2];
            R00[p] = rn00*inv_D;
            R01[p] = rn01*inv_D;
            R02[p] = rn02*inv_D;
            R11[p] = rn11*inv_D;
            R12[p] = rn12*inv_D;
            R22[p] = rn22*inv_D;
        }

        const double wx_l = wx(l);
        for (int p=0; p<np; p++) {
            const double wy_p = wy(p);
            for (int i=0; i<3; i++) {
                for (int j=i; j<3; j++) {
                    const int ij = Voigt_index[i][j];
                    const double rn_ij = R[i][j][p];
                    for (int k=0; k<3; k++) {
                        for (int ll=k; ll<3; ll++) {
                            const int kl = Voigt_index[k][ll];
                            const double H = rn_ij*X[ll][p]*X[k][p];      //calG
                            const double rg = H*wy_p+0.;                  //azimuthal weight
                            Gl[ij][kl] = Gl[ij][kl] + wx_l*rg;            //polar weight
                        }
                    }
                }
            }
        }
    }

    for (int i=0; i<6; i++)
        for (int j=0; j<6; j++)
            G(i,j) = Gl[i][j];
}
    
mat Eshelby(const mat &Lt, const double &a1, const double &a2, const double &a3, const vec &x, const vec &wx, const vec &y, const vec &wy, const int &mp, const int &np)
{
    mat G = zeros(6,6);
	Mat<int> Id(6,6);
    Id.zeros();

    mat S = zeros(6,6);
	int ij=0;
	int mn=0;
	int pq=0;
	int ip=0;
	int jq=0;
	int jp=0;
	int iq=0;
		
    Id(0,0) = 0;
    Id(0,1) = 3;
    Id(0,2) = 4;
    Id(1,0) = 3;
    Id(1,1) = 1;
    Id(1,2) = 5;
    Id(2,0) = 4;
    Id(2,1) = 5;
    Id(2,2) = 2;
    
	Gauss(Id, Lt, G, a1, a2, a3, x, wx, y, wy, mp, np);
	    
	for (int i=0; i<3; i++) {
		for (int j=i; j<3; j++) {
			ij = Id(i,j);
			for (int m=0; m<3; m++) {
				for (int n=m; n<3; n++) {
					mn = Id(m,n);
					for (int p=0; p<3; p++) {
						for (int q=0; q<3; q++) {
							pq = Id(p,q);
							ip = Id(i,p);
							jq = Id(j,q);
							jp = Id(j,p);
							iq = Id(i,q);
							S(ij,mn) = S(ij,mn)+Lt(pq,mn)*(G(ip,jq)+G(jp,iq));
						}
					}
				}
			}
		}
	}
	
	for (int i=0; i<3; i++) {
		for (int j=0; j<6; j++) {
			S(i,j) = S(i,j)*(1./(8.*simcoon::pi));
		}
	}

	for (int i=3; i<6; i++) {
		for (int j=0; j<6; j++) {
			S(i,j) = S(i,j)*(1./(4.*simcoon::pi));
		}
	}   
    return S;
}

mat Eshelby(const mat &Lt, const double &a1, const double &a2, const double &a3, const int &mp, const int &np) {
    
    const eshelby_quadrature &q = eshelby_quadrature_points(mp, np);
    return Eshelby(Lt, a1, a2, a3, q.x, q.wx, q.y, q.wy, mp, np);
}
    

mat T_II(const mat &Lt, const double &a1, const double &a2, const double &a3, const vec &x, const vec &wx, const vec &y, const vec &wy, const int &mp, const int &np)
{   
	mat G = zeros(6,6);
	Mat<int> Id(6,6);
    Id.zeros();
    
    mat T_II = zeros(6,6);
    mat I = eye(6,6);
    
	int ij=0;
	int mn=0;
	int pq=0;
	int ip=0;
	int jq=0;
	int jp=0;
	int iq=0;
    
    Id(0,0) = 0;
    Id(0,1) = 3;
    Id(0,2) = 4;
    Id(1,0) = 3;
    Id(1,1) = 1;
    Id(1,2) = 5;
    Id(2,0) = 4;
    Id(2,1) = 5;
    Id(2,2) = 2;
        
	Gauss(Id, Lt, G, a1, a2, a3, x, wx, y, wy, mp, np);
		
	for (int i=0; i<3; i++) {
		for (int j=i; j<3; j++) {
			ij = Id(i,j);
			for (int m=0; m<3; m++) {
				for (int n=m; n<3; n++) {
					mn = Id(m,n);
					for (int p=0; p<3; p++) {
						for (int q=0; q<3; q++) {
							pq = Id(p,q);
							ip = Id(i,p);
							jq = Id(j,q);
							jp = Id(j,p);
							iq = Id(i,q);
							T_II(ij,mn) = T_II(ij,mn)+I(pq,mn)*(G(ip,jq)+G(jp,iq));
						}
					}
				}
			}
		}
	}
	
	for (int i=0; i<3; i++) {
		for (int j=0; j<6; j++) {
			T_II(i,j) = T_II(i,j)*(1./(8.*simcoon::pi));
		}
	}
    
	for (int i=3; i<6; i++) {
		for (int j=0; j<6; j++) {
			T_II(i,j) = T_II(i,j)*(1./(4.*simcoon::pi));
		}
	}    
    return T_II;
}

mat T_II(const mat &Lt, const double &a1, const double &a2, const double &a3, const int &mp, const int &np) {
    
    const eshelby_quadrature &q = eshelby_quadrature_points(mp, np);
    return T_II(Lt, a1, a2, a3, q.x, q.wx, q.y, q.wy, mp, np);
}
    
    
void points(vec &x, vec &wx, vec &y, vec &wy, const int &mp, const int &np)
{
    
	double x1=0.;
	double x2=2.*simcoon::pi;
	int n=0;
	int m=0;
	double p1=0.;
	double p2=0.;
	double p3=0.;
	double pp=0.;
	double xl=0.;
	double xm=0.;
	double z=0.;
	double z1=0.;

	for (int II=1; II<=2; II++) {
		if (II == 2) {
			x1=-1.;
			x2=1.;
			n=mp;
		}
		else {
			n=np;
		}

		//The roots are symmetric in the interval, so we only have to find half of them
		m = (n+1)/2;
		xm = 0.5*(x2+x1);
		xl = 0.5*(x2-x1);
		
		//Loop over the desired roots
		for (int i=1; i<=m; i++) {
			z=cos(simcoon::pi*(i-0.25)/(n+0.5));
			z1=0.;

			while(fabs(z-z1)>simcoon::iota) {
				
				p1 = 1.;
				p2 = 0.;
				for (int j=1; j<=n; j++) {
					p3=p2;
					p2=p1;
					p1=((2.*j-1.)*z*p2-(j-1.)*p3)/j;
				}
				pp=n*(z*p1-p2)/(z*z-1.);
				z1=z;
				z=z1-p1/pp;
			}
			if (II == 2) {
				x(i-1)=xm-xl*z;
				x(n-i)=xm+xl*z;
				wx(i-1)=2.*xl/((1.-z*z)*pp*pp);
				wx(n-i)=wx(i-1);
			}
			else {
				y(i-1)=xm-xl*z;
				y(n-i)=xm+xl*z;
				wy(i-1)=2.*xl/((1.-z*z)*pp*pp);
				wy(n-i)=wy(i-1);
			}

		}

	}

}

eshelby_quadrature::eshelby_quadrature(const int &m, const int &n, const int &r, const bool &cf)
    : mp(r == 1 ? (m+1)/2 : m), np(n), rule(r), closed_form(cf), mp_asked(m)
{
    if (rule == 0) {
        x.set_size(mp);
        wx.set_size(mp);
        y.set_size(np);
        wy.set_size(np);
        points(x, wx, y, wy, mp, np);
    }
    else {
        //2*mp Gauss-Legendre points on [-1, 1], symmetric: the upper half, weights doubled
        vec x_full(2*mp), wx_full(2*mp), y_unused(1), wy_unused(1);
        points(x_full, wx_full, y_unused, wy_unused, 2*mp, 1);
        x = x_full.subvec(mp, 2*mp-1);
        wx = 2.*wx_full.subvec(mp, 2*mp-1);
        //trapezoidal rule over one period of phi
        y.set_size(np);
        wy.set_size(np);
        for (int i=0; i<np; i++) {
            y(i) = 2.*simcoon::pi*i/np;
            wy(i) = 2.*simcoon::pi/np;
        }
    }
}

const eshelby_quadrature &eshelby_quadrature_points(const int &mp, const int &np, const int &rule, const bool &closed_form)
{
    if (mp < 1 || np < 1) {
        throw std::invalid_argument("Eshelby integrals: mp = " + std::to_string(mp) + " and np = " + std::to_string(np)
                                    + " integration points, at least 1 of each is needed");
    }
    if (rule != 0 && rule != 1) {
        throw std::invalid_argument("Eshelby integrals: rule = " + std::to_string(rule) + " is neither 0 (Gauss-Legendre "
                                    "on the sphere) nor 1 (half sphere, trapezoidal rule in phi)");
    }
    // Last tables of this thread: no lock while a thread keeps integrating with one setting
    thread_local const eshelby_quadrature *last = nullptr;
    if (last != nullptr && last->mp_asked == mp && last->np == np && last->rule == rule && last->closed_form == closed_form)
        return *last;

    // never destroyed: no arma destructor at DLL unload (Windows)
    typedef std::tuple<int, int, int, bool> key_type;
    static std::mutex &tables_mutex = *new std::mutex();
    static std::map<key_type, const eshelby_quadrature *> &tables = *new std::map<key_type, const eshelby_quadrature *>();
    std::lock_guard<std::mutex> lock(tables_mutex);
    const key_type key(mp, np, rule, closed_form);
    auto it = tables.find(key);
    if (it == tables.end()) {
        it = tables.emplace(key, new eshelby_quadrature(mp, np, rule, closed_form)).first;
    }
    last = it->second;
    return *last;
}

bool Eshelby_closed_form(const mat &L, const double &a1, const double &a2, const double &a3, mat &S)
{
    if (a2 != a3)
        return false;
    const double ar = a1/a2;
    if (ar != 1. && fabs(ar - 1.) < 1.E-2)
        return false;
    const double lambda = L(0,1);
    const double mu = L(3,3);
    const double tol = 1.E-10*abs(L).max();
    for (int j=0; j<6; j++) {
        for (int i=0; i<6; i++) {
            double L_iso_ij = 0.;
            if (i < 3 && j < 3)
                L_iso_ij = (i == j) ? lambda + 2.*mu : lambda;
            else if (i == j)
                L_iso_ij = mu;
            if (!(fabs(L(i,j) - L_iso_ij) <= tol))
                return false;
        }
    }
    const double nu = lambda/(2.*(lambda + mu));
    if (ar == 1.)
        S = Eshelby_sphere(nu);
    else if (ar > 1.)
        S = Eshelby_prolate(nu, ar);
    else
        S = Eshelby_oblate(nu, ar);
    return true;
}

} //namespace simcoon
