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

///@file usermat_layout.cpp
///@brief A test USERMAT (Ansys interface) that exposes what it receives: a non-symmetric tangent
///       dsdde(i,j) = 100 i + j, rotateM and nStatev copied into ustatev, fixed energies
///@version 1.0

extern "C" void usermat_(int *, int *, int *, int *, int *, int *, int *, int *keycut,
                         int *, int *, int *ncomp, int *nStatev, int *,
                         double *, double *, double *, double *,
                         double *stress, double *ustatev, double *dsdde,
                         double *sedEl, double *sedPl, double *,
                         double *, double *, double *,
                         double *, double *, double *rotateM,
                         double *, double *, double *, double *,
                         int *, int *, int *, int *, int *)
{
    const int n = *ncomp;
    for (int j = 0; j < n; j++)
        for (int i = 0; i < n; i++) dsdde[i + n*j] = 100.*(i + 1) + (j + 1);   // Fortran dsdde(i,j)
    for (int i = 0; i < n; i++) stress[i] = 10.*(i + 1);
    for (int k = 0; k < 9 && k < *nStatev; k++) ustatev[k] = rotateM[k];
    if (*nStatev > 9) ustatev[9] = *nStatev;
    *sedEl = 1.5;
    *sedPl = 2.5;
    *keycut = 0;
}
