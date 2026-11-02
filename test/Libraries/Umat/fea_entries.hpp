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

///@file fea_entries.hpp
///@brief The FE entry points of software/ (umat_singleM/T: umat_, usermat_singleM: usermat_), as
///       the FE codes call them; each test executable links the entries it declares
///@version 1.0

#pragma once

extern "C" void umat_(double *stress, double *statev, double *ddsdde, double &sse, double &spd, double &scd, double &rpl, double *ddsddt, double *drplde, double &drpldt, const double *stran, const double *dstran, const double *time, const double &dtime, const double &temperature, const double &Dtemperature, const double &predef, const double &dpred, char *cmname, const int &ndi, const int &nshr, const int &ntens, const int &nstatev, const double *props, const int &nprops, const double &coords, const double *drot, double &pnewdt, const double &celent, const double *dfgrd0, const double *dfgrd1, const int &noel, const int &npt, const double &layer, const int &kspt, const int &kstep, const int &kinc);

extern "C" void usermat_(int *matId, int *elemId, int *kDomIntPt, int *kLayer, int *kSectPt, int *ldStep, int *iSubst, int *keycut, int *nDirect, int *nShear, int *ncomp, int *nStatev, int *nProp, int *nPt, double *Time, double *dTime, double *Temp, double *dTemp, double *stress, double *ustatev, double *dsdePl, double *sedEl, double *sedPl, double *epseq, double *strain, double *dstrain, double *epsPl, double *prop, double *coords, double *rotateM, double *defGrad_t, double *defGrad, double *tsstif);
