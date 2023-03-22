#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace simpy{

//Check the material symetries and the type of elastic response for a given stiffness tensor
pybind11::dict check_symetries(const pybind11::array_t<double> &);

//Return a list of elastic properties for the isotropic case (E,nu) from a stiffness tensor
pybind11::array_t<double> L_iso_props(const pybind11::array_t<double> &);

//Return a list of elastic properties for the isotropic case (E,nu) from a compliance tensor
pybind11::array_t<double> M_iso_props(const pybind11::array_t<double> &);
    
//Return a list of elastic properties for the transversely isotropic case (EL,ET,nuTL,nuTT,GLT) from a stiffness tensor
pybind11::array_t<double> L_isotrans_props(const pybind11::array_t<double> &, const int &);

//Return a list of elastic properties for the transversely isotropic case (EL,ET,nuTL,nuTT,GLT) from a compliance tensor
pybind11::array_t<double> M_isotrans_props(const pybind11::array_t<double> &, const int &);

//Return a list of elastic properties for the cubic case (E,nu,G) from a stiffness tensor
pybind11::array_t<double> L_cubic_props(const pybind11::array_t<double> &);

//Return a list of elastic properties for the cubic case (E,nu,G) from a compliance tensor
pybind11::array_t<double> M_cubic_props(const pybind11::array_t<double> &);

//Return a list of elastic properties for the orthtropic case (E1,E2,E3,nu12,nu13,nu23,G12,G13,G23) from a stiffness tensor
pybind11::array_t<double> L_ortho_props(const pybind11::array_t<double> &);

//Return a list of elastic properties for the orthtropic case (E1,E2,E3,nu12,nu13,nu23,G12,G13,G23) from a compliance tensor
pybind11::array_t<double> M_ortho_props(const pybind11::array_t<double> &);

//Return a list of elastic properties for the anisotropic case (E1,E2,E3,nu12,nu13,nu23,G12,G13,G23,deviations) from a compliance tensor
pybind11::array_t<double> M_aniso_props(const pybind11::array_t<double> &);
    
} //namespace simpy