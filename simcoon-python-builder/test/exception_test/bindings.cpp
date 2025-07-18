#include <pybind11/pybind11.h>
// include numpy header for usage of array_t
#include <pybind11/numpy.h>

#include <carma>
#include <armadillo>
#include <simcoon/exception.hpp>
#include "Teig_sym.hpp"

#include <string>

namespace py = pybind11;
using namespace testexception;

using namespace pybind11::literals;

PYBIND11_MODULE(test_exception, m) {

    //Register the test_eig_sym_val_vec
    m.def("test_eig_sym_val_modify", &test_eig_sym_val_modify, "matrix_fail_eigen"_a, "copy"_a=true, "This function is used to test the eig_sym exception (modify eigval)");
    m.def("test_eig_sym_val_affect", &test_eig_sym_val_affect, "matrix_fail_eigen"_a, "copy"_a=true, "This function is used to test the eig_sym exception (affectation)");
    m.def("test_eig_sym_val_vec", &test_eig_sym_val_vec, "matrix_fail_eigen"_a, "copy"_a=true, "This function is used to test the eig_sym exception (modify eigval and eigvec)");        


}
