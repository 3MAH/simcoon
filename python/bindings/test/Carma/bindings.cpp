#include <pybind11/pybind11.h>
// include numpy header for usage of array_t
#include <pybind11/numpy.h>

#include <carma>
#include <armadillo>
#include "Troundtrip.hpp"
#include "Tarr_to_mat.hpp"

#include <string>

namespace py = pybind11;
using namespace testcarma;

PYBIND11_MODULE(test_carma, m) {

    //Troundtrip
    m.def("mat_roundtrip", &testcarma::test_mat_roundtrip, "Test mat_roundtrip");
    m.def("row_roundtrip", &testcarma::test_row_roundtrip, "Test row_roundtrip");    
    m.def("col_roundtrip", &testcarma::test_col_roundtrip, "Test col_roundtrip");    
    m.def("cube_roundtrip", &testcarma::test_cube_roundtrip, "Test cube_roundtrip");    

    //arr_to_mat
    m.def("arr_to_mat_double", &testcarma::test_arr_to_mat_double, "Test arr_to_mat_double");
    m.def("arr_to_mat_long", &testcarma::test_arr_to_mat_long, "Test arr_to_mat_long");
    m.def("arr_to_mat_double_copy", &testcarma::test_arr_to_mat_double_copy, "Test arr_to_mat_double_copy");
    m.def("arr_to_mat_1d", &testcarma::test_arr_to_mat_1d, "Test arr_to_mat_1d");
    m.def("arr_to_col", &testcarma::test_arr_to_col, "Test arr_to_col");
    m.def("arr_to_row", &testcarma::test_arr_to_row, "Test arr_to_row");
    m.def("arr_to_cube", &testcarma::test_arr_to_cube, "Test arr_to_cube");
    m.def("to_arma_mat", &testcarma::test_to_arma_mat, "Test to_arma");
    m.def("to_arma_cube", &testcarma::test_to_arma_cube, "Test to_arma");
    m.def("to_arma_col", &testcarma::test_to_arma_col, "Test to_arma");
    m.def("to_arma_row", &testcarma::test_to_arma_row, "Test to_arma");

}
