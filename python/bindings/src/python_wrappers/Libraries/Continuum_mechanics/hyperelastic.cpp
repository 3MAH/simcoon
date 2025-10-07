
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <carma>
#include <armadillo>


#include <simcoon/Continuum_mechanics/Functions/hyperelastic.hpp>
#include <simcoon/Continuum_mechanics/Functions/transfer.hpp>
#include <simcoon/python_wrappers/Libraries/Continuum_mechanics/hyperelastic.hpp>

using namespace std;
using namespace arma;
namespace py=pybind11;

namespace simpy {

//This function returns the isochoric invariants
py::array_t<double> isochoric_invariants(const py::array_t<double> &input, const double &J, const bool &copy) {

    if(input.ndim() == 1) {
        if(input.size() == 3) {
            vec lambdas = carma::arr_to_col(input);
            vec t = simcoon::isochoric_invariants(lambdas, J);
            return carma::col_to_arr(t, copy);            
        }
        else if(input.size() == 6) {
            vec vec_input = carma::arr_to_col(input);
            mat b = simcoon::v2t_strain(vec_input);
            vec t = simcoon::isochoric_invariants(b, J);
            return carma::col_to_arr(t, copy);            
        }
        else 
            throw std::invalid_argument("Invalid size of the one-dimensional array. Expected 3 (from principal stretches) or 6 (from left Cauchy-Green tensor in Voigt notation)");
        }                
    if(input.ndim() == 2) {        
        if((input.shape(0) == 3)&&(input.shape(1) == 3)) {
            mat b = carma::arr_to_mat(input);
            vec t = simcoon::isochoric_invariants(b, J);
            return carma::col_to_arr(t, copy);                         
        }        
        else if((input.shape(0) == 6)&&(input.shape(1) == 1)) {
            mat mat_input = carma::arr_to_mat(input);            
            mat sim_input = simcoon::v2t_strain(mat_input.as_col());
            vec t = simcoon::isochoric_invariants(sim_input, J);
            return carma::col_to_arr(t, copy);            
        }
        else {
            throw std::invalid_argument("Invalid size of the two-dimensional array. Expected a 3x3 array or a 6x1 array considering Voigt notation");
        }         
    }
}

py::array_t<double> isochoric_pstretch(const py::array_t<double> &input, const std::string &input_tensor, const double &J, const bool &copy) {

    if(input.ndim() == 1) {
        if(input.size() == 6) {
            vec vec_input = carma::arr_to_col(input);
            if(input_tensor == "b") {
                mat b = simcoon::v2t_strain(vec_input);
                vec t = simcoon::isochoric_pstretch_from_b(b);                
                return carma::col_to_arr(t, copy);                            
            }
            else if(input_tensor == "V" || input_tensor == "v") {
                mat V = simcoon::v2t_strain(vec_input);
                vec t = simcoon::isochoric_pstretch_from_V(V);
                return carma::col_to_arr(t, copy);                            
            }
            else {
                throw std::invalid_argument("Invalid input string to describe the input vector: it should be *b* for left Cauchy-Green tensor or *v* or *V* for Eulerian stretch tensor");                
            }
        }
        else {
            throw std::invalid_argument("Invalid size of the one-dimensional array. Expected 3 (from principal stretches) or 6 (from left Cauchy-Green tensor in Voigt notation)");
        }
    }                
    if(input.ndim() == 2) {        
        mat mat_input = carma::arr_to_mat(input);
        if((input.shape(0) == 3)&&(input.shape(1) == 3)) {
            if(input_tensor == "b") {
                vec t = simcoon::isochoric_pstretch_from_b(mat_input);                
                return carma::col_to_arr(t, copy);                            
            }
            else if(input_tensor == "V" || input_tensor == "v") {
                vec t = simcoon::isochoric_pstretch_from_V(mat_input);
                return carma::col_to_arr(t, copy);                            
            }            
        }        
        else if((input.shape(0) == 6)&&(input.shape(1) == 1)) {
            mat sim_input = simcoon::v2t_strain(mat_input.as_col());
            if(input_tensor == "b") {
                vec t = simcoon::isochoric_pstretch_from_b(sim_input);                
                return carma::col_to_arr(t, copy);                            
            }
            else if(input_tensor == "V" || input_tensor == "v") {
                vec t = simcoon::isochoric_pstretch_from_V(sim_input);
                return carma::col_to_arr(t, copy);                            
            }            
        }
        else {
            throw std::invalid_argument("Invalid size of the two-dimensional array. Expected a 3x3 array or a 6x1 array considering Voigt notation");
        }         
    }
}


} //namepsace simpy