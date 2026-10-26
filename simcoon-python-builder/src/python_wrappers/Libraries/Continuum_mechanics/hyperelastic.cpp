
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <simcoon/python_wrappers/arma_to_numpy.hpp>
#include <simcoon/python_wrappers/numpy_to_arma.hpp>
#include <simcoon/python_wrappers/parallel_nogil.hpp>
#include <armadillo>
#include <simcoon/exception.hpp>
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
            vec lambdas = simpy::numpy_to_arma::arr_to_col(input);
            vec t = simcoon::isochoric_invariants(lambdas, J);
            return simpy::arma_to_numpy::col_to_arr(t, copy);            
        }
        else if(input.size() == 6) {
            vec vec_input = simpy::numpy_to_arma::arr_to_col(input);
            mat b = simcoon::v2t_strain(vec_input);
            vec t = simcoon::isochoric_invariants(b, J);
            return simpy::arma_to_numpy::col_to_arr(t, copy);            
        }
        else 
            throw std::invalid_argument("Invalid size of the one-dimensional array. Expected 3 (from principal stretches) or 6 (from left Cauchy-Green tensor in Voigt notation)");
        }                
    if(input.ndim() == 2) {        
        if((input.shape(0) == 3)&&(input.shape(1) == 3)) {
            mat b = simpy::numpy_to_arma::arr_to_mat(input);
            vec t = simcoon::isochoric_invariants(b, J);
            return simpy::arma_to_numpy::col_to_arr(t, copy);                         
        }        
        else if((input.shape(0) == 6)&&(input.shape(1) == 1)) {
            mat mat_input = simpy::numpy_to_arma::arr_to_mat(input);
            mat sim_input = simcoon::v2t_strain(mat_input.as_col());
            vec t = simcoon::isochoric_invariants(sim_input, J);
            return simpy::arma_to_numpy::col_to_arr(t, copy);            
        }
        else {
            throw std::invalid_argument("Invalid size of the two-dimensional array. Expected a 3x3 array or a 6x1 array considering Voigt notation");
        }         
    }
    throw std::invalid_argument("Invalid array dimensions");
}

py::array_t<double> isochoric_pstretch(const py::array_t<double> &input, const std::string &input_tensor, const double &J, const bool &copy) {

    // Batch: (3,3,N) tensors, or (6,N) Voigt columns with N > 1 -> (3,N) stretches, per point in
    // parallel with the GIL released. J is shared by every point (0: from each tensor).
    const bool tensors = input.ndim() == 3;
    const bool voigt_cols = input.ndim() == 2 && input.shape(0) == 6 && input.shape(1) > 1;
    if (tensors || voigt_cols) {
        const bool is_b = (input_tensor == "b");
        if (!is_b && input_tensor != "V" && input_tensor != "v") {
            throw std::invalid_argument("Invalid input string to describe the input tensor: it should be *b* for left Cauchy-Green tensor or *v* or *V* for Eulerian stretch tensor");
        }
        cube X_t;
        mat X_v;
        if (tensors) {
            X_t = simpy::numpy_to_arma::arr_to_cube_view(input);
            if (X_t.n_rows != 3 || X_t.n_cols != 3) {
                throw std::invalid_argument("Invalid shape of the three-dimensional array. Expected (3,3,N)");
            }
        }
        else {
            X_v = simpy::numpy_to_arma::arr_to_mat_view(input);
        }
        const int nb_points = tensors ? X_t.n_slices : X_v.n_cols;
        mat out(3, nb_points);   // every column written below
        parallel_for_nogil(nb_points, [&](int pt) {
            const mat X = tensors ? mat(const_cast<double *>(X_t.slice_memptr(pt)), 3, 3, false, true)
                                  : simcoon::v2t_strain(X_v.col(pt));
            out.col(pt) = is_b ? simcoon::isochoric_pstretch_from_b(X, J) : simcoon::isochoric_pstretch_from_V(X, J);
        });
        return simpy::arma_to_numpy::mat_to_arr(out, false);
    }

    if(input.ndim() == 1) {
        if(input.size() == 6) {
            vec vec_input = simpy::numpy_to_arma::arr_to_col(input);
            if(input_tensor == "b") {
                mat b = simcoon::v2t_strain(vec_input);
                vec t = simcoon::isochoric_pstretch_from_b(b);                
                return simpy::arma_to_numpy::col_to_arr(t, copy);                            
            }
            else if(input_tensor == "V" || input_tensor == "v") {
                mat V = simcoon::v2t_strain(vec_input);
                vec t = simcoon::isochoric_pstretch_from_V(V);
                return simpy::arma_to_numpy::col_to_arr(t, copy);                            
            }
            else {
                throw std::invalid_argument("Invalid input string to describe the input vector: it should be *b* for left Cauchy-Green tensor or *v* or *V* for Eulerian stretch tensor");                
            }
        }
        else {
            throw std::invalid_argument("Invalid size of the one-dimensional array. Expected 6 (from left Cauchy-Green tensor or Eulerian stretch tensor in Voigt notation)");
        }
    }                
    if(input.ndim() == 2) {        
        mat mat_input = simpy::numpy_to_arma::arr_to_mat(input);
        if((input.shape(0) == 3)&&(input.shape(1) == 3)) {
            if(input_tensor == "b") {
                vec t = simcoon::isochoric_pstretch_from_b(mat_input);                
                return simpy::arma_to_numpy::col_to_arr(t, copy);                            
            }
            else if(input_tensor == "V" || input_tensor == "v") {
                vec t = simcoon::isochoric_pstretch_from_V(mat_input);
                return simpy::arma_to_numpy::col_to_arr(t, copy);                            
            }            
        }        
        else if((input.shape(0) == 6)&&(input.shape(1) == 1)) {
            mat sim_input = simcoon::v2t_strain(mat_input.as_col());
            if(input_tensor == "b") {
                vec t = simcoon::isochoric_pstretch_from_b(sim_input);                
                return simpy::arma_to_numpy::col_to_arr(t, copy);                            
            }
            else if(input_tensor == "V" || input_tensor == "v") {
                vec t = simcoon::isochoric_pstretch_from_V(sim_input);
                return simpy::arma_to_numpy::col_to_arr(t, copy);                            
            }            
        }
        else {
            throw std::invalid_argument("Invalid size of the two-dimensional array. Expected a 3x3 array or a 6x1 array considering Voigt notation");
        }         
    }
    throw std::invalid_argument("Invalid array dimensions");
}

py::array_t<double> tau_iso_hyper_invariants(const double &dWdI_1_bar, const double &dWdI_2_bar, const py::array_t<double> &input, const double &J, const bool &copy) {

    if(input.ndim() == 1) {
        if(input.size() == 6) {
            vec vec_input = simpy::numpy_to_arma::arr_to_col(input);
            mat b_arma = simcoon::v2t_strain(vec_input);
            mat m = simcoon::tau_iso_hyper_invariants(dWdI_1_bar, dWdI_2_bar, b_arma, J);
            return simpy::arma_to_numpy::mat_to_arr(m, copy);
        }
        else {
            throw std::invalid_argument("Invalid size of the one-dimensional array. Expected 6 from left Cauchy-Green tensor in Voigt notation");
        }
    }                
    if(input.ndim() == 2) {        
        mat mat_input = simpy::numpy_to_arma::arr_to_mat(input);
        if((input.shape(0) == 3)&&(input.shape(1) == 3)) {
            vec t = simcoon::tau_iso_hyper_invariants(dWdI_1_bar, dWdI_2_bar, mat_input, J);
            return simpy::arma_to_numpy::mat_to_arr(t, copy);
        }
        else if((input.shape(0) == 6)&&(input.shape(1) == 1)) {
            mat sim_input = simcoon::v2t_strain(mat_input.as_col());
            mat m = simcoon::tau_iso_hyper_invariants(dWdI_1_bar, dWdI_2_bar, sim_input, J);
            return simpy::arma_to_numpy::mat_to_arr(m, copy);
        }
        else {
            throw std::invalid_argument("Invalid size of the two-dimensional array. Expected a 3x3 array or a 6x1 array considering Voigt notation");
        }         
    }
    throw std::invalid_argument("Invalid array dimensions");
}

py::array_t<double> sigma_iso_hyper_invariants(const double &dWdI_1_bar, const double &dWdI_2_bar, const py::array_t<double> &input, const double &J, const bool &copy) {

    if(input.ndim() == 1) {
        if(input.size() == 6) {
            vec vec_input = simpy::numpy_to_arma::arr_to_col(input);
            mat b_arma = simcoon::v2t_strain(vec_input);
            mat m = simcoon::sigma_iso_hyper_invariants(dWdI_1_bar, dWdI_2_bar, b_arma, J);
            return simpy::arma_to_numpy::mat_to_arr(m, copy);
        }
        else {
            throw std::invalid_argument("Invalid size of the one-dimensional array. Expected 6 from left Cauchy-Green tensor in Voigt notation");
        }
    }                
    if(input.ndim() == 2) {        
        mat mat_input = simpy::numpy_to_arma::arr_to_mat(input);
        if((input.shape(0) == 3)&&(input.shape(1) == 3)) {
            mat m = simcoon::sigma_iso_hyper_invariants(dWdI_1_bar, dWdI_2_bar, mat_input, J);                
            return simpy::arma_to_numpy::mat_to_arr(m, copy);                            
        }        
        else if((input.shape(0) == 6)&&(input.shape(1) == 1)) {
            mat sim_input = simcoon::v2t_strain(mat_input.as_col());
            mat m = simcoon::sigma_iso_hyper_invariants(dWdI_1_bar, dWdI_2_bar, sim_input, J);
            return simpy::arma_to_numpy::mat_to_arr(m, copy);
        }
        else {
            throw std::invalid_argument("Invalid size of the two-dimensional array. Expected a 3x3 array or a 6x1 array considering Voigt notation");
        }         
    }
    throw std::invalid_argument("Invalid array dimensions");
}

py::array_t<double> tau_vol_hyper(const double &dUdJ, const py::array_t<double> &input, const double &J, const bool &copy) {

    if(input.ndim() == 1) {
        if(input.size() == 6) {
            vec vec_input = simpy::numpy_to_arma::arr_to_col(input);
            mat b_arma = simcoon::v2t_strain(vec_input);
            mat m = simcoon::tau_vol_hyper(dUdJ, b_arma, J);
            return simpy::arma_to_numpy::mat_to_arr(m, copy);
        }
        else {
            throw std::invalid_argument("Invalid size of the one-dimensional array. Expected 6 from left Cauchy-Green tensor in Voigt notation");
        }
    }                
    if(input.ndim() == 2) {        
        mat mat_input = simpy::numpy_to_arma::arr_to_mat(input);
        if((input.shape(0) == 3)&&(input.shape(1) == 3)) {
            mat m = simcoon::tau_vol_hyper(dUdJ, mat_input, J);
            return simpy::arma_to_numpy::mat_to_arr(m, copy);
        }
        else if((input.shape(0) == 6)&&(input.shape(1) == 1)) {
            mat sim_input = simcoon::v2t_strain(mat_input.as_col());
            mat m = simcoon::tau_vol_hyper(dUdJ, sim_input, J);
            return simpy::arma_to_numpy::mat_to_arr(m, copy);
        }
        else {
            throw std::invalid_argument("Invalid size of the two-dimensional array. Expected a 3x3 array or a 6x1 array considering Voigt notation");
        }         
    }
    throw std::invalid_argument("Invalid array dimensions");
}

py::array_t<double> sigma_vol_hyper(const double &dUdJ, const py::array_t<double> &input, const double &J, const bool &copy) {

    if(input.ndim() == 1) {
        if(input.size() == 6) {
            vec vec_input = simpy::numpy_to_arma::arr_to_col(input);
            mat b_arma = simcoon::v2t_strain(vec_input);
            mat m = simcoon::sigma_vol_hyper(dUdJ, b_arma, J);
            return simpy::arma_to_numpy::mat_to_arr(m, copy);
        }
        else {
            throw std::invalid_argument("Invalid size of the one-dimensional array. Expected 6 from left Cauchy-Green tensor in Voigt notation");
        }
    }                
    if(input.ndim() == 2) {        
        mat mat_input = simpy::numpy_to_arma::arr_to_mat(input);
        if((input.shape(0) == 3)&&(input.shape(1) == 3)) {
            mat m = simcoon::sigma_vol_hyper(dUdJ, mat_input, J);
            return simpy::arma_to_numpy::mat_to_arr(m, copy);
        }
        else if((input.shape(0) == 6)&&(input.shape(1) == 1)) {
            mat sim_input = simcoon::v2t_strain(mat_input.as_col());
            mat m = simcoon::sigma_vol_hyper(dUdJ, sim_input, J);
            return simpy::arma_to_numpy::mat_to_arr(m, copy);
        }
        else {
            throw std::invalid_argument("Invalid size of the two-dimensional array. Expected a 3x3 array or a 6x1 array considering Voigt notation");
        }         
    }
    throw std::invalid_argument("Invalid array dimensions");
}

} //namespace simpy