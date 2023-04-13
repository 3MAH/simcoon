
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <string>
#include <carma>
#include <armadillo>
#include <assert.h>

#include <simcoon/Continuum_mechanics/Functions/stress.hpp>
#include <simcoon/Continuum_mechanics/Functions/transfer.hpp>
#include <simcoon/python_wrappers/Libraries/Continuum_mechanics/stress.hpp>

using namespace std;
using namespace arma;
namespace py=pybind11;

namespace simpy {

py::array_t<double> Cauchy2PKI(const py::array_t<double> &sigma, const py::array_t<double> &F, const double &J, const bool & copy) {

    if(sigma.ndim() == 1) {
        if(sigma.size() == 6) {
            mat sigma_cpp = simcoon::v2t_stress(carma::arr_to_col(sigma));
            mat F_cpp = carma::arr_to_mat(F);            
            vec PKI = simcoon::t2v_stress(simcoon::Cauchy2PKI(sigma_cpp, F_cpp, J));
            return carma::col_to_arr(PKI, copy);
        }
        else {
            throw std::invalid_argument("Invalid size of the one-dimensional array. Expected 6");
        }         
    }
    else if (sigma.ndim() == 2) {
        if((sigma.shape(0) == 3)&&(sigma.shape(1) == 3)) {
            mat sigma_cpp = carma::arr_to_mat(sigma);
            mat F_cpp = carma::arr_to_mat(F);
            mat PKI = simcoon::Cauchy2PKI(sigma_cpp, F_cpp, J);
            return carma::mat_to_arr(PKI, copy);
        }
        else if(sigma.shape(0) == 6) {
            assert(sigma.shape(1) == F.shape(2));
            mat PKI = zeros(6,sigma.shape(1));
            mat sigma_cpp_list = carma::arr_to_mat(sigma);
            cube F_cpp_list = carma::arr_to_cube(F);            

            cout << sigma_cpp_list;
            cout << F_cpp_list;            
            cout << PKI;            

            for (unsigned int i=0; i < sigma_cpp_list.n_cols; i++) {
                vec sigma_cpp = sigma_cpp_list.col(i);
                mat F_cpp = F_cpp_list.slice(i);

                cout << sigma_cpp;
                cout << F_cpp;                    

                PKI.col(i) = simcoon::t2v_stress(simcoon::Cauchy2PKI(simcoon::v2t_stress(sigma_cpp), F_cpp, J));
                cout << PKI;
            }   
            cout << PKI;
            return carma::mat_to_arr(PKI, copy);            
        }
        else {
            throw std::invalid_argument("Invalid shape of the two-dimensional array. Expected n rows, 6 columns (1 per component of the symmetric stress tensor)");
        }
    }
    else {
        throw std::invalid_argument("Invalid number of dimensions for input 'sigma'. Expected 1 or 2.");
    }    
}

} //namepsace simpy