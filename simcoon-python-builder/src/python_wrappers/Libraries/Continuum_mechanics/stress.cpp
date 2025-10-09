
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <string>
#include <armadillo>
#include <simcoon/python_wrappers/conversion_helpers.hpp>
#include <assert.h>

#include <simcoon/Continuum_mechanics/Functions/stress.hpp>
#include <simcoon/Continuum_mechanics/Functions/transfer.hpp>
#include <simcoon/python_wrappers/Libraries/Continuum_mechanics/stress.hpp>

using namespace std;
using namespace arma;
namespace py=pybind11;

namespace simpy {

py::array_t<double> stress_convert(const py::array_t<double> &sigma, const py::array_t<double> &F, const string &converter_key, const double &J, const bool & copy) {

    std::map<string, int> list_stress_convert;    
    list_stress_convert = {{"Cauchy2PKI",0},{"Cauchy2PKII",1},{"Cauchy2Kirchoff",2},{"Kirchoff2Cauchy",3},{"Kirchoff2PKI",4},{"Kirchoff2PKII",5},{"PKI2Kirchoff",6},{"PKII2Kirchoff",7},{"PKI2Cauchy",8},{"PKII2Cauchy",9}};

    int select = list_stress_convert[converter_key];
    mat (*functor_stress_converter)(const mat &, const mat &, const double &);

    switch (select) {
        case 0: {
            functor_stress_converter = &simcoon::Cauchy2PKI;
            break;
        }
        case 1: {
            functor_stress_converter = &simcoon::Cauchy2PKII;
            break;
        }
        case 2: {
            functor_stress_converter = &simcoon::Cauchy2Kirchoff;
            break;
        }
        case 3: {
            functor_stress_converter = &simcoon::Kirchoff2Cauchy;
            break;
        }  
        case 4: {
            functor_stress_converter = &simcoon::Kirchoff2PKI;
            break;
        }
        case 5: {
            functor_stress_converter = &simcoon::Kirchoff2PKII;
            break;
        }
        case 6: {
            functor_stress_converter = &simcoon::PKI2Kirchoff;
            break;
        }
        case 7: {
            functor_stress_converter = &simcoon::PKII2Kirchoff;
            break;
        }
        case 8: {
            functor_stress_converter = &simcoon::PKI2Cauchy;
            break;
        }
        case 9: {
            functor_stress_converter = &simcoon::PKII2Cauchy;
            break;
        }            
        default: {
            throw std::invalid_argument("Invalid input of the stress converter key: Cauchy2PKI, Cauchy2PKII, Cauchy2Kirchoff, Kirchoff2Cauchy, Kirchoff2PKI, Kirchoff2PKII, PKI2Kirchoff, PKII2Kirchoff, PKI2Cauchy or PKII2Cauchy");                
        }
    }

    /*auto functor_stress_converter = [](const arma::mat &sigma_cpp, const arma::mat &F_cpp, const double &J, const int &select) {
        switch (select) {

            case 0: {
                return simcoon::Cauchy2PKI(sigma_cpp, F_cpp, J);
            }
            case 1: {
                return simcoon::Cauchy2PKII(sigma_cpp, F_cpp, J);
            }
            case 2: {
                return simcoon::Cauchy2Kirchoff(sigma_cpp, F_cpp, J);
            }
            case 3: {
                return simcoon::Kirchoff2Cauchy(sigma_cpp, F_cpp, J);
            }  
            case 4: {
                return simcoon::Kirchoff2PKI(sigma_cpp, F_cpp, J);
            }
            case 5: {
                return simcoon::Kirchoff2PKII(sigma_cpp, F_cpp, J);
            }
            case 6: {
                return simcoon::PKI2Kirchoff(sigma_cpp, F_cpp, J);
            }
            case 7: {
                return simcoon::PKII2Kirchoff(sigma_cpp, F_cpp, J);
            }
            case 8: {
                return simcoon::PKI2Cauchy(sigma_cpp, F_cpp, J);
            }
            case 9: {
                return simcoon::PKII2Cauchy(sigma_cpp, F_cpp, J);
            }            
            default: {
                throw std::invalid_argument("Invalid input of the stress converter key: Cauchy2PKI, Cauchy2PKII, Cauchy2Kirchoff, Kirchoff2Cauchy, Kirchoff2PKI, Kirchoff2PKII, PKI2Kirchoff, PKII2Kirchoff, PKI2Cauchy or PKII2Cauchy");                
            }
        }
    };*/

    if(sigma.ndim() == 1) {
        if(sigma.size() == 6) {
            mat sigma_cpp = simcoon::v2t_stress(simpy::arr_to_col(sigma));
            mat F_cpp = simpy::arr_to_mat(F);            
            vec stress = simcoon::t2v_stress(functor_stress_converter(sigma_cpp, F_cpp, J));
            return simpy::col_to_arr(stress, copy);
        }
        else {
            throw std::invalid_argument("Invalid size of the one-dimensional array. Expected 6");
        }         
    }
    else if (sigma.ndim() == 2) {
        if((sigma.shape(0) == 3)&&(sigma.shape(1) == 3)) {
            mat sigma_cpp = simpy::arr_to_mat(sigma);
            mat F_cpp = simpy::arr_to_mat(F);
            mat stress = functor_stress_converter(sigma_cpp, F_cpp, J);
            return simpy::mat_to_arr(stress, copy);
        }
        else if(sigma.shape(0) == 6) {
            assert(sigma.shape(1) == F.shape(2));
            mat stress = zeros(6,sigma.shape(1));
            mat sigma_cpp_list = simpy::arr_to_mat_view(sigma);
            cube F_cpp_list = simpy::arr_to_cube_view(F);            

            for (unsigned int i=0; i < sigma_cpp_list.n_cols; i++) {
                vec sigma_cpp = sigma_cpp_list.unsafe_col(i);
                mat F_cpp = F_cpp_list.slice(i);

                //cout << sigma_cpp;
                //cout << F_cpp;                    

                stress.col(i) = simcoon::t2v_stress(functor_stress_converter(simcoon::v2t_stress(sigma_cpp), F_cpp, J));
            }   

            return simpy::mat_to_arr(stress, copy);            
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