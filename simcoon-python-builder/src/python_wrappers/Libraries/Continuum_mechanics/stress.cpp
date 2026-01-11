
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

py::array_t<double> stress_convert(const py::array_t<double> &sigma, const py::array_t<double> &F, const string &converter_key, const double &J, const bool & copy) {

    std::map<string, int> list_stress_convert;    
    // NOTE: Be strict about keys. Using operator[] would insert unknown keys
    // with a default value (0), silently mapping to Cauchy2PKI.
    list_stress_convert = {
        {"Cauchy2PKI",0},
        {"Cauchy2PKII",1},
        // Historically misspelled as Kirchoff in simcoon
        {"Cauchy2Kirchoff",2},
        // Common spelling aliases
        {"Cauchy2Kirchhoff",2},
        {"Kirchoff2Cauchy",3},
        {"Kirchhoff2Cauchy",3},
        {"Kirchoff2PKI",4},
        {"Kirchhoff2PKI",4},
        {"Kirchoff2PKII",5},
        {"Kirchhoff2PKII",5},
        {"PKI2Kirchoff",6},
        {"PKI2Kirchhoff",6},
        {"PKII2Kirchoff",7},
        {"PKII2Kirchhoff",7},
        {"PKI2Cauchy",8},
        {"PKII2Cauchy",9}
    };

    const auto it = list_stress_convert.find(converter_key);
    if (it == list_stress_convert.end()) {
        throw std::invalid_argument(
            "Invalid input of the stress converter key: "
            "Cauchy2PKI, Cauchy2PKII, Cauchy2Kirchoff (or Cauchy2Kirchhoff), "
            "Kirchoff2Cauchy (or Kirchhoff2Cauchy), "
            "Kirchoff2PKI (or Kirchhoff2PKI), "
            "Kirchoff2PKII (or Kirchhoff2PKII), "
            "PKI2Kirchoff (or PKI2Kirchhoff), "
            "PKII2Kirchoff (or PKII2Kirchhoff), "
            "PKI2Cauchy or PKII2Cauchy"
        );
    }
    const int select = it->second;
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
            mat sigma_cpp = simcoon::v2t_stress(carma::arr_to_col(sigma));
            // Normalize F to a single 3x3 matrix. Accept several shapes from Python:
            // (3,3), (3,3,1), (1,3,3), or (1,3,3) etc.
            py::array_t<double> F_in = F;
            if ((int)F.ndim() == 3) {
                // If shape is (3,3,1) or (3,3,...) with last dim 1 -> squeeze
                if (F.shape(0) == 3 && F.shape(1) == 3 && F.shape(2) == 1) {
                    F_in = F.attr("squeeze")().cast<py::array_t<double>>();
                }
                // If shape is (1,3,3) -> squeeze to (3,3)
                else if (F.shape(0) == 1 && F.shape(1) == 3 && F.shape(2) == 3) {
                    F_in = F.attr("squeeze")().cast<py::array_t<double>>();
                }
                // If shape is (N,3,3) with N>1 that's ambiguous for single-sigma input
                else if (F.shape(0) != 3 || F.shape(1) != 3) {
                    throw std::invalid_argument("For single sigma input, F must be a single 3x3 matrix or a (3,3,1)/(1,3,3) shaped array");
                }
            }
            else if ((int)F.ndim() != 2) {
                throw std::invalid_argument("For single sigma input, F must have shape (3,3)");
            }
            mat F_cpp = carma::arr_to_mat(F_in);
            vec stress = simcoon::t2v_stress(functor_stress_converter(sigma_cpp, F_cpp, J));
            return carma::col_to_arr(stress, copy);
        }
        else {
            throw std::invalid_argument("Invalid size of the one-dimensional array. Expected 6");
        }         
    }
    else if (sigma.ndim() == 2) {
        if((sigma.shape(0) == 3)&&(sigma.shape(1) == 3)) {
            mat sigma_cpp = carma::arr_to_mat(sigma);
            // Normalize single F similar to single-sigma case
            py::array_t<double> F_in = F;
            if ((int)F.ndim() == 3) {
                if (F.shape(0) == 3 && F.shape(1) == 3 && F.shape(2) == 1) {
                    F_in = F.attr("squeeze")().cast<py::array_t<double>>();
                } else if (F.shape(0) == 1 && F.shape(1) == 3 && F.shape(2) == 3) {
                    F_in = F.attr("squeeze")().cast<py::array_t<double>>();
                } else {
                    throw std::invalid_argument("For single 3x3 sigma, F must be (3,3) or a squeezed equivalent");
                }
            }
            mat F_cpp = carma::arr_to_mat(F_in);
            mat stress = functor_stress_converter(sigma_cpp, F_cpp, J);
            return carma::mat_to_arr(stress, copy);
        }
        else if(sigma.shape(0) == 6) {
            // Batch case: accept F shaped (3,3,N) or (N,3,3). Normalize to (3,3,N)
            const unsigned int N = sigma.shape(1);
            py::array_t<double> F_in = F;
            if ((int)F.ndim() == 3) {
                // If F is (N,3,3) transpose to (3,3,N)
                if ((unsigned int)F.shape(0) == N && F.shape(1) == 3 && F.shape(2) == 3) {
                    // transpose axes (0,1,2) -> (1,2,0)
                    F_in = F.attr("transpose")(py::make_tuple(1,2,0)).cast<py::array_t<double>>();
                }
                // If F is (3,3,N) it's already correct
                else if ((unsigned int)F.shape(0) == 3 && F.shape(1) == 3 && (unsigned int)F.shape(2) == N) {
                    // ok
                }
                else {
                    throw std::invalid_argument("For batch sigma (6,N), F must be shaped (3,3,N) or (N,3,3) with matching N");
                }
            }
            else {
                throw std::invalid_argument("For batch sigma (6,N), F must be a 3D array shaped (3,3,N) or (N,3,3)");
            }

            assert(N == (unsigned int)F_in.shape(2));
            mat stress = zeros(6,N);
            mat sigma_cpp_list = carma::arr_to_mat_view(sigma);
            cube F_cpp_list = carma::arr_to_cube_view(F_in);

            for (unsigned int i=0; i < sigma_cpp_list.n_cols; i++) {
                vec sigma_cpp = sigma_cpp_list.unsafe_col(i);
                mat F_cpp = F_cpp_list.slice(i);
                stress.col(i) = simcoon::t2v_stress(functor_stress_converter(simcoon::v2t_stress(sigma_cpp), F_cpp, J));
            }

            return carma::mat_to_arr(stress, copy);
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