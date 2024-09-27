
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <string>
#include <carma>
#include <armadillo>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <simcoon/Simulation/Maths/rotation.hpp>
#include <simcoon/Continuum_mechanics/Functions/objective_rates.hpp>
#include <simcoon/python_wrappers/Libraries/Continuum_mechanics/objective_rates.hpp>
#include <simcoon/Continuum_mechanics/Functions/transfer.hpp>

using namespace std;
using namespace arma;
namespace py=pybind11;

namespace simpy {

//This function computes the logarithmic strain velocity and the logarithmic spin, along with the correct rotation increment
py::tuple logarithmic(const py::array_t<double> &F0, const py::array_t<double> &F1, const double &DTime, const bool &copy) {
    mat F0_cpp = carma::arr_to_mat(F0);
    mat F1_cpp = carma::arr_to_mat(F1);
    mat DR = zeros(3,3);
    mat D = zeros(3,3);
    mat Omega = zeros(3,3);
    simcoon::logarithmic(DR, D, Omega, DTime, F0_cpp, F1_cpp);
    return py::make_tuple(carma::mat_to_arr(D, copy), carma::mat_to_arr(DR, copy), carma::mat_to_arr(Omega, copy));
}

//This function computes the logarithmic strain velocity and the logarithmic spin, along with the correct rotation increment
py::tuple logarithmic_R(const py::array_t<double> &F0, const py::array_t<double> &F1, const double &DTime, const bool &copy) {
    mat F0_cpp = carma::arr_to_mat(F0);
    mat F1_cpp = carma::arr_to_mat(F1);
    mat DR = zeros(3,3);
    mat D = zeros(3,3);
    mat N_1 = zeros(3,3);
    mat N_2 = zeros(3,3);    
    mat Omega = zeros(3,3);
    simcoon::logarithmic_R(DR, D, N_1, N_2, Omega, DTime, F0_cpp, F1_cpp);
    return py::make_tuple(carma::mat_to_arr(D, copy), carma::mat_to_arr(DR, copy), carma::mat_to_arr(Omega, copy), carma::mat_to_arr(N_1, copy), carma::mat_to_arr(N_2, copy));
}


//This function computes the logarithmic strain velocity and the logarithmic spin, along with the correct rotation increment
py::tuple objective_rate(const std::string& corate_name, const py::array_t<double> &F0, const py::array_t<double> &F1, const double &DTime, const bool &return_de, const unsigned int &n_threads) {
    std::map<string, int> list_corate;
    list_corate = { {"jaumann",0},{"green_naghdi",1},{"logarithmic",2},{"logarithmic_R",3}, {"gn",1},{"log",2},{"log_R",3}};
	int corate = list_corate[corate_name];

    void (*corate_function)(mat &, mat &, mat &, const double &, const mat &, const mat &); 
    void (*corate_function_2)(mat &, mat &, mat &, mat &, mat &, const double &, const mat &, const mat &);     
    switch (corate) {

        case 0: {
            corate_function = &simcoon::Jaumann;
            break;
        }
        case 1: {
            corate_function = &simcoon::Green_Naghdi;
            break;
        }
        case 2: {
            corate_function = &simcoon::logarithmic;
            break;
        }
        case 3: {
            corate_function_2 = &simcoon::logarithmic_R;
            break;
        }        
    }
    
    if (F1.ndim() == 2) {            
        if (F0.ndim() != 2) {
            throw std::invalid_argument("the number of dim of F1 should be the same as F0");
        }

        mat F0_cpp = carma::arr_to_mat_view(F0);
        mat F1_cpp = carma::arr_to_mat_view(F1);
        mat DR = zeros(3,3);
        mat D = zeros(3,3);
        mat Omega = zeros(3,3); 

        mat N_1 = zeros(3,3);
        mat N_2 = zeros(3,3);

		switch (corate) {

            case 0: case 1: case 2: {
                corate_function(DR, D, Omega, DTime, F0_cpp, F1_cpp);
                break;
            }
            case 3: {
                corate_function_2(DR, N_1, N_2, D, Omega, DTime, F0_cpp, F1_cpp);
                break;
            }
        }

        if (return_de) {
            //also return the strain increment
            vec de = zeros(6);

    		switch (corate) {
                case 0: case 1: case 2: {
                //could use simcoon::Delta_log_strain(D, Omega, DTime) but it would recompute DR (waste of time).
                //vec de = simcoon::t2v_strain(simcoon::Delta_log_strain(D, Omega, DTime));                 
                    de = (0.5*DTime)*simcoon::t2v_strain((D+(DR*D*DR.t())));
                    break;
                }
                case 3: {
                    mat I = eye(3,3);
                    mat DR_N = (inv(I-0.5*DTime*(N_1-N_2)))*(I+0.5*DTime*(N_1-N_2));
                    de = (0.5*DTime)*simcoon::t2v_strain((D+(DR*D*DR.t())));
                    de = simcoon::rotate_strain(de, DR_N);
                    break;
                }
            }
            return py::make_tuple(carma::col_to_arr(de,false), carma::mat_to_arr(D, false), carma::mat_to_arr(DR, false), carma::mat_to_arr(Omega, false));
        }
        else{
            return py::make_tuple(carma::mat_to_arr(D, false), carma::mat_to_arr(DR, false), carma::mat_to_arr(Omega, false));
        }
        
    }
    else if (F1.ndim() == 3) {
        cube F1_cpp = carma::arr_to_cube_view(F1);            
        int nb_points = F1_cpp.n_slices;
        cube DR(3,3,nb_points);
        cube D(3,3,nb_points);            
        cube N_1(3,3,nb_points);
        cube N_2(3,3,nb_points);                            
        cube Omega = zeros(3,3, nb_points); 	
        mat de;
        if(return_de) de.set_size(6,nb_points);
        mat DR_N = zeros(3,3);
        mat I = eye(3,3);        

        if (F0.ndim() == 2) {
            mat vec_F0 = carma::arr_to_mat_view(F0);
            for (int pt = 0; pt < nb_points; pt++) {

        		switch (corate) {
                    case 0: case 1: case 2: {
                        corate_function(DR.slice(pt), D.slice(pt), Omega.slice(pt), DTime, vec_F0, F1_cpp.slice(pt));
                        if (return_de) {
                            vec de_col = de.unsafe_col(pt);
                            de_col = (0.5*DTime) * simcoon::t2v_strain(D.slice(pt)+(DR.slice(pt)*D.slice(pt)*DR.slice(pt).t()));                                                          
                        }
                        break;
                    }
                    case 3: {
                        corate_function_2(DR.slice(pt), N_1.slice(pt), N_2.slice(pt), D.slice(pt), Omega.slice(pt), DTime, vec_F0, F1_cpp.slice(pt));
                        if (return_de) {
                            vec de_col = de.unsafe_col(pt);                                
                            DR_N = (inv(I-0.5*DTime*(N_1.slice(pt)-N_2.slice(pt))))*(I+0.5*DTime*(N_1.slice(pt)-N_2.slice(pt)));
                            de_col = (0.5*DTime) * simcoon::t2v_strain(D.slice(pt)+(DR.slice(pt)*D.slice(pt)*DR.slice(pt).t()));
                            de_col = simcoon::rotate_strain(de_col, DR_N);
                        }
                        break;
                    }
                }
            }
        }
        else if (F0.ndim() == 3) {
            cube F0_cpp = carma::arr_to_cube_view(F0); 
            if (F0_cpp.n_slices==1) {
                mat vec_F0 = F0_cpp.slice(0);

                #ifdef _OPENMP
                int max_threads = omp_get_max_threads();
                omp_set_num_threads(4);
                py::gil_scoped_release release;

                omp_set_max_active_levels(3);
                #pragma omp parallel for shared(DR, D, Omega, F1_cpp)    
    			#endif
                for (int pt = 0; pt < nb_points; pt++) {

            		switch (corate) {
                        case 0: case 1: case 2: {
                            corate_function(DR.slice(pt), D.slice(pt), Omega.slice(pt), DTime, vec_F0, F1_cpp.slice(pt));
                            if (return_de) {
                                vec de_col = de.unsafe_col(pt);
                                de_col = (0.5*DTime) * simcoon::t2v_strain(D.slice(pt)+(DR.slice(pt)*D.slice(pt)*DR.slice(pt).t()));                                                          
                            }
                            break;
                        }
                        case 3: {
                            corate_function_2(DR.slice(pt), N_1.slice(pt), N_2.slice(pt), D.slice(pt), Omega.slice(pt), DTime, vec_F0, F1_cpp.slice(pt));
                            if (return_de) {
                                vec de_col = de.unsafe_col(pt);                                
                                DR_N = (inv(I-0.5*DTime*(N_1.slice(pt)-N_2.slice(pt))))*(I+0.5*DTime*(N_1.slice(pt)-N_2.slice(pt)));
                                de_col = (0.5*DTime) * simcoon::t2v_strain(D.slice(pt)+(DR.slice(pt)*D.slice(pt)*DR.slice(pt).t()));
                                de_col = simcoon::rotate_strain(de_col, DR_N);
                            }
                            break;
                        }
                    }
                }
                #ifdef _OPENMP
                py::gil_scoped_acquire acquire;					
                omp_set_num_threads(max_threads);			
    			#endif                                
            }
            else {
                #ifdef _OPENMP                
                int max_threads = omp_get_max_threads();
                omp_set_num_threads(4);
                py::gil_scoped_release release;

                omp_set_max_active_levels(3);
                #pragma omp parallel for shared(DR, D, Omega, F0_cpp, F1_cpp)      
    			#endif
                for (int pt = 0; pt < nb_points; pt++) {

            		switch (corate) {                    
                        case 0: case 1: case 2: {
                            corate_function(DR.slice(pt), D.slice(pt), Omega.slice(pt), DTime, F0_cpp.slice(pt), F1_cpp.slice(pt));
                            if (return_de) {
                                vec de_col = de.unsafe_col(pt);
                                de_col = (0.5*DTime) * simcoon::t2v_strain(D.slice(pt)+(DR.slice(pt)*D.slice(pt)*DR.slice(pt).t()));                                                          
                            }
                            break;
                        }
                        case 3: {
                            corate_function_2(DR.slice(pt), N_1.slice(pt), N_2.slice(pt), D.slice(pt), Omega.slice(pt), DTime, F0_cpp.slice(pt), F1_cpp.slice(pt));
                            if (return_de) {
                                vec de_col = de.unsafe_col(pt);                                
                                DR_N = (inv(I-0.5*DTime*(N_1.slice(pt)-N_2.slice(pt))))*(I+0.5*DTime*(N_1.slice(pt)-N_2.slice(pt)));
                                de_col = (0.5*DTime) * simcoon::t2v_strain(D.slice(pt)+(DR.slice(pt)*D.slice(pt)*DR.slice(pt).t()));
                                de_col = simcoon::rotate_strain(de.col(pt), DR_N);
                            }
                            break;
                        }
                    }
                }
                #ifdef _OPENMP                                
                py::gil_scoped_acquire acquire;					
                omp_set_num_threads(max_threads);	
    			#endif                		
            }
        }
        if (return_de){	                     
            return py::make_tuple(carma::mat_to_arr(de, false), carma::cube_to_arr(D, false), carma::cube_to_arr(DR, false), carma::cube_to_arr(Omega, false));
        }
        else{
            return py::make_tuple(carma::cube_to_arr(D, false), carma::cube_to_arr(DR, false), carma::cube_to_arr(Omega, false));
        }        
    }
}

//This function computes the gradient of displacement (Eulerian) from the deformation gradient 
py::array_t<double> Delta_log_strain(const py::array_t<double> &D, const py::array_t<double> &Omega, const double &DTime, const bool &copy) {
    mat D_cpp = carma::arr_to_mat(D);
    mat Omega_cpp = carma::arr_to_mat(Omega);
    mat Delta_log_strain = simcoon::Delta_log_strain(D_cpp, Omega_cpp, DTime);
    return carma::mat_to_arr(Delta_log_strain, copy);
}

//This function computes the logarithmic strain velocity and the logarithmic spin, along with the correct rotation increment
py::array_t<double> Lt_convert(const py::array_t<double> &Lt, const py::array_t<double> &F, const py::array_t<double> &stress, const std::string &converter_key) {
    std::map<string, int> list_Lt_convert;
    list_Lt_convert = { {"DsigmaDe_2_DSDE",0},{"DsigmaDe_JaumannDD_2_DSDE",1}, {"DtauDe_JaumannDD_2_DSDE",2}, {"DsigmaDe_2_DtauDe",3}};
	int select = list_Lt_convert [converter_key];
    mat (*convert_function)(const mat &, const mat &, const mat &); 
    mat (*convert_function2)(const mat &, const double &);     

    switch (select) {
        case 0: {
            convert_function = &simcoon::DsigmaDe_2_DSDE;
            break;
        }
        case 1: {
            convert_function = &simcoon::DsigmaDe_JaumannDD_2_DSDE;
            break;
        }
        case 2: {
            convert_function = &simcoon::DtauDe_JaumannDD_2_DSDE;
            break;
        }
        case 3: {
            convert_function2 = &simcoon::DsigmaDe_2_DtauDe;            
            break;
        }
    }

    if (Lt.ndim() == 2) {            
        if ((F.ndim() != 2) or (stress.ndim() != 1))  {
            throw std::invalid_argument("the number of dim of Lt, F and stress are not consistent");
        }

        mat F_cpp = carma::arr_to_mat_view(F);
        mat Lt_cpp = carma::arr_to_mat_view(Lt);
        vec stress_cpp = carma::arr_to_col_view(stress);
        mat Lt_converted(6,6);

        switch (select) {             
            case 0: case 1: case 2: {
                Lt_converted = convert_function(Lt_cpp, F_cpp, stress_cpp);
                break;
            }
            case 3: {
                Lt_converted = convert_function2(F_cpp, det(F_cpp));
                break;
            }      
        }          
        return carma::mat_to_arr(Lt_converted,false);
    }
    else if (Lt.ndim() == 3) {
        cube F_cpp = carma::arr_to_cube_view(F);
        cube Lt_cpp = carma::arr_to_cube_view(Lt);
        mat stress_cpp = carma::arr_to_mat_view(stress);
        int nb_points = Lt_cpp.n_slices;
        cube Lt_converted(6,6,nb_points);

        mat stress_pt;

        #ifdef _OPENMP
        int max_threads = omp_get_max_threads();
        omp_set_num_threads(4);
        py::gil_scoped_release release;

        omp_set_max_active_levels(3);
        #pragma omp parallel for shared(Lt_converted)  
        #endif

        for (int pt = 0; pt < nb_points; pt++) {
            //vec stress_pt = stress_cpp.unsafe_col(pt); 
            stress_pt = simcoon::v2t_stress(stress_cpp.col(pt));
            switch (select) {             
                case 0: case 1: case 2: {
                    Lt_converted.slice(pt) = convert_function(Lt_cpp.slice(pt), F_cpp.slice(pt), stress_pt);
                    break;
                }
                case 3: {
                    Lt_converted.slice(pt) = convert_function2(Lt_cpp.slice(pt), det(F_cpp.slice(pt)));
                    break;
                }      
            }          
        }        
        #ifdef _OPENMP
        py::gil_scoped_acquire acquire;					
        omp_set_num_threads(max_threads);			                     
        #endif
        return carma::cube_to_arr(Lt_converted,false);
    }
}


} //namepsace simpy
