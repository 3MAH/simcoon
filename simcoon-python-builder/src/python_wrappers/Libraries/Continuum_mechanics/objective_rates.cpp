
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <string>
#include <carma>
#include <armadillo>

#include <simcoon/parallel.hpp>
#include <simcoon/exception.hpp>
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

//Log-strain concentration tensors A^R (rotated / log_R) and A^F (convected / log_F), 6x6 Voigt
py::array_t<double> A_R(const py::array_t<double> &F, const bool &copy) {
    mat F_cpp = carma::arr_to_mat(F);
    mat AR = simcoon::A_R(F_cpp);
    return carma::mat_to_arr(AR, copy);
}
py::array_t<double> A_F(const py::array_t<double> &F, const bool &copy) {
    mat F_cpp = carma::arr_to_mat(F);
    mat AF = simcoon::A_F(F_cpp);
    return carma::mat_to_arr(AF, copy);
}


//This function computes the logarithmic strain velocity and the logarithmic spin, along with the correct rotation increment
py::tuple objective_rate(const std::string& corate_name, const py::array_t<double> &F0, const py::array_t<double> &F1, const double &DTime, const bool &return_de, const unsigned int &n_threads) {
    std::map<string, int> list_corate;
    list_corate = { {"jaumann",0},{"green_naghdi",1},{"logarithmic",2},{"logarithmic_R",3},{"truesdell",4},{"logarithmic_F",5}, {"gn",1},{"log",2},{"log_R",3},{"log_F",5}};
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
        case 4: {
            corate_function = &simcoon::Truesdell;
            break;
        }
        case 5: {
            corate_function_2 = &simcoon::logarithmic_F;
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
        mat DR_N = zeros(3,3);
        mat N_1 = zeros(3,3);
        mat N_2 = zeros(3,3);

		switch (corate) {

            case 0: case 1: case 2: case 4: {
                corate_function(DR, D, Omega, DTime, F0_cpp, F1_cpp);
                break;
            }
            case 3: case 5: {
                corate_function_2(DR, N_1, N_2, D, Omega, DTime, F0_cpp, F1_cpp);
                break;
            }
        }

        if (return_de) {
            //also return the strain increment
            vec de = zeros(6);

    		switch (corate) {
                case 0: case 1: case 2: case 4: {
                    de = (0.5*DTime)*simcoon::t2v_strain((D+(DR*D*DR.t())));
                    break;
                }
                case 3: {
                    DR_N = simcoon::Hughes_Winget(N_1-N_2, DTime);
                    de = (0.5*DTime)*simcoon::t2v_strain((D+(DR*D*DR.t())));
                    de = simcoon::rotate_strain(de, DR_N);
                    break;
                }
                case 5: {
                    mat De_mat = (0.5*DTime)*(D+(DR*D*DR.t()));
                    DR_N = simcoon::Hughes_Winget(N_1-D, DTime);
                    mat inv_DR_N;
                    try {
                        inv_DR_N = inv(DR_N);
                    } catch (const std::runtime_error &e) {
                        cerr << "Error in inv: " << e.what() << endl;
                        throw simcoon::exception_inv("Error in inv function inside objective_rate (inv_DR_N).");
                    }
                    de = simcoon::t2v_strain(DR_N*De_mat*inv_DR_N);
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
        if(return_de) {
            de.set_size(6,nb_points);
        }
        mat I = eye(3,3);        

        if (F0.ndim() == 2) {
            mat vec_F0 = carma::arr_to_mat_view(F0);
            for (int pt = 0; pt < nb_points; pt++) {

        		switch (corate) {
                    case 0: case 1: case 2: case 4: {
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
                            mat DR_N = simcoon::Hughes_Winget(N_1.slice(pt)-N_2.slice(pt), DTime);
                            de_col = (0.5*DTime) * simcoon::t2v_strain(D.slice(pt)+(DR.slice(pt)*D.slice(pt)*DR.slice(pt).t()));
                            de_col = simcoon::rotate_strain(de_col, DR_N);
                        }
                        break;
                    }
                    case 5: {
                        corate_function_2(DR.slice(pt), N_1.slice(pt), N_2.slice(pt), D.slice(pt), Omega.slice(pt), DTime, vec_F0, F1_cpp.slice(pt));
                        if (return_de) {
                            vec de_col = de.unsafe_col(pt);
                            mat De_mat = (0.5*DTime)*(D.slice(pt)+(DR.slice(pt)*D.slice(pt)*DR.slice(pt).t()));
                            mat DR_N = simcoon::Hughes_Winget(N_1.slice(pt)-D.slice(pt), DTime);
                            mat inv_DR_N = inv(DR_N);
                            de_col = simcoon::t2v_strain(DR_N*De_mat*inv_DR_N);
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

                simcoon_parallel_for(nb_points, [&](int pt) {
                    switch (corate) {
                        case 0: case 1: case 2: case 4: {
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
                                mat DR_N = simcoon::Hughes_Winget(N_1.slice(pt)-N_2.slice(pt), DTime);
                                de_col = (0.5*DTime) * simcoon::t2v_strain(D.slice(pt)+(DR.slice(pt)*D.slice(pt)*DR.slice(pt).t()));
                                de_col = simcoon::rotate_strain(de_col, DR_N);
                            }
                            break;
                        }
                        case 5: {
                            corate_function_2(DR.slice(pt), N_1.slice(pt), N_2.slice(pt), D.slice(pt), Omega.slice(pt), DTime, vec_F0, F1_cpp.slice(pt));
                            if (return_de) {
                                vec de_col = de.unsafe_col(pt);
                                mat De_mat = (0.5*DTime)*(D.slice(pt)+(DR.slice(pt)*D.slice(pt)*DR.slice(pt).t()));
                                mat DR_N = simcoon::Hughes_Winget(N_1.slice(pt)-D.slice(pt), DTime);
                                mat inv_DR_N = inv(DR_N);
                                de_col = simcoon::t2v_strain(DR_N*De_mat*inv_DR_N);
                            }
                            break;
                        }
                    }
                });
            }
            else {
                simcoon_parallel_for(nb_points, [&](int pt) {
                    switch (corate) {
                        case 0: case 1: case 2: case 4: {
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
                                mat DR_N = simcoon::Hughes_Winget(N_1.slice(pt)-N_2.slice(pt), DTime);
                                de_col = (0.5*DTime) * simcoon::t2v_strain(D.slice(pt)+(DR.slice(pt)*D.slice(pt)*DR.slice(pt).t()));
                                de_col = simcoon::rotate_strain(de_col, DR_N);
                            }
                            break;
                        }
                        case 5: {
                            corate_function_2(DR.slice(pt), N_1.slice(pt), N_2.slice(pt), D.slice(pt), Omega.slice(pt), DTime, F0_cpp.slice(pt), F1_cpp.slice(pt));
                            if (return_de) {
                                vec de_col = de.unsafe_col(pt);
                                mat De_mat = (0.5*DTime)*(D.slice(pt)+(DR.slice(pt)*D.slice(pt)*DR.slice(pt).t()));
                                mat DR_N = simcoon::Hughes_Winget(N_1.slice(pt)-D.slice(pt), DTime);
                                mat inv_DR_N = inv(DR_N);
                                de_col = simcoon::t2v_strain(DR_N*De_mat*inv_DR_N);
                            }
                            break;
                        }
                    }
                });
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
    list_Lt_convert = { {"Dsigma_LieDD_2_DSDE",0}, {"DsigmaDe_2_DSDE",1},{"DsigmaDe_JaumannDD_2_DSDE",2}, {"Dsigma_LieDD_Dsigma_JaumannDD",3}, {"Dsigma_LieDD_Dsigma_GreenNaghdiDD",4}, {"Dsigma_LieDD_Dsigma_logarithmicDD",5}, {"DsigmaDe_GreenNaghdiDD_2_DSDE",6}, {"DSDE_2_Dsigma_GreenNaghdiDD",7}, {"DSDE_2_Dsigma_JaumannDD",8}, {"DSDE_2_Dsigma_LieDD",9}, {"DSDE_2_Dsigma_logarithmicDD",10}};
	int select = list_Lt_convert [converter_key];

    // The box tangent convention is Lt = d(tau_hat)/d(De): the Kirchhoff, log/Hencky-rate,
    // no-J corotational tangent that every simcoon UMAT now returns. Split by role:
    //  * INVERSE keys (consume the box Lt -> material dS/dE): use the no-J Dtau_* family so
    //    the box's Kirchhoff tangent is NOT spuriously multiplied by J (the J double-count
    //    fix). Stress is promoted to Kirchhoff (tau = J*sigma) internally.
    //  * FORWARD / CROSS keys (produce dsigma/dCorate, the Cauchy spatial tangent fedoo
    //    consumes -- "dsigma_dD"): KEEP the Cauchy (1/J) Dsigma_* family, taking the Cauchy
    //    sigma directly. (fedoo: dS/dE = DsigmaDe_2_DSDE(box_Lt); dsigma/dD = DSDE_2_Dsigma_*.)
    auto convert_pt = [select](const mat &Lt_pt, const mat &F_pt, const mat &sig_pt) -> mat {
        double Jdet = det(F_pt);
        mat tau = Jdet*sig_pt;                 // Kirchhoff = J * Cauchy (both spatial)
        switch (select) {
            // inverse: box d(tau_hat)/d(De) (Kirchhoff, no-J) -> material dS/dE
            case 0:  return simcoon::Dtau_LieDD_2_DSDE(Lt_pt, F_pt);
            case 1:  return simcoon::DtauDe_2_DSDE(Lt_pt, simcoon::get_BBBB(F_pt), F_pt, tau);
            case 2:  return simcoon::DtauDe_JaumannDD_2_DSDE(Lt_pt, F_pt, tau);
            case 6:  return simcoon::DtauDe_GreenNaghdiDD_2_DSDE(Lt_pt, F_pt, tau);
            // Cauchy spatial-tangent menu (dsigma/dD, 1/J) -- fedoo's material Jacobian
            case 3:  return simcoon::Dsigma_LieDD_Dsigma_JaumannDD(Lt_pt, sig_pt);
            case 4:  return simcoon::Dsigma_LieDD_Dsigma_GreenNaghdiDD(Lt_pt, F_pt, sig_pt);
            case 5:  return simcoon::Dsigma_LieDD_Dsigma_logarithmicDD(Lt_pt, F_pt, sig_pt);
            case 7:  return simcoon::DSDE_2_Dsigma_GreenNaghdiDD(Lt_pt, F_pt, sig_pt);
            case 8:  return simcoon::DSDE_2_Dsigma_JaumannDD(Lt_pt, F_pt, sig_pt);
            case 9:  return simcoon::DSDE_2_Dsigma_LieDD(Lt_pt, F_pt);
            case 10: return simcoon::DSDE_2_Dsigma_logarithmicDD(Lt_pt, F_pt, sig_pt);
        }
        return Lt_pt;
    };

    if (Lt.ndim() == 2) {
        if ((F.ndim() != 2) || (stress.ndim() != 1))  {
            throw std::invalid_argument("the number of dim of Lt, F and stress are not consistent");
        }
        mat F_cpp = carma::arr_to_mat_view(F);
        mat Lt_cpp = carma::arr_to_mat_view(Lt);
        vec stress_v = carma::arr_to_col_view(stress);
        mat sig_cpp = simcoon::v2t_stress(stress_v);
        mat Lt_converted = convert_pt(Lt_cpp, F_cpp, sig_cpp);
        return carma::mat_to_arr(Lt_converted, false);
    }
    else if (Lt.ndim() == 3) {
        cube F_cpp = carma::arr_to_cube_view(F);
        cube Lt_cpp = carma::arr_to_cube_view(Lt);
        mat stress_cpp = carma::arr_to_mat_view(stress);
        int nb_points = Lt_cpp.n_slices;
        cube Lt_converted = zeros(6,6,nb_points);
        for (int pt = 0; pt < nb_points; pt++) {
            mat sig_pt = simcoon::v2t_stress(stress_cpp.unsafe_col(pt));
            Lt_converted.slice(pt) = convert_pt(Lt_cpp.slice(pt), F_cpp.slice(pt), sig_pt);
        }
        return carma::cube_to_arr(Lt_converted, false);
    }
    throw std::invalid_argument("Lt.ndim() must be 2 or 3");
}


} //namepsace simpy
