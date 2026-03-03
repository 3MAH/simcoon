#include <pybind11/embed.h>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <optional>

#include <carma>
#include <armadillo>

#include <simcoon/parameter.hpp>
#include <simcoon/omp_compat.hpp>

#include <simcoon/python_wrappers/Libraries/Continuum_mechanics/umat.hpp>

#include <simcoon/Continuum_mechanics/Umat/Mechanical/External/external_umat.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Elasticity/elastic_isotropic.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Elasticity/elastic_transverse_isotropic.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Elasticity/elastic_orthotropic.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Plasticity/plastic_isotropic_ccp.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Plasticity/plastic_kin_iso_ccp.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Plasticity/Hill_chaboche_ccp.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Plasticity/Ani_chaboche_ccp.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Plasticity/DFA_chaboche_ccp.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Plasticity/plastic_chaboche_ccp.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Plasticity/Generic_chaboche_ccp.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Plasticity/Hill_isoh.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Plasticity/Hill_isoh_Nfast.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/SMA/unified_T.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/SMA/SMA_mono.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Damage/damage_LLD_0.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Viscoelasticity/Zener_fast.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Viscoelasticity/Zener_Nfast.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Viscoelasticity/Prony_Nfast.hpp>

#include <simcoon/Continuum_mechanics/Umat/Finite/generic_hyper_invariants.hpp>

#include <simcoon/Simulation/Maths/rotation.hpp> //for rotate_strain

#include <simcoon/Continuum_mechanics/Umat/Thermomechanical/External/external_umat.hpp>
#include <simcoon/Continuum_mechanics/Umat/Thermomechanical/Elasticity/elastic_isotropic.hpp>
#include <simcoon/Continuum_mechanics/Umat/Thermomechanical/Elasticity/elastic_isotropic.hpp>
#include <simcoon/Continuum_mechanics/Umat/Thermomechanical/Elasticity/elastic_transverse_isotropic.hpp>
#include <simcoon/Continuum_mechanics/Umat/Thermomechanical/Elasticity/elastic_orthotropic.hpp>
#include <simcoon/Continuum_mechanics/Umat/Thermomechanical/Plasticity/plastic_isotropic_ccp.hpp>
#include <simcoon/Continuum_mechanics/Umat/Thermomechanical/Plasticity/plastic_kin_iso_ccp.hpp>
#include <simcoon/Continuum_mechanics/Umat/Thermomechanical/Viscoelasticity/Zener_fast.hpp>
#include <simcoon/Continuum_mechanics/Umat/Thermomechanical/Viscoelasticity/Zener_Nfast.hpp>
#include <simcoon/Continuum_mechanics/Umat/Thermomechanical/Viscoelasticity/Prony_Nfast.hpp>

/*#include <simcoon/Continuum_mechanics/Functions/objective_rates.hpp>
#include <simcoon/Continuum_mechanics/Functions/transfer.hpp>
#include <simcoon/Continuum_mechanics/Functions/kinematics.hpp>
#include <simcoon/python_wrappers/Libraries/Continuum_mechanics/kinematics.hpp>
#include <simcoon/Continuum_mechanics/Functions/stress.hpp>
*/

#include <simcoon/Continuum_mechanics/Umat/Thermomechanical/SMA/unified_T.hpp>

using namespace std;
using namespace arma;
namespace py=pybind11;

namespace simpy {
	
	py::tuple launch_umat(const std::string &umat_name_py, const py::array_t<double> &etot_py, const py::array_t<double> &Detot_py, const py::array_t<double> &F0_py, const py::array_t<double> &F1_py, const py::array_t<double> &sigma_py, const py::array_t<double> &DR_py, const py::array_t<double> &props_py, const py::array_t<double> &statev_py, const float Time, const float DTime, const py::array_t<double> &Wm_py, const std::optional<py::array_t<double>> &T_py, const int &ndi, const unsigned int &n_threads){
		//Get the id of umat

		std::map<string, int> list_umat;
		list_umat = { {"UMEXT",0},{"UMABA",1},{"ELISO",2},{"ELIST",3},{"ELORT",4},{"EPICP",5},{"EPKCP",6},{"EPCHA",7},{"EPHIL",8},{"EPHAC",9},{"EPANI",10},{"EPDFA",11},{"EPHIN",12},{"SMAUT",13},{"SMANI",13},{"SMADI",13},{"SMADC",13},{"SMAAI",13},{"SMAAC",13},{"LLDM0",15},{"ZENER",16},{"ZENNK",17},{"PRONK",18},{"SMAMO",19},{"SMAMC",20},{"NEOHC",21},{"MOORI",22},{"YEOHH",23},{"ISHAH",24},{"GETHH",25},{"SWANH",26},{"EPCHG",27},{"MIHEN",100},{"MIMTN",101},{"MISCN",103},{"MIPLN",104} }; // TODO_2.0 SMAUT and SMANI compatibility to be removed in release 2.0 
		int id_umat = list_umat[umat_name_py];
		int arguments_type; //depends on the argument used in the umat

		// Unified small-strain function pointer: (umat_name, Etot, DEtot, sigma, Lt, L, DR, nprops, props, nstatev, statev, T, DT, Time, DTime, Wm, Wm_r, Wm_ir, Wm_d, ndi, nshr, start, tnew_dt)
		void (*umat_function)(const std::string &, const arma::vec &, const arma::vec &, arma::vec &, arma::mat &, arma::mat &, const arma::mat &, const int &, const arma::vec &, const int &, arma::vec &, const double &, const double &, const double &, const double &, double &, double &, double &, double &, const int &, const int &, const bool &, double &);
		// Unified finite-strain function pointer: (umat_name, etot, Detot, F0, F1, sigma, Lt, L, DR, nprops, props, nstatev, statev, T, DT, Time, DTime, Wm, Wm_r, Wm_ir, Wm_d, ndi, nshr, start, tnew_dt)
		void (*umat_function_finite)(const std::string &, const arma::vec &, const arma::vec &, const arma::mat &, const arma::mat &, arma::vec &, arma::mat &, arma::mat &, const arma::mat &, const int &, const arma::vec &, const int &, arma::vec &, const double &, const double &, const double &, const double &, double &, double &, double &, double &, const int &, const int &, const bool &, double &);
		const int ncomp=6;
		int nshr;
		if (ndi==3) {
			nshr=3;
		} else if (ndi==1) {
			nshr=0;
		} else if (ndi==2) {
			nshr=1;
		} else {
			throw std::invalid_argument( "ndi should be 1, 2 or 3 dimenions" );
		}

		bool start = true;

		if (Time > simcoon::limit) {
			start = false;
		}

		double tnew_dt = 0;//usefull ?		
		//bool use_temp;
		//if (T.n_elem == 0.) use_temp = false; 
		//else use_temp = true;

		bool use_temp = false;
		vec vec_T;
		if (T_py.has_value()) {
			vec_T = carma::arr_to_col_view(T_py.value());
			use_temp = true;
		}
		else {
			use_temp = false; 
		}

		mat list_etot = carma::arr_to_mat_view(etot_py);		
		int nb_points = list_etot.n_cols; //number of material points
		mat list_Detot = carma::arr_to_mat_view(Detot_py); 
		mat list_sigma = carma::arr_to_mat(std::move(sigma_py)); //copy data because values are changed by the umat and returned to python
		cube DR = carma::arr_to_cube_view(DR_py); 
		cube F0, F1;

		vec props;
		mat list_props = carma::arr_to_mat_view(props_py);
		auto shape = props_py.shape();

		bool unique_props = false;
		if (shape[1] == 1) {
			props = list_props.col(0);
			unique_props = true;
		}		

		mat list_statev = carma::arr_to_mat(std::move(statev_py)); //copy data because values are changed by the umat and returned to python
		mat list_Wm = carma::arr_to_mat(std::move(Wm_py)); //copy data because values are changed by the umat and returned to python
		cube L(ncomp, ncomp, nb_points);
		cube Lt(ncomp, ncomp, nb_points);
		int nprops = list_props.n_rows;
		int nstatev = list_statev.n_rows;

		switch (id_umat) {
			case 2: {
				umat_function = &simcoon::umat_elasticity_iso;
				arguments_type = 1;
				break;
			}
			case 3: {
				umat_function = &simcoon::umat_elasticity_trans_iso;
				arguments_type = 1;
				break;
			}
			case 4: {
				umat_function = &simcoon::umat_elasticity_ortho;
				arguments_type = 1;
				break;
			}
			case 5: {
				umat_function = &simcoon::umat_plasticity_iso_CCP;
				arguments_type = 1;
				break;
			}
			case 6: {
				umat_function = &simcoon::umat_plasticity_kin_iso_CCP;
				arguments_type = 1;
				break;
			}
			case 7: {
				umat_function = &simcoon::umat_plasticity_chaboche_CCP;
				arguments_type = 1;
				break;
			}
			case 8: {
				umat_function = &simcoon::umat_plasticity_hill_isoh_CCP;
				arguments_type = 1;
				break;
			}
			case 9: {
				umat_function = &simcoon::umat_hill_chaboche_CCP;
				arguments_type = 1;
				break;
			}
			case 10: {
				umat_function = &simcoon::umat_ani_chaboche_CCP;
				arguments_type = 1;
				break;
			}
			case 11: {
				umat_function = &simcoon::umat_dfa_chaboche_CCP;
				arguments_type = 1;
				break;
			}
			case 12: {
				umat_function = &simcoon::umat_plasticity_hill_isoh_CCP_N;
				arguments_type = 1;
				break;
			}
			case 13: {
				umat_function = &simcoon::umat_sma_unified_T;
				arguments_type = 1;
				break;
			}
			case 15: {
				umat_function = &simcoon::umat_damage_LLD_0;
				arguments_type = 1;
				break;
			}
			case 16: {
				umat_function = &simcoon::umat_zener_fast;
				arguments_type = 1;
				break;
			}
			case 17: {
				umat_function = &simcoon::umat_zener_Nfast;
				arguments_type = 1;
				break;
			}
			case 18: {
				umat_function = &simcoon::umat_prony_Nfast;
				arguments_type = 1;
				break;
			}
			case 19: case 20: {
				umat_function = &simcoon::umat_sma_mono;
				arguments_type = 1;
				break;
			}
			case 21: case 22: case 23: case 24: case 26: {
				F0 = carma::arr_to_cube_view(F0_py);
				F1 = carma::arr_to_cube_view(F1_py);
				umat_function_finite = &simcoon::umat_generic_hyper_invariants;
				arguments_type = 2;
				break;
			}
			case 27: {
				umat_function = &simcoon::umat_generic_chaboche_CCP;
				arguments_type = 1;
				break;
			}				
			default: {
				throw std::invalid_argument( "The choice of Umat could not be found in the umat library." );
			}
		}

		#ifdef _OPENMP
			int max_threads = omp_get_max_threads();
			omp_set_num_threads(n_threads);
			omp_set_active_levels(3);
			#pragma omp parallel for shared(Lt, L, DR)
		#endif
		for (int pt = 0; pt < nb_points; pt++) {

			//if (use_temp) T = list_T(pt);
			if (unique_props == false) {
				props = list_props.col(pt); //if list_props has only one element, we keep only this one (assuming homogeneous material)		
			} 	
			vec statev = list_statev.unsafe_col(pt);
			vec sigma = list_sigma.unsafe_col(pt); 

			vec etot = list_etot.unsafe_col(pt);
			vec Detot = list_Detot.unsafe_col(pt);
			vec Wm = list_Wm.unsafe_col(pt);				
			
			double T = 0.0, DT = 0.0;
			if (use_temp && pt < vec_T.n_elem) {
				T = vec_T(pt);
			}

			switch (arguments_type) {

				case 1: {
					umat_function(umat_name_py, etot, Detot, sigma, Lt.slice(pt), L.slice(pt), DR.slice(pt), nprops, props, nstatev, statev, T, DT, Time, DTime, Wm(0), Wm(1), Wm(2), Wm(3), ndi, nshr, start, tnew_dt);
					break;
				}
				case 2: {
					umat_function_finite(umat_name_py, etot, Detot, F0.slice(pt), F1.slice(pt), sigma, Lt.slice(pt), L.slice(pt), DR.slice(pt), nprops, props, nstatev, statev, T, DT, Time, DTime, Wm(0), Wm(1), Wm(2), Wm(3), ndi, nshr, start, tnew_dt);
					break;
				}
			}
		}
		#ifdef _OPENMP
		omp_set_num_threads(max_threads);		
		#endif		
		return py::make_tuple(carma::mat_to_arr(list_sigma, false), carma::mat_to_arr(list_statev, false), carma::mat_to_arr(list_Wm, false), carma::cube_to_arr(Lt, false));

	}
}

/* py::tuple launch_umat_T(const std::string& umat_name_py, const py::array_t<double> &etot_py, const py::array_t<double> &Detot_py, const py::array_t<double> &sigma_py, const py::array_t<double> &DR_py, const py::array_t<double> &props_py, const py::array_t<double> &statev_py, const float Time, const float DTime, const py::array_t<double> &Wm_py, py::array_t<double> &T){
		//Get the id of umat
		std::map<string, int> list_umat;
		list_umat = { {"UMEXT",0},{"UMABA",1},{"ELISO",2},{"ELIST",3},{"ELORT",4},{"EPICP",5},{"EPKCP",6},{"EPCHA",7},{"SMADI",8},{"SMADC",8},{"SMAAI",8},{"SMAAC",8},{"LLDM0",9},{"ZENER",10},{"ZENNK",11},{"PRONK",12},{"EPHIC",17},{"EPHIN",18},{"SMAMO",19},{"SMAMC",20},{"MIHEN",100},{"MIMTN",101},{"MISCN",103},{"MIPLN",104} };
		int id_umat = list_umat[umat_name_py];
		int arguments_type; //depends on the argument used in the umat
				
		void (*umat_function)(const arma::vec &, const arma::vec &, arma::vec &, double &, arma::mat &, arma::mat &, arma::mat &, arma::mat &, const arma::mat &, const int &, const arma::vec &, const int &, arma::vec &, const double &, const double &,const double &,const double &, double &, double &, double &, double &, double &, double &, double &, const int &, const int &, const bool &, double &);

		switch (id_umat) {

				case 2: {
					umat_function = &simcoon::umat_elasticity_iso_T;
					arguments_type = 1;
					//umat_elasticity_iso_T(umat_T->Etot, umat_T->DEtot, umat_T->sigma, umat_T->r, umat_T->dSdE, umat_T->dSdT, umat_T->drdE, umat_T->drdT, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_T->nstatev, umat_T->statev, umat_T->T, umat_T->DT, Time, DTime, umat_T->Wm(0), umat_T->Wm(1), umat_T->Wm(2), umat_T->Wm(3), umat_T->Wt(0), umat_T->Wt(1), umat_T->Wt(2), ndi, nshr, start, tnew_dt);										
					break;
				}
				case 3: {
					umat_function = &simcoon::umat_elasticity_trans_iso_T;
					arguments_type = 4;
					break;
				}
				case 4: {
					umat_function = &simcoon::umat_elasticity_ortho_T;
					arguments_type = 4;
					break;
				}
				case 5: {
					umat_function = &simcoon::umat_plasticity_iso_CCP;
					arguments_type = 1;
					//simcoon::umat_plasticity_iso_CCP(etot, Detot, sigma, Lt, L, sigma_in, DR, nprops, props, nstatev, statev, T, DT, Time, DTime, Wm, Wm_r, Wm_ir, Wm_d, ndi, nshr, start, solver_type, tnew_dt);
					break;
				}
				case 6: {
					if (thermomechancial){					
						arguments_type = 0;
					}
					else {
						umat_function = &simcoon::umat_plasticity_kin_iso_CCP;
						arguments_type = 1;
						//simcoon::umat_plasticity_kin_iso_CCP(etot, Detot, sigma, Lt, L, sigma_in, DR, nprops, props, nstatev, statev, T, DT, Time, DTime, Wm, Wm_r, Wm_ir, Wm_d, ndi, nshr, start, solver_type, tnew_dt);
					}
					break;
				}
				case 7: {
					if (thermomechancial){					
						arguments_type = 0;
					}
					else {
						umat_function = &simcoon::umat_plasticity_chaboche_CCP;
						arguments_type = 1;
						//simcoon::umat_plasticity_chaboche_CCP(etot, Detot, sigma, Lt, L, sigma_in, DR, nprops, props, nstatev, statev, T, DT, Time, DTime, Wm, Wm_r, Wm_ir, Wm_d, ndi, nshr, start, solver_type, tnew_dt);
					}
					break;
				}				
				case 8: {		
					if (thermomechancial){					
						arguments_type = 0;
					}
					else {
						umat_function_2 = &simcoon::umat_sma_unified_T;
						arguments_type = 2;				
						//simcoon::umat_sma_unified_T(etot, Detot, sigma, Lt, DR, nprops, props, nstatev, statev, T, DT, Time, DTime, Wm, Wm_r, Wm_ir, Wm_d, ndi, nshr, start, tnew_dt);
					}
					break;
				}
				case 9: {
					if (thermomechancial){					
						arguments_type = 0;
					}
					else {
						umat_function = &simcoon::umat_damage_LLD_0;
						arguments_type = 1;
						//simcoon::umat_damage_LLD_0(etot, Detot, sigma, Lt, L, sigma_in, DR, nprops, props, nstatev, statev, T, DT, Time, DTime, Wm, Wm_r, Wm_ir, Wm_d, ndi, nshr, start, solver_type, tnew_dt);
					}
					break;
				}				
				case 10: {
					if (thermomechancial){					
						arguments_type = 0;
					}
					else {
						umat_function_2 = &simcoon::umat_zener_fast;
						arguments_type = 2;		
						//simcoon::umat_zener_fast(etot, Detot, sigma, Lt, DR, nprops, props, nstatev, statev, T, DT, Time, DTime, Wm, Wm_r, Wm_ir, Wm_d, ndi, nshr, start, tnew_dt);
					}
					break;
				}				
				case 11: {
					if (thermomechancial){					
						arguments_type = 0;
					}
					else {
						umat_function_2 = &simcoon::umat_zener_Nfast;
						arguments_type = 2;
						//simcoon::umat_zener_Nfast(etot, Detot, sigma, Lt, DR, nprops, props, nstatev, statev, T, DT, Time, DTime, Wm, Wm_r, Wm_ir, Wm_d, ndi, nshr, start, tnew_dt);
					}						
					break;
				}
				case 12: {
					if (thermomechancial){					
						arguments_type = 0;
					}
					else {
						umat_function_2 = &simcoon::umat_prony_Nfast;
						arguments_type = 2;
						//simcoon::umat_prony_Nfast(etot, Detot, sigma, Lt, DR, nprops, props, nstatev, statev, T, DT, Time, DTime, Wm, Wm_r, Wm_ir, Wm_d, ndi, nshr, start, tnew_dt);
					}
					break;
				}
				case 17: {
					if (thermomechancial){					
						arguments_type = 0;
					}
					else {						
						umat_function_2 = &simcoon::umat_plasticity_hill_isoh_CCP;
						arguments_type = 2;
						//simcoon::umat_plasticity_hill_isoh_CCP(etot, Detot, sigma, Lt, DR, nprops, props, nstatev, statev, T, DT, Time, DTime, Wm, Wm_r, Wm_ir, Wm_d, ndi, nshr, start, tnew_dt);
					}
					break;
				}
				case 18: {
					if (thermomechancial){					
						arguments_type = 0;
					}
					else {
						umat_function_2 = &simcoon::umat_plasticity_hill_isoh_CCP_N;
						arguments_type = 2;
						//simcoon::umat_plasticity_hill_isoh_CCP_N(etot, Detot, sigma, Lt, DR, nprops, props, nstatev, statev, T, DT, Time, DTime, Wm, Wm_r, Wm_ir, Wm_d, ndi, nshr, start, tnew_dt);
					}
					break;
				}
				case 19: {
					if (thermomechancial){					
						arguments_type = 0;
					}
					else {
						umat_function_3 = &simcoon::umat_sma_mono;
						arguments_type = 3;
						//simcoon::umat_sma_mono(etot, Detot, sigma, Lt, L, DR, nprops, props, nstatev, statev, T, DT, Time, DTime, Wm, Wm_r, Wm_ir, Wm_d, ndi, nshr, start, tnew_dt);
					}
					break;
				}
				case 20: {
					if (thermomechancial){					
						arguments_type = 0;
					}
					else {
						umat_function_3 = &simcoon::umat_sma_mono_cubic;
						arguments_type = 3;
						//simcoon::umat_sma_mono_cubic(etot, Detot, sigma, Lt, L, DR, nprops, props, nstatev, statev, T, DT, Time, DTime, Wm, Wm_r, Wm_ir, Wm_d, ndi, nshr, start, tnew_dt);
					}
					break;
				}

				default: {
					//py::print("Error: The choice of Umat could not be found in the umat library \n");
					throw std::invalid_argument( "The choice of Umat could not be found in the umat library." );
					//exit(0);
				}
			}
		
		if (arguments_type == 0) throw std::invalid_argument( "The choice of Umat could not be found in the umat library. Check the thermomecanical argument." );

		//py::print("id_umat ok");

		//scalar needed to launch umat
		const int solver_type = 0;
		const int ndi=3;
		const int nshr=3;
		const int ncomp=ndi+nshr;
		const bool start = false;
		double tnew_dt = 0;//usefull ?		
		double T = 0; double DT = 0;  //modify to let the program set the actual temperature
		//bool use_temp;
		//if (T.n_elem == 0.) use_temp = false; 
		//else use_temp = true;

		//py::print("ndim = ", etot_py.ndim());
		if (etot_py.ndim() == 1) {			
			vec etot = carma::arr_to_col_view(etot_py);
			vec Detot = carma::arr_to_col_view(Detot_py); 
			vec sigma = carma::arr_to_col(sigma_py); //copy data because values are changed by the umat and returned to python
			mat DR = carma::arr_to_mat_view(DR_py); 
			vec props = carma::arr_to_col_view(props_py);
			vec statev = carma::arr_to_col(statev_py); //copy data because values are changed by the umat and returned to python
			vec Wm = carma::arr_to_col(Wm_py); //copy data because values are changed by the umat and returned to python
			mat L(ncomp, ncomp);
			mat Lt(ncomp, ncomp);
			vec sigma_in = zeros(1); //not used
			int nprops = props.n_elem;
			int nstatev = statev.n_elem;
			
			switch (arguments_type) {

				case 1: {
					umat_function(etot, Detot, sigma, Lt, L, sigma_in, DR, nprops, props, nstatev, statev, T, DT, Time, DTime, Wm(0), Wm(1), Wm(2), Wm(3), ndi, nshr, start, solver_type, tnew_dt);
					break;
				}
				case 2: {
					umat_function_2(etot, Detot, sigma, Lt, DR, nprops, props, nstatev, statev, T, DT, Time, DTime, Wm(0), Wm(1), Wm(2), Wm(3), ndi, nshr, start, tnew_dt);
					break;
				}
				case 3: {
					umat_function_3(etot, Detot, sigma, Lt, L, DR, nprops, props, nstatev, statev, T, DT, Time, DTime, Wm(0), Wm(1), Wm(2), Wm(3), ndi, nshr, start, tnew_dt);
					break;
				}
				case 4: {
					umat_elasticity_iso_T(etot, Detot, sigma, umat_T->r, umat_T->dSdE, umat_T->dSdT, umat_T->drdE, umat_T->drdT, DR, nprops, props, nstatev, statev, T, DT, Time, DTime, Wm(0), Wm(1), Wm(2), Wm(3), umat_T->Wt(0), umat_T->Wt(1), umat_T->Wt(2), ndi, nshr, start, tnew_dt);
					break;
				}
			}
			//py::print("umat_done");
			return py::make_tuple(carma::col_to_arr(sigma, false), carma::col_to_arr(statev, false), carma::col_to_arr(Wm, false), carma::mat_to_arr(Lt, false));
		}
		else if (etot_py.ndim() == 2) {
			
			mat etot = carma::arr_to_mat_view(etot_py);
			int nb_points = etot.n_cols; //number of material points
			mat Detot = carma::arr_to_mat_view(Detot_py); 
			mat list_sigma = carma::arr_to_mat(sigma_py); //copy data because values are changed by the umat and returned to python
			cube DR = carma::arr_to_cube_view(DR_py); 
			mat list_props = carma::arr_to_mat_view(props_py);
			mat list_statev = carma::arr_to_mat(statev_py); //copy data because values are changed by the umat and returned to python
			mat Wm = carma::arr_to_mat(Wm_py); //copy data because values are changed by the umat and returned to python
			cube L(ncomp, ncomp, nb_points);
			cube Lt(ncomp, ncomp, nb_points);
			vec sigma_in = zeros(1); //not used
			int nprops = list_props.n_rows;
			int nstatev = list_statev.n_rows;

			vec props(ncomp);
	
			for (int pt = 0; pt < nb_points; pt++) {
				//if (use_temp) T = list_T(pt);
				if (pt < list_props.n_cols) props = list_props.col(pt); //if list_props has only one element, we keep only this one (assuming homogeneous material)			
				vec statev = list_statev.unsafe_col(pt);
				vec sigma = list_sigma.unsafe_col(pt); 

				switch (arguments_type) {

					case 1: {
						umat_function(etot.col(pt), Detot.col(pt), sigma, Lt.slice(pt), L.slice(pt), sigma_in, DR.slice(pt), nprops, props, nstatev, statev, T, DT, Time, DTime, Wm(0,pt), Wm(1,pt), Wm(2,pt), Wm(3,pt), ndi, nshr, start, solver_type, tnew_dt);
						break;
					}
					case 2: {
						umat_function_2(etot.col(pt), Detot.col(pt), sigma, Lt.slice(pt), DR.slice(pt), nprops, props, nstatev, statev, T, DT, Time, DTime, Wm(0,pt), Wm(1,pt), Wm(2,pt), Wm(3,pt), ndi, nshr, start, tnew_dt);
						break;
					}
					case 3: {
						umat_function_3(etot.col(pt), Detot.col(pt), sigma, Lt.slice(pt), L.slice(pt), DR.slice(pt), nprops, props, nstatev, statev, T, DT, Time, DTime, Wm(0,pt), Wm(1,pt), Wm(2,pt), Wm(3,pt), ndi, nshr, start, tnew_dt);
						break;
					}
				}
			}
			return py::make_tuple(carma::mat_to_arr(list_sigma, false), carma::mat_to_arr(list_statev, false), carma::mat_to_arr(Wm, false), carma::cube_to_arr(Lt, false));
		}
	}
}

*/
