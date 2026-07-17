#include <pybind11/embed.h>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <optional>

#include <carma>
#include <armadillo>

#include <simcoon/parameter.hpp>
#include <simcoon/parallel.hpp>

#include <simcoon/python_wrappers/Libraries/Continuum_mechanics/umat.hpp>

#include <simcoon/Continuum_mechanics/Umat/Mechanical/External/external_umat.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Plasticity/plastic_isotropic_ccp.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Plasticity/plastic_chaboche_ccp.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/SMA/unified_T.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/SMA/unified_TR.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/SMA/SMA_mono.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Damage/damage_LLD_0.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Viscoelasticity/Zener_fast.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Viscoelasticity/Zener_Nfast.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Viscoelasticity/Prony_Nfast.hpp>

#include <simcoon/Continuum_mechanics/Umat/Finite/generic_hyper_invariants.hpp>
#include <simcoon/Continuum_mechanics/Umat/Finite/saint_venant.hpp>
#include <simcoon/Continuum_mechanics/Umat/Finite/neo_hookean_incomp.hpp>

#include <simcoon/Continuum_mechanics/Umat/Modular/modular_umat.hpp>
#include <simcoon/Continuum_mechanics/Umat/Modular/legacy_adapters.hpp>

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
	
	py::tuple launch_umat(const std::string &umat_name_py, const py::array_t<double> &etot_py, const py::array_t<double> &Detot_py, const py::array_t<double> &F0_py, const py::array_t<double> &F1_py, const py::array_t<double> &sigma_py, const py::array_t<double> &DR_py, const py::array_t<double> &props_py, const py::array_t<double> &statev_py, const float Time, const float DTime, const py::array_t<double> &Wm_py, const std::optional<py::array_t<double>> &T_py, const int &ndi, const unsigned int &n_threads, const int &tangent_mode){
		// tangent_mode: 0 = none (explicit integration, Lt = elastic L),
		//               1 = continuum, 2 = algorithmic (Simo-Hughes, DEFAULT),
		//               3 = closest-point (reserved). See parameter.hpp tangent_* constants.

		// Validate up front, in serial context: the per-point dispatch below
		// runs inside a non-exception-safe parallel region (GCD/OpenMP) where a
		// throw would std::terminate the host process.
		if (tangent_mode < simcoon::tangent_none || tangent_mode > simcoon::tangent_algorithmic) {
			throw std::invalid_argument("tangent_mode must be 0 (none), 1 (continuum) or 2 (algorithmic); got "
			                            + std::to_string(tangent_mode) + " (3 = closest-point is reserved)");
		}
		std::map<string, int> list_umat;
		list_umat = { {"UMEXT",0},{"UMABA",1},{"ELISO",2},{"ELIST",3},{"ELORT",4},{"EPICP",5},{"EPKCP",6},{"EPCHA",7},{"EPHIL",8},{"EPHAC",9},{"EPANI",10},{"EPDFA",11},{"EPHIN",12},{"SMAUT",13},{"SMANI",13},{"SMADI",13},{"SMADC",13},{"SMAAI",13},{"SMAAC",13},{"LLDM0",15},{"ZENER",16},{"ZENNK",17},{"PRONK",18},{"SMAMO",19},{"SMAMC",20},{"NEOHC",21},{"MOORI",22},{"YEOHH",23},{"ISHAH",24},{"GETHH",25},{"SWANH",26},{"EPCHG",27},{"SMRDI",28},{"SMRDC",28},{"SMRAI",28},{"SMRAC",28},{"SNTVE",29},{"NEOHI",30},{"MODUL",200},{"MIHEN",100},{"MIMTN",101},{"MISCN",103},{"MIPLN",104} }; // TODO_2.0 SMAUT and SMANI compatibility to be removed in release 2.0
		int id_umat = list_umat[umat_name_py];
		int arguments_type; //depends on the argument used in the umat

		// Unified small-strain function pointer: (umat_name, Etot, DEtot, sigma, Lt, L, DR, nprops, props, nstatev, statev, T, DT, Time, DTime, Wm, Wm_r, Wm_ir, Wm_d, ndi, nshr, start, tnew_dt, tangent_mode)
		void (*umat_function)(const std::string &, const arma::vec &, const arma::vec &, arma::vec &, arma::mat &, arma::mat &, const arma::mat &, const int &, const arma::vec &, const int &, arma::vec &, const double &, const double &, const double &, const double &, double &, double &, double &, double &, const int &, const int &, const bool &, double &, const int &);
		// Unified finite-strain function pointer: (umat_name, etot, Detot, F0, F1, sigma, Lt, L, DR, nprops, props, nstatev, statev, T, DT, Time, DTime, Wm, Wm_r, Wm_ir, Wm_d, ndi, nshr, start, tnew_dt, tangent_mode)
		void (*umat_function_finite)(const std::string &, const arma::vec &, const arma::vec &, const arma::mat &, const arma::mat &, arma::vec &, arma::mat &, arma::mat &, const arma::mat &, const int &, const arma::vec &, const int &, arma::vec &, const double &, const double &, const double &, const double &, double &, double &, double &, double &, const int &, const int &, const bool &, double &, const int &);
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
				umat_function = &simcoon::umat_legacy_modular; // legacy name -> modular adapter
				arguments_type = 1;
				break;
			}
			case 3: {
				umat_function = &simcoon::umat_legacy_modular; // legacy name -> modular adapter
				arguments_type = 1;
				break;
			}
			case 4: {
				umat_function = &simcoon::umat_legacy_modular; // legacy name -> modular adapter
				arguments_type = 1;
				break;
			}
			case 5: {
				umat_function = &simcoon::umat_plasticity_iso_CCP;
				arguments_type = 1;
				break;
			}
			case 6: {
				umat_function = &simcoon::umat_legacy_modular; // legacy name -> modular adapter
				arguments_type = 1;
				break;
			}
			case 7: {
				umat_function = &simcoon::umat_plasticity_chaboche_CCP;
				arguments_type = 1;
				break;
			}
			case 8: {
				umat_function = &simcoon::umat_legacy_modular; // legacy name -> modular adapter
				arguments_type = 1;
				break;
			}
			case 9: {
				umat_function = &simcoon::umat_legacy_modular; // legacy name -> modular adapter
				arguments_type = 1;
				break;
			}
			case 10: {
				umat_function = &simcoon::umat_legacy_modular; // legacy name -> modular adapter
				arguments_type = 1;
				break;
			}
			case 11: {
				umat_function = &simcoon::umat_legacy_modular; // legacy name -> modular adapter
				arguments_type = 1;
				break;
			}
			case 12: {
				umat_function = &simcoon::umat_legacy_modular; // legacy name -> modular adapter
				arguments_type = 1;
				break;
			}
			case 13: {
				umat_function = &simcoon::umat_sma_unified_T;
				arguments_type = 1;
				break;
			}
			case 28: { // SMA_TR (unified reduced-tangent model)
				umat_function = &simcoon::umat_sma_unified_TR;
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
			case 200: { // MODUL (modular UMAT, small-strain)
				umat_function = &simcoon::umat_modular;
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
				umat_function = &simcoon::umat_legacy_modular; // legacy name -> modular adapter
				arguments_type = 1;
				break;
			}
			case 29: { // SNTVE (Saint-Venant-Kirchhoff, finite)
				F0 = carma::arr_to_cube_view(F0_py);
				F1 = carma::arr_to_cube_view(F1_py);
				umat_function_finite = &simcoon::umat_saint_venant;
				arguments_type = 2;
				break;
			}
			case 30: { // NEOHI (Neo-Hookean incompressible, finite)
				F0 = carma::arr_to_cube_view(F0_py);
				F1 = carma::arr_to_cube_view(F1_py);
				umat_function_finite = &simcoon::umat_neo_hookean_incomp;
				arguments_type = 2;
				break;
			}
			default: {
				throw std::invalid_argument( "The choice of Umat could not be found in the umat library." );
			}
		}

		simcoon_parallel_for(nb_points, [&](int pt) {
			// Alias the props column without copying so the parallel region makes no
			// NumPy-backed (carma) allocation: GCD/OpenMP workers then never call
			// PyDataMem_NEW (which needs the GIL) -> no GIL deadlock, no GIL handling.
			// props (unique) / list_props (per-point) outlive the lambda and are read-only.
			const double* _props_ptr = unique_props ? props.memptr() : list_props.colptr(pt);
			const vec local_props(const_cast<double*>(_props_ptr), nprops, false, true);
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
					umat_function(umat_name_py, etot, Detot, sigma, Lt.slice(pt), L.slice(pt), DR.slice(pt), nprops, local_props, nstatev, statev, T, DT, Time, DTime, Wm(0), Wm(1), Wm(2), Wm(3), ndi, nshr, start, tnew_dt, tangent_mode);
					break;
				}
				case 2: {
					umat_function_finite(umat_name_py, etot, Detot, F0.slice(pt), F1.slice(pt), sigma, Lt.slice(pt), L.slice(pt), DR.slice(pt), nprops, local_props, nstatev, statev, T, DT, Time, DTime, Wm(0), Wm(1), Wm(2), Wm(3), ndi, nshr, start, tnew_dt, tangent_mode);
					break;
				}
			}
		});
		return py::make_tuple(carma::mat_to_arr(list_sigma, false), carma::mat_to_arr(list_statev, false), carma::mat_to_arr(list_Wm, false), carma::cube_to_arr(Lt, false));

	}

	py::tuple launch_umat_T(const std::string &umat_name_py, const py::array_t<double> &etot_py, const py::array_t<double> &Detot_py, const py::array_t<double> &sigma_py, const py::array_t<double> &DR_py, const py::array_t<double> &props_py, const py::array_t<double> &statev_py, const float Time, const float DTime, const py::array_t<double> &Wm_py, const py::array_t<double> &Wt_py, const py::array_t<double> &T_py, const py::array_t<double> &DT_py, const int &ndi, const unsigned int &n_threads, const int &tangent_mode){
		// Point-wise thermomechanical UMAT batch entry (small strain), mirroring launch_umat.
		// Dispatch follows the select_umat_T table (umat_smart.cpp).
		// Returns (sigma, statev, Wm, Wt, r, dSdE, dSdT, drdE, drdT).

		if (tangent_mode < simcoon::tangent_none || tangent_mode > simcoon::tangent_algorithmic) {
			throw std::invalid_argument("tangent_mode must be 0 (none), 1 (continuum) or 2 (algorithmic); got "
			                            + std::to_string(tangent_mode) + " (3 = closest-point is reserved)");
		}
		std::map<string, int> list_umat;
		list_umat = { {"ELISO",1},{"ELIST",2},{"ELORT",3},{"EPICP",4},{"EPKCP",5},{"ZENER",6},{"ZENNK",7},{"PRONK",8},{"SMAUT",9},{"SMANI",9},{"SMADI",9},{"SMADC",9},{"SMAAI",9},{"SMAAC",9} }; // TODO_2.0 SMAUT and SMANI compatibility to be removed in release 2.0
		if (list_umat.count(umat_name_py) == 0) {
			throw std::invalid_argument("The choice of thermomechanical Umat could not be found in the umat library: " + umat_name_py);
		}
		int id_umat = list_umat[umat_name_py];
		int arguments_type; //depends on the argument used in the umat

		// Unified thermomechanical function pointer: (Etot, DEtot, sigma, r, dSdE, dSdT, drdE, drdT, DR, nprops, props, nstatev, statev, T, DT, Time, DTime, Wm, Wm_r, Wm_ir, Wm_d, Wt, Wt_r, Wt_ir, ndi, nshr, start, tnew_dt, tangent_mode)
		void (*umat_function)(const arma::vec &, const arma::vec &, arma::vec &, double &, arma::mat &, arma::mat &, arma::mat &, arma::mat &, const arma::mat &, const int &, const arma::vec &, const int &, arma::vec &, const double &, const double &, const double &, const double &, double &, double &, double &, double &, double &, double &, double &, const int &, const int &, const bool &, double &, const int &);
		// SMA family variant carrying the umat_name as leading argument
		void (*umat_function_named)(const std::string &, const arma::vec &, const arma::vec &, arma::vec &, double &, arma::mat &, arma::mat &, arma::mat &, arma::mat &, const arma::mat &, const int &, const arma::vec &, const int &, arma::vec &, const double &, const double &, const double &, const double &, double &, double &, double &, double &, double &, double &, double &, const int &, const int &, const bool &, double &, const int &);

		const int ncomp = 6;
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
		double tnew_dt = 0;

		mat list_etot = carma::arr_to_mat_view(etot_py);
		unsigned int nb_points = list_etot.n_cols; //number of material points
		mat list_Detot = carma::arr_to_mat_view(Detot_py);
		mat list_sigma = carma::arr_to_mat(std::move(sigma_py)); //copy: modified by the umat and returned
		cube DR = carma::arr_to_cube_view(DR_py);
		vec vec_T = carma::arr_to_col_view(T_py);
		vec vec_DT = carma::arr_to_col_view(DT_py);

		vec props;
		//n_cols (not the raw numpy shape) so a 1-D (nprops,) array is a valid single-props input
		mat list_props = carma::arr_to_mat_view(props_py);
		bool unique_props = false;
		if (list_props.n_cols == 1) {
			props = list_props.col(0);
			unique_props = true;
		}
		else if (list_props.n_cols != nb_points) {
			throw std::invalid_argument("umat_T: props must have 1 column (shared) or one column per material point; got "
			                            + std::to_string(list_props.n_cols) + " columns for " + std::to_string(nb_points) + " points");
		}

		mat list_statev = carma::arr_to_mat(std::move(statev_py)); //copy: modified by the umat and returned
		mat list_Wm = carma::arr_to_mat(std::move(Wm_py)); //copy: modified by the umat and returned
		mat list_Wt = carma::arr_to_mat(std::move(Wt_py)); //copy: modified by the umat and returned

		//Validate every batch dimension here, in serial context: an out-of-range
		//access inside the non-exception-safe parallel region would terminate the process
		if (list_Detot.n_cols != nb_points || list_sigma.n_cols != nb_points || list_statev.n_cols != nb_points
		    || list_Wm.n_cols != nb_points || list_Wt.n_cols != nb_points || DR.n_slices != nb_points) {
			throw std::invalid_argument("umat_T: Detot, sigma, statev, Wm, Wt and DR must have one column (resp. slice) per material point (" + std::to_string(nb_points) + ")");
		}
		if (vec_T.n_elem != nb_points || vec_DT.n_elem != nb_points) {
			throw std::invalid_argument("umat_T: T and DT must have one entry per material point (" + std::to_string(nb_points) + ")");
		}
		if (list_etot.n_rows != 6 || list_Detot.n_rows != 6 || list_sigma.n_rows != 6 || list_Wm.n_rows != 4 || list_Wt.n_rows != 3) {
			throw std::invalid_argument("umat_T: expected shapes (6,N) for etot/Detot/sigma, (4,N) for Wm and (3,N) for Wt");
		}
		vec list_r(nb_points, fill::zeros);
		cube dSdE(ncomp, ncomp, nb_points);
		cube dSdT(ncomp, 1, nb_points);
		cube drdE(ncomp, 1, nb_points); //T UMATs write drdE as a (6,1) column (e.g. drdE = zeros(6))
		cube drdT(1, 1, nb_points);
		int nprops = list_props.n_rows;
		int nstatev = list_statev.n_rows;

		switch (id_umat) {
			case 1: {
				umat_function = &simcoon::umat_elasticity_iso_T;
				arguments_type = 1;
				break;
			}
			case 2: {
				umat_function = &simcoon::umat_elasticity_trans_iso_T;
				arguments_type = 1;
				break;
			}
			case 3: {
				umat_function = &simcoon::umat_elasticity_ortho_T;
				arguments_type = 1;
				break;
			}
			case 4: {
				umat_function = &simcoon::umat_plasticity_iso_CCP_T;
				arguments_type = 1;
				break;
			}
			case 5: {
				umat_function = &simcoon::umat_plasticity_kin_iso_CCP_T;
				arguments_type = 1;
				break;
			}
			case 6: {
				umat_function = &simcoon::umat_zener_fast_T;
				arguments_type = 1;
				break;
			}
			case 7: {
				umat_function = &simcoon::umat_zener_Nfast_T;
				arguments_type = 1;
				break;
			}
			case 8: {
				umat_function = &simcoon::umat_prony_Nfast_T;
				arguments_type = 1;
				break;
			}
			case 9: {
				umat_function_named = &simcoon::umat_sma_unified_T_T;
				arguments_type = 2;
				break;
			}
			default: {
				throw std::invalid_argument( "The choice of thermomechanical Umat could not be found in the umat library." );
			}
		}

		simcoon_parallel_for(nb_points, [&](int pt) {
			// props aliased without copying: no NumPy-backed allocation in the
			// parallel region (same GIL-safety pattern as launch_umat)
			const double* _props_ptr = unique_props ? props.memptr() : list_props.colptr(pt);
			const vec local_props(const_cast<double*>(_props_ptr), nprops, false, true);
			vec statev = list_statev.unsafe_col(pt);
			vec sigma = list_sigma.unsafe_col(pt);

			vec etot = list_etot.unsafe_col(pt);
			vec Detot = list_Detot.unsafe_col(pt);
			vec Wm = list_Wm.unsafe_col(pt);
			vec Wt = list_Wt.unsafe_col(pt);

			double T = vec_T(pt);
			double DT = vec_DT(pt);
			double tnew_dt_pt = tnew_dt;

			switch (arguments_type) {
				case 1: {
					umat_function(etot, Detot, sigma, list_r(pt), dSdE.slice(pt), dSdT.slice(pt), drdE.slice(pt), drdT.slice(pt), DR.slice(pt), nprops, local_props, nstatev, statev, T, DT, Time, DTime, Wm(0), Wm(1), Wm(2), Wm(3), Wt(0), Wt(1), Wt(2), ndi, nshr, start, tnew_dt_pt, tangent_mode);
					break;
				}
				case 2: {
					umat_function_named(umat_name_py, etot, Detot, sigma, list_r(pt), dSdE.slice(pt), dSdT.slice(pt), drdE.slice(pt), drdT.slice(pt), DR.slice(pt), nprops, local_props, nstatev, statev, T, DT, Time, DTime, Wm(0), Wm(1), Wm(2), Wm(3), Wt(0), Wt(1), Wt(2), ndi, nshr, start, tnew_dt_pt, tangent_mode);
					break;
				}
			}
		});

		// post-loop repacking (serial): dSdT (6,1,N) -> (6,N), drdE (1,6,N) -> (6,N), drdT (1,1,N) -> (N)
		mat dSdT_out(ncomp, nb_points);
		mat drdE_out(ncomp, nb_points);
		vec drdT_out(nb_points);
		for (unsigned int pt = 0; pt < nb_points; pt++) {
			dSdT_out.col(pt) = dSdT.slice(pt);
			drdE_out.col(pt) = drdE.slice(pt);
			drdT_out(pt) = drdT(0, 0, pt);
		}
		// copy=true throughout: with few points these arrays fit armadillo's internal
		// (pre-allocated) buffer and a zero-copy steal would hand numpy a dangling pointer
		return py::make_tuple(carma::mat_to_arr(list_sigma, true), carma::mat_to_arr(list_statev, true), carma::mat_to_arr(list_Wm, true), carma::mat_to_arr(list_Wt, true), carma::col_to_arr(list_r, true), carma::cube_to_arr(dSdE, true), carma::mat_to_arr(dSdT_out, true), carma::mat_to_arr(drdE_out, true), carma::col_to_arr(drdT_out, true));
	}
}
