
#include <armadillo>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <simcoon/arma2numpy/numpy_arma.hpp>

#include <simcoon/python_wrappers/Libraries/Continuum_mechanics/RunUmat.hpp>

#include <simcoon/Continuum_mechanics/Umat/Mechanical/External/external_umat.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Elasticity/elastic_isotropic.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Elasticity/elastic_transverse_isotropic.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Elasticity/elastic_orthotropic.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Plasticity/plastic_isotropic_ccp.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Plasticity/plastic_kin_iso_ccp.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Plasticity/plastic_chaboche_ccp.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Plasticity/Hill_isoh.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Plasticity/Hill_isoh_Nfast.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/SMA/unified_T.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/SMA/SMA_mono.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/SMA/SMA_mono_cubic.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Damage/damage_LLD_0.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Viscoelasticity/Zener_fast.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Viscoelasticity/Zener_Nfast.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Viscoelasticity/Prony_Nfast.hpp>
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

#include <simcoon/Continuum_mechanics/Functions/objective_rates.hpp>
#include <simcoon/Continuum_mechanics/Functions/transfer.hpp>
#include <simcoon/Continuum_mechanics/Functions/kinematics.hpp>
#include <simcoon/python_wrappers/Libraries/Continuum_mechanics/kinematics.hpp>
#include <simcoon/Continuum_mechanics/Functions/stress.hpp>

#include <simcoon/Continuum_mechanics/Umat/Thermomechanical/SMA/unified_T.hpp>

namespace bn = boost::python::numpy;
namespace bp = boost::python;
using namespace std;
using namespace arma;
using namespace arma2numpy;

namespace simpy {

	bp::tuple RunUmat_fedoo(const std::string& umat_name_py, bn::ndarray& etot_py, const bn::ndarray& F0_py, const bn::ndarray& F1_py, bn::ndarray& sigma_py, bn::ndarray& Lt_py, bn::ndarray& L_py, const bn::ndarray& props_py, bn::ndarray& statev_py, const double T, const double DT, const double Time, const double DTime, const int ndi, const int corate) {
		
		mat list_etot = array2mat_inplace(etot_py); //total strain cumulated
		cube listF0 = array2cube_inplace(F0_py);
		cube listF1 = array2cube_inplace(F1_py);
		cube list_Lt = array2cube_inplace(Lt_py);
		cube list_L = array2cube_inplace(L_py);

		const int solver_type = 0;
		const bool start = false;

		double tnew_dt = 0;//usefull ?
		//put Wm, Wm_r, Wm_ir, Wm_d in an array
		double Wm = 0.;
		double Wm_r = 0.;
		double Wm_ir = 0.;
		double Wm_d = 0.;
		const int& nshr = ndi;

		// convert array to armdillo vec
		const mat list_props = array2mat_inplace(props_py); //each col is prop for one pg
		const int nprops = list_props.n_rows;
		mat list_statev = array2mat_inplace(statev_py); //should be mat
		const int nstatev = list_statev.n_rows;

		mat list_sigma = array2mat_inplace(sigma_py); ////each col is sigma for one pg
		vec sigma = zeros(ndi + nshr);
		vec sigma_in = zeros(1); //not used
		
		vec statev; //vec sigma; //mat Lt; mat L; 
		mat F0; mat F1; vec props;
		//vec Etot; vec DEtot; //Lagrangian strain
		vec etot; vec Detot; //Eulerian strain
		mat DR = eye(3, 3);
		mat D = zeros(3, 3);
		mat Omega = zeros(3, 3);

		std::map<string, int> list_umat;
		list_umat = { {"UMEXT",0},{"UMABA",1},{"ELISO",2},{"ELIST",3},{"ELORT",4},{"EPICP",5},{"EPKCP",6},{"EPCHA",7},{"SMAUT",8},{"LLDM0",9},{"ZENER",10},{"ZENNK",11},{"PRONK",12},{"EPHIC",17},{"EPHIN",18},{"SMAMO",19},{"SMAMC",20},{"MIHEN",100},{"MIMTN",101},{"MISCN",103},{"MIPLN",104} };

		for (int pg = 0; pg < listF1.n_slices; pg++) {
			if (pg < list_props.n_cols) props = list_props.col(pg); //if list_props has only one element, we keep only this one			
			statev = list_statev.col(pg);			
			sigma = list_sigma.col(pg);			

			if (pg < listF0.n_slices) F0 = listF0.slice(pg); //if listF0 has only one element, we keep only this one			
			F1 = listF1.slice(pg);
			mat Lt = list_Lt.slice(pg);
			mat L = list_L.slice(pg);
			
			if (corate == -1) {
				//Green-Lagrange Strain
				etot = simcoon::t2v_strain(simcoon::Green_Lagrange(F0));
				Detot = simcoon::t2v_strain(simcoon::Green_Lagrange(F1)) - etot;
				}
			else if (corate == 0) {
				//Jaumann
				etot = list_etot.col(pg);
				simcoon::Jaumann(DR, D, Omega, DTime, F0, F1); //to compute D, W, Omega
				Detot = simcoon::t2v_strain(simcoon::Delta_log_strain(D, Omega, DTime)); //etot : VR decomposition, then ln(V) equals the logarithmic strain			
				}	
			else if (corate == 1) {
				//Green Naghdi
				etot = list_etot.col(pg);
				simcoon::Green_Naghdi(DR, D, Omega, DTime, F0, F1); //to compute D, W, Omega
				Detot = simcoon::t2v_strain(simcoon::Delta_log_strain(D, Omega, DTime)); //etot : VR decomposition, then ln(V) equals the logarithmic strain			
				}
			else if (corate == 2) {
				//Log Strain with etot cumulated
				etot = list_etot.col(pg);
				simcoon::logarithmic(DR, D, Omega, DTime, F0, F1); //to compute D, W, Omega
				Detot = simcoon::t2v_strain(simcoon::Delta_log_strain(D, Omega, DTime)); //etot : VR decomposition, then ln(V) equals the logarithmic strain			
				}
			else if (corate == 3) {
				//Log Strain with etot = ln(V)
				etot = simcoon::t2v_strain(0.5 * logmat_sympd(simcoon::L_Cauchy_Green(F0))); //etot : VR decomposition, then ln(V) equals the logarithmic strain
				simcoon::logarithmic(DR, D, Omega, DTime, F0, F1); //to compute D, W, Omega
				Detot = simcoon::t2v_strain(simcoon::Delta_log_strain(D, Omega, DTime)); //etot : VR decomposition, then ln(V) equals the logarithmic strain			
				}		

			//Launch the simcoon umat
			switch (list_umat[umat_name_py]) {

				/*case 0: {
					//                umat_external(umat_M->etot, umat_M->Detot, umat_M->sigma, umat_M->Lt, umat_M->L, umat_M->sigma_in, DR, rve.sptr_matprops->nprops, rve.sptr_matprops->props, umat_M->nstatev, umat_M->statev, umat_M->T, umat_M->DT, Time, DTime, umat_M->Wm(0), umat_M->Wm(1), umat_M->Wm(2), umat_M->Wm(3), ndi, nshr, start, solver_type, tnew_dt);

					boost::filesystem::path lib_path("external");          // argv[1] contains path to directory with our plugin library
					boost::shared_ptr<umat_plugin_ext_api> external_umat;            // variable to hold a pointer to plugin variable
					external_umat = dll::import<umat_plugin_ext_api>(          // type of imported symbol is located between `<` and `>`
						lib_path / "umat_plugin_ext",                     // path to the library and library name
						"external_umat",                                       // name of the symbol to import
						dll::load_mode::append_decorations              // makes `libmy_plugin_sum.so` or `my_plugin_sum.dll` from `my_plugin_sum`
						);

					external_umat->umat_external_M(etot, Detot, sigma, Lt, L, sigma_in, DR, nprops, props, nstatev, statev, T, DT, Time, DTime, Wm, Wm_r, Wm_ir, Wm_d, ndi, nshr, start, solver_type, tnew_dt);

					break;
				}*/
				/*case 1: {
					//
					boost::filesystem::path lib_path("external");
					boost::shared_ptr<umat_plugin_aba_api> abaqus_umat;            // variable to hold a pointer to plugin variable

					abaqus_umat = dll::import<umat_plugin_aba_api>(lib_path / "umat_plugin_aba", "abaqus_umat", dll::load_mode::append_decorations);
					abaqus_umat->umat_abaqus(rve, DR, Time, DTime, ndi, nshr, start, solver_type, tnew_dt);
					break;
				}*/
				case 2: {
					simcoon::umat_elasticity_iso(etot, Detot, sigma, Lt, L, sigma_in, DR, nprops, props, nstatev, statev, T, DT, Time, DTime, Wm, Wm_r, Wm_ir, Wm_d, ndi, nshr, start, solver_type, tnew_dt);
					break;
				}
				case 3: {
					simcoon::umat_elasticity_trans_iso(etot, Detot, sigma, Lt, L, sigma_in, DR, nprops, props, nstatev, statev, T, DT, Time, DTime, Wm, Wm_r, Wm_ir, Wm_d, ndi, nshr, start, solver_type, tnew_dt);
					break;
				}
				case 4: {
					simcoon::umat_elasticity_ortho(etot, Detot, sigma, Lt, L, sigma_in, DR, nprops, props, nstatev, statev, T, DT, Time, DTime, Wm, Wm_r, Wm_ir, Wm_d, ndi, nshr, start, solver_type, tnew_dt);
					break;
				}
				case 5: {
					simcoon::umat_plasticity_iso_CCP(etot, Detot, sigma, Lt, L, sigma_in, DR, nprops, props, nstatev, statev, T, DT, Time, DTime, Wm, Wm_r, Wm_ir, Wm_d, ndi, nshr, start, solver_type, tnew_dt);
					break;
				}
				case 6: {
					simcoon::umat_plasticity_kin_iso_CCP(etot, Detot, sigma, Lt, L, sigma_in, DR, nprops, props, nstatev, statev, T, DT, Time, DTime, Wm, Wm_r, Wm_ir, Wm_d, ndi, nshr, start, solver_type, tnew_dt);
					break;
				}
				case 7: {
					simcoon::umat_plasticity_chaboche_CCP(etot, Detot, sigma, Lt, L, sigma_in, DR, nprops, props, nstatev, statev, T, DT, Time, DTime, Wm, Wm_r, Wm_ir, Wm_d, ndi, nshr, start, solver_type, tnew_dt);
					break;
				}
				case 8: {
					simcoon::umat_sma_unified_T(etot, Detot, sigma, Lt, DR, nprops, props, nstatev, statev, T, DT, Time, DTime, Wm, Wm_r, Wm_ir, Wm_d, ndi, nshr, start, tnew_dt);
					break;
				}
				case 9: {
					simcoon::umat_damage_LLD_0(etot, Detot, sigma, Lt, L, sigma_in, DR, nprops, props, nstatev, statev, T, DT, Time, DTime, Wm, Wm_r, Wm_ir, Wm_d, ndi, nshr, start, solver_type, tnew_dt);
					break;
				}
				case 10: {
					simcoon::umat_zener_fast(etot, Detot, sigma, Lt, DR, nprops, props, nstatev, statev, T, DT, Time, DTime, Wm, Wm_r, Wm_ir, Wm_d, ndi, nshr, start, tnew_dt);
					break;
				}
				case 11: {
					simcoon::umat_zener_Nfast(etot, Detot, sigma, Lt, DR, nprops, props, nstatev, statev, T, DT, Time, DTime, Wm, Wm_r, Wm_ir, Wm_d, ndi, nshr, start, tnew_dt);
					break;
				}
				case 12: {
					simcoon::umat_prony_Nfast(etot, Detot, sigma, Lt, DR, nprops, props, nstatev, statev, T, DT, Time, DTime, Wm, Wm_r, Wm_ir, Wm_d, ndi, nshr, start, tnew_dt);
					break;
				}
				case 17: {
					simcoon::umat_plasticity_hill_isoh_CCP(etot, Detot, sigma, Lt, DR, nprops, props, nstatev, statev, T, DT, Time, DTime, Wm, Wm_r, Wm_ir, Wm_d, ndi, nshr, start, tnew_dt);
					break;
				}
				case 18: {
					simcoon::umat_plasticity_hill_isoh_CCP_N(etot, Detot, sigma, Lt, DR, nprops, props, nstatev, statev, T, DT, Time, DTime, Wm, Wm_r, Wm_ir, Wm_d, ndi, nshr, start, tnew_dt);
					break;
				}
				case 19: {
					simcoon::umat_sma_mono(etot, Detot, sigma, Lt, L, DR, nprops, props, nstatev, statev, T, DT, Time, DTime, Wm, Wm_r, Wm_ir, Wm_d, ndi, nshr, start, tnew_dt);

					break;
				}
				case 20: {
					simcoon::umat_sma_mono_cubic(etot, Detot, sigma, Lt, L, DR, nprops, props, nstatev, statev, T, DT, Time, DTime, Wm, Wm_r, Wm_ir, Wm_d, ndi, nshr, start, tnew_dt);
					break;
				}
							/*case 100: case 101: case 103: case 104: {
								umat_multi(rve, DR, Time, DTime, ndi, nshr, start, solver_type, tnew_dt, list_umat[rve.sptr_matprops->umat_name]);
								break;
							}*/
				default: {
					cout << "Error: The choice of Umat could not be found in the umat library :" << umat_name_py << "\n";
					exit(0);
				}
			}			

			list_statev.col(pg) = statev;

			if (corate == -1) {
				//without conversion: constitutive eq assumed expressed in GL/PK2
				list_sigma.col(pg) = sigma;
				list_Lt.slice(pg) = Lt;
				list_L.slice(pg) = L;
				list_etot.col(pg) = etot;
				}
			else if ((corate == 0) || (corate == 1) || (corate == 2)) {
				//Constitutive eq assumed expressed in Kirchoff/Logstrain
				//Convert kirkoff/Logstrain to PKII/GLstrain
				const mat sigma_t = simcoon::v2t_stress(sigma);
				list_Lt.slice(pg) = simcoon::DtauDe_2_DSDE(Lt, simcoon::get_BBBB(F1), F1, sigma_t); //transform the tangeant matrix into pkII/green lagrange
				list_L.slice(pg) = simcoon::DtauDe_2_DSDE(L, simcoon::get_BBBB(F1), F1, sigma_t); //transform the elastic matrix into pkII/green lagrange
				list_sigma.col(pg) = simcoon::t2v_stress(simcoon::Kirchoff2PKII(sigma_t, F1));
				list_etot.col(pg) = simcoon::rotate_strain(etot, DR) + Detot;
				}
			else if (corate == 3) {
				//Constitutive eq assumed expressed in Kirchoff/Logstrain
				//Convert kirkoff/Logstrain to PKII/GLstrain
				const mat sigma_t = simcoon::v2t_stress(sigma);
				list_Lt.slice(pg) = simcoon::DtauDe_2_DSDE(Lt, simcoon::get_BBBB(F1), F1, sigma_t); //transform the tangeant matrix into pkII/green lagrange
				list_L.slice(pg) = simcoon::DtauDe_2_DSDE(L, simcoon::get_BBBB(F1), F1, sigma_t); //transform the elastic matrix into pkII/green lagrange
				list_sigma.col(pg) = simcoon::t2v_stress(simcoon::Kirchoff2PKII(sigma_t, F1));
				list_etot.col(pg) = simcoon::t2v_strain(0.5 * logmat_sympd(simcoon::L_Cauchy_Green(F1)));
				//list_etot.col(pg) = simcoon::rotate_strain(etot, DR) + Detot;
			}

			/*else if ((corate == 0) || (corate == 1) || (corate == 2)) {
				//Constitutive eq assumed expressed in Cauchy/Logstrain. 
				//Convert Cauchy/Logstrain to PKII/GLstrain the tangent matrix is not well defined (possibility of convergence problem)
				const mat sigma_t = simcoon::v2t_stress(sigma);
				list_Lt.slice(pg) = simcoon::DtauDe_2_DSDE(Lt, simcoon::get_BBBB(F1), F1, sigma_t); //transform the tangeant matrix into pkII/green lagrange
				list_L.slice(pg) = simcoon::DtauDe_2_DSDE(L, simcoon::get_BBBB(F1), F1, sigma_t); //transform the elastic matrix into pkII/green lagrange
				list_sigma.col(pg) = simcoon::t2v_stress(simcoon::Cauchy2PKII(sigma_t, F1, 0));
				}*/
		}

		//convert modified vec 2 numpy array
		//sigma_py = vec2array(sigma); //not required in principle
		//Lt_py = mat2array(Lt);
		//L_py = mat2array(L);
		//sigma_in_py = vec2array(sigma_in);
		//statev_py = vec2array(statev);

		//return bp::make_tuple(cube2array_inplace(etot), cube2array_inplace(Detot));
		return bp::make_tuple(mat2array(F0), mat2array(F1), vec2array(etot), vec2array(Detot), cube2array(list_Lt), cube2array(list_L), vec2array(statev), vec2array(props), mat2array_inplace(list_sigma));
	}



}