
#include <armadillo>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <simcoon/arma2numpy/numpy_arma.hpp>

#include <simcoon/python_wrappers/Libraries/Continuum_mechanics/Umat_fedoo.hpp>

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
	

//=====Private methods for phase_multi===================================

//=====Public methods for phase_multi====================================


  //Constructor
	//-------------------------------------------------------------
	Umat_fedoo::Umat_fedoo(const std::string& umat_name_py, const bn::ndarray& props_py, const int& corate_py, const int& ndi_py, const int& nshr_py) {
	//-------------------------------------------------------------			
		//Get the id of umat
		std::map<string, int> list_umat;
		list_umat = { {"UMEXT",0},{"UMABA",1},{"ELISO",2},{"ELIST",3},{"ELORT",4},{"EPICP",5},{"EPKCP",6},{"EPCHA",7},{"SMAUT",8},{"LLDM0",9},{"ZENER",10},{"ZENNK",11},{"PRONK",12},{"EPHIC",17},{"EPHIN",18},{"SMAMO",19},{"SMAMC",20},{"MIHEN",100},{"MIMTN",101},{"MISCN",103},{"MIPLN",104} };
		id_umat = list_umat[umat_name_py];

		//Get the list of propertie
		//list_props = arrayT2mat_inplace(props_py);  //each col is prop for one pg
		list_props = array2mat(props_py, false);  //each col is prop for one pg

		nprops = list_props.n_rows;
		nstatev = list_statev.n_rows;

		corate = corate_py;
		ndi = ndi_py; nshr = nshr_py; ncomp = ndi + nshr;
		
		nb_points = 0; //call init_statev to initialize list_statev, nb_points and all the array defined for each pg.
		//nb_points = 0 means than statev is not defined. Once defined, nb_points will be the total number of material points where a constitutive law is required (gauss point for standard FE analysis).
	}

	//-------------------------------------------------------------
	void Umat_fedoo::Initialize(const double &Time_py, bn::ndarray& statev_py, const int& nlgeom_py){
	//-------------------------------------------------------------
		//TODO: rajouter update_dt dans les arguments pour alleger le stockage dans le cas o� on n'a pas besoin de modifier le pas de temps

		Time = Time_py;
		nlgeom = nlgeom_py;
		
		//list_statev = array2mat_inplace(statev_py);		
		list_statev = array2mat(statev_py);
		nb_points = list_statev.n_cols;

		// initialization of Wm and etot to 0
		list_Wm.zeros(4, nb_points);
		list_etot.zeros(ncomp, nb_points);
		list_Detot.zeros(ncomp, nb_points);

		// initialization of F0 to identity matrix
		listF0.set_size(ndi, ndi, 1);
		listF0.slice(0).eye(ndi,ndi);
		
		listF1.set_size(ndi, ndi, nb_points);
		for (int pg = 0; pg < nb_points; pg++) listF1.slice(pg).eye(ndi, ndi);

		// allow memory for variable that will computed in the Umat		
		if (nlgeom == 1) {
			//total lagrangian method based on PKII
			list_cauchy.zeros(ncomp, nb_points);
			list_PKII.zeros(ncomp, nb_points);
			list_PKII_start.zeros(ncomp, nb_points);
		}
		else {
			// nlgeom = false or update lagrangian method based on cauchy
			list_cauchy.zeros(ncomp, nb_points);
			list_cauchy_start.zeros(ncomp, nb_points);
		}

		list_R.set_size(ndi, ndi, nb_points);
		for (int pg = 0; pg < nb_points; pg++) list_R.slice(pg).eye(ndi,ndi);		 
		list_DR.set_size(ndi, ndi, nb_points); 
		for (int pg = 0; pg < nb_points; pg++) list_DR.slice(pg).eye(ndi,ndi);

		list_Lt.set_size(ncomp, ncomp, nb_points); //tangent rigitidy matrix in PKII/GL
		list_L.set_size(ncomp, ncomp, nb_points); //elastic rigidity matrix in Cauchy/logstrain

		list_statev_start = list_statev;
		list_Wm_start = list_Wm;
	}

	//-------------------------------------------------------------
	void Umat_fedoo::set_start() {
	//-------------------------------------------------------------
		list_statev_start = list_statev;
		list_Wm_start = list_Wm;
		mat F1, sigma_t, tau_t;

		if (nlgeom==1) {
			list_PKII_start = list_PKII;
			for (int pg = 0; pg < nb_points; pg++) {
				list_Lt.slice(pg) = simcoon::rotate_stress(list_Lt.slice(pg), list_R.slice(pg));
				list_R.slice(pg) = list_DR.slice(pg)*list_R.slice(pg);			
				list_cauchy_start.col(pg) = simcoon::rotate_stress(list_cauchy.col(pg), list_DR.slice(pg));				

				if (corate == 2) {
					list_etot.col(pg) = simcoon::t2v_strain(simcoon::Log_strain(listF1.slice(pg)));
				}
				else {
					list_etot.col(pg) = simcoon::rotate_strain(list_etot.col(pg), list_DR.slice(pg)) + list_Detot.col(pg);
				}
				
				// Replace the tangent matrix with the elastic matrix for a prediction
				// Constitutive eq assumed expressed in Cauchy / Logstrain 		
				if (corate == 0) {
					//Convert cauchy/Logstrain to PKII/GLstrain
					F1 = listF1.slice(pg);
					sigma_t = simcoon::v2t_stress(list_cauchy.col(pg));
					tau_t = simcoon::Cauchy2Kirchoff(sigma_t, F1);
					list_Lt.slice(pg) = simcoon::DsigmaDe_JaumannDD_2_DSDE(list_L.slice(pg), F1, tau_t); //transform the tangent matrix into pkII/green lagrange				
				}
				else if (corate == 1) {
					//Convert Cauchy/Logstrain to PKII/GLstrain
					F1 = listF1.slice(pg);
					sigma_t = simcoon::v2t_stress(list_cauchy.col(pg));
					tau_t = simcoon::Cauchy2Kirchoff(sigma_t, F1);
					list_Lt.slice(pg) = simcoon::DsigmaDe_2_DSDE(list_L.slice(pg), simcoon::get_BBBB_GN(F1), F1, tau_t); //transform the tangent matrix into pkII/green lagrange
				}
				else if (corate == 2) {
					//Convert Cauchy/Logstrain to PKII/GLstrain
					F1 = listF1.slice(pg);
					sigma_t = simcoon::v2t_stress(list_cauchy.col(pg));
					tau_t = simcoon::Cauchy2Kirchoff(sigma_t, F1);
					list_Lt.slice(pg) = simcoon::DsigmaDe_2_DSDE(list_L.slice(pg), simcoon::get_BBBB(F1), F1, tau_t); //transform the tangent matrix into pkII/green lagrange
				}
			}
		}
		else if (nlgeom == 2) {
			//to use with update lagrangian mathod. Tangent matrix is expressed on the current configuration 
			for (int pg = 0; pg < nb_points; pg++) {
				list_Lt.slice(pg) = simcoon::rotate_stress(list_Lt.slice(pg), list_R.slice(pg));
				list_R.slice(pg) = list_DR.slice(pg)*list_R.slice(pg);			
				list_cauchy_start.col(pg) = simcoon::rotate_stress(list_cauchy.col(pg), list_DR.slice(pg));				
				if (corate == 2) {
					list_etot.col(pg) = simcoon::t2v_strain(simcoon::Log_strain(listF1.slice(pg)));
				}
				else {
					list_etot.col(pg) = simcoon::rotate_strain(list_etot.col(pg), list_DR.slice(pg)) + list_Detot.col(pg);				
				}

			}
		}
		else {
			//nlgeom == 0 
			list_cauchy_start = list_cauchy;
			list_etot = list_etot + list_Detot;
			list_Lt = list_L;
		}

		list_Lt_start = list_Lt;

		Time = Time + DTime;
		listF0 = listF1; //tester si &listF0 = listF1 est mieux ?
		//T += DT; to add if required				
	}

	//-------------------------------------------------------------
	void Umat_fedoo::to_start() {
	//-------------------------------------------------------------
		// Replace the tangent matrix with the elastic matrix for a prediction		
		list_Lt = list_Lt_start;
		if (nlgeom == 1) list_PKII = list_PKII_start;
		else list_cauchy = list_cauchy_start;
		list_statev = list_statev_start;
		list_Wm = list_Wm_start;
	}

	//-------------------------------------------------------------
	void Umat_fedoo::compute_Detot(const double& DTime, const bn::ndarray& F1_py)
	//-------------------------------------------------------------
	{
		//check if this work in 2D
		// Variable containing list for all pg values
		listF1 = array2cube(F1_py, false); //inplace (without copy).
		
		//mat listDetot(6, listF1.n_slices); //Strain increment (eulerian) 
		//cube listDR(3, 3, listF1.n_slices);
		
		//nb_points = listF1.n_slices;

		// Variable for only one pg 
		mat F0(ndi, ndi); mat F1(ndi, ndi); 
		mat DR(ndi, ndi);
		mat D(ndi, ndi);
		mat Omega(ndi, ndi);
		mat dF;

		if (!nlgeom) {
			DR.eye(ndi, ndi);
			dF.set_size(ndi, ndi);
		}

		for (int pg = 0; pg < nb_points; pg++) {
			if (pg < listF0.n_slices) F0 = listF0.slice(pg); //if listF0 has only one element, we keep only this one			
			F1 = listF1.slice(pg);

			if (!nlgeom) {
				//Green-Lagrange Strain
				//etot = simcoon::t2v_strain(simcoon::Green_Lagrange(F0));
				//list_Detot.col(pg) = simcoon::t2v_strain(simcoon::Green_Lagrange(F1)) - simcoon::t2v_strain(simcoon::Green_Lagrange(F0));				
				dF = F1 - F0;
				list_Detot.col(pg) = 0.5*simcoon::t2v_strain(dF + dF.t());
			}
			else if (corate == 0) {
				//Jaumann
				simcoon::Jaumann(DR, D, Omega, DTime, F0, F1); //to compute D, W, Omega
				list_Detot.col(pg) = simcoon::t2v_strain(simcoon::Delta_log_strain(D, Omega, DTime)); //etot : VR decomposition, then ln(V) equals the logarithmic strain			
			}
			else if (corate == 1) {
				//Green Naghdi
				simcoon::Green_Naghdi(DR, D, Omega, DTime, F0, F1); //to compute D, W, Omega
				list_Detot.col(pg) = simcoon::t2v_strain(simcoon::Delta_log_strain(D, Omega, DTime)); //etot : VR decomposition, then ln(V) equals the logarithmic strain			
			}
			else if (corate == 2) {
				//Log Strain 
				simcoon::logarithmic(DR, D, Omega, DTime, F0, F1); //to compute D, W, Omega
				list_Detot.col(pg) = simcoon::t2v_strain(simcoon::Log_strain(listF1.slice(pg)))-list_etot.col(pg);
				//list_Detot.col(pg) = simcoon::t2v_strain(simcoon::Delta_log_strain(D, Omega, DTime)); //etot : VR decomposition, then ln(V) equals the logarithmic strain			
			}
			list_DR.slice(pg) = DR;
		}
	}

	//-------------------------------------------------------------
	bp::tuple Umat_fedoo::Run(const double& DTime_py) {
	//-------------------------------------------------------------

		//DTime == 0. for the very first increment to comput the elastic matrix
		
		//scalar needed to launch umat
		const int solver_type = 0;
		const bool start = false;
		double tnew_dt = 0;//usefull ?		
		double T = 0; double DT = 0;  //modify to let the program set the actual temperature
		bool use_temp;
		if (list_T.n_elem == 0.) use_temp = false; 
		else use_temp = true;


		DTime = DTime_py; 
								 
		//initial temp local variable for umat
		vec sigma(ncomp);
		vec sigma_in = zeros(1); //not used
		vec statev; 
		vec props;
		vec etot; vec Detot; //Eulerian strain
		mat L(ncomp, ncomp);
		mat Lt(ncomp, ncomp);
		mat DR(ndi, ndi);
		mat F1(ndi, ndi);
		mat sigma_t(ndi, ndi);
		mat tau_t(ndi, ndi);

		double Wm;
		double Wm_r;
		double Wm_ir;
		double Wm_d;

		for (int pg = 0; pg < nb_points; pg++) {
			if (use_temp) T = list_T(pg);
			if (pg < list_props.n_cols) props = list_props.col(pg); //if list_props has only one element, we keep only this one (assuming homogeneous material)			
			statev = list_statev_start.col(pg);
			sigma = list_cauchy_start.col(pg);
			DR = list_DR.slice(pg);

			etot = list_etot.col(pg);
			Detot = list_Detot.col(pg);
			//L = list_L.slice(pg);
			Lt = list_Lt.slice(pg);

			Wm = list_Wm_start(0, pg);
			Wm_r = list_Wm_start(1, pg);
			Wm_ir = list_Wm_start(2, pg);
			Wm_d = list_Wm_start(3, pg);

			//Launch the simcoon umat
			switch (id_umat) {

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
					cout << "Error: The choice of Umat could not be found in the umat library \n";
					exit(0);
				}
			}

			list_statev.col(pg) = statev;
			list_Wm(0, pg) = Wm;
			list_Wm(1, pg) = Wm_r;
			list_Wm(2, pg) = Wm_ir;
			list_Wm(3, pg) = Wm_d;

			if (DTime == 0.) { // if DTime = 0., init the elastic matrix for future use in set_start method, sigma should be = 0)
				list_L.slice(pg) = Lt;
			}
			
			if (nlgeom== 0 || nlgeom == 2) {
				//without conversion of stress and rigidity matrix
				list_cauchy.col(pg) = sigma;
				list_Lt.slice(pg) = Lt;
			}
			else if (corate == 0) {
				//Constitutive eq assumed expressed in Cauchy/Logstrain
				//Convert kirkoff/Logstrain to PKII/GLstrain
				list_cauchy.col(pg) = sigma;
				F1 = listF1.slice(pg);
				sigma_t = simcoon::v2t_stress(sigma);
				tau_t = simcoon::Cauchy2Kirchoff(sigma_t, F1);
				list_Lt.slice(pg) = simcoon::DsigmaDe_JaumannDD_2_DSDE(Lt, F1, tau_t); //transform the tangeant matrix into pkII/green lagrange
				list_PKII.col(pg) = simcoon::t2v_stress(simcoon::Cauchy2PKII(sigma_t, F1));
			}
			else if (corate == 1) {
				//Constitutive eq assumed expressed in Cauchy/Logstrain
				//Convert cauchy/Logstrain to PKII/GLstrain
				list_cauchy.col(pg) = sigma;
				F1 = listF1.slice(pg);
				sigma_t = simcoon::v2t_stress(sigma);
				tau_t = simcoon::Cauchy2Kirchoff(sigma_t, F1);
				list_Lt.slice(pg) = simcoon::DsigmaDe_2_DSDE(Lt, simcoon::get_BBBB_GN(F1), F1, tau_t); //transform the tangeant matrix into pkII/green lagrange
				list_PKII.col(pg) = simcoon::t2v_stress(simcoon::Cauchy2PKII(sigma_t, F1));
			}
			else if (corate == 2) {
				//Constitutive eq assumed expressed in Cauchy/Logstrain
				//Convert cauchy/Logstrain to PKII/GLstrain
				list_cauchy.col(pg) = sigma;
				F1 = listF1.slice(pg);
				sigma_t = simcoon::v2t_stress(sigma);
                tau_t = simcoon::Cauchy2Kirchoff(sigma_t, F1);
				list_Lt.slice(pg) = simcoon::DsigmaDe_2_DSDE(Lt, simcoon::get_BBBB(F1), F1, tau_t); //transform the tangeant matrix into pkII/green lagrange
				list_PKII.col(pg) = simcoon::t2v_stress(simcoon::Cauchy2PKII(sigma_t, F1));
			}
		}
		if (DTime == 0.) { list_Lt_start = list_Lt; } //1st iteration only -> inint Lt_start

		//return bp::make_tuple(mat2array(DR), cube2array(list_DR, false), vec2array(Detot), vec2array(statev));
		return bp::make_tuple();
	}

	//-------------------------------------------------------------
	void Umat_fedoo::set_T(bn::ndarray& T_py) {
	//-------------------------------------------------------------		
		list_T = array2vec(T_py); //copy -> todo : check if we can do without copying ?
	}

	//-------------------------------------------------------------
	bn::ndarray Umat_fedoo::Get_T() {
	//-------------------------------------------------------------
		return vec2array(list_T, false); //inplace (without copy). Return the transpose by default (C_contiguous)
	}

	//-------------------------------------------------------------
	bn::ndarray Umat_fedoo::Get_props() {
	//-------------------------------------------------------------
		return mat2array(list_props, false); //inplace (without copy). Return the transpose by default (C_contiguous)
	}

	//-------------------------------------------------------------
	bn::ndarray Umat_fedoo::Get_PKII() {
	//-------------------------------------------------------------
		if (nlgeom==1) {
			return mat2array(list_PKII, false); //inplace (without copy). Return the transpose by default (C_contiguous)
		}
		else if (nlgeom==0) {
			return mat2array(list_cauchy, false); //inplace (without copy). Return the transpose by default (C_contiguous)
		}
		else {
			if (list_PKII.is_empty()) {
				list_PKII.set_size(ncomp, nb_points);
			}

			mat F1(ndi, ndi);
			mat sigma_t;

			for (int pg = 0; pg < nb_points; pg++) {
				F1 = listF1.slice(pg);
				sigma_t = simcoon::v2t_stress(list_cauchy.col(pg));
				list_PKII.col(pg) = simcoon::t2v_stress(simcoon::Cauchy2PKII(sigma_t, F1));
			}
			return mat2array(list_cauchy, false); //inplace (without copy). Return the transpose by default (C_contiguous)
		}
	}

	//-------------------------------------------------------------
	bn::ndarray Umat_fedoo::Get_Cauchy() {
		//-------------------------------------------------------------
		return mat2array(list_cauchy, false); //inplace (without copy). Return the transpose by default (C_contiguous)
	}

	//-------------------------------------------------------------
	bn::ndarray Umat_fedoo::Get_Kirchhoff() {
	//-------------------------------------------------------------
		if (!nlgeom) return mat2array(list_cauchy, false); //inplace (without copy). Return the transpose by default (C_contiguous)
		
		if (list_kirchoff.is_empty()) {
			list_kirchoff.set_size(ncomp, nb_points);
		}

		mat F1(ndi, ndi);
		vec cauchy;

		for (int pg = 0; pg < nb_points; pg++) {			
			F1 = listF1.slice(pg);
			cauchy = list_cauchy.col(pg);
			list_kirchoff.col(pg) = simcoon::Cauchy2Kirchoff(cauchy, F1);
		}
		return mat2array(list_cauchy, false); //inplace (without copy). Return the transpose by default (C_contiguous)
	}

	//-------------------------------------------------------------
	bn::ndarray Umat_fedoo::Get_etot() {
	//-------------------------------------------------------------
		return mat2array(list_etot, false); //inplace (without copy). Return the transpose by default (C_contiguous)
	}

	//-------------------------------------------------------------
	bn::ndarray Umat_fedoo::Get_Detot() {
	//-------------------------------------------------------------
		return mat2array(list_Detot, false); //inplace (without copy). Return the transpose by default (C_contiguous)
	}

	//-------------------------------------------------------------
	bn::ndarray Umat_fedoo::Get_statev() {
	//-------------------------------------------------------------
		return mat2array(list_statev, false); //inplace (without copy). Return the transpose by default (C_contiguous)
	}

	//-------------------------------------------------------------
	bn::ndarray Umat_fedoo::Get_L() {
	//-------------------------------------------------------------
		//L is expressed according to the choosen strain and stress mesure (mainly cauchy/Logstrain)
		return cube2array(list_L, false); //inplace (without copy).
	}

	//-------------------------------------------------------------
	bn::ndarray Umat_fedoo::Get_Lt() {
	//-------------------------------------------------------------
		// Lt is expressed PKII/GLstrain
		return cube2array(list_Lt, false); //inplace (without copy).
	}

	//-------------------------------------------------------------
	bn::ndarray Umat_fedoo::Get_R() {
	//-------------------------------------------------------------
		return cube2array(list_R, false); //inplace (without copy).
	}

	//-------------------------------------------------------------
	bn::ndarray Umat_fedoo::Get_DR() {
	//-------------------------------------------------------------
		return cube2array(list_DR, false); //inplace (without copy).
	}

	//-------------------------------------------------------------
	bn::ndarray Umat_fedoo::Get_Wm() {
	//-------------------------------------------------------------
		return mat2array(list_Wm, false); //inplace (without copy). Return the transpose by default (C_contiguous)
	}

	//-------------------------------------------------------------
	bn::ndarray Umat_fedoo::Get_F0() {
	//-------------------------------------------------------------
		return cube2array(listF0, false); //inplace (without copy).
	}

	//-------------------------------------------------------------
	bn::ndarray Umat_fedoo::Get_F1() {
	//-------------------------------------------------------------
		return cube2array(listF1, false);//inplace (without copy).
	}

	/*bp::tuple Log_strain_fedoo(const bn::ndarray& F_py) {
		cube listF = array2cube(F_py, false); //inplace (without copy).
		mat list_etot(6,listF.n_slices); //work for 2D problem ?

		// convert array to armdillo vec
		mat F;

		for (int pg = 0; pg < listF.n_slices; pg++) {
			F = listF.slice(pg);
			list_etot.col(pg) = simcoon::Log_strain(F);
		}

		return bp::make_tuple(mat2array(list_etot, false, "F"));
	}*/
}
