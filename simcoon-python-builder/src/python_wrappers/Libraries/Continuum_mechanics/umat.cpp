#include <pybind11/embed.h>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <optional>
#include <vector>

#include <simcoon/python_wrappers/arma_to_numpy.hpp>
#include <simcoon/python_wrappers/numpy_to_arma.hpp>
#include <armadillo>

#include <simcoon/parameter.hpp>
#include <simcoon/python_wrappers/parallel_nogil.hpp>

#include <simcoon/python_wrappers/Libraries/Continuum_mechanics/umat.hpp>

#include <simcoon/Continuum_mechanics/Umat/umat_smart.hpp>
#include <simcoon/Continuum_mechanics/Umat/condensed_state.hpp>
#include <simcoon/Continuum_mechanics/Functions/kinematics.hpp>
#include <simcoon/Continuum_mechanics/Umat/umat_callback.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/External/external_umat.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Plasticity/plastic_isotropic.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Plasticity/plastic_chaboche.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Plasticity/plastic_johnson_cook_ccp.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/SMA/unified_T.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/SMA/unified_TR.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/SMA/SMA_mono.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Damage/damage_LLD_0.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Viscoelasticity/Zener_fast.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Viscoelasticity/Zener_Nfast.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Viscoelasticity/Prony_Nfast.hpp>

#include <simcoon/Continuum_mechanics/Umat/Finite/generic_hyper_invariants.hpp>
#include <simcoon/Continuum_mechanics/Umat/Finite/generic_hyper_pstretch.hpp>
#include <simcoon/Continuum_mechanics/Umat/Finite/saint_venant.hpp>
#include <simcoon/Continuum_mechanics/Umat/Finite/hypoelastic_orthotropic.hpp>
#include <simcoon/Continuum_mechanics/Umat/Finite/neo_hookean_incomp.hpp>

#include <simcoon/Continuum_mechanics/Umat/Modular/modular_umat.hpp>
#include <simcoon/Continuum_mechanics/Umat/Modular/legacy_adapters.hpp>

#include <simcoon/Simulation/Maths/rotation.hpp> //for rotate_strain
#include <simcoon/Continuum_mechanics/Functions/transfer.hpp>
#include <simcoon/Continuum_mechanics/Functions/objective_rates.hpp>

#include <simcoon/Continuum_mechanics/Umat/Thermomechanical/External/external_umat.hpp>
#include <simcoon/Continuum_mechanics/Umat/Thermomechanical/Elasticity/elastic_isotropic.hpp>
#include <simcoon/Continuum_mechanics/Umat/Thermomechanical/Elasticity/elastic_isotropic.hpp>
#include <simcoon/Continuum_mechanics/Umat/Thermomechanical/Elasticity/elastic_transverse_isotropic.hpp>
#include <simcoon/Continuum_mechanics/Umat/Thermomechanical/Elasticity/elastic_orthotropic.hpp>
#include <simcoon/Continuum_mechanics/Umat/Thermomechanical/Plasticity/plastic_isotropic.hpp>
#include <simcoon/Continuum_mechanics/Umat/Thermomechanical/Plasticity/plastic_kin_iso.hpp>
#include <simcoon/Continuum_mechanics/Umat/Thermomechanical/Plasticity/plastic_johnson_cook_ccp.hpp>
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

namespace {

// Lab start stress sym(DR^-1 tau_start_tr DR) from the one passed (transported by the caller),
// for Delta_work_conjugacy. Fixed size and the bool inv(): no heap (numpy allocator, GIL), no
// throw in the parallel region; tau_start_tr itself on a singular DR.
arma::vec::fixed<6> lab_start_stress(const arma::vec::fixed<6> &tau_start_tr, const arma::mat &DR) {
	arma::mat::fixed<3,3> DR_inv;
	if (!arma::inv(DR_inv, arma::mat::fixed<3,3>(DR)))
		return tau_start_tr;
	const arma::mat::fixed<3,3> X = DR_inv*simcoon::v2t_stress(tau_start_tr)*DR;
	return simcoon::t2v_stress(0.5*(X + X.t()));
}

// Raise simcoon.StepCut for a batch whose kernels asked for a smaller increment (tnew_dt < 1),
// with the smallest ratio asked for. Serial context, GIL held. The inputs are copied before the
// kernels run, so the caller's arrays are untouched and the call can simply be retried.
[[noreturn]] void raise_step_cut(const std::string &entry, const std::vector<double> &tnew_dt) {
    std::vector<size_t> points;
    double ratio = 1.;
    for (size_t pt = 0; pt < tnew_dt.size(); ++pt) {
        if (tnew_dt[pt] < 1.) {
            points.push_back(pt);
            ratio = std::min(ratio, tnew_dt[pt]);
        }
    }
    const std::string msg = entry + ": the law requested a step cut at "
        + std::to_string(points.size()) + " material point(s), the first being point "
        + std::to_string(points.front()) + ". The batch entry cannot subdivide the increment: "
        "discard this call and retry with a smaller increment (the input arrays are untouched).";
    py::object exc;
    try {
        exc = py::module_::import("simcoon.pyumat").attr("StepCut")(py::arg("ratio") = ratio,
                                                                    py::arg("msg") = msg);
    } catch (py::error_already_set &) {   // bare _core use, without the python package
        exc = py::module_::import("simcoon._core").attr("StepCut")(msg);
        exc.attr("ratio") = ratio;
    }
    PyErr_SetObject(reinterpret_cast<PyObject *>(Py_TYPE(exc.ptr())), exc.ptr());
    throw py::error_already_set();
}

// Code of a kinematic form given as "total", "rate" (log_R), or an objective rate by name or
// code; -1 when not given. The names are those of simcoon.solver.maps.
int kinematic_code(const py::object &value) {
    if (value.is_none()) return -1;
    if (py::isinstance<py::str>(value)) {
        static const std::map<std::string, int> names = {
            {"total", simcoon::corate_type_total}, {"rate", 3}, {"jaumann", 0}, {"green_naghdi", 1},
            {"logarithmic", 2}, {"logarithmic_R", 3}, {"truesdell", 4}, {"logarithmic_F", 5}};
        const auto it = names.find(value.cast<std::string>());
        if (it == names.end()) {
            throw std::invalid_argument("kinematic_form (or corate) must be 'total', 'rate' or an objective rate "
                                        "('jaumann', 'green_naghdi', 'logarithmic', 'logarithmic_R', 'truesdell', "
                                        "'logarithmic_F'); got '" + value.cast<std::string>() + "'");
        }
        return it->second;
    }
    const int code = value.cast<int>();
    if (code < 0 || code > 5) {
        throw std::invalid_argument("corate must be 0 (Jaumann), 1 (Green-Naghdi), 2 (logarithmic), "
                                    "3 (logarithmic_R), 4 (Truesdell) or 5 (logarithmic_F); got "
                                    + std::to_string(code));
    }
    return code;
}

// (ndi, nshr) of a batch: nshr defaults to Abaqus' for ndi; ndi < 3 is a condensed state,
// solved per point around the 3D law (condensed_state.hpp)
simcoon::condensed_state batch_stress_state(const int &ndi, const std::optional<int> &nshr_py, int &nshr) {
	if (ndi < 1 || ndi > 3)
		throw std::invalid_argument("ndi should be 1, 2 or 3");
	if (nshr_py && (*nshr_py < 0 || *nshr_py > 3))
		throw std::invalid_argument("nshr should be 0, 1, 2 or 3");
	nshr = nshr_py.value_or(simcoon::default_nshr(ndi));
	return simcoon::make_condensed_state(ndi, nshr);
}

// Deformation gradient of a condensed finite-strain point. The caller's F holds the rotation and
// the kept stretches only: F = exp(ln V) R with ln V on the free components replaced by the strain
// state, etot at the start of the increment, etot + Detot_pt at its end
struct condensed_F {
	arma::vec lnV;
	arma::mat R;
	explicit condensed_F(const arma::mat &F) {
		arma::mat V;
		simcoon::VR_decomposition(V, R, F);
		lnV = simcoon::t2v_strain(arma::logmat_sympd(V));
	}
	arma::mat operator()(const simcoon::condensed_state &cs, const arma::vec &e) const {
		arma::vec l = lnV;
		l(cs.C) = e(cs.C);
		return simcoon::eR_to_F(simcoon::v2t_strain(l), R);
	}
	// the total form: the strain state is the whole H = ln U, F = R exp(H)
	arma::mat total(const arma::vec &H) const { return R*arma::expmat_sym(simcoon::v2t_strain(H)); }
};

// condense = 'approximate': the stress the law returns vanishes on the free set (tau_hat_C under the
// rate form, T_C under the total one), Newton with the law's own tangent. 'exact': the second
// Piola-Kirchhoff stress S_C vanishes, i.e. the traction on the free material faces (P N = 0), with
// the Green-Lagrange E_C as unknowns and dS/dE as Jacobian; it needs F0, F1. Both coincide at small
// strain and when the free directions are principal (condensed_state.hpp)
bool condense_exact(const std::string &condense) {
	if (condense == "approximate") return false;
	if (condense == "exact") return true;
	throw std::invalid_argument("condense must be 'approximate' or 'exact'; got '" + condense + "'");
}

// ndi < 3: the solved strain goes back into the caller's Detot (rows: the free components, or all
// six under condense='exact', where Detot is the increment of the law's own strain), so that
// etot + Detot carries it to the next increment. Written once the batch has converged, through the
// array's own strides (a copy was taken when it is not F-contiguous)
void write_back_free_strain(const py::array_t<double> &Detot_py, const arma::uvec &rows, const arma::mat &Detot_C) {
	py::array_t<double> out = Detot_py;   // same buffer
	char *base = reinterpret_cast<char *>(out.mutable_data());
	const py::ssize_t s0 = out.strides(0), s1 = (out.ndim() > 1) ? out.strides(1) : 0;
	for (arma::uword pt = 0; pt < Detot_C.n_cols; pt++)
		for (arma::uword j = 0; j < rows.n_elem; j++) {
			const double v = Detot_C(j, pt);
			std::memcpy(base + rows(j)*s0 + pt*s1, &v, sizeof(double));
		}
}

void require_writable_Detot(const py::array_t<double> &Detot_py, const simcoon::condensed_state &cs, const std::string &fname) {
	if (cs.active() && !Detot_py.writeable())
		throw std::invalid_argument(fname + ": with ndi < 3 the solved stress-free components are written back into Detot, which must be writable");
}

}  // namespace
	
	py::tuple launch_umat(const std::string &umat_name_py, const py::array_t<double> &etot_py, const py::array_t<double> &Detot_py, const py::array_t<double> &F0_py, const py::array_t<double> &F1_py, const py::array_t<double> &sigma_py, const py::array_t<double> &DR_py, const py::array_t<double> &props_py, const py::array_t<double> &statev_py, const double Time, const double DTime, const py::array_t<double> &Wm_py, const std::optional<py::array_t<double>> &T_py, const int &ndi, const std::optional<int> &nshr_py, const unsigned int &n_threads, const int &tangent_mode, const py::object &corate_py, const bool &work_correction_on, const std::string &tangent_output, const std::optional<bool> &start_py, const py::object &kinematic_form_py, const std::optional<py::array_t<double>> &DT_py, const std::string &condense){
		// tangent_mode: 0 = none (explicit integration, Lt = elastic L),
		//               1 = continuum, 2 = algorithmic (Simo-Hughes, DEFAULT),
		//               3 = closest-point (reserved). See parameter.hpp tangent_* constants.

		// Validate up front, in serial context: the per-point dispatch below
		// runs inside a non-exception-safe parallel region (GCD/OpenMP) where a
		// throw would std::terminate the host process.
		if (tangent_mode < simcoon::tangent_none || tangent_mode > simcoon::tangent_closest_point) {
			throw std::invalid_argument("tangent_mode must be 0 (none), 1 (continuum), 2 (algorithmic) or 3 (closest-point); got "
			                            + std::to_string(tangent_mode));
		}
		// kinematic_form and corate (its former name) are one value; the dispatch carries the
		// total formulation in the corate slot
		const int form = kinematic_code(kinematic_form_py), former = kinematic_code(corate_py);
		int corate_type = (form >= 0) ? form : ((former >= 0) ? former : 3);
		if (form >= 0 && former >= 0 && form != former) {
			const bool default_rate = py::isinstance<py::str>(kinematic_form_py) && kinematic_form_py.cast<std::string>() == "rate";
			if (!default_rate || former == simcoon::corate_type_total) {
				throw std::invalid_argument("umat: kinematic_form and corate (its former name) are one value; got two that differ");
			}
			corate_type = former;
		}
		const bool total = (corate_type == simcoon::corate_type_total);
		// tangent_output: the box d(tau_hat)/dDe as the kernels return it, or the material dS/dE
		// (TL) / spatial Lie dsigma/dD (UL) tangent, converted per point inside the parallel loop
		enum { tangent_out_box, tangent_out_material, tangent_out_spatial };
		static const std::map<string, int> list_tangent_output = { {"box", tangent_out_box}, {"material", tangent_out_material}, {"spatial", tangent_out_spatial} };
		const auto it_tangent_out = list_tangent_output.find(tangent_output);
		if (it_tangent_out == list_tangent_output.end()) {
			throw std::invalid_argument("tangent_output must be 'box', 'material' or 'spatial'; got '" + tangent_output + "'");
		}
		const int tangent_out = it_tangent_out->second;
		static const std::map<string, int> list_umat = { {"UMEXT",0},{"UMABA",1},{"ELISO",201},{"ELIST",201},{"ELORT",201},{"EPICP",5},{"EPKCP",201},{"EPCHA",7},{"EPHIL",201},{"EPTRI",201},{"EPHAC",201},{"EPANI",201},{"EPDFA",201},{"EPHIN",201},{"SMADI",13},{"SMADC",13},{"SMAAI",13},{"SMAAC",13},{"LLDM0",15},{"ZENER",16},{"ZENNK",17},{"PRONK",18},{"SMAMO",19},{"SMAMC",20},{"NEOHC",21},{"MOORI",22},{"YEOHH",23},{"ISHAH",24},{"GETHH",25},{"SWANH",26},{"HOLZA",27},{"EPCHG",201},{"SMRDI",28},{"SMRDC",28},{"SMRAI",28},{"SMRAC",28},{"SNTVE",29},{"NEOHI",30},{"OGDEN",31},{"HYPOO",32},{"EPJCK",33},{"MODUL",200},{"MIHEN",100},{"MIMTN",101},{"MISCN",103},{"MIPLN",104},{"PYEXT",300} };
		// guarded lookup (serial context): operator[] would default-insert
		// 0 = UMEXT, silently routing typos to the external-plugin path
		const auto it_umat = list_umat.find(umat_name_py);
		if (it_umat == list_umat.end()) {
			throw std::invalid_argument("Unknown umat name: " + umat_name_py);
		}
		const int id_umat = it_umat->second;
		int arguments_type; //depends on the argument used in the umat
		// PYEXT (Python callback law) re-enters the interpreter: it must run on the calling
		// thread, never inside the GCD/OpenMP region (parallel.hpp contract).
		const bool serial = (id_umat == 300);

		// Unified small-strain function pointer: (umat_name, Etot, DEtot, sigma, Lt, L, DR, nprops, props, nstatev, statev, T, DT, Time, DTime, Wm, Wm_r, Wm_ir, Wm_d, ndi, nshr, start, tnew_dt, tangent_mode)
		void (*umat_function)(const std::string &, const arma::vec &, const arma::vec &, arma::vec &, arma::mat &, arma::mat &, const arma::mat &, const int &, const arma::vec &, const int &, arma::vec &, const double &, const double &, const double &, const double &, double &, double &, double &, double &, const int &, const int &, const bool &, double &, const int &);
		// Unified finite-strain function pointer: (umat_name, etot, Detot, F0, F1, sigma, Lt, L, DR, nprops, props, nstatev, statev, T, DT, Time, DTime, Wm, Wm_r, Wm_ir, Wm_d, ndi, nshr, start, tnew_dt, corate_type, tangent_mode)
		void (*umat_function_finite)(const std::string &, const arma::vec &, const arma::vec &, const arma::mat &, const arma::mat &, arma::vec &, arma::mat &, arma::mat &, const arma::mat &, const int &, const arma::vec &, const int &, arma::vec &, const double &, const double &, const double &, const double &, double &, double &, double &, double &, const int &, const int &, const bool &, double &, const int &, const int &);
		const int ncomp=6;
		int nshr;
		const simcoon::condensed_state cs = batch_stress_state(ndi, nshr_py, nshr);

		// start re-initialises the point (T_init, stress, internal variables, Wm):
		// the caller's choice when given, otherwise inferred from Time.
		const bool start = start_py.value_or(Time <= simcoon::limit);

		//bool use_temp;
		//if (T.n_elem == 0.) use_temp = false; 
		//else use_temp = true;

		bool use_temp = false;
		vec vec_T;
		if (T_py.has_value()) {
			vec_T = simpy::numpy_to_arma::arr_to_col_view(T_py.value());
			use_temp = true;
		}
		vec vec_DT;
		if (DT_py.has_value()) {
			if (!use_temp)
				throw std::invalid_argument("umat: DT is the increment of temp, which must be given too");
			vec_DT = simpy::numpy_to_arma::arr_to_col_view(DT_py.value());
		}

		mat list_etot = simpy::numpy_to_arma::arr_to_mat_view(etot_py);
		int nb_points = list_etot.n_cols; //number of material points
		// temp and DT: one value shared by every point, or one per point
		if (use_temp && vec_T.n_elem != 1 && vec_T.n_elem != static_cast<arma::uword>(nb_points))
			throw std::invalid_argument("umat: temp must have 1 entry (shared) or one per material point ("
			                            + std::to_string(nb_points) + "); got " + std::to_string(vec_T.n_elem));
		if (DT_py.has_value() && vec_DT.n_elem != 1 && vec_DT.n_elem != static_cast<arma::uword>(nb_points))
			throw std::invalid_argument("umat: DT must have 1 entry (shared) or one per material point ("
			                            + std::to_string(nb_points) + "); got " + std::to_string(vec_DT.n_elem));
		mat list_Detot = simpy::numpy_to_arma::arr_to_mat_view(Detot_py);
		require_writable_Detot(Detot_py, cs, "umat");
		mat list_sigma = simpy::numpy_to_arma::arr_to_mat(sigma_py); // copy: the umat updates it and it is returned to python
		cube DR = simpy::numpy_to_arma::arr_to_cube_view(DR_py);
		cube F0, F1;

		// props: (nprops, 1) or a 1-D (nprops,) vector = shared by all points; (nprops, n) = per point.
		// (a 1-D array has no shape[1]: reading it was undefined behaviour and fed garbage props
		// to every point after the first)
		vec props;
		mat list_props;
		bool unique_props = false;
		if (props_py.ndim() == 1) {
			props = simpy::numpy_to_arma::arr_to_col(props_py);
			list_props = mat(props.memptr(), props.n_elem, 1, false, true);
			unique_props = true;
		} else {
			if (props_py.ndim() != 2) {
				throw std::invalid_argument("umat: props must be a (nprops,) vector or a (nprops, 1 | n_points) array");
			}
			list_props = simpy::numpy_to_arma::arr_to_mat_view(props_py);
			if (props_py.shape(1) == 1) {
				props = list_props.col(0);
				unique_props = true;
			} else if (props_py.shape(1) != nb_points) {
				throw std::invalid_argument("umat: props must have one column, or one column per material point");
			}
		}

		mat list_statev = simpy::numpy_to_arma::arr_to_mat(statev_py); // copy: the umat updates it and it is returned to python
		mat list_Wm = simpy::numpy_to_arma::arr_to_mat(Wm_py); // copy: the umat updates it and it is returned to python
		cube L(ncomp, ncomp, nb_points);
		cube Lt(ncomp, ncomp, nb_points);
		int nprops = list_props.n_rows;
		int nstatev = list_statev.n_rows;

		switch (id_umat) {
			case 201: { // legacy names served by the modular adapter
				// one shared id (mirrors umat_smart): the adapter dispatches on
				// umat_name, so every legacy-modular name takes this case
				umat_function = &simcoon::umat_legacy_modular;
				arguments_type = 1;
				break;
			}
			case 5: {
				umat_function = &simcoon::umat_plasticity_iso;
				arguments_type = 1;
				break;
			}
			case 33: {
				umat_function = &simcoon::umat_plasticity_johnson_cook_CCP;
				arguments_type = 1;
				break;
			}
			case 7: {
				umat_function = &simcoon::umat_plasticity_chaboche;
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
				arguments_type = 3;   // called by name: it takes the kinematic form
				break;
			}
			case 21: case 22: case 23: case 24: case 25: case 26: case 27: {
				umat_function_finite = &simcoon::umat_generic_hyper_invariants;
				arguments_type = 2;
				break;
			}
			case 32: { // HYPOO (hypoelastic orthotropic, finite): corotational Kirchhoff rate
				umat_function_finite = &simcoon::umat_hypoelasticity_ortho;
				arguments_type = 2;
				break;
			}
			case 29: { // SNTVE (Saint-Venant-Kirchhoff, finite)
				umat_function_finite = &simcoon::umat_saint_venant;
				arguments_type = 2;
				break;
			}
			case 30: { // NEOHI (Neo-Hookean incompressible, finite)
				umat_function_finite = &simcoon::umat_neo_hookean_incomp;
				arguments_type = 2;
				break;
			}
			case 31: { // OGDEN (isochoric principal stretches, finite)
				umat_function_finite = &simcoon::umat_generic_hyper_pstretch;
				arguments_type = 2;
				break;
			}
			case 300: {
				// PYEXT: process-wide callback UMAT (umat_callback.hpp), i.e. a Python law.
				// It re-enters the interpreter, so it runs serially on the calling thread.
				umat_function = &simcoon::umat_callback_M;
				arguments_type = 1;
				break;
			}
			default: {
				throw std::invalid_argument( "The choice of Umat could not be found in the umat library." );
			}
		}

		// Kirchhoff-convention normalization at the python boundary
		// (stress_output_is_kirchhoff = the shared name set). The python
		// contract is: sigma in/out = Cauchy, Lt = the KIRCHHOFF box
		// d(tau_hat)/dDe with no J — the measure the genuine finite kernels
		// already emit and the one Lt_convert's exact box->DSDE map consumes
		// (rescaling Lt by 1/J here would break that consumer by exactly J).
		// So with deformation gradients provided: stress converted on the
		// way in (x J0) and out (/ J1), Lt passed through untouched.
		const bool kirchhoff_normalize = simcoon::stress_output_is_kirchhoff(umat_name_py)
				&& F0_py.size() > 0 && F1_py.size() > 0;
		// Total formulation: the caller passes etot = H(F0), Detot = H(F1) - H(F0), H = ln U, and
		// DR = I. Same admission and stress contract as select_umat_M_finite.
		const bool total_box = total && simcoon::total_formulation_serves(umat_name_py)
				&& simcoon::output_convention_of(umat_name_py).material_frame;
		if (total) {
			if (!simcoon::total_formulation_serves(umat_name_py)) {
				throw std::invalid_argument("umat: the total formulation (strain space) does not serve '"
				                            + umat_name_py + "'; use kinematic_form='rate'");
			}
			if (!kirchhoff_normalize) {
				throw std::invalid_argument("umat: the total formulation (strain space) needs F0 and F1 to map the stress to and from the Cauchy stress");
			}
		}
		// condense='exact' needs the deformation gradients (without them both conditions are sigma_C = 0)
		const bool exact = cs.active() && condense_exact(condense) && kirchhoff_normalize;
		if (exact && !total && corate_type != 3) {
			throw std::invalid_argument("umat: condense='exact' runs under kinematic_form='rate' (log_R) or 'total'");
		}
		if (cs.active() && !exact && tangent_out != tangent_out_box) {
			throw std::invalid_argument("umat: condensed stress states (ndi < 3) return the box tangent (tangent_output='box'), "
			                            "or dS/dE under condense='exact' (tangent_output='material')");
		}
		if (exact && tangent_out != tangent_out_material) {
			throw std::invalid_argument("umat: condense='exact' returns dS/dE condensed on the kept components: "
			                            "tangent_output='material'");
		}
		// the solved strain written back into Detot: the free components, or the whole law increment
		const uvec wb_rows = exact ? arma::regspace<uvec>(0, 5) : cs.C;
		mat Detot_W = cs.active() ? mat(list_Detot.rows(wb_rows)) : mat();
		if (tangent_out != tangent_out_box) {
			if (F1_py.ndim() != 3 || F1_py.shape(2) != nb_points) {
				throw std::invalid_argument("umat: tangent_output='" + tangent_output
				                            + "' needs F1 with one 3x3 slice per material point");
			}
		}
		// F0/F1 are converted once: a strict alias of the numpy buffers (see
		// numpy_to_arma.hpp), never re-assigned. F0 is optional for tangent_output alone.
		const bool need_F0 = arguments_type == 2 || kirchhoff_normalize;
		if (need_F0 || tangent_out != tangent_out_box) {
			if (need_F0 || F0_py.size() > 0) {
				F0 = simpy::numpy_to_arma::arr_to_cube_view(F0_py);
			}
			F1 = simpy::numpy_to_arma::arr_to_cube_view(F1_py);
		}
		if (kirchhoff_normalize) {
			// loud on a shape mismatch: silently skipping would return
			// Kirchhoff under the documented Cauchy contract
			if (F0.n_slices != (arma::uword)nb_points
					|| F1.n_slices != (arma::uword)nb_points) {
				throw std::invalid_argument(
					"umat: F0/F1 must carry one 3x3 slice per material point when provided");
			}
		}

		// Step-cut request of each point (tnew_dt < 1). One slot per point, sized here in serial
		// context: no shared write, and no NumPy-backed allocation, in the parallel region.
		std::vector<double> tnew_dt(nb_points, 1.);
		// log-corate work correction (finite-strain calls only): decided once, not per point
		const bool work_correction = work_correction_on && kirchhoff_normalize
		                             && simcoon::work_correction_applies(corate_type);
		const mat I3 = eye(3, 3);
		auto point_kernel = [&](int pt) {
			// Alias the props column without copying: the parallel region makes no
			// allocation for it, and never touches Python (no GIL needed by workers).
			// props (unique) / list_props (per-point) outlive the lambda and are read-only.
			const double* _props_ptr = unique_props ? props.memptr() : list_props.colptr(pt);
			const vec local_props(const_cast<double*>(_props_ptr), nprops, false, true);
			vec statev = list_statev.unsafe_col(pt);
			vec sigma = list_sigma.unsafe_col(pt);

			vec etot = list_etot.unsafe_col(pt);
			vec Detot = list_Detot.unsafe_col(pt);
			vec Wm = list_Wm.unsafe_col(pt);

			const double T = use_temp ? vec_T(vec_T.n_elem == 1 ? 0 : pt) : 0.0;
			const double DT = vec_DT.n_elem ? vec_DT(vec_DT.n_elem == 1 ? 0 : pt) : 0.0;
			double tnew_dt_pt = 1.;   // per-point: no shared write inside the parallel region

			// F0_pt, F1_pt alias the caller's slices unless condensed: a condensed finite point
			// takes its free stretches from the strain state (condensed_F); F1 follows the solve
			const bool have_F0 = F0.n_slices == static_cast<arma::uword>(nb_points);
			const bool have_F1 = F1.n_slices == static_cast<arma::uword>(nb_points);
			mat F0_pt = have_F0 ? mat(F0.slice(pt).memptr(), 3, 3, cs.active(), true) : mat();
			mat F1_pt = have_F1 ? mat(F1.slice(pt).memptr(), 3, 3, cs.active(), true) : mat();
			const bool rebuild_F = cs.active() && have_F0 && have_F1
			                       && (total || arguments_type == 2 || !arma::approx_equal(F0.slice(pt), F1.slice(pt), "absdiff", 0.));
			std::optional<condensed_F> F1_of;
			auto F_at = [&](const condensed_F &Fc, const vec &e) { return total ? Fc.total(e) : Fc(cs, e); };
			if (rebuild_F) {
				try {   // no throw in the parallel region: a degenerate F is a step-cut request
					F0_pt = F_at(condensed_F(F0.slice(pt)), etot);
					F1_of.emplace(F1.slice(pt));
				}
				catch (const std::exception &) {
					tnew_dt[pt] = 0.5;
					return;
				}
			}

			bool total_cut = false;
			if (kirchhoff_normalize) {
				// python contract stress (Cauchy) -> kernel internal
				// (Kirchhoff). Degenerate F (e.g. legacy zero-filled
				// placeholders for the positional F arguments) means no
				// meaningful finite-strain state: pass through unscaled
				// (previous behavior) instead of producing 0/NaN — no throw
				// here, this runs inside the parallel region.
				const double J0 = arma::det(F0_pt);
				if (J0 > simcoon::iota) sigma *= J0;
				if (total) {
					// Kirchhoff -> dual stress T. No throw in the parallel region: a degenerate F
					// becomes a step-cut request, raised as StepCut in serial context below.
					try { sigma = simcoon::t2v_stress(simcoon::Kirchoff2LogStressU(simcoon::v2t_stress(sigma), F0_pt)); }
					catch (const std::exception &) { total_cut = true; }
				}
			}
			// the total formulation has no frame increment, whatever DR the caller holds
			mat DR_pt = total ? I3 : mat(DR.slice(pt));
			double tau_start[6];
			for (int k = 0; k < 6; k++) tau_start[k] = sigma(k);
			// the law integrates Detot_pt: Detot (the caller's view) except on the stress-free
			// components, solved below; plain aliases when nothing is condensed
			vec Detot_pt(Detot.memptr(), Detot.n_elem, cs.active(), true);
			const int ndi_law = cs.active() ? 3 : ndi, nshr_law = cs.active() ? 3 : nshr;
			auto run_law = [&](double &tnew) {
				switch (arguments_type) {
					case 1: {
						umat_function(umat_name_py, etot, Detot_pt, sigma, Lt.slice(pt), L.slice(pt), DR_pt, nprops, local_props, nstatev, statev, T, DT, Time, DTime, Wm(0), Wm(1), Wm(2), Wm(3), ndi_law, nshr_law, start, tnew, tangent_mode);
						break;
					}
					case 3: {
						simcoon::umat_modular(umat_name_py, etot, Detot_pt, sigma, Lt.slice(pt), L.slice(pt), DR_pt, nprops, local_props, nstatev, statev, T, DT, Time, DTime, Wm(0), Wm(1), Wm(2), Wm(3), ndi_law, nshr_law, start, tnew, tangent_mode, total ? simcoon::KinematicForm::total : simcoon::KinematicForm::rate);
						break;
					}
					case 2: {
						umat_function_finite(umat_name_py, etot, Detot_pt, F0_pt, F1_pt, sigma, Lt.slice(pt), L.slice(pt), DR_pt, nprops, local_props, nstatev, statev, T, DT, Time, DTime, Wm(0), Wm(1), Wm(2), Wm(3), ndi_law, nshr_law, start, tnew, corate_type, tangent_mode);
						break;
					}
				}
			};
			if (!cs.active()) {
				run_law(tnew_dt_pt);
			}
			else {
				// Newton on the free strain components (solve_condensed), the law restored to its start
				// state before every evaluation; no throw out of the parallel region: a failure is a
				// step-cut request
				const vec sigma0 = sigma, statev0 = statev, Wm0 = Wm;
				auto restore = [&]() { sigma = sigma0; statev = statev0; Wm = Wm0; };
				try {
					if (!exact) {
						// the stress the law returns vanishes on C, unknowns Detot_pt(C). Under finite
						// strain F1 follows: an F-based law reads it at every evaluation (and under the
						// total form returns tau, with dT/dH: its residual is T at that F1), a box law
						// only through J, rebuilt once at the end
						const bool T_residual = total && arguments_type == 2;
						vec res_T;
						auto run = [&](double &tnew) {
							if (rebuild_F && arguments_type == 2) F1_pt = F_at(*F1_of, etot + Detot_pt);
							run_law(tnew);
							if (T_residual) res_T = simcoon::t2v_stress(simcoon::Kirchoff2LogStressU(simcoon::v2t_stress(sigma), F1_pt));
						};
						if (simcoon::solve_condensed(cs, Detot_pt, T_residual ? res_T : sigma, Lt.slice(pt), run, restore, tnew_dt_pt)) {
							simcoon::condense_tangent(cs, Lt.slice(pt));
							Detot_W.col(pt) = Detot_pt(wb_rows);
							if (rebuild_F) F1_pt = F_at(*F1_of, etot + Detot_pt);
						}
					}
					else {
						// S_C = 0, unknowns the Green-Lagrange E1(C): E1(K) from the caller's F1,
						// F1 = R1 U1(E1), and the law fed the increment of its own strain from F0 to
						// that F1 (H1 - etot, or the log_R increment and its rotation)
						vec E1 = simcoon::t2v_strain(simcoon::Green_Lagrange(F_at(*F1_of, etot + Detot_pt)));
						E1(cs.K) = simcoon::t2v_strain(simcoon::Green_Lagrange(F1.slice(pt)))(cs.K);
						vec S;
						mat M;
						const double dt = (DTime > simcoon::iota) ? DTime : 1.;
						auto run = [&](double &tnew) {
							const mat C1 = arma::eye(3, 3) + 2.*simcoon::v2t_strain(E1);
							F1_pt = F1_of->R*arma::sqrtmat_sympd(C1);
							if (total) Detot_pt = simcoon::t2v_strain(0.5*arma::logmat_sympd(C1)) - etot;
							else {
								mat D, Omega;
								simcoon::corate_kinematics(3, DR_pt, D, Omega, F0_pt, F1_pt, dt);
								Detot_pt = simcoon::t2v_strain(simcoon::Delta_log_strain_corate(F0_pt, F1_pt, DR_pt, D, Omega, dt, 3));
							}
							run_law(tnew);
							mat tau = simcoon::v2t_stress(sigma);
							if (total_box) tau = simcoon::LogStressU2Kirchoff(tau, F1_pt);
							const mat Fi = arma::inv(F1_pt);
							S = simcoon::t2v_stress(Fi*tau*Fi.t());
							M = simcoon::DtauDe_corate_2_DSDE(Lt.slice(pt), corate_type, F1_pt, tau);
						};
						if (simcoon::solve_condensed(cs, E1, S, M, run, restore, tnew_dt_pt)) {
							simcoon::condense_tangent(cs, M);
							Lt.slice(pt) = M;
							Detot_W.col(pt) = Detot_pt;
						}
					}
				}
				catch (const std::exception &) {
					restore();
					tnew_dt_pt = std::min(tnew_dt_pt, 0.5);
				}
			}
			if (total) {   // back to Kirchhoff, before the /J below
				try {
					double dW_T;
					sigma = simcoon::total_formulation_tau(sigma, Detot_pt, F1_pt, total_box, dW_T);
					Wm(0) += dW_T;
					Wm(1) += dW_T;
				}
				catch (const std::exception &) { total_cut = true; }
				if (total_cut) tnew_dt_pt = std::min(tnew_dt_pt, 0.5);
			}
			tnew_dt[pt] = tnew_dt_pt;   // own slot: no shared write in the parallel region
			if (work_correction && !arma::approx_equal(F0_pt, F1_pt, "absdiff", 0.)) {
				// true work under the log corates, as select_umat_M_finite. Only when F0 -> F1 is the
				// increment: identical F (small-strain use with placeholder F) carries no D, and the
				// correction would cancel the kernel's work.
				const arma::vec::fixed<6> ts(tau_start);
				const double dW = simcoon::Delta_work_conjugacy(lab_start_stress(ts, DR_pt), ts, sigma,
				                                                Detot_pt, F0_pt, F1_pt, corate_type);
				Wm(0) += dW;
				Wm(1) += dW;
			}
			if (kirchhoff_normalize) {
				// kernel internal (Kirchhoff) -> python contract (Cauchy);
				// Lt is deliberately NOT rescaled (see the block above).
				// Same degenerate-F passthrough as the input side.
				const double J1 = arma::det(F1_pt);
				if (J1 > simcoon::iota) sigma /= J1;
			}
			if (tangent_out != tangent_out_box && !exact) {   // exact: already dS/dE
				// same maps as Lt_convert on the returned (Cauchy) stress, with the corate the law ran with
				const mat dSdE = simcoon::DtauDe_corate_2_DSDE(Lt.slice(pt), corate_type, F1.slice(pt),
				                                               arma::det(F1.slice(pt))*simcoon::v2t_stress(sigma));
				Lt.slice(pt) = (tangent_out == tangent_out_material) ? dSdE : simcoon::DSDE_2_Dsigma_LieDD(dSdE, F1.slice(pt));
			}
		};
		if (serial) {
			// PYEXT calls back into Python: it must stay on the calling thread (which holds the
			// GIL, so the gil_scoped_acquire in the bridge is a no-op). A plain loop, not
			// simcoon_parallel_for_safe with a large cutoff: the OpenMP build of the helper
			// captures the exception of a failing point and keeps calling the law for all the
			// remaining ones, which the GCD build does not — here the first error must stop
			// the batch on every platform.
			for (int pt = 0; pt < nb_points; pt++) {
				point_kernel(pt);
				if (tnew_dt[pt] < 1.) raise_step_cut("umat", tnew_dt);
			}
		} else {
			parallel_for_nogil(nb_points, point_kernel, n_threads);
		}
		// A built-in kernel asks for a smaller increment through tnew_dt (e.g. the modular
		// engine on a non-finite or runaway return mapping, statev left untouched): surface it.
		if (std::any_of(tnew_dt.begin(), tnew_dt.end(), [](double r) { return r < 1.; }))
			raise_step_cut("umat", tnew_dt);
		if (cs.active()) write_back_free_strain(Detot_py, wb_rows, Detot_W);
		return py::make_tuple(simpy::arma_to_numpy::mat_to_arr(list_sigma, false), simpy::arma_to_numpy::mat_to_arr(list_statev, false), simpy::arma_to_numpy::mat_to_arr(list_Wm, false), simpy::arma_to_numpy::cube_to_arr(Lt, false));

	}

	py::tuple launch_umat_T(const std::string &umat_name_py, const py::array_t<double> &etot_py, const py::array_t<double> &Detot_py, const py::array_t<double> &sigma_py, const py::array_t<double> &DR_py, const py::array_t<double> &props_py, const py::array_t<double> &statev_py, const double Time, const double DTime, const py::array_t<double> &Wm_py, const py::array_t<double> &Wt_py, const py::array_t<double> &T_py, const py::array_t<double> &DT_py, const int &ndi, const std::optional<int> &nshr_py, const unsigned int &n_threads, const int &tangent_mode, const std::optional<bool> &start_py, const std::optional<py::array_t<double>> &F0_py, const std::optional<py::array_t<double>> &F1_py, const py::object &kinematic_form_py, const std::string &condense, const std::string &tangent_output){
		// Point-wise thermomechanical UMAT batch entry, mirroring launch_umat.
		// Dispatch follows the select_umat_T table (umat_smart.cpp).
		// Returns (sigma, statev, Wm, Wt, r, dSdE, dSdT, drdE, drdT).

		if (tangent_mode < simcoon::tangent_none || tangent_mode > simcoon::tangent_closest_point) {
			throw std::invalid_argument("tangent_mode must be 0 (none), 1 (continuum), 2 (algorithmic) or 3 (closest-point); got "
			                            + std::to_string(tangent_mode));
		}
		static const std::map<string, int> list_umat = { {"ELISO",1},{"ELIST",2},{"ELORT",3},{"EPICP",4},{"EPKCP",5},{"ZENER",6},{"ZENNK",7},{"PRONK",8},{"SMADI",9},{"SMADC",9},{"SMAAI",9},{"SMAAC",9},{"EPJCK",10} };
		auto it_umat = list_umat.find(umat_name_py);
		if (it_umat == list_umat.end()) {
			throw std::invalid_argument("The choice of thermomechanical Umat could not be found in the umat library: " + umat_name_py);
		}
		int id_umat = it_umat->second;
		int arguments_type; //depends on the argument used in the umat

		// Unified thermomechanical function pointer: (Etot, DEtot, sigma, r, dSdE, dSdT, drdE, drdT, DR, nprops, props, nstatev, statev, T, DT, Time, DTime, Wm, Wm_r, Wm_ir, Wm_d, Wt, Wt_r, Wt_ir, ndi, nshr, start, tnew_dt, tangent_mode)
		void (*umat_function)(const arma::vec &, const arma::vec &, arma::vec &, double &, arma::mat &, arma::mat &, arma::mat &, arma::mat &, const arma::mat &, const int &, const arma::vec &, const int &, arma::vec &, const double &, const double &, const double &, const double &, double &, double &, double &, double &, double &, double &, double &, const int &, const int &, const bool &, double &, const int &);
		// SMA family variant carrying the umat_name as leading argument
		void (*umat_function_named)(const std::string &, const arma::vec &, const arma::vec &, arma::vec &, double &, arma::mat &, arma::mat &, arma::mat &, arma::mat &, const arma::mat &, const int &, const arma::vec &, const int &, arma::vec &, const double &, const double &, const double &, const double &, double &, double &, double &, double &, double &, double &, double &, const int &, const int &, const bool &, double &, const int &);

		const int ncomp = 6;
		int nshr;
		const simcoon::condensed_state cs = batch_stress_state(ndi, nshr_py, nshr);

		// start re-initialises the point (T_init, stress, internal variables, Wm):
		// the caller's choice when given, otherwise inferred from Time.
		const bool start = start_py.value_or(Time <= simcoon::limit);
		mat list_etot = simpy::numpy_to_arma::arr_to_mat_view(etot_py);
		unsigned int nb_points = list_etot.n_cols; //number of material points

		// Finite strain (F0, F1): every thermomechanical law, fed the log_R increment (rate) or
		// ln U (total); as launch_umat, the kernels run in the lab frame with DR
		if (F0_py.has_value() != F1_py.has_value())
			throw std::invalid_argument("umat_T: give both F0 and F1 (finite strain) or neither (small strain)");
		const bool finite = F0_py.has_value();
		const int form = kinematic_code(kinematic_form_py);
		if (!finite && form >= 0)
			throw std::invalid_argument("umat_T: kinematic_form needs F0 and F1");
		const bool total = (form == simcoon::corate_type_total);
		cube F0, F1;
		if (finite) {
			if (form >= 0 && form != 3 && !total)
				throw std::invalid_argument("umat_T: finite strain serves kinematic_form='rate' (log_R, corate 3) "
				                            "or 'total'; got code " + std::to_string(form));
			F0 = simpy::numpy_to_arma::arr_to_cube_view(F0_py.value());
			F1 = simpy::numpy_to_arma::arr_to_cube_view(F1_py.value());
			if (F0.n_rows != 3 || F0.n_cols != 3 || F0.n_slices != nb_points
			    || F1.n_rows != 3 || F1.n_cols != 3 || F1.n_slices != nb_points)
				throw std::invalid_argument("umat_T: F0/F1 must carry one 3x3 slice per material point");
		}
		// as launch_umat; under 'exact' dSdE, dSdT and drdE come back as the material derivatives
		// dS/dE, dS/dtheta, dr/dE condensed on the kept components
		const bool exact = cs.active() && condense_exact(condense) && finite;
		// tangent_output, as launch_umat: the box derivatives as the kernels return them, or the
		// material (TL) / spatial (UL) ones of the stress and the heat, converted per point
		enum { tangent_out_box, tangent_out_material, tangent_out_spatial };
		static const std::map<string, int> list_tangent_output = { {"box", tangent_out_box}, {"material", tangent_out_material}, {"spatial", tangent_out_spatial} };
		const auto it_tangent_out = list_tangent_output.find(tangent_output);
		if (it_tangent_out == list_tangent_output.end())
			throw std::invalid_argument("tangent_output must be 'box', 'material' or 'spatial'; got '" + tangent_output + "'");
		const int tangent_out = it_tangent_out->second;
		if (tangent_out != tangent_out_box && !finite)
			throw std::invalid_argument("umat_T: tangent_output='" + tangent_output + "' needs F0 and F1");
		if (cs.active() && !exact && tangent_out != tangent_out_box)
			throw std::invalid_argument("umat_T: condensed stress states (ndi < 3) return the box derivatives (tangent_output='box'), "
			                            "or the material ones under condense='exact' (tangent_output='material')");
		if (exact && tangent_out != tangent_out_material)
			throw std::invalid_argument("umat_T: condense='exact' returns dS/dE, dS/dtheta and dr/dE condensed on the kept "
			                            "components: tangent_output='material'");
		std::vector<double> tnew_dt(nb_points, 1.);   // step-cut request, one slot per point
		mat list_Detot = simpy::numpy_to_arma::arr_to_mat_view(Detot_py);
		require_writable_Detot(Detot_py, cs, "umat_T");
		const uvec wb_rows = exact ? arma::regspace<uvec>(0, 5) : cs.C;   // as launch_umat
		mat Detot_W = cs.active() ? mat(list_Detot.rows(wb_rows)) : mat();
		mat list_sigma = simpy::numpy_to_arma::arr_to_mat(std::move(sigma_py)); //copy: modified by the umat and returned
		cube DR = simpy::numpy_to_arma::arr_to_cube_view(DR_py);
		vec vec_T = simpy::numpy_to_arma::arr_to_col_view(T_py);
		vec vec_DT = simpy::numpy_to_arma::arr_to_col_view(DT_py);

		vec props;
		//n_cols (not the raw numpy shape) so a 1-D (nprops,) array is a valid single-props input
		mat list_props = simpy::numpy_to_arma::arr_to_mat_view(props_py);
		bool unique_props = false;
		if (list_props.n_cols == 1) {
			props = list_props.col(0);
			unique_props = true;
		}
		else if (list_props.n_cols != nb_points) {
			throw std::invalid_argument("umat_T: props must have 1 column (shared) or one column per material point; got "
			                            + std::to_string(list_props.n_cols) + " columns for " + std::to_string(nb_points) + " points");
		}

		mat list_statev = simpy::numpy_to_arma::arr_to_mat(std::move(statev_py)); //copy: modified by the umat and returned
		mat list_Wm = simpy::numpy_to_arma::arr_to_mat(std::move(Wm_py)); //copy: modified by the umat and returned
		mat list_Wt = simpy::numpy_to_arma::arr_to_mat(std::move(Wt_py)); //copy: modified by the umat and returned

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
				umat_function = &simcoon::umat_plasticity_iso_T;
				arguments_type = 1;
				break;
			}
			case 10: {
				umat_function = &simcoon::umat_plasticity_johnson_cook_CCP_T;
				arguments_type = 1;
				break;
			}
			case 5: {
				umat_function = &simcoon::umat_plasticity_kin_iso_T;
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

		const mat I3 = eye(3, 3);
		auto thermal_point_kernel = [&](int pt) {
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
			double tnew_dt_pt = 1.;

			// as launch_umat: F0_pt, F1_pt alias the caller's slices unless condensed
			mat F0_pt = finite ? mat(F0.slice(pt).memptr(), 3, 3, cs.active(), true) : mat();
			mat F1_pt = finite ? mat(F1.slice(pt).memptr(), 3, 3, cs.active(), true) : mat();
			std::optional<condensed_F> F1_of;
			auto F_at = [&](const condensed_F &Fc, const vec &e) { return total ? Fc.total(e) : Fc(cs, e); };
			if (finite && cs.active()) {
				try {
					F0_pt = F_at(condensed_F(F0.slice(pt)), etot);
					F1_of.emplace(F1.slice(pt));
				}
				catch (const std::exception &) {
					tnew_dt[pt] = 0.5;
					return;
				}
			}

			// finite strain, as launch_umat: Cauchy -> Kirchhoff (x J0) -> T (total form); a
			// degenerate F is a step-cut request (no throw in the parallel region)
			bool finite_cut = false;
			double tau_start[6];
			if (finite) {
				const double J0 = arma::det(F0_pt);
				if (J0 > simcoon::iota) sigma *= J0;
				else finite_cut = true;
				for (int k = 0; k < 6; k++) tau_start[k] = sigma(k);
				if (total) {
					try { sigma = simcoon::t2v_stress(simcoon::Kirchoff2LogStressU(simcoon::v2t_stress(sigma), F0_pt)); }
					catch (const std::exception &) { finite_cut = true; }
				}
			}
			mat DR_pt = total ? I3 : mat(DR.slice(pt));

			// as launch_umat: Detot_pt and F1_pt alias the caller's arrays unless condensed
			vec Detot_pt(Detot.memptr(), Detot.n_elem, cs.active(), true);
			double dlnJ_dT = 0.;   // condensed finite point: J follows theta through the free stretches
			const int ndi_law = cs.active() ? 3 : ndi, nshr_law = cs.active() ? 3 : nshr;
			auto run_law = [&](double &tnew) {
				switch (arguments_type) {
					case 1: {
						umat_function(etot, Detot_pt, sigma, list_r(pt), dSdE.slice(pt), dSdT.slice(pt), drdE.slice(pt), drdT.slice(pt), DR_pt, nprops, local_props, nstatev, statev, T, DT, Time, DTime, Wm(0), Wm(1), Wm(2), Wm(3), Wt(0), Wt(1), Wt(2), ndi_law, nshr_law, start, tnew, tangent_mode);
						break;
					}
					case 2: {
						umat_function_named(umat_name_py, etot, Detot_pt, sigma, list_r(pt), dSdE.slice(pt), dSdT.slice(pt), drdE.slice(pt), drdT.slice(pt), DR_pt, nprops, local_props, nstatev, statev, T, DT, Time, DTime, Wm(0), Wm(1), Wm(2), Wm(3), Wt(0), Wt(1), Wt(2), ndi_law, nshr_law, start, tnew, tangent_mode);
						break;
					}
				}
			};
			if (!cs.active()) {
				run_law(tnew_dt_pt);
			}
			else {
				// as launch_umat; no law here reads F1, rebuilt once at the end for J
				const vec sigma0 = sigma, statev0 = statev, Wm0 = Wm, Wt0 = Wt;
				auto restore = [&]() { sigma = sigma0; statev = statev0; Wm = Wm0; Wt = Wt0; list_r(pt) = 0.; };
				try {   // no throw out of the parallel region: a failure is a step-cut request
					if (!exact) {
						if (simcoon::solve_condensed(cs, Detot_pt, sigma, dSdE.slice(pt), run_law, restore, tnew_dt_pt)) {
							if (finite) {   // d(strain_C)/dtheta = -L_CC^-1 dstress_C/dtheta, its direct part is dlnJ/dtheta
								const vec dC = -arma::solve(mat(dSdE.slice(pt)(cs.C, cs.C)), vec(mat(dSdT.slice(pt)).rows(cs.C)));
								for (arma::uword j = 0; j < cs.C.n_elem; j++)
									if (cs.C(j) < 3) dlnJ_dT += dC(j);
							}
							simcoon::condense_tangent_T(cs, dSdE.slice(pt), dSdT.slice(pt), drdE.slice(pt), drdT.slice(pt));
							Detot_W.col(pt) = Detot_pt(wb_rows);
							if (finite) F1_pt = F_at(*F1_of, etot + Detot_pt);
						}
					}
					else {
						// as launch_umat (S_C = 0 on the Green-Lagrange E1(C)); the box laws here return
						// tau (rate) or T (total)
						vec E1 = simcoon::t2v_strain(simcoon::Green_Lagrange(F_at(*F1_of, etot + Detot_pt)));
						E1(cs.K) = simcoon::t2v_strain(simcoon::Green_Lagrange(F1.slice(pt)))(cs.K);
						vec S;
						mat M;
						const double dt = (DTime > simcoon::iota) ? DTime : 1.;
						auto kirchhoff = [&](const mat &law_stress) { return total ? simcoon::LogStressU2Kirchoff(law_stress, F1_pt) : law_stress; };
						auto run = [&](double &tnew) {
							const mat C1 = arma::eye(3, 3) + 2.*simcoon::v2t_strain(E1);
							F1_pt = F1_of->R*arma::sqrtmat_sympd(C1);
							if (total) Detot_pt = simcoon::t2v_strain(0.5*arma::logmat_sympd(C1)) - etot;
							else {
								mat D, Omega;
								simcoon::corate_kinematics(3, DR_pt, D, Omega, F0_pt, F1_pt, dt);
								Detot_pt = simcoon::t2v_strain(simcoon::Delta_log_strain_corate(F0_pt, F1_pt, DR_pt, D, Omega, dt, 3));
							}
							run_law(tnew);
							const mat tau = kirchhoff(simcoon::v2t_stress(sigma));
							const mat Fi = arma::inv(F1_pt);
							S = simcoon::t2v_stress(Fi*tau*Fi.t());
							M = simcoon::DtauDe_corate_2_DSDE(dSdE.slice(pt), total ? simcoon::corate_type_total : 3, F1_pt, tau);
						};
						if (simcoon::solve_condensed(cs, E1, S, M, run, restore, tnew_dt_pt)) {
							// material derivatives at fixed F1: dS/dtheta = F^-1 dtau/dtheta F^-T,
							// dr/dE = dr/dH dH/dE (rate: dH rotated by R1 into the box)
							const mat Fi = arma::inv(F1_pt);
							mat dSdT_S = simcoon::t2v_stress(Fi*kirchhoff(simcoon::v2t_stress(vec(dSdT.slice(pt))))*Fi.t());
							mat drdE_E = simcoon::box_covector_2_dE(arma::vectorise(drdE.slice(pt)), F1_pt, total ? simcoon::corate_type_total : 3);   // drdE is stored (6, 1)
							simcoon::condense_tangent_T(cs, M, dSdT_S, drdE_E, drdT.slice(pt));
							dSdE.slice(pt) = M;
							dSdT.slice(pt) = dSdT_S;
							drdE.slice(pt) = drdE_E;
							Detot_W.col(pt) = Detot_pt;
						}
					}
				}
				catch (const std::exception &) {
					restore();
					tnew_dt_pt = std::min(tnew_dt_pt, 0.5);
				}
			}
			if (finite) {
				// back to Kirchhoff (total: T and dT/dtheta mapped at fixed F1), then Cauchy;
				// dSdE and drdE box derivatives unless tangent_output says otherwise; r per reference
				// volume except under tangent_output='spatial'
				if (total) {
					try {
						double dW_T;
						sigma = simcoon::total_formulation_tau(sigma, Detot_pt, F1_pt, true, dW_T);
						if (!exact)
							dSdT.slice(pt) = simcoon::t2v_stress(simcoon::LogStressU2Kirchoff(simcoon::v2t_stress(vec(dSdT.slice(pt))), F1_pt));
					}
					catch (const std::exception &) { finite_cut = true; }
				}
				else if (!arma::approx_equal(F0_pt, F1_pt, "absdiff", 0.)) {
					// true work under log_R, as launch_umat
					const arma::vec::fixed<6> ts(tau_start);
					const double dW = simcoon::Delta_work_conjugacy(lab_start_stress(ts, DR_pt), ts, sigma,
					                                                Detot_pt, F0_pt, F1_pt, 3);
					Wm(0) += dW;
					Wm(1) += dW;
				}
				const double J1 = arma::det(F1_pt);
				if (J1 > simcoon::iota) {
					sigma /= J1;
					if (!exact) {   // the Cauchy derivative; exact returns dS/dtheta
						dSdT.slice(pt) /= J1;
						dSdT.slice(pt) -= sigma * dlnJ_dT;
					}
					if (tangent_out != tangent_out_box && !exact) {   // exact: already material
						// at fixed F1: dS/dE by the form's box map, dS/dtheta = F^-1 dtau/dtheta F^-T,
						// dr/dE = P^T dr/dh; spatial: the Lie tangent, dsigma/dtheta and the heat per current
						// volume, r/J, d(r/J)/dd = (F dr/dE F^T - r I)/J, drdT/J
						try {
							const int code = total ? simcoon::corate_type_total : 3;
							const mat dSdE_S = simcoon::DtauDe_corate_2_DSDE(dSdE.slice(pt), code, F1_pt, J1*simcoon::v2t_stress(sigma));
							const vec drdE_E = simcoon::box_covector_2_dE(arma::vectorise(drdE.slice(pt)), F1_pt, code);
							if (tangent_out == tangent_out_material) {
								const mat Fi = arma::inv(F1_pt);
								dSdE.slice(pt) = dSdE_S;
								dSdT.slice(pt) = simcoon::t2v_stress(Fi*(J1*simcoon::v2t_stress(vec(dSdT.slice(pt))))*Fi.t());
								drdE.slice(pt) = drdE_E;
							}
							else {
								dSdE.slice(pt) = simcoon::DSDE_2_Dsigma_LieDD(dSdE_S, F1_pt);
								const vec I = {1., 1., 1., 0., 0., 0.};
								drdE.slice(pt) = (simcoon::t2v_stress(F1_pt*simcoon::v2t_stress(drdE_E)*F1_pt.t()) - list_r(pt)*I)/J1;
								list_r(pt) /= J1;
								drdT.slice(pt) /= J1;
							}
						}
						catch (const std::exception &) { finite_cut = true; }
					}
				}
				else finite_cut = true;
				if (finite_cut) tnew_dt_pt = std::min(tnew_dt_pt, 0.5);
			}
			tnew_dt[pt] = tnew_dt_pt;
		};
		parallel_for_nogil(nb_points, thermal_point_kernel, n_threads);
		if (std::any_of(tnew_dt.begin(), tnew_dt.end(), [](double r) { return r < 1.; }))
			raise_step_cut("umat_T", tnew_dt);
		if (cs.active()) write_back_free_strain(Detot_py, wb_rows, Detot_W);

		// post-loop repacking (serial): dSdT (6,1,N) -> (6,N), drdE (1,6,N) -> (6,N), drdT (1,1,N) -> (N)
		mat dSdT_out(ncomp, nb_points);
		mat drdE_out(ncomp, nb_points);
		vec drdT_out(nb_points);
		for (unsigned int pt = 0; pt < nb_points; pt++) {
			dSdT_out.col(pt) = dSdT.slice(pt);
			drdE_out.col(pt) = drdE.slice(pt);
			drdT_out(pt) = drdT(0, 0, pt);
		}
		return py::make_tuple(simpy::arma_to_numpy::mat_to_arr(list_sigma, true), simpy::arma_to_numpy::mat_to_arr(list_statev, true), simpy::arma_to_numpy::mat_to_arr(list_Wm, true), simpy::arma_to_numpy::mat_to_arr(list_Wt, true), simpy::arma_to_numpy::col_to_arr(list_r, true), simpy::arma_to_numpy::cube_to_arr(dSdE, true), simpy::arma_to_numpy::mat_to_arr(dSdT_out, true), simpy::arma_to_numpy::mat_to_arr(drdE_out, true), simpy::arma_to_numpy::col_to_arr(drdT_out, true));
	}
}
