#include <string>
#include <vector>
#include <memory>
#include <cstring>
#include <stdexcept>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <carma>
#include <armadillo>

#include <simcoon/parameter.hpp>
#include <simcoon/exception.hpp>
#include <simcoon/Simulation/Solver/block.hpp>
#include <simcoon/Simulation/Solver/step.hpp>
#include <simcoon/Simulation/Solver/step_meca.hpp>
#include <simcoon/Simulation/Solver/step_thermomeca.hpp>
#include <simcoon/Simulation/Solver/read.hpp>
#include <simcoon/Simulation/Solver/output.hpp>
#include <simcoon/Simulation/Solver/solver_sink.hpp>
#include <simcoon/python_wrappers/Libraries/Solver/solver_run.hpp>

using namespace std;
using namespace arma;
namespace py = pybind11;

namespace simpy {

namespace {

template <typename T>
T dget(const py::dict &d, const char *key, const T &dflt) {
    return d.contains(key) ? d[key].cast<T>() : dflt;
}

//Fill the fields shared by step_meca and step_thermomeca (identical member names)
template <typename StepPtr>
void fill_step_common(StepPtr &sptr, const py::dict &sd, const unsigned int &control_type,
                      const unsigned int &size_meca, const int &number, const int &block_type,
                      const double &T_init) {

    sptr->number = number;
    sptr->control_type = control_type;
    sptr->mode = sd["mode"].cast<int>();
    if ((sptr->mode < 1) || (sptr->mode > 3)) {
        throw std::invalid_argument("solver_run: step mode must be 1 (linear), 2 (sinusoidal) or 3 (tabular); got " + std::to_string(sptr->mode));
    }

    sptr->Dn_init = dget<double>(sd, "Dn_init", 1.);
    sptr->Dn_mini = dget<double>(sd, "Dn_mini", 1.);
    if (sptr->mode < 3) {
        sptr->Dn_inc = sd["Dn_inc"].cast<double>();
        sptr->BC_Time = sd["time"].cast<double>();
    }

    //Mechanical control flags and targets, in internal Voigt order [11,22,33,12,13,23] (or the 9 F components row-major for control_type 5/6)
    py::sequence cbc = sd["cBC_meca"].cast<py::sequence>();
    if (static_cast<unsigned int>(py::len(cbc)) != size_meca) {
        throw std::invalid_argument("solver_run: cBC_meca must have " + std::to_string(size_meca) + " components for control_type " + std::to_string(control_type));
    }
    for (unsigned int k = 0; k < size_meca; k++) {
        int flag = cbc[k].cast<int>();
        if ((flag < 0) || (flag > 2)) {
            throw std::invalid_argument("solver_run: cBC_meca flags must be 0 (strain), 1 (stress) or 2 (mode-3 zero-hold); got " + std::to_string(flag));
        }
        if ((flag == 2) && (sptr->mode != 3)) {
            throw std::invalid_argument("solver_run: cBC_meca flag 2 (zero-hold) is only valid for tabular steps (mode 3)");
        }
        if ((flag == 1) && (control_type >= 5)) {
            throw std::invalid_argument("solver_run: control types 5 and 6 are fully kinematic; stress control (cBC_meca == 1) is not allowed");
        }
        sptr->cBC_meca(k) = flag;
    }
    if (sd.contains("BC_meca")) {
        vec BC_meca = carma::arr_to_col(sd["BC_meca"].cast<py::array_t<double>>());
        if (BC_meca.n_elem != size_meca) {
            throw std::invalid_argument("solver_run: BC_meca must have " + std::to_string(size_meca) + " components");
        }
        sptr->BC_meca = BC_meca;
    }
    else if (sptr->mode < 3) {
        throw std::invalid_argument("solver_run: BC_meca is required for steps with mode 1 or 2");
    }

    //Rotation rate for the mixed finite-strain control types
    if (sd.contains("BC_w")) {
        mat BC_w = carma::arr_to_mat(sd["BC_w"].cast<py::array_t<double>>());
        if ((BC_w.n_rows != 3) || (BC_w.n_cols != 3)) {
            throw std::invalid_argument("solver_run: BC_w must be a 3x3 matrix");
        }
        sptr->BC_w = BC_w;
    }

    //Thermal boundary condition; the default holds the temperature at T_init
    sptr->cBC_T = dget<int>(sd, "cBC_T", (sptr->mode == 3) ? 2 : 0);
    sptr->BC_T = dget<double>(sd, "BC_T", T_init);
    if (block_type == 1) {
        if ((sptr->cBC_T != 0) && (sptr->cBC_T != 2)) {
            throw std::invalid_argument("solver_run: mechanical blocks only accept cBC_T = 0 (temperature) or 2 (constant, mode 3); got " + std::to_string(sptr->cBC_T));
        }
    }
    else {
        if ((sptr->cBC_T < 0) || (sptr->cBC_T > 3)) {
            throw std::invalid_argument("solver_run: thermomechanical cBC_T must be 0 (temperature), 1 (heat flux), 2 (constant, mode 3) or 3 (convection); got " + std::to_string(sptr->cBC_T));
        }
        if ((sptr->cBC_T == 2) && (sptr->mode != 3)) {
            throw std::invalid_argument("solver_run: cBC_T = 2 (constant temperature) is only valid for tabular steps (mode 3)");
        }
    }

    //In-memory tabular data (mode 3)
    if (sd.contains("tab_data")) {
        if (sptr->mode != 3) {
            throw std::invalid_argument("solver_run: tab_data is only valid for tabular steps (mode 3)");
        }
        sptr->tab_data = carma::arr_to_mat(sd["tab_data"].cast<py::array_t<double>>());
    }
    if ((sptr->mode == 3) && (sptr->tab_data.n_rows == 0)) {
        throw std::invalid_argument("solver_run: tabular steps (mode 3) require tab_data");
    }
}

//one row per record, vecs/mats flattened row-major (matrices reshape to
//(N, n_rows, n_cols) C-order on the Python side). copy=true: small histories can
//fit armadillo's internal pre-allocated buffer, which must never be handed to
//numpy as a zero-copy steal.
template <typename ArmaT>
py::array_t<double> rows_to_arr(const std::vector<ArmaT> &v) {
    arma::mat M(v.size(), v.empty() ? 0 : v[0].n_elem);
    for (size_t i = 0; i < v.size(); i++) {
        M.row(i) = arma::vectorise(v[i].t()).t();
    }
    return carma::mat_to_arr(M, true);
}

template <typename T>
py::array_t<T> scalars_to_arr(const std::vector<T> &v) {
    py::array_t<T> a(v.size());
    std::memcpy(a.mutable_data(), v.data(), v.size()*sizeof(T));
    return a;
}

} //anonymous namespace

py::dict solver_run(const py::list &blocks_py, const double &T_init,
                    const std::string &umat_name, const py::array_t<double> &props_py,
                    const int &nstatev, const double &psi_rve, const double &theta_rve, const double &phi_rve,
                    const int &solver_type, const int &corate_type,
                    const py::dict &params_py, const bool &record_tangent) {

    vec props = carma::arr_to_col(props_py);

    //Numeric solver controls
    simcoon::solver_params ctrl;
    ctrl.div_tnew_dt = dget<double>(params_py, "div_tnew_dt", ctrl.div_tnew_dt);
    ctrl.mul_tnew_dt = dget<double>(params_py, "mul_tnew_dt", ctrl.mul_tnew_dt);
    ctrl.miniter = dget<int>(params_py, "miniter", ctrl.miniter);
    ctrl.maxiter = dget<int>(params_py, "maxiter", ctrl.maxiter);
    ctrl.inforce = dget<int>(params_py, "inforce", ctrl.inforce);
    ctrl.precision = dget<double>(params_py, "precision", ctrl.precision);
    //"lambda_solver" is the documented spelling (lambda is a Python keyword); "lambda" kept as alias
    ctrl.lambda = dget<double>(params_py, "lambda_solver", dget<double>(params_py, "lambda", ctrl.lambda));
    ctrl.tangent_mode = dget<int>(params_py, "tangent_mode", ctrl.tangent_mode);

    //Loading blocks
    unsigned int nblocks = static_cast<unsigned int>(py::len(blocks_py));
    if (nblocks == 0) {
        throw std::invalid_argument("solver_run: at least one loading block is required");
    }
    std::vector<simcoon::block> blocks;
    blocks.resize(nblocks);    //filled in place: the block copy constructor does not carry control_type

    for (unsigned int i = 0; i < nblocks; i++) {
        py::dict bd = blocks_py[i].cast<py::dict>();

        blocks[i].number = i+1;
        blocks[i].type = bd["type"].cast<unsigned int>();
        blocks[i].control_type = bd["control_type"].cast<unsigned int>();
        blocks[i].ncycle = dget<unsigned int>(bd, "ncycle", 1);

        if ((blocks[i].type < 1) || (blocks[i].type > 2)) {
            throw std::invalid_argument("solver_run: block type must be 1 (mechanical) or 2 (thermomechanical); got " + std::to_string(blocks[i].type));
        }
        if (blocks[i].type != blocks[0].type) {
            //the engine constructs the RVE state variables once (first block): a later
            //block of the other type would dereference a null cast and crash
            throw std::invalid_argument("solver_run: all blocks must share the same type (mechanical or thermomechanical); block " + std::to_string(i+1) + " differs from block 1");
        }
        if ((blocks[i].control_type < 1) || (blocks[i].control_type > 6)) {
            throw std::invalid_argument("solver_run: control_type must be in [1, 6]; got " + std::to_string(blocks[i].control_type));
        }
        if ((blocks[i].type == 2) && (blocks[i].control_type != 1)) {
            throw std::invalid_argument("solver_run: thermomechanical blocks only support control_type 1 (small strain)");
        }

        py::list steps_py = bd["steps"].cast<py::list>();
        blocks[i].nstep = static_cast<unsigned int>(py::len(steps_py));
        if (blocks[i].nstep == 0) {
            throw std::invalid_argument("solver_run: block " + std::to_string(i+1) + " has no steps");
        }
        blocks[i].generate();  //allocates step_meca / step_thermomeca with the proper 6/9 sizing

        unsigned int size_meca = (blocks[i].control_type >= 5) ? 9 : 6;

        for (unsigned int j = 0; j < blocks[i].nstep; j++) {
            py::dict sd = steps_py[j].cast<py::dict>();
            if (blocks[i].type == 1) {
                auto sptr = std::dynamic_pointer_cast<simcoon::step_meca>(blocks[i].steps[j]);
                fill_step_common(sptr, sd, blocks[i].control_type, size_meca, j+1, blocks[i].type, T_init);
            }
            else {
                auto sptr = std::dynamic_pointer_cast<simcoon::step_thermomeca>(blocks[i].steps[j]);
                fill_step_common(sptr, sd, blocks[i].control_type, size_meca, j+1, blocks[i].type, T_init);
            }
        }
    }

    //Output cadence: record every increment
    simcoon::solver_output so(nblocks);
    for (unsigned int i = 0; i < nblocks; i++) {
        so.o_type(i) = 1;
        so.o_nfreq(i) = 1;
    }
    simcoon::check_path_output(blocks, so);

    simcoon::solver_memory_sink sink;
    sink.record_tangent = record_tangent;

    int status = 0;
    {
        py::gil_scoped_release nogil;
        status = simcoon::solver_run(blocks, T_init, so, umat_name, props, nstatev,
                                     psi_rve, theta_rve, phi_rve, solver_type, corate_type, ctrl, sink);
    }

    //Assemble the results (canonical state; measure transforms happen in Python)
    py::dict res;
    res["status"] = status;
    res["sv_type"] = sink.sv_type;
    res["block"] = scalars_to_arr(sink.blocks_i);
    res["cycle"] = scalars_to_arr(sink.cycles_i);
    res["step"] = scalars_to_arr(sink.steps_i);
    res["inc"] = scalars_to_arr(sink.incs_i);
    res["time"] = scalars_to_arr(sink.time);
    res["Etot"] = rows_to_arr(sink.Etot);
    res["etot"] = rows_to_arr(sink.etot);
    res["PKII"] = rows_to_arr(sink.PKII);
    res["tau"] = rows_to_arr(sink.tau);
    res["sigma"] = rows_to_arr(sink.sigma);
    res["F1"] = rows_to_arr(sink.F1);
    res["R"] = rows_to_arr(sink.R);
    res["DR"] = rows_to_arr(sink.DR);
    res["T"] = scalars_to_arr(sink.T);
    res["Wm"] = rows_to_arr(sink.Wm);
    res["statev"] = rows_to_arr(sink.statev);
    if (sink.sv_type == 2) {
        res["Q"] = scalars_to_arr(sink.Q);
        res["r"] = scalars_to_arr(sink.r);
        res["Wt"] = rows_to_arr(sink.Wt);
        if (record_tangent) {
            res["dSdE"] = rows_to_arr(sink.dSdE);
            res["dSdT"] = rows_to_arr(sink.dSdT);
            res["drdE"] = rows_to_arr(sink.drdE);
            res["drdT"] = rows_to_arr(sink.drdT);
        }
    }
    else if (record_tangent) {
        res["Lt"] = rows_to_arr(sink.Lt);
    }
    return res;
}

} //namespace simpy
