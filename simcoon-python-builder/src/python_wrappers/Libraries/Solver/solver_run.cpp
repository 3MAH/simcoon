#include <algorithm>
#include <cmath>
#include <string>
#include <vector>
#include <memory>
#include <utility>
#include <cstring>
#include <stdexcept>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <simcoon/python_wrappers/arma_to_numpy.hpp>
#include <simcoon/python_wrappers/numpy_to_arma.hpp>
#include <armadillo>

#include <simcoon/parameter.hpp>
#include <simcoon/exception.hpp>
#include <simcoon/Simulation/Solver/block.hpp>
#include <simcoon/Simulation/Solver/step.hpp>
#include <simcoon/Simulation/Solver/step_meca.hpp>
#include <simcoon/Simulation/Solver/step_thermomeca.hpp>
#include <simcoon/Simulation/Solver/solver_assembly.hpp>
#include <simcoon/Simulation/Solver/output.hpp>
#include <simcoon/Simulation/Solver/solver_sink.hpp>
#include <simcoon/Continuum_mechanics/Umat/umat_smart.hpp>
#include <simcoon/Simulation/Phase/state_variables_M.hpp>
#include <simcoon/Simulation/Phase/state_variables_T.hpp>
#include <simcoon/python_wrappers/dict_get.hpp>
#include <simcoon/python_wrappers/Libraries/Phase/phases.hpp>
#include <simcoon/python_wrappers/Libraries/Solver/solver_run.hpp>

using namespace std;
using namespace arma;
namespace py = pybind11;

namespace simpy {

namespace {

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
        //vec_of copies: the caller's array is never borrowed nor modified
        vec BC_meca = vec_of(sd["BC_meca"], "solver_run: BC_meca");
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
        const py::array_t<double> BC_w_py = sd["BC_w"].cast<py::array_t<double>>();
        mat BC_w = simpy::numpy_to_arma::arr_to_mat(BC_w_py);   //const-ref overload: a copy
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
        const py::array_t<double> tab_py = sd["tab_data"].cast<py::array_t<double>>();
        sptr->tab_data = simpy::numpy_to_arma::arr_to_mat(tab_py);   //const-ref overload: a copy
    }
    if ((sptr->mode == 3) && (sptr->tab_data.n_rows == 0)) {
        throw std::invalid_argument("solver_run: tabular steps (mode 3) require tab_data");
    }
}

template <typename T>
py::array_t<T> scalars_to_arr(const std::vector<T> &v) {
    py::array_t<T> a(v.size());
    std::memcpy(a.mutable_data(), v.data(), v.size()*sizeof(T));
    return a;
}

//History of one recorded vector or matrix as the (N, n_comp) Fortran-ordered array it is returned
//as (one column per component, a matrix flattened row by row), filled in place. Sized at once for the
//expected count of records (untouched pages cost address space only), doubling from a smaller start
//beyond 1 GiB per field; grown if need be, cut to the records of an aborted run.
class flat_history {
public:
    void set_capacity(const std::size_t &n_records) { capacity_ = n_records; }

    void push(const arma::mat &m) {
        if (n_records_ == 0) {
            n_comp_ = m.n_elem;
            const std::size_t bounded = std::max<std::size_t>(first_rows, max_doubles/std::max<std::size_t>(n_comp_, 1));
            buf_.set_size(std::max<std::size_t>(std::min(capacity_, bounded), 1), n_comp_);
        }
        else if (m.n_elem != n_comp_) {
            throw std::length_error("solver_run: a recorded quantity changed size during the run");
        }
        if (n_records_ == buf_.n_rows) {
            grow();
        }
        const arma::uword ld = buf_.n_rows;
        double *out = buf_.memptr() + n_records_;
        if (m.n_rows == 1 || m.n_cols == 1) {
            const double *in = m.memptr();
            for (arma::uword k = 0; k < n_comp_; k++)
                out[k*ld] = in[k];
        }
        else {
            for (arma::uword a = 0; a < m.n_rows; a++)
                for (arma::uword b = 0; b < m.n_cols; b++, out += ld)
                    *out = m(a, b);
        }
        n_records_++;
    }

    //The (N, n_comp) array (Fortran order): the Python side takes .T of a vector history and
    //reshapes a matrix one to (N, n_rows, n_cols), then (n_rows, n_cols, N), both views
    py::array_t<double> columns() {
        if (n_records_ == 0 || n_comp_ == 0) {
            return simpy::arma_to_numpy::mat_to_arr(arma::mat(n_records_, n_comp_));
        }
        if (n_records_ < buf_.n_rows) {
            buf_.shed_rows(n_records_, buf_.n_rows - 1);
        }
        return simpy::arma_to_numpy::mat_to_arr(std::move(buf_));   // handed over, not copied
    }

private:
    //doubles, but stops at the expected count when it lies in between
    void grow() {
        std::size_t rows = 2*buf_.n_rows;
        if (capacity_ > buf_.n_rows && capacity_ < rows) rows = capacity_;
        arma::mat bigger(rows, n_comp_, arma::fill::none);
        for (arma::uword k = 0; k < n_comp_; k++)
            std::memcpy(bigger.colptr(k), buf_.colptr(k), n_records_*sizeof(double));
        buf_ = std::move(bigger);
    }

    static constexpr std::size_t first_rows = 1024;                ///< smallest first allocation (records)
    static constexpr std::size_t max_doubles = std::size_t(1) << 27;   ///< 1 GiB: larger expected histories double from here
    arma::mat buf_;
    std::size_t capacity_ = 0;
    std::size_t n_comp_ = 0;
    std::size_t n_records_ = 0;
};

//The memory sink of the bindings: the converged state of every output point, as
//simcoon::solver_memory_sink captures it, in flat histories that become the arrays of the
//result dict.
class history_sink : public simcoon::solver_results_sink {
public:
    bool record_tangent = true;
    int sv_type = 0;   ///< 1: mechanical, 2: thermomechanical (set by init())

    explicit history_sink(const std::size_t &n_hint) : n_hint_(n_hint) {}

    void init(simcoon::phase_characteristics &rve) override {
        auto *sv = rve.sptr_sv_global.get();
        if (dynamic_cast<simcoon::state_variables_T *>(sv) != nullptr)
            sv_type = 2;
        else if (dynamic_cast<simcoon::state_variables_M *>(sv) != nullptr)
            sv_type = 1;
        for (auto *v : {&blocks_i, &cycles_i, &steps_i, &incs_i}) v->reserve(n_hint_);
        for (auto *v : {&time, &T, &Q, &r}) v->reserve(n_hint_);
        for (auto *h : {&Etot, &etot, &PKII, &tau, &sigma, &F1, &R, &DR, &statev, &Wm, &Wt,
                        &Lt, &dSdE, &dSdT, &drdE, &drdT}) h->set_capacity(n_hint_);
    }

    void record(simcoon::phase_characteristics &rve, const simcoon::solver_output &, const int &kblock,
                const int &kcycle, const int &kstep, const int &kinc, const double &Time) override {
        simcoon::state_variables *sv = rve.sptr_sv_global.get();
        auto *sv_T = (sv_type == 2) ? static_cast<simcoon::state_variables_T *>(sv) : nullptr;
        auto *sv_M = (sv_type == 1) ? static_cast<simcoon::state_variables_M *>(sv) : nullptr;

        blocks_i.push_back(kblock);
        cycles_i.push_back(kcycle);
        steps_i.push_back(kstep);
        incs_i.push_back(kinc);
        time.push_back(Time);

        Etot.push(sv->Etot);
        etot.push(sv->etot);
        PKII.push(sv->PKII);
        tau.push(sv->tau);
        // Cauchy stress, as solver_memory_sink (o_stress_type 4 of phase_characteristics::output)
        sigma.push(sv->tau / det(sv->F1));
        F1.push(sv->F1);
        R.push(sv->R);
        DR.push(sv->DR);
        T.push_back(sv->T);
        statev.push(sv->statev);

        if (sv_T != nullptr) {
            Wm.push(sv_T->Wm);
            Wt.push(sv_T->Wt);
            Q.push_back(sv_T->Q);
            r.push_back(sv_T->r);
            if (record_tangent) {
                dSdE.push(sv_T->dSdE);
                dSdT.push(sv_T->dSdT);
                drdE.push(sv_T->drdE);
                drdT.push(arma::vec{sv_T->drdT});
            }
        }
        else if (sv_M != nullptr) {
            Wm.push(sv_M->Wm);
            if (record_tangent) {
                Lt.push(sv_M->Lt);
            }
        }
    }

    std::vector<int> blocks_i, cycles_i, steps_i, incs_i;
    std::vector<double> time, T, Q, r;
    flat_history Etot, etot, PKII, tau, sigma, F1, R, DR, statev, Wm, Wt, Lt, dSdE, dSdT, drdE, drdT;

private:
    std::size_t n_hint_;
};

//Records of a complete run, every increment being an output point here (as step::generate()
//counts the increments): the capacity of the histories
std::size_t expected_records(const std::vector<simcoon::block> &blocks) {
    std::size_t n = 0;
    for (const auto &bl : blocks) {
        std::size_t n_inc = 0;
        for (const auto &sptr : bl.steps) {
            if (sptr->mode == 3) n_inc += sptr->tab_data.n_rows;
            else if (sptr->Dn_inc > 0.) n_inc += static_cast<std::size_t>(std::round(1./sptr->Dn_inc));
        }
        n += n_inc*bl.ncycle;
    }
    return n;
}

} //anonymous namespace

py::dict solver_run(const py::list &blocks_py, const double &T_init,
                    const std::string &umat_name, const py::array_t<double> &props_py,
                    const int &nstatev, const py::object &orientation,
                    const int &solver_type, const int &corate_type,
                    const py::dict &params_py, const bool &record_tangent,
                    const py::object &phases) {

    if (nstatev < 0) {   //the library takes an unsigned: -1 would ask for 4 billion state variables
        throw std::invalid_argument("solver_run: nstatev = " + std::to_string(nstatev) + " is negative");
    }
    vec props = simpy::numpy_to_arma::arr_to_col(props_py);
    double psi_rve, theta_rve, phi_rve;
    angles_of(orientation, "solver_run, orientation", psi_rve, theta_rve, phi_rve);

    //Numeric solver controls
    simcoon::solver_params ctrl;
    ctrl.div_tnew_dt = dget<double>(params_py, "div_tnew_dt", ctrl.div_tnew_dt);
    ctrl.mul_tnew_dt = dget<double>(params_py, "mul_tnew_dt", ctrl.mul_tnew_dt);
    ctrl.miniter = dget<int>(params_py, "miniter", ctrl.miniter);
    ctrl.maxiter = dget<int>(params_py, "maxiter", ctrl.maxiter);
    ctrl.inforce = dget<int>(params_py, "inforce", ctrl.inforce);
    ctrl.precision = dget<double>(params_py, "precision", ctrl.precision);
    ctrl.precision_rel = dget<double>(params_py, "precision_rel", ctrl.precision_rel);
    if (ctrl.precision_rel < 0.) {
        throw std::invalid_argument("solver_run: precision_rel must be >= 0 (0: absolute tolerance only)");
    }
    //"lambda_solver" is the documented spelling (lambda is a Python keyword); "lambda" kept as alias
    ctrl.lambda = dget<double>(params_py, "lambda_solver", dget<double>(params_py, "lambda", ctrl.lambda));
    // integrator (cutting plane / closest point) and tangent operator, folded into the kernel code
    ctrl.tangent_mode = simcoon::kernel_tangent_mode(dget<int>(params_py, "integrator", simcoon::integrator_cutting_plane),
                                                     dget<int>(params_py, "tangent_mode", ctrl.tangent_mode));
    ctrl.ncp = dget<int>(params_py, "ncp", ctrl.ncp);

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

    history_sink sink(expected_records(blocks));
    sink.record_tangent = record_tangent;

    //Sub-phases of a mean-field model, built HERE, while the GIL is still held: they are read
    //from Python objects. umat_multi checks them against the props of the model (check_sub_phases).
    const std::vector<simcoon::phase_characteristics> sub_phases =
        make_sub_phases(phases, umat_name, T_init);

    int status = 0;
    {
        py::gil_scoped_release nogil;
        status = simcoon::solver_run(blocks, T_init, so, umat_name, props, nstatev,
                                     psi_rve, theta_rve, phi_rve, solver_type, corate_type, ctrl, sink,
                                     sub_phases);
    }

    //Assemble the results (canonical state; measure transforms happen in Python): one row per
    //record, (N, n) in Fortran order, a matrix flattened row by row
    py::dict res;
    res["status"] = status;
    res["sv_type"] = sink.sv_type;
    res["block"] = scalars_to_arr(sink.blocks_i);
    res["cycle"] = scalars_to_arr(sink.cycles_i);
    res["step"] = scalars_to_arr(sink.steps_i);
    res["inc"] = scalars_to_arr(sink.incs_i);
    res["time"] = scalars_to_arr(sink.time);
    res["Etot"] = sink.Etot.columns();
    res["etot"] = sink.etot.columns();
    res["PKII"] = sink.PKII.columns();
    res["tau"] = sink.tau.columns();
    res["sigma"] = sink.sigma.columns();
    res["F1"] = sink.F1.columns();
    res["R"] = sink.R.columns();
    res["DR"] = sink.DR.columns();
    res["T"] = scalars_to_arr(sink.T);
    res["Wm"] = sink.Wm.columns();
    res["statev"] = sink.statev.columns();
    if (sink.sv_type == 2) {
        res["Q"] = scalars_to_arr(sink.Q);
        res["r"] = scalars_to_arr(sink.r);
        res["Wt"] = sink.Wt.columns();
        if (record_tangent) {
            res["dSdE"] = sink.dSdE.columns();
            res["dSdT"] = sink.dSdT.columns();
            res["drdE"] = sink.drdE.columns();
            res["drdT"] = sink.drdT.columns();
        }
    }
    else if (record_tangent) {
        res["Lt"] = sink.Lt.columns();
    }
    return res;
}

} //namespace simpy
