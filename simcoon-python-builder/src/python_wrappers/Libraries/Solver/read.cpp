
#include <armadillo>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <simcoon/arma2numpy/numpy_arma.hpp>

#include <simcoon/Simulation/Solver/read.hpp>
#include <simcoon/python_wrappers/Libraries/Solver/read.hpp>
#include <simcoon/python_wrappers/Libraries/Solver/step_meca.hpp>
#include <simcoon/python_wrappers/Libraries/Solver/step_thermomeca.hpp>

namespace bp = boost::python;
namespace bn = boost::python::numpy;
using namespace std;
using namespace arma;
using namespace arma2numpy;

namespace simpy {

//This function reads material properties to prepare a simulation
bp::tuple read_matprops(const std::string &path_data_py, const std::string &materialfile_py) {
    unsigned int nprops;
    unsigned int nstatev;
    double psi_rve;
    double theta_rve;
    double phi_rve;
    vec v;
    string umat_name;
//    string path_data = bp::extract<std::string>(path_data_py);
//    string materialfile = bp::extract<std::string>(materialfile_py);
    simcoon::read_matprops(umat_name, nprops, v, nstatev, psi_rve, theta_rve, phi_rve, path_data_py, materialfile_py);
    return bp::make_tuple(nprops, nstatev, psi_rve, theta_rve, phi_rve, vec2array(v));
}

bp::tuple read_path(const std::string &path_data_py, const std::string &pathfile_py) {

    double T;
    bp::list blocks_py;
    bp::list cycles_per_blocks_py;
    std::vector<simcoon::block> blocks;
    
//    string path_data = bp::extract<std::string>(path_data_py);
//    string materialfile = bp::extract<std::string>(materialfile_py);
    simcoon::read_path(blocks, T, path_data_py, pathfile_py);

    //blocks loop
    for(unsigned int i = 0 ; i < blocks.size() ; i++) {
        bp::list ith_block_py;
        switch(blocks[i].type) {
            case 1: { //Mechanical

                //cycle loop
                for(unsigned int n = 0; n < blocks[i].ncycle; n++) {
                    // Step loop
                    for(unsigned int j = 0; j < blocks[i].nstep; j++) {
                        shared_ptr<simcoon::step_meca> sptr_meca = std::dynamic_pointer_cast<simcoon::step_meca>(blocks[i].steps[j]);
                        step_meca_py stm_py(*sptr_meca);
                        ith_block_py.append(stm_py);
                    }
                }
                blocks_py.append(ith_block_py);
                cycles_per_blocks_py.append(blocks[i].ncycle);
                break;
            }
            case 2: { //Thermomechanical

                //cycle loop
                for(unsigned int n = 0; n < blocks[i].ncycle; n++) {
                    // Step loop
                    for(unsigned int j = 0; j < blocks[i].nstep; j++) {
                        shared_ptr<simcoon::step_thermomeca> sptr_thermomeca = std::dynamic_pointer_cast<simcoon::step_thermomeca>(blocks[i].steps[j]);
                        step_thermomeca_py sttm_py(*sptr_thermomeca);
                        ith_block_py.append(sttm_py);
                    }
                }
                blocks_py.append(ith_block_py);
                cycles_per_blocks_py.append(blocks[i].ncycle);
                break;
            }
            default: {
                cout << "the block type is not defined!\n";
                break;
            }
        }

    }
    return bp::make_tuple(T, cycles_per_blocks_py, blocks_py);
}

        
} //namepsace simpy
