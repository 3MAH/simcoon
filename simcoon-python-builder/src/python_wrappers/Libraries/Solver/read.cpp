
#include <armadillo>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <simcoon/arma2numpy/numpy_arma.hpp>

#include <simcoon/Simulation/Solver/read.hpp>
#include <simcoon/python_wrappers/Libraries/Solver/read.hpp>

namespace bp = boost::python;
namespace bn = boost::python::numpy;
using namespace std;
using namespace arma;
using namespace arma2numpy;

namespace simpy {

//This function reads material properties to prepare a simulation
void read_matprops(int &nprops, bn::ndarray &props, int &nstatev, double &psi_rve, double &theta_rve, double &phi_rve, const bp::str &path_data_py, const bp::str &materialfile_py) {
    vec v = array2vec(props);
    string umat_name;
    string path_data = bp::extract<std::string>(path_data_py);
    string materialfile = bp::extract<std::string>(materialfile_py);
    simcoon::read_matprops(umat_name, nprops, v, nstatev, psi_rve, theta_rve, phi_rve, path_data, materialfile);
    props = vec2array(v);
}
        
} //namepsace simpy