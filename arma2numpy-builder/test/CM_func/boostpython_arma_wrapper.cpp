
#include <armadillo>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <simcoon/arma2numpy/numpy_arma.hpp>
#include <simcoon/arma2numpy/list_vector.hpp>
#include "TCM_func.hpp"

namespace bp = boost::python;
namespace bn = boost::python::numpy;
using namespace std;
using namespace arma;
using namespace arma2numpy;

BOOST_PYTHON_MODULE(TCM_func) {

    Py_Initialize();
    bn::initialize();
    
    // Register the from-python converters
    bp::def("L_iso", L_iso);
    
//    Py_Finalize();
    
}
