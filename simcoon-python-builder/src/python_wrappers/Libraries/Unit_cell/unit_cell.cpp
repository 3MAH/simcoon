
#include <armadillo>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <simcoon/arma2numpy/numpy_arma.hpp>
#include <simcoon/arma2numpy/list_vector.hpp>

#include <simcoon/Continuum_Mechanics/Unit_cell/write.hpp>
#include <simcoon/python_wrappers/Libraries/Unit_cell/unit_cell.hpp>

using namespace std;
using namespace arma;
using namespace arma2numpy;
namespace bp = boost::python;
namespace bn = boost::python::numpy;

namespace simpy {
    
//This function reads the essentials information for the loading of a Unit Cell
/*void unit_cell_essentials(int &loading_type, int &BC_type, int &max_temp, const bp::str &path_py, const bp::str &filename_py) {

    string path_data = bp::extract<std::string>(path_py);
    string filename = bp::extract<std::string>(filename_py);
    
    simcoon::unit_cell_essentials(loading_type, BC_type, max_temp, path, filename);
}

void apply_PBC(const bp::str &path_py, const bp::str &inputfile_py) {
    
    std::vector<Node> nodes;
    read_mesh(nodes, path_data, pointsfile);
    nb_nodes = nodes.size();
    unsigned int nb_nodes_init = nb_nodes;
    cubic_mesh cm(nodes, pointsfile);
    cm.get_domain();
    //cm.find_pairs();
    cm.construct_lists();
    
    cubic_mesh cm_perio = perio_RVE(cm, nb_nodes);
    cm_perio.construct_lists();
    
    
    std::vector<Node> nodes;
    string path_data = bp::extract<std::string>(path_py);
    string inputfile = bp::extract<std::string>(inputfile_py);
    
    simcoon::read_mesh(nodes, path_data, inputfile)
    return std_vector_to_py_list_node(nodes);
}*/
    
} //namepsace simpy