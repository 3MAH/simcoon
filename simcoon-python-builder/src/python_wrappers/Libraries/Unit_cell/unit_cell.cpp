#include <assert.h>
#include <armadillo>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <simcoon/arma2numpy/numpy_arma.hpp>
#include <simcoon/arma2numpy/numpy_cgal.hpp>
#include <simcoon/arma2numpy/list_vector.hpp>
#include <iostream>
#include <fstream>
#include <CGAL/Simple_cartesian.h>
#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/node.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/section_characteristics.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/read.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/write.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/geom_functions.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/cubic_mesh.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/step_meca.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/step_thermomeca.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/interpolate.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/materials.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/interpolate.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/equation.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/component.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/cubic_equation.hpp>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;

using namespace std;
using namespace arma;
using namespace arma2numpy;
using namespace cgal2numpy;
namespace bp = boost::python;
namespace bn = boost::python::numpy;

namespace simpy {

//-------------------------------------------------------------
simcoon::Node build_node(const int &pynumber, const bn::ndarray &pycoords)
//-------------------------------------------------------------
{
    simcoon::Node n;
    
    n.number = pynumber;
    n.coords = array2Point(pycoords);
    return n;
}

//------------------------------------------------------
bn::ndarray Node_get_input_coords(simcoon::Node &n) {
    return Point2array(n.coords);
}
//------------------------------------------------------

//------------------------------------------------------
void Node_set_input_coords(simcoon::Node &self, const bn::ndarray &mcoords_py) {
    self.coords = array2Point(mcoords_py);
}
//------------------------------------------------------

//-------------------------------------------------------------
simcoon::cubic_mesh build_cubic_mesh(const bp::str &Node_list_name_py, const bp::list &pynodes)
//-------------------------------------------------------------
{
    std::vector<simcoon::Node> Node_list;
    Node_list = py_list_to_std_vector_Node(pynodes);
    std::string Node_list_name = bp::extract<std::string>(Node_list_name_py);
    simcoon::cubic_mesh cm(Node_list, Node_list_name);
    return cm;
}

//-------------------------------------------------------------
void get_domain(simcoon::cubic_mesh &self)
//-------------------------------------------------------------
{
    self.get_domain();
}

//-------------------------------------------------------------
void construct_lists(simcoon::cubic_mesh &self)
//-------------------------------------------------------------
{
    self.construct_lists();
}

bp::list read_mesh(const bp::str &path_data_py, const bp::str &inputfile_py)
{
  std::vector<simcoon::Node> nd;
  std::string path_data= bp::extract<std::string>(path_data_py);
  std::string inputfile= bp::extract<std::string>(inputfile_py);
  simcoon::read_mesh(nd,path_data,inputfile);
  return std_vector_to_py_list_Node(nd);
}

bp::list read_sections(const int &loading_type, const bp::str &path_data_py, const bp::str &inputfile_py)
{
  std::vector<simcoon::section_characteristics> sections;
  std::string path_data= bp::extract<std::string>(path_data_py);
  std::string inputfile= bp::extract<std::string>(inputfile_py);
  simcoon::read_sections(sections,loading_type,path_data,inputfile);
  return std_vector_to_py_list_section_characteristics(sections);
}

/*bp::list read_path(double &T, const bp::str &path_data_py, const bp::str &pathfile_py)
{
  std::vector<simcoon::block> blocks;
  std::string path_data= bp::extract<std::string>(path_data_py);
  std::string pathfile= bp::extract<std::string>(pathfile_py);
  simcoon::read_path(blocks,T,path_data,pathfile);
  return std_vector_to_py_list_block(blocks);
}*/

/*bp::cubic_mesh perio_RVE(bp::cubic_mesh &RVE_py, bp::int &nb_nodes_py)
{
  simcoon::cubic_mesh RVE= bp::extract<simcoon::cubic_mesh>(RVE_py);
  unsigned int nb_nodes= bp::extract<unsigned int>(nb_nodes_py);
  return simcoon::perio_RVE(RVE,nb_nodes);
}

bp::list update_steps(const bp::list &blocks_py, const bool &nlgeom, const int &loading_type, const int &max_temp)
{
  std::vector<std::shared_ptr<simcoon::step> > aba_steps;
  std::vector<simcoon::block> blocks= py_list_to_std_vector_block(blocks_py);
  simcoon::update_steps(aba_steps,blocks,nlgeom,loading_type,max_temp);
  return std_vector_to_py_list_shptr_step(aba_steps);
}

void write_steps(bp::list &aba_steps_py, const int &ldg_type, const double &temp_ini, const bp::str &path_data_py, const bp::str &outputfile_py)
{
  std::vector<std::shared_ptr<simcoon::step> > aba_steps= py_list_to_std_vector_shptr_step(aba_steps_py);
  std::string path_data= bp::extract<std::string>(path_data_py);
  std::string outputfile= bp::extract<std::string>(outputfile_py);
  simcoon::write_steps(aba_steps,ldg_type,temp_ini,path_data,outputfile);
}

void write_sections(bp::section_characteristics &section_rve_py, const bp::int &loading_type_py, const bp::str &path_data_py, const bp::str &outputfile_py)
{
  simcoon::section_characteristics section_rve= bp::extract<simcoon::section_characteristics>(section_rve_py);
  std::string path_data= bp::extract<std::string>(path_data_py);
  std::string outputfile= bp::extract<std::string>(outputfile_py);
  unsigned int loading_type= bp::extract<unsigned int>(loading_type_py);
  simcoon::write_sections(section_rve,loading_type,path_data,outputfile);
}

void write_PBC(const bp::cubic_mesh &cm_py, const bp::int &nb_nodes_py, const bp::str &path_data_py, const bp::str &outputfile_py)
{
  simcoon::cubic_mesh cm= bp::extract<simcoon::cubic_mesh>(cm_py);
  std::string path_data= bp::extract<std::string>(path_data_py);
  std::string outputfile= bp::extract<std::string>(outputfile_py);
  unsigned int nb_nodes= bp::extract<unsigned int>(nb_nodes_py);
  simcoon::write_PBC(cm,nb_nodes,path_data,outputfile);
}

void write_NonPerio_CDN(const bp::cubic_mesh &cm_y, const bp::cubic_mesh &cm_perio_py, const bp::int &loading_type_py, const bp::str &path_data_py, const bp::str &outputfile_py)
{
  simcoon::cubic_mesh cm= bp::extract<simcoon::cubic_mesh>(cm_py);
  simcoon::cubic_mesh cm_perio= bp::extract<simcoon::cubic_mesh>(cm_perio_py);
  std::string path_data= bp::extract<std::string>(path_data_py);
  std::string outputfile= bp::extract<std::string>(outputfile_py);
  unsigned int loading_type= bp::extract<unsigned int>(loading_type_py);
  simcoon::write_NonPerio_CDN(cm,cm_perio,loading_type,path_data,outputfile);
}

void write_TIE(const bp::cubic_mesh &cm_py, const bp::cubic_mesh &cm_perio_py, const bp::str &path_data_py, const bp::str &outputfile_py)
{
  simcoon::cubic_mesh cm= bp::extract<simcoon::cubic_mesh>(cm_py);
  simcoon::cubic_mesh cm_perio= bp::extract<simcoon::cubic_mesh>(cm_perio_py);
  std::string path_data= bp::extract<std::string>(path_data_py);
  std::string outputfile= bp::extract<std::string>(outputfile_py);
  simcoon::write_TIE(cm,cm_perio,loading_type,path_data,outputfile);
}

void write_CDN(const bp::cubic_mesh &cm_py, const bp::str &path_data_py, const bp::str &outputfile_py)
{
simcoon::cubic_mesh cm= bp::extract<simcoon::cubic_mesh>(cm_py);
std::string path_data= bp::extract<std::string>(path_data_py);
std::string outputfile= bp::extract<std::string>(outputfile_py);
simcoon::write_CDN(cm,path_data,outputfile);
}
*/

}
