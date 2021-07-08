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

#include <simcoon/python_wrappers/Libraries/Unit_cell/unit_cell.hpp>

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
bp::tuple test_mesh(const bn::ndarray &nodes_coords_py, const double &min_dist, const double &dist_replace)
//-------------------------------------------------------------
{
    mat nodes_coords = array2mat(nodes_coords_py, false);
    std::vector<simcoon::Node> nodes_full;

    for (unsigned int i=0; i <nodes_coords.n_cols; i++) {
        vec coords = nodes_coords.col(i);
        simcoon::Node node_temp(i+1,Point(coords(0),coords(1),coords(2)));
        nodes_full.push_back(node_temp);
    }
    unsigned int nb_nodes_full = nodes_full.size();
    
    simcoon::cubic_mesh cm(nodes_full, "nodes_full");
    cm.get_domain();
    cm.construct_lists();
    cm.find_pairs(1.E-6, 1.E-4);
    
    std::vector<simcoon::Node> nodes_all_modified = *cm.Face_listXm;
    nodes_all_modified.insert(nodes_all_modified.end(), cm.Face_listXp->begin(), cm.Face_listXp->end());
    nodes_all_modified.insert(nodes_all_modified.end(), cm.Face_listYm->begin(), cm.Face_listYm->end());
    nodes_all_modified.insert(nodes_all_modified.end(), cm.Face_listYp->begin(), cm.Face_listYp->end());
    nodes_all_modified.insert(nodes_all_modified.end(), cm.Face_listZm->begin(), cm.Face_listZm->end());
    nodes_all_modified.insert(nodes_all_modified.end(), cm.Face_listZp->begin(), cm.Face_listZp->end());

    nodes_all_modified.insert(nodes_all_modified.end(), cm.Edge_listXmYm->begin(), cm.Edge_listXmYm->end());
    nodes_all_modified.insert(nodes_all_modified.end(), cm.Edge_listXpYm->begin(), cm.Edge_listXpYm->end());
    nodes_all_modified.insert(nodes_all_modified.end(), cm.Edge_listXpYp->begin(), cm.Edge_listXpYp->end());
    nodes_all_modified.insert(nodes_all_modified.end(), cm.Edge_listXmYp->begin(), cm.Edge_listXmYp->end());

    nodes_all_modified.insert(nodes_all_modified.end(), cm.Edge_listXmZm->begin(), cm.Edge_listXmZm->end());
    nodes_all_modified.insert(nodes_all_modified.end(), cm.Edge_listXpZm->begin(), cm.Edge_listXpZm->end());
    nodes_all_modified.insert(nodes_all_modified.end(), cm.Edge_listXpZp->begin(), cm.Edge_listXpZp->end());
    nodes_all_modified.insert(nodes_all_modified.end(), cm.Edge_listXmZp->begin(), cm.Edge_listXmZp->end());

    nodes_all_modified.insert(nodes_all_modified.end(), cm.Edge_listYmZm->begin(), cm.Edge_listYmZm->end());
    nodes_all_modified.insert(nodes_all_modified.end(), cm.Edge_listYpZm->begin(), cm.Edge_listYpZm->end());
    nodes_all_modified.insert(nodes_all_modified.end(), cm.Edge_listYpZp->begin(), cm.Edge_listYpZp->end());
    nodes_all_modified.insert(nodes_all_modified.end(), cm.Edge_listYmZp->begin(), cm.Edge_listYmZp->end());
    
    nodes_all_modified.push_back(*cm.Corner_listXmYmZm);
    nodes_all_modified.push_back(*cm.Corner_listXpYmZm);
    nodes_all_modified.push_back(*cm.Corner_listXmYpZm);
    nodes_all_modified.push_back(*cm.Corner_listXpYpZm);
    nodes_all_modified.push_back(*cm.Corner_listXmYmZp);
    nodes_all_modified.push_back(*cm.Corner_listXpYmZp);
    nodes_all_modified.push_back(*cm.Corner_listXmYpZp);
    nodes_all_modified.push_back(*cm.Corner_listXpYpZp);
    
    std::sort(nodes_all_modified.begin(), nodes_all_modified.end());

    mat nodes_coords_modified(nodes_all_modified.size(), 3);
    int i = 0;
    for (auto n:nodes_all_modified) {
        nodes_coords_modified(i,0) = nodes_all_modified[i].coords.x();
        nodes_coords_modified(i,1) = nodes_all_modified[i].coords.y();
        nodes_coords_modified(i,2) = nodes_all_modified[i].coords.z();
        ++i;
    }
    return bp::make_tuple(cm.is_perio, mat2array(nodes_coords_modified, true, "C"));
}

//-------------------------------------------------------------
bp::list build_MPC_from_cubic_mesh(const bn::ndarray &nodes_coords_py)
//-------------------------------------------------------------
{
    
    mat nodes_coords = array2mat(nodes_coords_py, false);
    std::vector<simcoon::Node> nodes_full;

    for (unsigned int i=0; i <nodes_coords.n_cols; i++) {
        vec coords = nodes_coords.col(i);
        simcoon::Node node_temp(i+1,Point(coords(0),coords(1),coords(2)));
        nodes_full.push_back(node_temp);
    }
    unsigned int nb_nodes_full = nodes_full.size();
    
    simcoon::cubic_mesh cm(nodes_full, "nodes_full");
    cm.get_domain();
    cm.construct_lists();
    
    simcoon::cubic_mesh cm_perio = perio_RVE(cm, nb_nodes_full);
    cm_perio.construct_lists();

    //MEchanical only for now
    unsigned int loading_type = 1;
    unsigned int control_type = 1;
    
    simcoon::cubic_equation cubic_eq(cm, cm_perio, loading_type, control_type);
    std::vector<simcoon::equation> MPC_equations = simcoon::MPC_equations_non_perio(cm, cm_perio, cubic_eq, loading_type, control_type);
    
    bp::list MPC_equations_list;
    for(auto eq_it : MPC_equations) {
        cout << eq_it << endl;
        bp::list MPC_equation_temp;
        for(auto comp_it : eq_it.components) {
            simcoon::Node node_comp_it = comp_it.node;
            MPC_equation_temp.append(node_comp_it.number);
            MPC_equation_temp.append(comp_it.dof);
            MPC_equation_temp.append(comp_it.coef);
        }
        MPC_equations_list.append(MPC_equation_temp);
    }
    
    
//    = py_list_to_std_vector<simcoon::equation>(MPC_equations);
    
    return MPC_equations_list;
}

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
simcoon::cubic_mesh build_cubic_mesh(const std::string &Node_list_name_py, const bp::list &pynodes)
//-------------------------------------------------------------
{
    std::vector<simcoon::Node> Node_list;
    Node_list = py_list_to_std_vector_Node(pynodes);
//    std::string Node_list_name = bp::extract<std::string>(Node_list_name_py);
    simcoon::cubic_mesh cm(Node_list, Node_list_name_py);
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

bp::list read_nodes_file(const std::string &path_data_py, const std::string &inputfile_py)
{
  std::vector<simcoon::Node> nd;
//  std::string path_data= bp::extract<std::string>(path_data_py);
//  std::string inputfile= bp::extract<std::string>(inputfile_py);
  simcoon::read_nodes_file(nd,path_data_py,inputfile_py);
  return std_vector_to_py_list_Node(nd);
}

bp::list read_sections(const int &loading_type, const std::string &path_data_py, const std::string &inputfile_py)
{
  std::vector<simcoon::section_characteristics> sections;
//  std::string path_data= bp::extract<std::string>(path_data_py);
//  std::string inputfile= bp::extract<std::string>(inputfile_py);
  simcoon::read_sections(sections,loading_type,path_data_py,inputfile_py);
  return std_vector_to_py_list_section_characteristics(sections);
}

/*bp::list read_path(double &T, const std::string &path_data_py, const std::string &pathfile_py)
{
  std::vector<simcoon::block> blocks;
//  std::string path_data= bp::extract<std::string>(path_data_py);
//  std::string pathfile= bp::extract<std::string>(pathfile_py);
  simcoon::read_path(blocks,T,path_data_py,pathfile_py);
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

void write_steps(bp::list &aba_steps_py, const int &ldg_type, const double &temp_ini, const std::string &path_data_py, const std::string &outputfile_py)
{
  std::vector<std::shared_ptr<simcoon::step> > aba_steps= py_list_to_std_vector_shptr_step(aba_steps_py);
//  std::string path_data= bp::extract<std::string>(path_data_py);
//  std::string outputfile= bp::extract<std::string>(outputfile_py);
  simcoon::write_steps(aba_steps,ldg_type,temp_ini,path_data_py,outputfile_py);
}

void write_sections(bp::section_characteristics &section_rve_py, const bp::int &loading_type_py, const std::string &path_data_py, const std::string &outputfile_py)
{
  simcoon::section_characteristics section_rve= bp::extract<simcoon::section_characteristics>(section_rve_py);
//  std::string path_data= bp::extract<std::string>(path_data_py);
//  std::string outputfile= bp::extract<std::string>(outputfile_py);
  unsigned int loading_type= bp::extract<unsigned int>(loading_type_py);
  simcoon::write_sections(section_rve,loading_type,path_data_py,outputfile_py);
}

void write_PBC(const bp::cubic_mesh &cm_py, const bp::int &nb_nodes_py, const std::string &path_data_py, const std::string &outputfile_py)
{
  simcoon::cubic_mesh cm= bp::extract<simcoon::cubic_mesh>(cm_py);
//  std::string path_data= bp::extract<std::string>(path_data_py);
//  std::string outputfile= bp::extract<std::string>(outputfile_py);
  unsigned int nb_nodes= bp::extract<unsigned int>(nb_nodes_py);
  simcoon::write_PBC(cm,nb_nodes,path_data_py,outputfile_py);
}

void write_NonPerio_CDN(const bp::cubic_mesh &cm_y, const bp::cubic_mesh &cm_perio_py, const bp::int &loading_type_py, const std::string &path_data_py, const std::string &outputfile_py)
{
  simcoon::cubic_mesh cm= bp::extract<simcoon::cubic_mesh>(cm_py);
  simcoon::cubic_mesh cm_perio= bp::extract<simcoon::cubic_mesh>(cm_perio_py);
//  std::string path_data= bp::extract<std::string>(path_data_py);
//  std::string outputfile= bp::extract<std::string>(outputfile_py);
  unsigned int loading_type= bp::extract<unsigned int>(loading_type_py);
  simcoon::write_NonPerio_CDN(cm,cm_perio,loading_type,path_data_py,outputfile_py);
}

void write_TIE(const bp::cubic_mesh &cm_py, const bp::cubic_mesh &cm_perio_py, const std::string &path_data_py, const std::string &outputfile_py)
{
  simcoon::cubic_mesh cm= bp::extract<simcoon::cubic_mesh>(cm_py);
  simcoon::cubic_mesh cm_perio= bp::extract<simcoon::cubic_mesh>(cm_perio_py);
//  std::string path_data= bp::extract<std::string>(path_data_py);
//  std::string outputfile= bp::extract<std::string>(outputfile_py);
  simcoon::write_TIE(cm,cm_perio,loading_type,path_data_py,outputfile_py);
}

void write_CDN(const bp::cubic_mesh &cm_py, const std::string &path_data_py, const std::string &outputfile_py)
{
simcoon::cubic_mesh cm= bp::extract<simcoon::cubic_mesh>(cm_py);
//std::string path_data= bp::extract<std::string>(path_data_py);
//std::string outputfile= bp::extract<std::string>(outputfile_py);
simcoon::write_CDN(cm,path_data_py,outputfile_py);
}
*/

}
