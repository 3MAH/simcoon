///@file Abaqus_apply_all.cpp
///@brief Abaqus_apply_all : Apply Materials, Steps and periodic Boundary Conditions to a periodic box
///@version 1.0

#include <fstream>
#include <iterator>
#include <armadillo>
#include <CGAL/Simple_cartesian.h>
#include <boost/filesystem.hpp>
#include <simcoon/parameter.hpp>
#include <simcoon/Simulation/Solver/block.hpp>
#include <simcoon/Simulation/Solver/read.hpp>
#include <simcoon/Continuum_Mechanics/Unit_cell/node.hpp>
#include <simcoon/Continuum_Mechanics/Unit_cell/steps.hpp>
#include <simcoon/Continuum_Mechanics/Unit_cell/read.hpp>
#include <simcoon/Continuum_Mechanics/Unit_cell/write.hpp>
#include <simcoon/Continuum_Mechanics/Unit_cell/geom_functions.hpp>
#include <simcoon/Continuum_Mechanics/Unit_cell/cubic_mesh.hpp>
#include <simcoon/Continuum_Mechanics/Unit_cell/interpolate.hpp>

using namespace std;
using namespace arma;
using namespace simcoon;

int main() {
    
    string path_data = "data_aba";
    string path_run = "run_aba";
    string pointsfile = "node_nperio.inp";
    string points_out = "node_perio0.inp";
    string el_file = "element.inp";
    
    //copy the original node file
    string src_file = path_data + "/" + pointsfile;
    string dst_file = path_run + "/" + points_out;
    boost::filesystem::copy_file(src_file,dst_file,boost::filesystem::copy_option::overwrite_if_exists);

    //copy the original element file
    src_file = path_data + "/" + el_file;
    dst_file = path_run + "/" + el_file;
    boost::filesystem::copy_file(src_file,dst_file,boost::filesystem::copy_option::overwrite_if_exists);
    
    string buffer;
    
    //Action 0 : Defining the mesh and the periodic mesh
    
    unsigned int nb_nodes = 0;
    unsigned int loading_type = 2;
    std::vector<Node> nodes;
    read_mesh(nodes, path_data, pointsfile);
    nb_nodes = nodes.size();
    cubic_mesh cm(nodes, pointsfile);
    cm.get_domain();
    //cm.find_pairs();
    cm.construct_lists();
    
    cubic_mesh cm_perio = perio_RVE(cm, nb_nodes);
    cm_perio.construct_lists();
    
    //Action 1 : The materials and sections
    string umat_name;
    string sections_file = "Nsections.dat";
    string sections_out = "mat_sec.inp";
    
    //Sections
    std::vector<section_characteristics> sections;
    read_sections(sections, path_data, sections_file);
    write_sections(sections, path_run, sections_out);
    
    //Action 2 : the Periocidc Boundary Conditions : PBC & Constraint drivers : CDN
    
    string PBC_file_name = "PBC_file_name.inp";
    string TIE_file_name = "TIE_file_name.inp";
    string CDN_file_name = "CDN_file_name.inp";
    write_PBC(cm, nb_nodes, path_run, PBC_file_name);
//    write_TIE(cm, cm_perio, path_run, TIE_file_name);
    write_NonPerio2_CDN(cm, cm_perio, loading_type, path_run, CDN_file_name);
//    write_CDN(cm, path_run, CDN_file_name);
    
    //Action 4 : The steps
    string pathfile = "path.txt";
    bool nlgeom = false;
    
    ///Usefull UMAT variables
    std::vector<block> blocks;  //loading blocks
    
    double T_init = 0.;
    //Read the loading path
    read_path(blocks, T_init, path_data, pathfile);
    
    //Phases
    string steps_file_name = "steps_file_name.inp";
    std::vector<aba_step_meca> aba_steps;
    update_steps(aba_steps, blocks, nlgeom);
    write_steps(aba_steps, T_init, path_run, steps_file_name);
    
    //Finally
    string run_file = "run_aba.inp";
    string run_out = path_run + "/" + run_file;
    std::ofstream aba_head;
    
    aba_head.open(run_out, ios::out);

    aba_head << "*Heading" << "\n";
    aba_head << "**" << "\n";
    
    aba_head << "*INCLUDE, INPUT= " << points_out << "\n";
    aba_head << "*INCLUDE, INPUT= " << el_file << "\n";
    aba_head << "*INCLUDE, INPUT= " << sections_out << "\n";
    aba_head << "*INCLUDE, INPUT= " << PBC_file_name << "\n";
    aba_head << "*INCLUDE, INPUT= " << CDN_file_name << "\n";
    aba_head << "*INCLUDE, INPUT= " << steps_file_name << "\n";
        
    return 0;
}
