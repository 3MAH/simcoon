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
#include <simcoon/Continuum_mechanics/Unit_cell/node.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/step_meca.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/step_thermomeca.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/read.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/write.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/geom_functions.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/cubic_mesh.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/interpolate.hpp>

using namespace std;
using namespace arma;
using namespace simcoon;

int main() {
    
    string path_data = "data_aba";
    string path_run = "run_aba";
    string pointsfile = "node_nperio.inp";
    string points_out = "node_perio0.inp";
    string el_file = "element.inp";
    string uc_essentials = "unit_cell_essentials.inp";
    string postproc = "postproc_info.dat";
    
    //copy the original node file
    string src_file = path_data + "/" + pointsfile;
    string dst_file = path_run + "/" + points_out;
    boost::filesystem::copy_file(src_file,dst_file,boost::filesystem::copy_option::overwrite_if_exists);

    //copy the original element file
    src_file = path_data + "/" + el_file;
    dst_file = path_run + "/" + el_file;
    boost::filesystem::copy_file(src_file,dst_file,boost::filesystem::copy_option::overwrite_if_exists);
    
    string postproc_file = path_run + "/" + postproc;
    
    string buffer;
    
    //Action 0 : Defining the mesh and the periodic mesh
    
    int nb_nodes = 0;
    int max_temp = 0;
    
    //Action 1 : The materials and sections
    string umat_name;
    string sections_file = "Nsections.dat";
    string sections_out = "mat_sec.inp";
    
    //Action 2 : The steps
    string pathfile = "path.txt";
    bool nlgeom = false;
    double T_init = 0.;
    std::vector<block> blocks;  //loading blocks
    //Read the loading path
    read_path(blocks, T_init, path_data, pathfile);
    unsigned int loading_type = blocks[0].type;
    unsigned int control_type = blocks[0].control_type;

    
    //unit_cell_essentials(BC_type, max_temp, path_data, uc_essentials); //Unused
    
    std::vector<Node> nodes;
    read_nodes_file(nodes, path_data, pointsfile);
    nb_nodes = nodes.size();
    int nb_nodes_init = nb_nodes;
    
    std::vector<int> NodeCD;
    if((loading_type == 1) || (loading_type == 2)){
        if(control_type == 1){
            NodeCD = {nb_nodes+1, nb_nodes+1, nb_nodes+1, nb_nodes+2, nb_nodes+2, nb_nodes+2};
        }
        else if(control_type > 1) {
            NodeCD = {nb_nodes+1, nb_nodes+1, nb_nodes+1, nb_nodes+2, nb_nodes+2, nb_nodes+2, nb_nodes+3, nb_nodes+3, nb_nodes+3};
        }
    }
    //Thermal
    else if(loading_type == 3) {
        NodeCD = {nb_nodes+1, nb_nodes+2, nb_nodes+3};
    }
    else{
        cout << "Error in software/Salome_apply_inter.cpp : loading_type should take the following values : 1 for mechanical loading, 2 for thermomechanical loading and 3 for pure thermal (heat transfer) analysis" << endl;
    }
    
    cubic_mesh cm(nodes, pointsfile);
    cm.get_domain();
    //cm.find_pairs();
    cm.construct_lists();
    
    cubic_mesh cm_perio = perio_RVE(cm, nb_nodes);
    cm_perio.construct_lists();
        
    //Phases
    string steps_file_name = "steps_file_name.inp";
    std::vector<std::shared_ptr<step>> aba_steps;
    
    update_steps(aba_steps, blocks, nlgeom, loading_type, max_temp);
    write_steps(aba_steps, loading_type, T_init, path_run, steps_file_name);
    
    
    //Action 3 : Sections and materials
    std::vector<section_characteristics> sections;
    read_sections(sections, loading_type, path_data, sections_file);
    write_sections(sections, loading_type, path_run, sections_out);
    
    //Action 4 : the Periocidc Boundary Conditions : PBC & Constraint drivers : CDN
    
    string PBC_file_name = "PBC_file_name.inp";
    string TIE_file_name = "TIE_file_name.inp";
    string CDN_file_name = "CDN_file_name.inp";
//    write_PBC(cm, nb_nodes, path_run, PBC_file_name);
    write_PBC(cm, path_run, PBC_file_name);
//    write_TIE(cm, cm_perio, path_run, TIE_file_name);
    int n_neigh = 4;
    double pow_int = 1.0;
    write_NonPerio_CDN(cm, cm_perio, NodeCD, loading_type, control_type, n_neigh, pow_int, path_run, CDN_file_name);
//    write_CDN(cm, path_run, CDN_file_name);
    
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
    aba_head.close();

    aba_head.open(postproc_file, ios::out);
    aba_head << "volume\t" << cm.volume << "\n";
    aba_head.close();
    
    return 0;
}
