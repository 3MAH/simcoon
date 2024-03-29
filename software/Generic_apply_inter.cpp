///@file Salome_apply_inter.cpp
///@brief Salome_apply_inter : Apply Materials, Steps and periodic Boundary Conditions to a periodic box from a Salome model
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
    
    string path_data = "data";
    string path_run = "run_aba";
    string nodes_file = "nodes.inp";
    string elements_file = "elements.inp";
    string sets_file = "sets.inp";
    string ori_file = "orientations.inp";
    string Mesh_ori_file = "Mesh.ori";
    string uc_essentials = "unit_cell_essentials.inp";
    string steps_file_name = "run_strain_3D.inp";
    
    string buffer;
    
    //Read the loading path - to get the loading_path informations
    std::vector<block> blocks;  //loading blocks
    double T_init = 0.;
    string pathfile = "path.txt";
    bool nlgeom = false;
    read_path(blocks, T_init, path_data, pathfile);
    unsigned int loading_type = blocks[0].type;
    unsigned int control_type = blocks[0].control_type;

    int max_temp = 0;
    
    //std::vector<section_characteristics> sections;

    std::vector<Node> nodes_full;
    read_nodes_file(nodes_full, path_data, nodes_file);
    
    //copy the original node file
    string src_file = path_data + "/" + nodes_file;
    string dst_file = path_run + "/" + nodes_file;
    boost::filesystem::copy_file(src_file,dst_file,boost::filesystem::copy_option::overwrite_if_exists);

    //copy the original element file
    src_file = path_data + "/" + elements_file;
    dst_file = path_run + "/" + elements_file;
    boost::filesystem::copy_file(src_file,dst_file,boost::filesystem::copy_option::overwrite_if_exists);

    //copy the original sections file
    src_file = path_data + "/" + elements_file;
    dst_file = path_run + "/" + elements_file;
    boost::filesystem::copy_file(src_file,dst_file,boost::filesystem::copy_option::overwrite_if_exists);

    //copy the original orientations file
    src_file = path_data + "/" + ori_file;
    dst_file = path_run + "/" + ori_file;
    boost::filesystem::copy_file(src_file,dst_file,boost::filesystem::copy_option::overwrite_if_exists);

    //copy the original Mesh orientation file
    src_file = path_data + "/" + Mesh_ori_file;
    dst_file = path_run + "/" + Mesh_ori_file;
    boost::filesystem::copy_file(src_file,dst_file,boost::filesystem::copy_option::overwrite_if_exists);

    //copy the original Mesh orientation file
    src_file = path_data + "/" + steps_file_name;
    dst_file = path_run + "/" + steps_file_name;
    boost::filesystem::copy_file(src_file,dst_file,boost::filesystem::copy_option::overwrite_if_exists);

    
/*    std::vector<Element> elements_full;
    for(auto sc : sections) {
        nodes_full.insert(nodes_full.end(), sc.nodes.begin(), sc.nodes.end());
        elements_full.insert(elements_full.end(), sc.elements.begin(), sc.elements.end());
    }*/
    
    int nb_nodes_full = nodes_full.size();
    int nb_nodes_init = nb_nodes_full;
    
    std::vector<int> NodeCD;
    if((loading_type == 1) || (loading_type == 2)){
        if(control_type == 1){
            NodeCD = {nb_nodes_full+1, nb_nodes_full+1, nb_nodes_full+1, nb_nodes_full+2, nb_nodes_full+2, nb_nodes_full+2};
        }
        else if(control_type > 1) {
            NodeCD = {nb_nodes_full+1, nb_nodes_full+1, nb_nodes_full+1, nb_nodes_full+2, nb_nodes_full+2, nb_nodes_full+2, nb_nodes_full+3, nb_nodes_full+3, nb_nodes_full+3};
        }
    }
    //Thermal
    else if(loading_type == 3) {
        NodeCD = {nb_nodes_full+1, nb_nodes_full+2, nb_nodes_full+3};
    }
    else{
        cout << "Error in software/Salome_apply_inter.cpp : loading_type should take the following values : 1 for mechanical loading, 2 for thermomechanical loading and 3 for pure thermal (heat transfer) analysis" << endl;
    }
    
    cubic_mesh cm(nodes_full, "nodes_full");
    cm.get_domain();
    //cm.find_pairs();
    cm.construct_lists();
    
    cubic_mesh cm_perio = perio_RVE(cm, nb_nodes_full);
    cm_perio.construct_lists();
    
    //Action_0 : The node and element file
//    write_nodes_file(nodes_full, path_run, nodes_file);
//    write_elements_file(elements_full, path_run, elements_file);
//    write_sets_file(sections, nodes_full, path_run, sets_file);
    
    //Action 1 : The materials and sections
    string umat_name;
    string sections_out = "mat_sec.inp";
    
    //Phases
//    string steps_file_name = "steps_file_name.inp";
//    std::vector<std::shared_ptr<step>> aba_steps;
    
//    update_steps(aba_steps, blocks, nlgeom, loading_type, max_temp);
//    write_steps(aba_steps, loading_type, T_init, path_run, steps_file_name);
    
    //Action 3 : Sections and materials

//    read_sections(sections, loading_type, path_data, sections_file);
//    write_sections(sections, loading_type, path_run, sections_out);
    
    //Action 4 : the Periocidc Boundary Conditions : PBC & Constraint drivers : CDN
    
    string PBC_file_name = "PBC_file_name.inp";
    string CDN_file_name = "CDN_file_name.inp";
    write_PBC(cm, path_run, PBC_file_name);
    int n_neigh = 4;
    double pow_int = 1.0;
    write_NonPerio_CDN(cm, cm_perio, NodeCD, loading_type, control_type, n_neigh, pow_int, path_run, CDN_file_name);
    
    //Finally
    string run_file = "run_aba.inp";
    string run_out = path_run + "/" + run_file;
    std::ofstream aba_head;
    
    aba_head.open(run_out, ios::out);

    aba_head << "*Heading" << "\n";
    aba_head << "**" << "\n";
    
    aba_head << "*INCLUDE, INPUT= " << nodes_file << "\n";
    aba_head << "*INCLUDE, INPUT= " << elements_file << "\n";
    aba_head << "*INCLUDE, INPUT= " << sets_file << "\n";
    aba_head << "*INCLUDE, INPUT= " << sections_out << "\n";
    aba_head << "*INCLUDE, INPUT= " << PBC_file_name << "\n";
    aba_head << "*INCLUDE, INPUT= " << CDN_file_name << "\n";
    aba_head << "*INCLUDE, INPUT= " << steps_file_name << "\n";
    aba_head.close();

//    aba_head.open(postproc_file, ios::out);
//    aba_head << "volume\t" << cm.volume << "\n";
//    aba_head.close();
    
    return 0;
}
