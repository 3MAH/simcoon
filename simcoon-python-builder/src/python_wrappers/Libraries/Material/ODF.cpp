#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <carma>
#include <armadillo>

#include <simcoon/Simulation/Phase/phase_characteristics.hpp>
#include <simcoon/Simulation/Phase/read.hpp>
#include <simcoon/Simulation/Phase/write.hpp>
#include <simcoon/Continuum_mechanics/Material/ODF.hpp>
#include <simcoon/Continuum_mechanics/Material/ODF2Nphases.hpp>
#include <simcoon/Continuum_mechanics/Material/read.hpp>

#include <simcoon/python_wrappers/Libraries/Material/ODF.hpp>

using namespace std;
using namespace arma;
namespace py=pybind11;

namespace simpy{
    
py::array_t<double> get_densities_ODF(const py::array_t<double> &x_py, const std::string &path_data_py, const std::string &peak_file_py, const bool &radian) {
    
//    string path_data = bp::extract<std::string>(path_data_py);
//    string peak_file = bp::extract<std::string>(peak_file_py);
    
    //transform x in a vec
    vec x = carma::arr_to_col(x_py);
    //Get the densities
    vec y = simcoon::get_densities_ODF(x, path_data_py, peak_file_py, radian);
    //Get the densities
    return carma::col_to_arr(y);
}
    
void ODF_discretization(const int &nphases_rve, const int &num_phase_disc, const double &angle_min, const double &angle_max, const std::string &umat_name_py, const py::array_t<double> &props_py, const std::string &path_data_py, const std::string &peak_file_py, const std::string &rve_init_file_py, const std::string &rve_disc_file_py, const int &angle_mat) {

//    string umat_name = bp::extract<std::string>(umat_name_py);
//    string path_data = bp::extract<std::string>(path_data_py);
//    string peak_file = bp::extract<std::string>(peak_file_py);
//    string rve_init_file = bp::extract<std::string>(rve_init_file_py);
//    string rve_disc_file = bp::extract<std::string>(rve_disc_file_py);
    
    vec props = carma::arr_to_col(props_py);
        
    double psi_rve = 0.;
    double theta_rve = 0.;
    double phi_rve = 0.;
    
    simcoon::phase_characteristics rve_init;
    rve_init.sptr_matprops->update(0, umat_name_py, 1, psi_rve, theta_rve, phi_rve, props.n_elem, props);
    
    string inputfile;
    string outputfile;
    
    std::map<string, int> list_umat;
    list_umat = {{"MIHEN",100},{"MIMTN",101},{"MISCN",102},{"MIPCW",103},{"MIPLN",104}};
    
    //Here we read everything about the initial rve
    switch (list_umat[rve_init.sptr_matprops->umat_name]) {
            
        case 100: case 101: case 102: case 103: {
            rve_init.construct(2,1); //The rve is supposed to be mechanical only here
            simcoon::read_ellipsoid(rve_init, path_data_py, rve_init_file_py);
            break;
        }
        case 104: {
            rve_init.construct(1,1); //The rve is supposed to be mechanical only here
            simcoon::read_layer(rve_init, path_data_py, rve_init_file_py);
            break;
        }
    }
    
    simcoon::ODF odf_rve(0, false, angle_min, angle_max);
    simcoon::read_peak(odf_rve, path_data_py, peak_file_py);
    
    simcoon::phase_characteristics rve = discretize_ODF(rve_init, odf_rve, num_phase_disc, nphases_rve, angle_mat);
    
    if(rve.shape_type == 0) {
        simcoon::write_phase(rve, path_data_py, rve_disc_file_py);
    }
    if(rve.shape_type == 1) {
        simcoon::write_layer(rve, path_data_py, rve_disc_file_py);
    }
    else if(rve.shape_type == 2) {
        simcoon::write_ellipsoid(rve, path_data_py, rve_disc_file_py);
    }
    else if(rve.shape_type == 3) {
        simcoon::write_cylinder(rve, path_data_py, rve_disc_file_py);
    }
    
    return;
}
    
} //namespace simcoon_py
