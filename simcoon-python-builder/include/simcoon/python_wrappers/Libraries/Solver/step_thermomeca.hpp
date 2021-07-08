#pragma once
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <simcoon/Simulation/Solver/step.hpp>
#include <simcoon/Simulation/Solver/step_thermomeca.hpp>

namespace simpy{

//======================================
class step_thermomeca_py : public simcoon::step_thermomeca
//======================================
{
    private:

    protected:

    public :

    step_thermomeca_py();
    step_thermomeca_py(const int &, const double &, const double &, const double &, const int &, const unsigned int &, const boost::python::numpy::ndarray&, const boost::python::numpy::ndarray&, const boost::python::numpy::ndarray&, const double&, const int&, const boost::python::numpy::ndarray&, const boost::python::numpy::ndarray&, const boost::python::numpy::ndarray&); //Constructor with parameters
        
    step_thermomeca_py(const simcoon::step_thermomeca &);
    
    //using simcoon::step_meca::generate;
    virtual void generate(const double&, const boost::python::numpy::ndarray&, const boost::python::numpy::ndarray&, const double&);
    //using simcoon::step_meca::generate_kin;
    virtual void generate_kin(const double&, const boost::python::numpy::ndarray&, const double &);
//    using step_meca::assess_inc;
//    virtual void assess_inc(const double &, double &, const double &, simcoon::phase_characteristics &, double &, const double &, const arma::mat &, const int &);

    boost::python::numpy::ndarray Get_times();

    boost::python::numpy::ndarray Get_cBC_meca();
    void Set_cBC_meca(const boost::python::numpy::ndarray &cBC_meca);
    
    boost::python::numpy::ndarray Get_BC_meca();
    void Set_BC_meca(const boost::python::numpy::ndarray &BC_meca);
    
    boost::python::numpy::ndarray Get_mecas();
    boost::python::numpy::ndarray Get_BC_w();
    void Set_BC_w(const boost::python::numpy::ndarray &BC_w);
    
    boost::python::numpy::ndarray Get_BC_R();
    void Set_BC_R(const boost::python::numpy::ndarray &BC_R);
    
    boost::python::numpy::ndarray Get_Ts();
};

} //namespace simpy

