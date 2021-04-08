#pragma once
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <simcoon/Simulation/Phase/state_variables.hpp>
#include <simcoon/Simulation/Phase/state_variables_M.hpp>
#include <simcoon/Simulation/Phase/state_variables_T.hpp>

namespace simpy{

//======================================
class state_variables_T_py : public simcoon::state_variables_T
//======================================
{
    private:

    protected:

    public :

        state_variables_T_py();
        state_variables_T_py(const boost::python::numpy::ndarray&, const boost::python::numpy::ndarray&, const boost::python::numpy::ndarray&, const boost::python::numpy::ndarray&, const boost::python::numpy::ndarray&, const double&, const double&);

        boost::python::numpy::ndarray Get_F0();

        void Set_F0(const boost::python::numpy::ndarray &F0_py);

        boost::python::numpy::ndarray Get_F1();

        void Set_F1(const boost::python::numpy::ndarray &F1_py);

        boost::python::numpy::ndarray Get_etot();

        boost::python::numpy::ndarray Get_Detot();

        boost::python::numpy::ndarray Get_Etot();

        boost::python::numpy::ndarray Get_DEtot();

        boost::python::numpy::ndarray Get_statev();

        boost::python::numpy::ndarray Get_R();

        boost::python::numpy::ndarray Get_DR();
    
        boost::python::numpy::ndarray Get_Wm();
    
        boost::python::numpy::ndarray Get_Wt();
    
        boost::python::numpy::ndarray Get_dSdE();
    
        boost::python::numpy::ndarray Get_dSdEt();
    
        boost::python::numpy::ndarray Get_dSdT();
    
        boost::python::numpy::ndarray Get_drdE();
    
        boost::python::numpy::ndarray Get_drdT();
    
        using simcoon::state_variables_T::rotate_l2g;
        virtual state_variables_T_py& rotate_l2g(const state_variables_T_py&, const double&, const double&, const double&);
    
        using simcoon::state_variables_T::rotate_g2l;
        virtual state_variables_T_py& rotate_g2l(const state_variables_T_py&, const double&, const double&, const double&);

};

} //namespace simpy
