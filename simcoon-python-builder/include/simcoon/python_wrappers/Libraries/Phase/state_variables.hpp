#pragma once
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <simcoon/Simulation/Phase/state_variables.hpp>
#include <simcoon/Simulation/Phase/state_variables_M.hpp>
#include <simcoon/Simulation/Phase/state_variables_T.hpp>

namespace simpy{

//======================================
class state_variables_py : public simcoon::state_variables
//======================================
{
    private:

    protected:

    public :

        state_variables_py();
        state_variables_py(const boost::python::numpy::ndarray&, const boost::python::numpy::ndarray&, const boost::python::numpy::ndarray&, const boost::python::numpy::ndarray&, const boost::python::numpy::ndarray&, const double&, const double&);

        boost::python::numpy::ndarray Get_F0();

        void Set_F0(const boost::python::numpy::ndarray &);

        boost::python::numpy::ndarray Get_F1();

        void Set_F1(const boost::python::numpy::ndarray &);

        boost::python::numpy::ndarray Get_etot();

        boost::python::numpy::ndarray Get_Detot();

        boost::python::numpy::ndarray Get_Etot();

        boost::python::numpy::ndarray Get_DEtot();

        boost::python::numpy::ndarray Get_statev();

        boost::python::numpy::ndarray Get_R();

        boost::python::numpy::ndarray Get_DR();

        state_variables_py rotate_l2g_py(const state_variables_py&, const double&, const double&, const double&);
    
        state_variables_py rotate_g2l_py(const state_variables_py&, const double&, const double&, const double&);
};

} //namespace simpy
