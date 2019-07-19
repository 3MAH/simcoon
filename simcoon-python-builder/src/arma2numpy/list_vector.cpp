/* This file is part of arma2numpy.

 arma2numpy is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 arma2numpy is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with arma2numpy.  If not, see <http://www.gnu.org/licenses/>.

 */

#include <iostream>
#include <vector>
#include <memory>
#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/node.hpp>
#include <simcoon/Simulation/Solver/block.hpp>
#include <simcoon/Simulation/Solver/step.hpp>
#include <simcoon/Simulation/Identification/constants.hpp>
#include <simcoon/Simulation/Identification/parameters.hpp>
#include <simcoon/arma2numpy/list_vector.hpp>

namespace bp = boost::python;

namespace arma2numpy {

template<typename T> std::vector<T> py_list_to_std_vector(const bp::object& iterable)
{
    return std::vector<T>( bp::stl_input_iterator<T>( iterable ),
                            bp::stl_input_iterator<T>() );
}

std::vector<double> py_list_to_std_vector_double(const bp::object& iterable)
{
    return py_list_to_std_vector<double>(iterable);
}

std::vector<int> py_list_to_std_vector_int(const bp::object& iterable)
{
    return py_list_to_std_vector<int>(iterable);
}

std::vector<std::string> py_list_to_std_vector_string(const bp::object& iterable)
{
    return py_list_to_std_vector<std::string>(iterable);
}

std::vector<simcoon::constants> py_list_to_std_vector_constants(const bp::object& iterable)
{
    return py_list_to_std_vector<simcoon::constants>(iterable);
}

std::vector<simcoon::parameters> py_list_to_std_vector_parameters(const bp::object& iterable)
{
    return py_list_to_std_vector<simcoon::parameters>(iterable);
}

std::vector<simcoon::Node> py_list_to_std_vector_Node(const bp::object& iterable)
{
    return py_list_to_std_vector<simcoon::Node>(iterable);
}

std::vector<simcoon::block> py_list_to_std_vector_block(const bp::object& iterable)
{
    return py_list_to_std_vector<simcoon::block>(iterable);
}

std::vector<std::shared_ptr<simcoon::step>> py_list_to_std_vector_shptr_step(const bp::object& iterable)
{
    return py_list_to_std_vector<std::shared_ptr<simcoon::step>>(iterable);
}

/*std::vector<simcoon::section_characteristics> py_list_to_std_vector_section_characteristics(const bp::object& iterable)
{
    return py_list_to_std_vector<simcoon::section_characteristics>(iterable);
}*/

/// @brief Transfer ownership to a Python object.  If the transfer fails,
///        then object will be destroyed and an exception is thrown.
template <typename T> boost::python::object transfer_to_python(T* t)
{
    // Transfer ownership to a smart pointer, allowing for proper cleanup
    // incase Boost.Python throws.
    std::unique_ptr<T> ptr(t);

    // Create a functor with a call policy that will have Boost.Python
    // manage the new object, then invoke it.
    bp::object object = bp::make_function([t]() { return t; },
                                          bp::return_value_policy<bp::manage_new_object>(),
                                          boost::mpl::vector<T*>())();

    // As the Python object now has ownership, release ownership from
    // the smart pointer.
    ptr.release();
    return object;
}

template <class T> bp::list std_vector_to_py_list(const std::vector<T> &vector) {
    bp::list list;
    for (auto const & x : vector) {
        list.append(x);
    }
    return list;
}

template <class T> bp::list std_vector_to_py_list_class(const std::vector<T> &vector) {
    bp::list list;
    for (auto r : vector) {
        T *x = new T(r);
        list.append(transfer_to_python<T>(x));
    }
    return list;
}

boost::python::list std_vector_to_py_list_double(const std::vector<double> &vector)
{
    return std_vector_to_py_list<double>(vector);
}

boost::python::list std_vector_to_py_list_int(const std::vector<int> &vector)
{
    return std_vector_to_py_list<int>(vector);
}

boost::python::list std_vector_to_py_list_string(const std::vector<std::string> &vector)
{
    return std_vector_to_py_list<std::string>(vector);
}

boost::python::list std_vector_to_py_list_constants(const std::vector<simcoon::constants> &vector)
{
    return std_vector_to_py_list_class<simcoon::constants>(vector);
}

boost::python::list std_vector_to_py_list_parameters(const std::vector<simcoon::parameters> &vector)
{
    return std_vector_to_py_list_class<simcoon::parameters>(vector);
}

boost::python::list std_vector_to_py_list_Node(const std::vector<simcoon::Node> &vector)
{
    return std_vector_to_py_list_class<simcoon::Node>(vector);
}

/*boost::python::list std_vector_to_py_list_shptr_step(const std::vector<std::shared_ptr<simcoon::step> > &vector)
{
    return std_vector_to_py_list_class<std::shared_ptr<simcoon::step>>(vector);
}*/

/*boost::python::list std_vector_to_py_list_section_characteristics(const std::vector<simcoon::section_characteristics> &vector)
{
    return std_vector_to_py_list_class<simcoon::section_characteristics>(vector);
}*/

} //end of namespace arma2numpy
