The Function_N Library
========================

.. default-domain:: cpp

.. function:: func_N(const vec &params, const vec &variables, const string& N_file, const string& outputfile, const string& path_data, const string& path_results)

    This function computes the result of a single function math:`y=f(x)`, providing the vector of input values :math:`x`.
    The list of values shall be stored in a file named *N_file*, in the data folder *path_data*

    .. warning:: This is a temproary function and necessitate to modify the code with your own function. This will be deprecated in a future release
    You shall define your own function along with the definition of the vector, for example by adding

    .. code-block:: cpp

        vec y = p_cumulative(N, variables(0), variables(1), params); Insert here the fonction you want


    in the file func_N.cpp. You should then reinstall the library

    The x and y values are written in a file named *outputfile*, in the data folder *path_results*

    Example:

    .. code-block:: cpp

    string outputfile = "results.txt";
    string N_file = "list_inputs.txt";
    string path_data = "data";
    string path_data = "results";

    vec props = {1., 2.} //A vector utilized to define the parameters
    vec sigma = randu(6);
    double sigma_eq = Mises_stress(sigma);
    vec variables = {sigma_eq} //A vector utilized to define the variables

    func_N(props, variables, N_file, outputfile, path_data, path_results);



