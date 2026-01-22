
1. Compile the plugin in the external directory
clang++ -w -c -std=c++14 -fPIC my_plugin_sum.cpp
clang++ -w -std=c++14 -shared my_plugin_sum.o -o my_plugin_sum -lsimcoon -larmadillo

2. Compile the software test
clang++ -w -std=c++14 -o test test.cpp -lsimcoon -larmadillo


For the simcoon plugin

clang++ -w -c -std=c++14 -fPIC umat_plugin_ext.cpp
clang++ -w -std=c++14 -shared umat_plugin_ext.o -o umat_plugin_ext -lsimcoon -larmadillo


to compile Abaqus subroutines:
export LIBRARY_PATH=/usr/local/lib/gcc/9
gfortran -w -c UMAT_ABAQUS_ELASTIC.for

clang++ -w -c -std=c++14 -fPIC umat_plugin_aba.cpp
clang++ -w -std=c++14 -shared umat_plugin_aba.o UMAT_ABAQUS_ELASTIC.o -o umat_plugin_aba -lsimcoon -larmadillo -lgfortran


