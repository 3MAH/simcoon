External
========

Compile an external plugin
--------------------------

Simcoon offers the possibility to develop your own user material plugin. You can develop such an external plugin using the Simcoon API, but this isn't mandatory. You can also rely on other libraries or APIs.
Two examples are provided in the folder 'External'.

There is two formats to define your User Material plugin. The first one is based on Armadillo vectors and matrices, and correspond to the classical parameters


.. code-block:: none

    void umat_(double *stress, double *statev, double *ddsdde, double &sse, double &spd, double &scd, double &rpl, double *ddsddt, double *drplde, double &drpldt, const double *stran, const double *dstran, const double *time, const double &dtime, const double &temperature, const double &Dtemperature, const double &predef, const double &dpred, char *cmname, const int &ndi, const int &nshr, const int &ntens, const int &nstatev, const double *props, const int &nprops, const double &coords, const double *drot, double &pnewdt, const double &celent, const double *dfgrd0, const double *dfgrd1, const int &noel, const int &npt, const double &layer, const int &kspt, const int &kstep, const int &kinc)



The format of the external User Material plugin is similar to that of the Finite Element Analysis software Abaqus. It accepts the following arguments:

.. code-block:: none

    void umat_(double *stress, double *statev, double *ddsdde, double &sse, double &spd, double &scd, double &rpl, double *ddsddt, double *drplde, double &drpldt, const double *stran, const double *dstran, const double *time, const double &dtime, const double &temperature, const double &Dtemperature, const double &predef, const double &dpred, char *cmname, const int &ndi, const int &nshr, const int &ntens, const int &nstatev, const double *props, const int &nprops, const double &coords, const double *drot, double &pnewdt, const double &celent, const double *dfgrd0, const double *dfgrd1, const int &noel, const int &npt, const double &layer, const int &kspt, const int &kstep, const int &kinc)

==============  ==========
parameter       definition
==============  ==========
stress          array containing the components of the stress tensor (dimension ntens)
statev          array containing the evolution variables (dimension nstatev)
ddsdde          array containing the mechanical tangent operator (dimension ntens*ntens)
sse             unused
spd             unused
scd             unused
rpl             unused
ddsddt          array containing the thermal tangent operator
drple           unused
drpldt          unused
stran           array containing total strain component (dimension ntens) at the beginning of increment
dstran          array containing the component of total strain increment (dimension ntens)
time            two compoenent array : first component is the value of step time at the beginning of the current increment and second component is the value of total time at the beginning of the current increment
dtime           time increment
temperature     temperature avlue at the beginning of increment
Dtemperature    temperature increment
predef          unused
dpred           unused
cmname          user-defined material name
ndi             number of direct stress components
nshr            number of shear stress components
ntens           number stress and strain components
nstatev         number of evolution variables
props           array containing material properties
nprops          number of material properties
coords          coordinates of the considered point
drot            rotation increment matrix (dimension 3*3)
pnewdt          ratio of suggested new time increment
celent          characteristic element length
dfgrd0          array containing the deformation gradient at the beginning of increment (dimension 3*3)
dfgrd1          array containing the deformation gradient at the end of increment (dimension 3*3)
noel            element number
npt             integration point number
layer           layer number - not used
kspt            section point number within the current layer - not used
kstep           step number
kinc            increment number
==============  ==========

To compile them you can utilize the following commands (I have used clang as a compiler)

.. code-block:: none

    clang++ -c -fPIC -std=c++14 umat_plugin_aba.cpp
    clang++ -std=c++14 -shared -lsimcoon -larmadillo -o libumat_plugin_aba.dylib umat_plugin_aba.o

The first example is a replica of the ``UMAT_ABAQUS_ELASTIC.for`` interface.
