External
========

Simcoon offers the possibility to develop your own user material plugin that can be loaded at runtime by the solver.
You can develop such an external plugin using the Simcoon API. Two plugin formats are available:

- **UMEXT** (``umat_plugin_ext_api``): Native simcoon format using Armadillo vectors and matrices
- **UMABA** (``umat_plugin_aba_api``): Abaqus UMAT-compatible format for existing Fortran/C++ UMATs

Both examples are provided in the ``external/`` folder.

UMEXT Plugin Format
-------------------

The UMEXT format is the native simcoon plugin format that uses Armadillo vectors and matrices directly.
Your plugin must inherit from ``umat_plugin_ext_api`` and implement the following interface:

.. code-block:: cpp

    #include <simcoon/Continuum_mechanics/Umat/umat_plugin_api.hpp>
    
    class umat_plugin_ext : public umat_plugin_ext_api {
    public:
        std::string name() const override {
            return "umext";  // Material name used in material.dat
        }
    
        void umat_external_M(
            const arma::vec &Etot,      // Total strain at end of increment
            const arma::vec &DEtot,     // Strain increment
            arma::vec &sigma,           // Stress tensor (in/out)
            arma::mat &Lt,              // Tangent modulus (output)
            arma::mat &L,               // Elastic stiffness (output)
            arma::vec &sigma_in,        // Internal stress (output)
            const arma::mat &DR,        // Rotation increment matrix
            const int &nprops,          // Number of material properties
            const arma::vec &props,     // Material properties
            const int &nstatev,         // Number of state variables
            arma::vec &statev,          // State variables (in/out)
            const double &T,            // Temperature at start
            const double &DT,           // Temperature increment
            const double &Time,         // Current time
            const double &DTime,        // Time increment
            double &Wm,                 // Mechanical work (in/out)
            double &Wm_r,               // Reversible work (in/out)
            double &Wm_ir,              // Irreversible work (in/out)
            double &Wm_d,               // Dissipated work (in/out)
            const int &ndi,             // Number of direct stress components
            const int &nshr,            // Number of shear stress components
            const bool &start,          // True if first increment
            const int &solver_type,     // Solver type (0=Newton, 1=RNL)
            double &tnew_dt             // Suggested time step ratio
        ) override;
    };
    
    // Required factory functions
    extern "C" umat_plugin_ext_api* create_api() {
        return new umat_plugin_ext();
    }
    
    extern "C" void destroy_api(umat_plugin_ext_api* p) {
        delete p;
    }

The file ``external/umat_plugin_ext.cpp`` provides a complete example implementing a thermoelastic material.

UMABA Plugin Format (Abaqus Compatibility)
------------------------------------------

The UMABA format provides a compatibility layer for existing Abaqus UMAT subroutines.
Your plugin wraps an external ``umat_`` function with the standard Abaqus signature.

.. code-block:: cpp

    #include <simcoon/Continuum_mechanics/Umat/umat_plugin_api.hpp>
    
    // Declare the external UMAT subroutine (Fortran or C++)
    extern "C" {
        void umat_(double *stress, double *statev, double *ddsdde, ...);
    }
    
    class umat_plugin_aba : public umat_plugin_aba_api {
    public:
        std::string name() const override {
            return "umaba";  // Material name used in material.dat
        }
    
        void umat_abaqus(
            simcoon::phase_characteristics &rve,  // RVE state (contains material props, state vars)
            const arma::mat &DR,                   // Rotation increment
            const double &Time,                    // Current time
            const double &DTime,                   // Time increment
            const int &ndi,                        // Number of direct stress components
            const int &nshr,                       // Number of shear stress components
            bool &start,                           // True if first increment
            const int &solver_type,                // Solver type
            double &tnew_dt                        // Suggested time step ratio
        ) override;
    };

The Abaqus UMAT function signature:

.. code-block:: none

    void umat_(double *stress, double *statev, double *ddsdde, double &sse, double &spd, double &scd, double &rpl, double *ddsddt, double *drplde, double &drpldt, const double *stran, const double *dstran, const double *time, const double &dtime, const double &temperature, const double &Dtemperature, const double &predef, const double &dpred, char *cmname, const int &ndi, const int &nshr, const int &ntens, const int &nstatev, const double *props, const int &nprops, const double &coords, const double *drot, double &pnewdt, const double &celent, const double *dfgrd0, const double *dfgrd1, const int &noel, const int &npt, const double &layer, const int &kspt, const int &kstep, const int &kinc)

==============  ==========
parameter       definition
==============  ==========
stress          array containing the components of the stress tensor (dimension ntens)
statev          array containing the evolution variables (dimension nstatev)
ddsdde          array containing the mechanical tangent operator (dimension ntens*ntens)
sse             specific elastic strain energy
spd             plastic dissipation
scd             creep dissipation
rpl             volumetric heat generation (coupled analysis only)
ddsddt          array containing the thermal tangent operator
drplde          variation of RPL with strain increments
drpldt          variation of RPL with temperature
stran           array containing total strain component (dimension ntens) at the beginning of increment
dstran          array containing the component of total strain increment (dimension ntens)
time            two component array: step time and total time at beginning of increment
dtime           time increment
temperature     temperature value at the beginning of increment
Dtemperature    temperature increment
predef          predefined field variables
dpred           predefined field increments
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

The file ``external/umat_plugin_aba.cpp`` provides a complete example, and ``external/UMAT_ABAQUS_ELASTIC.for`` is a sample Fortran UMAT.

UMANS Plugin Format (Ansys Compatibility)
-----------------------------------------

The UMANS format provides a compatibility layer for existing Ansys USERMAT subroutines.
Your plugin wraps an external ``usermat_`` function with the standard Ansys signature.

.. code-block:: cpp

    #include <simcoon/Continuum_mechanics/Umat/umat_plugin_api.hpp>
    
    // Declare the external USERMAT subroutine (Fortran)
    extern "C" {
        void usermat_(int *matId, int *elemId, int *kDomIntPt, ...);
    }
    
    class umat_plugin_ans : public umat_plugin_ans_api {
    public:
        std::string name() const override {
            return "umans";  // Material name used in material.dat
        }
    
        void umat_ansys(
            simcoon::phase_characteristics &rve,
            const arma::mat &DR,
            const double &Time, const double &DTime,
            const int &ndi, const int &nshr,
            bool &start, const int &solver_type,
            double &tnew_dt
        ) override;
    };

The Ansys USERMAT function signature:

.. code-block:: fortran

    subroutine usermat(
         matId, elemId, kDomIntPt, kLayer, kSectPt,
         ldstep, isubst, keycut,
         nDirect, nShear, ncomp, nStatev, nProp,
         Time, dTime, Temp, dTemp,
         stress, ustatev, dsdde,
         sedEl, sedPl, epseq,
         Strain, dStrain, epsPl,
         prop, coords, rotateM,
         defGrad_t, defGrad,
         tsstif, epsZZ,
         var1, var2, var3, var4, var5)

==============  ==========
parameter       definition
==============  ==========
matId           Material ID number
elemId          Element number
kDomIntPt       Current integration point in element
kLayer          Current layer number
kSectPt         Current section point within layer
ldstep          Current load step
isubst          Current substep
keycut          Cutback flag (output: set to 1 to request smaller time step)
nDirect         Number of direct stress components (typically 3)
nShear          Number of shear stress components (typically 3)
ncomp           Number of stress/strain components (nDirect + nShear)
nStatev         Number of state variables
nProp           Number of material properties
Time            Current time
dTime           Time increment
Temp            Current temperature
dTemp           Temperature increment
stress          Stress tensor (in: at t, out: at t+dt), dimension ncomp
ustatev         State variables (in/out), dimension nStatev
dsdde           Material tangent modulus (output), dimension ncomp×ncomp
sedEl           Elastic strain energy density (in/out)
sedPl           Plastic strain energy density (in/out)
epseq           Equivalent plastic strain (in/out)
Strain          Total strain at t, dimension ncomp
dStrain         Strain increment, dimension ncomp
epsPl           Plastic strain components, dimension ncomp
prop            Material property array, dimension nProp
coords          Integration point coordinates (x, y, z)
rotateM         Rotation matrix (3×3)
defGrad_t       Deformation gradient at t (3×3)
defGrad         Deformation gradient at t+dt (3×3)
tsstif          Transverse shear stiffness (for shells)
epsZZ           Out-of-plane strain (for plane stress)
var1-var5       Reserved for future use
==============  ==========

.. note::

    **Voigt notation difference**: Ansys uses the ordering (11, 22, 33, 12, 23, 13) for stress/strain 
    components, while simcoon uses (11, 22, 33, 12, 13, 23). The plugin handles this conversion 
    automatically.

The file ``external/umat_plugin_ans.cpp`` provides a complete example, and ``external/USERMAT_ANSYS_ELASTIC.for`` is a sample Fortran USERMAT.

Compiling External Plugins
--------------------------

To compile plugins you can utilize the following commands (using clang as compiler):

For UMEXT plugins:

.. code-block:: bash

    clang++ -c -fPIC -std=c++14 umat_plugin_ext.cpp -I/path/to/simcoon/include
    clang++ -std=c++14 -shared -lsimcoon -larmadillo -o libumat_plugin_ext.dylib umat_plugin_ext.o

For UMABA plugins (with Fortran UMAT):

.. code-block:: bash

    gfortran -c -fPIC UMAT_ABAQUS_ELASTIC.for -o umat_fortran.o
    clang++ -c -fPIC -std=c++14 umat_plugin_aba.cpp -I/path/to/simcoon/include
    clang++ -std=c++14 -shared -lsimcoon -larmadillo -lgfortran -o libumat_plugin_aba.dylib umat_plugin_aba.o umat_fortran.o

For UMANS plugins (with Ansys Fortran USERMAT):

.. code-block:: bash

    gfortran -c -fPIC USERMAT_ANSYS_ELASTIC.for -o usermat_fortran.o
    clang++ -c -fPIC -std=c++14 umat_plugin_ans.cpp -I/path/to/simcoon/include
    clang++ -std=c++14 -shared -lsimcoon -larmadillo -lgfortran -o libumat_plugin_ans.dylib umat_plugin_ans.o usermat_fortran.o

On Linux, replace ``.dylib`` with ``.so``. On Windows, use ``.dll``.

The compiled plugin library must be placed in a location where it can be found at runtime (e.g., the working directory or a path in ``LD_LIBRARY_PATH``/``DYLD_LIBRARY_PATH``).

Setting Up Input Files
----------------------

To use an external UMAT, configure your input files as follows:

**material.dat**

Specify the material name matching the plugin's ``name()`` return value:

.. code-block:: none

    Material
    Name    UMEXT
    Number_of_material_parameters   3
    Number_of_internal_variables    1
    
    #Orientation
    psi     0
    theta   0
    phi     0
    
    #Mechanical
    E       70000
    nu      0.3
    alpha   0.

For UMABA format, use ``Name  UMABA`` instead. For Ansys USERMAT, use ``Name  UMANS``.

**solver_essentials.inp**

.. code-block:: none

    Solver_type_0_Newton_tangent_1_RNL
    0
    Rate_type
    2

**path.txt**

Standard loading path definition (see :doc:`simulation/solver` for details):

.. code-block:: none

    #Initial_temperature
    290
    #Number_of_blocks
    1
    
    #Block
    1
    #Loading_type
    1
    #Control_type(NLGEOM)
    1
    #Repeat
    1
    #Steps
    2
    
    #Mode
    1
    #Dn_init 1.
    #Dn_mini 1.
    #Dn_inc 0.01
    #time
    1
    #prescribed_mechanical_state
    E 0.02
    S 0 S 0
    S 0 S 0 S 0
    #prescribed_temperature_state
    T 290

Testing External Plugins
------------------------

External UMAT plugins can be tested using the Python API:

.. code-block:: python

    import numpy as np
    from simcoon.solver import Solver, Block, StepMeca

    # Material properties for your external UMAT
    props = np.array([...])  # UMAT-specific properties

    # Define loading step (uniaxial tension to 2% strain)
    step = StepMeca(
        DEtot_end=np.array([0.02, 0, 0, 0, 0, 0]),
        control=['strain', 'stress', 'stress', 'stress', 'stress', 'stress'],
        Dn_init=1,
        Dn_inc=100
    )

    # Create block with external UMAT
    block = Block(
        steps=[step],
        umat_name='UMEXT',  # or 'UMABA' for Abaqus-compatible
        props=props,
        nstatev=1
    )

    # Run simulation
    solver = Solver(blocks=[block])
    history = solver.solve()

    # Access results
    for state in history:
        print(f"Strain: {state.Etot[0]:.4f}, Stress: {state.sigma[0]:.2f}")

**TUMEXT Test**

The C++ test located in ``test/Umats/UMEXT/`` validates the UMEXT plugin format.
Test data is located in ``testBin/Umats/UMEXT/data/``.

**TUMABA Test**

Located in ``test_extern/Umats/UMABA/``, this test validates the Abaqus UMAT compatibility layer.
Test data is located in ``testBin/Umats/UMABA/data/``.

.. note::

    The TUMABA test is located in ``test_extern/`` because it requires linking against an external 
    Fortran UMAT library, which is compiled separately from the main simcoon library.

Running the Tests
-----------------

To run the external UMAT tests:

.. code-block:: bash

    # Build simcoon with tests
    mkdir build && cd build
    cmake .. -DBUILD_TESTING=ON
    make
    
    # Run UMEXT test
    cd testBin/Umats/UMEXT
    ./TUMEXT
    
    # Run UMABA test (requires external library)
    cd ../../../testBin/Umats/UMABA
    ./TUMABA

Both tests apply a uniaxial tensile loading to 2% strain on a thermoelastic material (E=70000 MPa, ν=0.3) 
and verify that the computed stress-strain response matches the reference solution.
