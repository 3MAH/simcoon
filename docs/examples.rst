Examples
========

Run a simulation
----------------

Here is a simple example to run your first simulation using SMART+.
Open the file Path.txt in the folder Control
You should find something that look like that

.. code-block:: none

    #Initial_temperature
    293.5
    #Number_of_blocks
    1

    #Block
    1
    #Loading_type
    1
    #Repeat
    1
    #Steps
    1

    #Mode
    1
    #Dn_init 1.
    #Dn_mini 0.1
    #Dn_inc 0.01
    #time
    30.
    #mechanical_state
    E 0.3 
    S 0 S 0
    S 0 S 0 S 0
    #temperature_state
    T 293.5

The first part of the file describe the initial temperature conditions, under the tag #Initial_temperature.

Just below, under the tag #Number_of_blocks, you define the number of blocks. Here we will start with a single block, so this value is set to 1.
The next part is to define the first block:
#Block defines the block number
#Loading_type defines the physical problem to solve, which is:
1 – mechanical; 2 – thermomechanical
#Repeat is the number of time the block is repeated
#Steps is the number of steps of the block

The next part of the file defines the steps of the first block. It always starts with the mode of the step (#mode), which is:
1 – linear; 2 – sinusoidal; 3 – tabular (from a file)

In this example we will consider that the step mode is linear. We therefore need to set up the following

#. The mode of the step (under #mode)
#. The initial size of the first increment (usually 1.), under #Dn_init
#. The minimal size of an increment (usually less than 1.), under #Dn_mini
#. The size of the increment as a fraction of the step $\delta n$, under #Dn_inc. (#Dn_inc 0.01 means that 100 increments will be utilized to simulate the step)
#. The time $\Delta t$ of the step (under #time). Note that the increment of time for any increment is $\delta t = \Delta t \delta n$
#. The mechanical loading stage at the end of the step (#mechanical_state)
The elements are organized such that either stress or strain components are defined in the following order:
11
12 22
13 23 33
The letter 'S' in front of any component means that a stress control is considered in that direction, and the letter 'E' stands for a strain control. Note that those values indicate the state at the end of the step
#. The thermal loading stage at the end of the step (#temperature_state)
For mechanical loading, the letter 'T' is followed by the temperature at the end of the step. For thermomechanical loading, either the final temperature can be considered (with the letter 'T'), or the thermal flux (with the letter 'Q'). This last quantity is defined as the rate of heat that flows to the material representative volume element considered.

If the #mode is set to 2 - sinusoidal, a sinusoidal evolution is considered automatically between the state of the previous step and the final values indicated in the current step.

If the #mode is set to 3 - tabular, a prescribed evolution is considered, and a file that contains such prescribed evolution must be indicated. In that case, the step block has to be defined like:

.. code-block:: none

	#Mode
	3
	#File
	tabular_file.txt
	#Dn_init 1.
	#Dn_mini 0.01
	#Consigne
	S
	0  S
	0  0  0
	#T_is_set
	0

In the following example, a biaxial test in the directions 11 and 22 is considered, with a stress control. The temperature is not set, which means that it is constant throughout the step and keep its value from the previous step (or the intial temperature if this is the first step). Note that the time is always indicated in the tabular_file.txt.
The struture of the tabular file will be the following:

.. code-block:: none

	0	0.0	10	10		
	1	0.01	20	20
	2	0.02	30	30
	3	0.03	30	30
	...

The columns define the quantities in the followin order : #ninc, #time, #S11, #S22.
The order of the mechanical quantities is always 11,12,22,13,23,33, and if the temperature is set (with the letter 'T' instead of '0'), the following order is always considered: #ninc, #time, #T, #S11, #S22 in the case of the biaxial loading.

.. code-block:: none

	0	0.0	293.15	10	10		
	1	0.01	294.15	20	20
	2	0.02	295.15	30	30
	3	0.03	296.15	30	30
	...



Set up a micro mechanical model
-------------------------------

The first thing you want to do when setting up a micro mechanical model is to define the microstructure. At a certain scale, you should inform the model about the phases, their volume fraction, geometry and their properties.

First, in the file data/material.dat, you need to enter the material properties corresponding to the micro mechanical model you selected:

For Mori-Tanaka and Self-Consistent: 4 material parameters (and a consequent number of state_variables)



#. props(0) : Number of phases
#. props(1) : File number that stores the microstructure properties
#. props(2) : Number of integration points in the 1 direction
#. props(3) : Number of integration points in the 2 direction

For Periodic layers: 2 material parameters (and a consequent number of state_variables)


#. props(0) : Number of phases
#. props(1) : File number that stores the microstructure properties

The file data/material.dat should look like this for a 2-phase material using a Mor-Tanaka model:

.. code-block:: none

	Material
	Name    MIMTN
	Number_of_material_parameters   4
	Number_of_internal_variables    10000

	#Thermal
	density 1.12
	c_p   1.64

	#Mechancial
	nphases 2
	file_number 0
	nItg1 20
	nItg2 20

The density and specific heat capacity c_p are utilized only if you want to solve a thermomechanical boundary-value problem.

The file number represents the number of the Nphases[i].dat file, where [i] is replaced by the number value. In this case we should fill the file Nphases0.dat, which looks like this:

.. code-block:: none

    Number  Coatingof  umat   c    phi_mat  theta_mat  psi_mat  a1  a2  a3  phi_geom  theta_geom  psi_geom  nprops  nstatev  props
    0       0          ELISO  0.8  0        0          0        1   1   1   0.        0.          0.        3       1        3000    0.4   1.E-5
    1       0          ELISO  0.2  0        0          0        1   1   1   0.        0.          0.        3       1        70000   0.4   1.E-5

Note that for Mori-Tanaka the first phase in the file should always be the matrix.
The characteristics of the phases are described below:

#. Number : The number of the phase
#. Coatingof : If the model is a coating of an other phase. 0 if the phase is not a coating
#. umat : Constitutive model considered
#. c : Volume fraction of the phase
#. phi_mat: First Euler angle corresponding to the material orientation
#. theta_mat: Second Euler angle corresponding to the material orientation
#. psi_mat: Third Euler angle corresponding to the material orientation
#. a1:
#. a2:
#. a3:
#. phi_geom: First Euler angle corresponding to the ellipsoid orientation
#. theta_geom: Second Euler angle corresponding to the ellipsoid orientation
#. psi_geom: Third Euler angle corresponding to the ellipsoid orientation
#. npros: Number of material properties
#. nstatev: Number of scalar internal variables
#. props: The list of material properties

For a wide majority of composites, the orientation of the material coincides with the orientation of the reinforcement (For instance transversely isotropic carbon fibers).
However, for metallic polycristals, the two materials systems have to be considered to separate the orientation of the lattice with the orientation of the ellipsoid that represent a grain.
This version of SMART+ currently does not support coated inclusions, but the files Nphase[i].dat is prepared so that you can easily add this to a custom micromechancial model.

Note that the Euler system reference utilised (3-1-3 for the most common) is defined in the parameter.hpp file. For instance this system is defined by default in the parameter.hpp:

.. code-block:: none

    #ifndef axis_psi
    #define axis_psi 3
    #endif

    #ifndef axis_theta
    #define axis_theta 1
    #endif

    #ifndef axis_phi
    #define axis_phi 3
    #endif

In the example here we are defining a 2-phase composite, with spherical reinforcements, considering two phases:

#. An epoxy matrix, 80% volume, with E=3000MPa and nu=0.4, and alpha=1.E-5
#. Aluminium reinforcements: 20% volume, with E=70000MPa and nu=0.3, and alpha=5.E-5

Once these files have been set up, you can run a simulation using the classical solver.
