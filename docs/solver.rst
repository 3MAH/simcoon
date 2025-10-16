Use the solver
================================

Elastic tensile test
--------------------

Probably the first thing you would like to do with Simcoon is to simulate the mechanical response corresponding of a simple tension test, considering an elastic isotropic material:

We first import *simmit* (the python simulation module of simcoon) and *numpy* 

.. code-block:: python

	import numpy as np
	from simcoon import simmit as sim

Next we shall define the material constitutive law to be utilized and the associated material properties. We will pass them as a numpy array:

.. code-block:: python

	umat_name = 'ELISO' #This is the 5 character code for the elastic-isotropic subroutine
	nstatev = 1 #The number of scalar variables required, only the initial temperature is stored here to consider the thermal expansion if temperature changes.

	E = 700000. #The Young modulus
	nu = 0.2 #The Poisson coefficient
	alpha = 1.E-5 #The coefficient of thermal expansion

	#Three Euler angles to represent the material orientation with respect to the reference basis (in which the loading is expressed)
	psi_rve = 0.
	theta_rve = 0.
	phi_rve = 0.

	#Solver_type define the solver strategy (only a classical newton scheme is actually implemented for now), and the corate_type define the type of corotational spin rate (0 for Jauman, 1 for Green-Naghdi, 2 for logarithmic)
	solver_type = 0
	corate_type = 2

	props = np.array([E, nu, alpha])

We shall then define, the location of the data input files (i.e., to define the loading path), and the results outut file and location:

.. code-block:: python

	path_data = 'data'
	path_results = 'results'
	pathfile = 'path.txt'
	outputfile = 'results_ELISO.txt'


The last part is to define the loading path. Further details about this file is given in the example just below but for this first example you could just create a folder 'data' and create a text file named 'path.txt' with the following inside:

.. code-block:: none

	#Initial_temperature
	293.5
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
	1

	#Mode
	1
	#Dn_init 1.
	#Dn_mini 0.1
	#Dn_inc 0.01
	#time
	30.
	#mechanical_state
	E 0.01 
	S 0 S 0
	S 0 S 0 S 0
	#temperature_state
	T 293.5

The latter correspond to a pure strain-controlled tension test in the direction *1* up to 1% strain, at the temperature 293.5K.

And finally we can call the solver function:

.. code-block:: python

	sim.solver(
    umat_name,
    props,
    nstatev,
    psi_rve,
    theta_rve,
    phi_rve,
    solver_type,
    corate_type,
    path_data,
    path_results,
    pathfile,
    outputfile,
)

You can now run your just created python file (you could also create a jupyter notebook, or run the code for linear isotropic materials  that you can download in the examples). You will now find in the 'results' folder a file named *results_ELISO.txt*. Have a look at the existing notebook or in the documentation to know how to analyse the result file.

Define the loading path
-----------------------------------------

If you navigate into the 'data' folder, you shall find the following files:

#. path.txt, which is structured as this:

.. code-block:: none

    #Initial_temperature
    293.5
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

The first part of the file describe the initial temperature conditions, under the tag **#Initial_temperature**.

Just below, under the tag **#Number_of_blocks**, you define the number of blocks. Here we will start with a single block, so this value is set to 1.
The next part is to define the first block:
**#Block** defines the block number
**#Loading_type** defines the physical problem to solve, which is:
1 for mechanical; 2 for thermomechanical
**#Control_type** defines if the mechanical part of the problem to solve is controlled from infinitesimal strains/stress, or if a finite deformation framework is utilized.
1 for infinitesimal strains/stress; 2 for finite deformation using Lagrangian control (Green-Lagrange strain / Piola-Kirchoff II stress); 3 for finite deformation using logarithmic strain / Kirchoff stress
**#Repeat** is the number of time the block is repeated
**#Steps** is the number of steps of the block

The next part of the file defines the steps of the first block. It always starts with the mode of the step (#mode), which is:
1 for linear; 2 for sinusoidal; 3 for tabular (from a file)

In this example we will consider that the step mode is linear. We therefore need to set up the following

#. The mode of the step (under **#mode**)
#. The initial size of the first increment (usually 1.), under **#Dn_init**
#. The minimal size of an increment (usually less than 1.), under **#Dn_mini**
#. The size of the increment as a fraction of the step $\delta n$, under **#Dn_inc**. (**#Dn_inc 0.01** means that 100 increments will be utilized to simulate the step)
#. The time :math:`\Delta t` of the step (under **#time**). Note that the increment of time for any increment is :math:`\delta t = \Delta t \delta n`
#. The mechanical loading stage at the end of the step (**#mechanical_state**)

If Control_type=1, the elements are organized such that either stress or strain components are defined in the following order:
11
12 22
13 23 33
The letter 'S' in front of any component means that a stress control is considered in that direction, and the letter 'E' stands for a strain control. Note that those values indicate the state at the end of the step
#. The thermal loading stage at the end of the step (**#temperature_state**)
If Control_type=2, the elements are organized such that either the first Piola stress $\Sigma$ or displacement gradient :math:`\nabla u` components are defined in the following order:
11 12 13
21 22 23
31 32 33
The letter **'S'** in front of any component means that a stress control is considered in that direction, and the letter **'E'** stands for a kinematic control. Note that those values indicate the state at the end of the step

#. The thermal loading stage at the end of the step (**#temperature_state**)
For mechanical loading, the letter 'T' is followed by the temperature at the end of the step. For thermomechanical loading, either the final temperature can be considered (with the letter 'T'), or the thermal flux (with the letter 'Q'). This last quantity is defined as the rate of heat that flows to the material representative volume element considered.

If the **#mode** is set to 2 - sinusoidal, a sinusoidal evolution is considered automatically between the state of the previous step and the final values indicated in the current step.

If the **#mode** is set to 3 - tabular, a prescribed evolution is considered, and a file that contains such prescribed evolution must be indicated. In that case, the step block has to be defined like:

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

	0	0.0		10	10		
	1	0.01	20	20
	2	0.02	30	30
	3	0.03	30	30
	...

The columns define the quantities in the following order : **ninc**, **time**, **S11**, **S22**.
The order of the mechanical quantities is always 11,12,22,13,23,33, and if the temperature is set (with the letter 'T' instead of '0'), the following order is always considered: **ninc**, **time**, **T**, **S11**, **S22** in the case of the biaxial loading.

.. code-block:: none

	0	0.0		293.15	10	10		
	1	0.01	294.15	20	20
	2	0.02	295.15	30	30
	3	0.03	296.15	30	30
	...
``