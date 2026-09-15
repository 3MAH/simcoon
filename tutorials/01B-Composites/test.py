import numpy as np
import matplotlib.pyplot as plt
import simcoon as sim
import os
from simcoon.solver.micromechanics import Ellipsoid, to_phase_dicts

dir = os.path.dirname(os.path.realpath('__file__'))

nstatev = 0

nphases = 2 #The number of phases
int1 = 50
int2 = 50
n_matrix = 0

# props[1] used to be the number of the Nellipsoids<N>.dat file to read. The phases are
# now handed to L_eff directly, so that slot is kept only for the layout of the vector.
props = np.array([nphases, 0, int1, int2, n_matrix],  dtype='float')

path_data = dir + '/data'

# The volume fractions used to be swept by substituting @0p / @1p into a copy of
# keys/Nellipsoids0.dat; they are set on the phase objects instead.
matrix = Ellipsoid(number=0, umat_name='ELISO', save=1, concentration=0.8, nstatev=1,
                   props=np.array([2250., 0.19, 8.8e-5]))
reinforcement = Ellipsoid(number=1, umat_name='ELISO', save=1, concentration=0.2, nstatev=1,
                          props=np.array([73000., 0.19, 0.5e-6]))

psi_rve = 0.
theta_rve = 0.
phi_rve = 0.

concentration = np.arange(0.,0.51,0.01)

E_MT = np.zeros(len(concentration))
umat_name = 'MIMTN'
for i, x in enumerate (concentration):

    reinforcement.concentration = x
    matrix.concentration = 1.-x

    L = sim.L_eff(umat_name, props, nstatev, psi_rve, theta_rve, phi_rve,
                  to_phase_dicts([matrix, reinforcement]))
    p = sim.L_iso_props(L).flatten()
    print(p)
    E_MT[i] = p[0]


E_SC = np.zeros(len(concentration))
umat_name = 'MISCN'
for i, x in enumerate (concentration):

    reinforcement.concentration = x
    matrix.concentration = 1.-x

    L = sim.L_eff(umat_name, props, nstatev, psi_rve, theta_rve, phi_rve,
                  to_phase_dicts([matrix, reinforcement]))
    p = sim.L_iso_props(L).flatten()
    E_SC[i] = p[0]

print(props)

fig = plt.figure()
plt.plot(concentration,E_MT, c='blue')
plt.plot(concentration,E_SC, c='red')
expfile = path_data + '/' + 'E_exp.txt'
c,E = np.loadtxt(expfile, usecols=(0,1), unpack=True)
plt.plot(c,E,linestyle='None', marker='x', color='black', markersize=10)
plt.show()
