import numpy as np
import matplotlib.pyplot as plt
import simcoon as sim
import os
from simcoon.solver.micromechanics import Ellipsoid

try:
    script_dir = os.path.dirname(os.path.abspath(__file__))
except NameError:   # executed by the docs gallery, from the example's directory
    script_dir = os.getcwd()

nstatev = 0

int1 = 50  # Integration points of the Eshelby integrals, first direction
int2 = 50  # Integration points of the Eshelby integrals, second direction
n_matrix = 0  # Index of the matrix phase in the list of phases

# Mori-Tanaka props: [int1, int2, n_matrix]; the phases themselves are passed as objects
props = np.array([int1, int2, n_matrix], dtype="float")

path_data = os.path.join(script_dir, 'data')

# The two phases; their volume fractions are swept below, on the objects
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

    L = sim.L_eff(umat_name, props, nstatev, orientation=(psi_rve, theta_rve, phi_rve),
                  phases=[matrix, reinforcement])
    p = sim.L_iso_props(L).flatten()
    print(p)
    E_MT[i] = p[0]


E_SC = np.zeros(len(concentration))
umat_name = 'MISCN'
for i, x in enumerate (concentration):

    reinforcement.concentration = x
    matrix.concentration = 1.-x

    L = sim.L_eff(umat_name, props, nstatev, orientation=(psi_rve, theta_rve, phi_rve),
                  phases=[matrix, reinforcement])
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
