import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import simcoon as sim
import os
import itertools

dir = os.path.dirname(os.path.realpath('__file__'))

nstatev = 0

nphases = 2 #The number of phases
num_file = 0 #The num of the file that contains the subphases
int1 = 50
int2 = 50
n_matrix = 0

props = np.array([nphases, num_file, int1, int2, n_matrix],  dtype='float')

NPhases_file = dir + '/keys/Nellipsoids0.dat'
NPhases = pd.read_csv(NPhases_file, delimiter=r'\s+', index_col=False, engine='python')
#NPhases[::]

path_data = dir + '/data'
path_keys = dir + '/keys'
pathfile = 'path.txt'

nparams = 4
param_list = sim.read_parameters(nparams)

psi_rve = 0.
theta_rve = 0.
phi_rve = 0.

concentration = np.arange(0.,0.51,0.01)

E_MT = np.zeros(len(concentration))
umat_name = 'MIMTN'
for i, x in enumerate (concentration):
 
    param_list[1].value = x
    param_list[0].value = 1.-x
    
    sim.copy_parameters(param_list, path_keys, path_data)
    sim.apply_parameters(param_list, path_data)

    L = sim.L_eff(umat_name, props, nstatev, psi_rve, theta_rve, phi_rve)
    p = sim.L_iso_props(L)
    print(p)
    E_MT[i] = p[0]

    
E_SC = np.zeros(len(concentration))
umat_name = 'MISCN'
for i, x in enumerate (concentration):
 
    param_list[1].value = x
    param_list[0].value = 1.-x
    
    sim.copy_parameters(param_list, path_keys, path_data)
    sim.apply_parameters(param_list, path_data)

    L = sim.L_eff(umat_name, props, nstatev, psi_rve, theta_rve, phi_rve)
    p = sim.L_iso_props(L)
    E_SC[i] = p[0]

print(props)

fig = plt.figure()
plt.plot(concentration,E_MT, c='blue')
plt.plot(concentration,E_SC, c='red')
expfile = path_data + '/' + 'E_exp.txt'
c,E = np.loadtxt(expfile, usecols=(0,1), unpack=True)
plt.plot(c,E,linestyle='None', marker='x', color='black', markersize=10)
plt.show()


