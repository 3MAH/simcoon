#!/usr/bin/python

import numpy as np
#import Tarma2numpy
from simcoon import simmit as sim

const1 = sim.constants(2,1)
const2 = sim.constants(3,2)
const_list = [const1, const2]
print(const_list[0].number)
print(const_list[1].number)
#y = Tarma2numpy.test_vector_list_constants(const_list)
#print(y[0].number)
#print(y[1].number)

param1 = sim.parameters(1, 0., 4.)
#param2 = sim.parameters(2, 0., 4., '@2p', 1, ['source.exp'])
print(param1.number)

F0 = np.zeros((3,3))
F1 = np.zeros((3,3))
sigma = np.zeros(6)
etot = np.zeros(6)
statev = np.zeros(1)
T = 273.15
DT = 2.

sv_1 = sim.state_variables()
sv_2 = sim.state_variables(etot, F0, F1, sigma, statev, T, DT)

print(sv_2.etot)
