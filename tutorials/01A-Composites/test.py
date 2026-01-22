import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import simcoon as sim
import os
import itertools

dir = os.path.dirname(os.path.realpath("__file__"))

nstatev = 0

nphases = 2  # The number of phases
num_file = 0  # The num of the file that contains the subphases
int1 = 50
int2 = 50
n_matrix = 0

props = np.array([nphases, num_file, int1, int2, n_matrix], dtype="float")

NPhases_file = dir + "/keys/Nellipsoids0.dat"
NPhases = pd.read_csv(NPhases_file, delimiter=r"\s+", index_col=False, engine="python")

psi_rve = 0.0
theta_rve = 0.0
phi_rve = 0.0

umat_name = "MIMTN"

L = sim.L_eff(umat_name, props, nstatev, psi_rve, theta_rve, phi_rve)
p = sim.L_iso_props(L)
print(p)
