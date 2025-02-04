# -*- coding: utf-8 -*-

# This is a very simple matplotlib script to plot results from 'results_job.txt'
# from matplotlib import cm
# from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np
import os

fig = plt.figure()
path = os.path.dirname(os.path.realpath("__file__")) + "/results/"

PL = path + "results_job_global-0.txt"

# valid = data_path + 'valid.txt'

e11, e22, e33, e12, e13, e23, s11, s22, s33, s12, s13, s23 = np.loadtxt(
    PL, usecols=(8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19), unpack=True
)
# time, T, Q, r, Wm, Wm_r, Wm_ir, Wm_d, Wt, Wt_r, Wt_ir = np.loadtxt(PL, usecols=(4,5,6,7,20,21,22,23,24,25,26), unpack=True)
time, T, Q, r = np.loadtxt(PL, usecols=(4, 5, 6, 7), unpack=True)

Wm, Wm_r, Wm_ir, Wm_d = np.loadtxt(PL, usecols=(20, 21, 22, 23), unpack=True)

ax = fig.add_subplot(2, 2, 1)
plt.grid(True)
plt.tick_params(axis="both", which="major", labelsize=15)
plt.xlabel("Strain", size=15)
plt.ylabel("Stress (MPa)", size=15)
plt.plot(e11, s11, c="black", label="direction 1")
# plt.xlim(-0.051,0.051)
plt.legend(loc=2)

ax = fig.add_subplot(2, 2, 2)
plt.grid(True)
plt.tick_params(axis="both", which="major", labelsize=15)
plt.xlabel("time (s)", size=15)
plt.ylabel("T", size=15)
plt.plot(time, T, c="black", label="Temperature")
plt.legend(loc=2)


ax = fig.add_subplot(2, 2, 4)
plt.grid(True)
plt.tick_params(axis="both", which="major", labelsize=15)
plt.xlabel("time (s)", size=15)
plt.ylabel("Wm", size=15)
plt.plot(time, Wm, c="black", label="Wm")
plt.plot(time, Wm_r, c="green", label="Wm_r")
plt.plot(time, Wm_ir, c="blue", label="Wm_ir")
plt.plot(time, Wm_d, c="red", label="Wm_d")
plt.legend(loc=2)

plt.show()
