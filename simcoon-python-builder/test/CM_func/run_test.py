import numpy as np
import TCM_func
import unittest

print(np.__version__)
print(np.__path__)

E = 700000.
nu = 0.2

L = TCM_func.L_iso(E,nu,"Enu")
print(L)
