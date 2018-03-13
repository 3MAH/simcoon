#!/usr/bin/python

import numpy as np
import Tarma2numpy
from simcoon import identify as iden

const1 = iden.constants(2,1)
const2 = iden.constants(3,2)
const_list = [const1, const2]
print(const_list[0].number)
print(const_list[1].number)
y = Tarma2numpy.test_vector_list_constants(const_list)
print(y[0].number)
print(y[1].number)