#!/usr/bin/python

import numpy as np
import Tarma2numpy

x = np.random.randint(20, size=6)
y = Tarma2numpy.test_vec_int(x)
np.testing.assert_array_equal(x,y)

x = np.random.randint(20, size=(3, 2))
y = Tarma2numpy.test_mat_int(x)
np.testing.assert_array_equal(x,y)

x = np.random.rand(6)
y = Tarma2numpy.test_vec_double(x)
np.testing.assert_array_equal(x,y)

x = np.random.rand(3,2)
y = Tarma2numpy.test_mat_double(x)
np.testing.assert_array_equal(x,y)

words = ["Hello", "World", "!"]
y = Tarma2numpy.test_vector_list_string(words)
np.testing.assert_array_equal(words,y)

words = [0, 1, 2]
y = Tarma2numpy.test_vector_list_int(words)
np.testing.assert_array_equal(words,y)

words = [4.6, 1.1, 2.2]
y = Tarma2numpy.test_vector_list_double(words)
np.testing.assert_array_equal(words,y)