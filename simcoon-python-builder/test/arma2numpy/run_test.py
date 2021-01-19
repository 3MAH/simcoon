import numpy as np
import Tarma2numpy
import unittest

print(np.__version__)
print(np.__path__)

x = np.random.randint(20, size=6)
y = Tarma2numpy.test_vec_int(x)
a = np.testing.assert_array_equal(x,y)

x = np.random.randint(20, size=(3, 2))
y = Tarma2numpy.test_mat_int(x)
np.testing.assert_array_equal(x,y)

x = np.random.rand(6)
y = Tarma2numpy.test_vec_double(x)
np.testing.assert_array_equal(x,y)

x = np.random.rand(6,6)
y = Tarma2numpy.test_mat_double(x)
np.testing.assert_array_equal(x,y)

words = ["Hello", "World", "!"]
y = Tarma2numpy.test_vector_list_string(words)
if y == words:
    print('test_vector_list_string: Success')
else:
    print('test_vector_list_string: Fail')

words = [0, 1, 2]
y = Tarma2numpy.test_vector_list_int(words)
if y == words:
    print('test_vector_list_int: Success')
else:
    print('test_vector_list_int: Fail')

words = [4.6, 1.1, 2.2]
y = Tarma2numpy.test_vector_list_double(words)
if y == words:
    print('test_vector_list_double: Success')
else:
    print('test_vector_list_double: Fail')

