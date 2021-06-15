import numpy as np
import Tarma2numpy
import unittest

print(np.__version__)
print(np.__path__)

#np.set_printoptions(precision=2, suppress=True)

x = np.random.randint(20, size=6, dtype='int32')
y = Tarma2numpy.test_vec_int(x)
a = np.testing.assert_array_equal(x,y)

x = np.random.randint(20, size=(3, 2), dtype='int32')
y = Tarma2numpy.test_mat_int(x)
np.testing.assert_array_equal(x,y)

x = np.random.rand(6)
y = Tarma2numpy.test_vec_double(x)
np.testing.assert_array_equal(x,y)

print("Test 1")
x = np.random.rand(3,6)
print(x)
print(x.flags)
print(x.strides)
y = Tarma2numpy.test_mat_double(x)
print(y)
print(y.flags)
print(y.strides)
np.testing.assert_array_equal(x,y)

print("Test 2")
x2 = np.asfortranarray(x)
print(x2)
print(x2.flags)
print(x2.strides)
y2 = Tarma2numpy.test_mat_double(x2)
print(y2)
print(y2.flags)
print(y2.strides)
np.testing.assert_array_equal(x2,y2)

print("Test 3")
x3 = np.random.rand(4,6)
print(x3)
print(x3.flags)
print(x3.strides)
y3 = Tarma2numpy.test_mat_inplace_double(x3)
print(y3)
print(y3.flags)
print(y3.strides)
print(y3[0,0])
y3 = y3+0.
np.testing.assert_array_equal(x3,y3)

print("Test 4")
x4 = np.asfortranarray(x3)
print(x4)
print(x4.flags)
print(x4.strides)
y4 = Tarma2numpy.test_mat_inplace_double(x4)
print(y4)
print(y4.flags)
print(y4.strides)
y4 = y4+0.
print(y4.flags)
np.testing.assert_array_equal(x4,y4)

#x2 = np.asfortranarray(x)
#print(x2.flags)
#print(x2.strides)
#print('x2 = ', x2)
#y2 = Tarma2numpy.test_mat_inplace_double(x2)
#print('y2 = ', y2)
#print(y.flags)
#print(y2.strides)
#np.testing.assert_array_equal(x2,y2)

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





