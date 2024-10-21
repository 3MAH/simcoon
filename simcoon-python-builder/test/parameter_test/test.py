import pytest
import numpy as np
import numpy.typing as npt
from simcoon import simmit as sim
from simcoon import parameter as par

dir = os.path.dirname(os.path.realpath('__file__'))

nstatev = 0

nphases = 2 #The number of phases
num_file = 0 #The num of the file that contains the subphases
int1 = 50
int2 = 50
n_matrix = 0

props = np.array([nphases, num_file, int1, int2, n_matrix],  dtype='float')

path_data = dir + '/data'
path_keys = dir + '/keys'
pathfile = 'path.txt'

param_list = par.read_parameters()

psi_rve = 0.
theta_rve = 0.
phi_rve = 0.


from simcoon import parameter as par

dir = os.path.dirname(os.path.realpath('__file__'))
pylab.rcParams['figure.figsize'] = (18.0, 8.0) #configure the figure output size

nstatev = 0

nphases = 2 #The number of phases
num_file = 0 #The num of the file that contains the subphases
int1 = 50
int2 = 50
n_matrix = 0

props = np.array([nphases, num_file, int1, int2, n_matrix],  dtype='float')

#NPhases_file = dir + '/keys/Nellipsoids0.dat'
#NPhases = pd.read_csv(NPhases_file, delimiter=r'\s+', index_col=False, engine='python')
#NPhases[::]

path_data = dir + '/data'
path_keys = dir + '/keys'
pathfile = 'path.txt'

param_list = par.read_parameters()

psi_rve = 0.
theta_rve = 0.
phi_rve = 0.


param_list = par.read_parameters(ncjlknl;jw)
