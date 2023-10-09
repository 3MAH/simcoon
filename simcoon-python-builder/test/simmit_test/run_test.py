import pytest
import numpy as np
import numpy.typing as npt
from simcoon import simmit as sim

@pytest.fixture(scope='session')
def L_gt() -> npt.NDArray[np.float_]:
    L_gt = np.array([
        [1.34615385, 0.57692308, 0.57692308, 0., 0., 0.],
        [0.57692308, 1.34615385, 0.57692308, 0., 0., 0.],
        [0.57692308, 0.57692308, 1.34615385, 0., 0., 0.],
        [0., 0., 0., 0.38461538, 0., 0.],
        [0., 0., 0., 0., 0.38461538, 0.],
        [0., 0., 0., 0., 0., 0.38461538]
    ])
    return L_gt

def test_L_iso(L_gt):
    """Test L_iso with E=1., nu=0.3."""
    L = sim.L_iso([1., 0.3], 'Enu')
    assert np.allclose(L, L_gt)