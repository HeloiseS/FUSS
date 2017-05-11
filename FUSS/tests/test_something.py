import FUSS.datred as r
import numpy as np

def test_zero_angle():
    np.loadtxt(r.zero_angles, unpack = True, usecols=(0,1))
