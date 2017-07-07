import FUSS.statistics as s
import numpy as np

def test_chi2():
    x = np.arange(0,10,0.1)
    y = x*2+5
    y_r = np.array([0.1]*len(y))
    chi = s.chi2(y, y_r, y+0.1)
    assert int(chi) == 100
    
