import envellope as env
import os

def test_pickle_phot():
    env.pickle_phot(2, 100, [0.99, 0.99])
    os.remove("*.pkl")

def test_avg_pol():
    test_phot = env.Envellope(0.99, 0.99, 100, 100)
    test_phot.avg_pol()
    assert test_phot.avg_q is not None and test_phot.avg_u is not None
