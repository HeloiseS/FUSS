import FUSS.polmisc as misc
    
def test_getspctr():
    misc.get_spctr('tests/dc_11hs_ep1_clean.flx')

def test_getpol():
    pol = misc.get_pol('tests/dc_11hs_ep1.pol', wlmin=4500, wlmax=5600)
    
