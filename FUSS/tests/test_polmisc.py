import FUSS.polmisc as misc
import matplotlib.pyplot as plt
    
def test_getspctr():
    misc.get_spctr('dc_11hs_ep1_clean.flx')

def test_getpol():
    pol = misc.get_pol('dc_11hs_ep1.pol', wlmin=4500, wlmax=5600)
    
class TestPolData():
    def test_init(self):
        testObj = misc.PolData('test', 'dc_11hs_ep1.pol', wlmin=4500, wlmax=9000)
    
    def test_add_flux_data(self):
        # For some reason I need to create the object again here as not remembered 
        # from test_init()
        testObj = misc.PolData('test', 'dc_11hs_ep1.pol', wlmin=4500, wlmax=9000)
        testObj.add_flux_data('dc_11hs_ep1_clean.flx')
        
    def test_flu_n_pol(self):
        testObj = misc.PolData('test', 'dc_11hs_ep1.pol', wlmin=4500, wlmax=9000)
        testObj.flu_n_pol()
    
    def test_find_isp(self):
        testObj = misc.PolData('test', 'dc_11hs_ep1.pol', wlmin=4500, wlmax=9000)
        testObj.find_isp(wlmin=5000, wlmax=6000)
        
    def test_qu_plt(self):
        testObj = misc.PolData('test', 'dc_11hs_ep1.pol', wlmin=4500, wlmax=9000)
        testObj.qu_plt()
        plt.show()
        
