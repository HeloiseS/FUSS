import FUSS.polmisc as misc
import matplotlib.pyplot as plt
import numpy as np
    
def test_getspctr():
    misc.get_spctr('dc_11hs_ep1_clean.flx')

def test_getpol():
    pol = misc.get_pol('dc_11hs_ep1.pol', wlmin=4500, wlmax=5600)
    
class TestPolData():
    def test_init(self):
        testObj = misc.PolData('dc_11hs_ep1.pol', wlmin=4500, wlmax=9000)
    
    def test_add_flux_data(self):
        # For some reason I need to create the object again here as not remembered 
        # from test_init()
        testObj = misc.PolData('dc_11hs_ep1.pol', wlmin=4500, wlmax=9000)
        testObj.add_flux_data('dc_11hs_ep1_clean.flx')
        
    def test_flu_n_pol(self):
        testObj = misc.PolData('dc_11hs_ep1.pol', wlmin=4500, wlmax=9000)
        testObj.flu_n_pol()
    
    def test_find_isp(self):
        testObj = misc.PolData('dc_11hs_ep1.pol', wlmin=4500, wlmax=9000)
        testObj.find_isp(wlmin=5000, wlmax=6000)
        
    def test_qu_plt(self):
        testObj = misc.PolData('dc_11hs_ep1.pol', wlmin=4500, wlmax=9000)
        testObj.qu_plt()
        plt.show()
        
    def test_add_isp1(self):
        dat = misc.PolData('dc_11hs_ep1.pol' , wlmin=4400)
        dat.add_isp(linearisp_params=[[0.000263, 0.000019], [-2.68, 0.14], 
                                      [0.000361, 0.000017], [-3.94, 0.12],
                                      [-0.00000261,-0.000001934]])
        gradq, errq = dat.gradq
        gradu, erru = dat.gradu
        assert np.isclose(gradq, 0.000263)
        assert np.isclose(errq, 0.000019)
        assert np.isclose(gradu, 0.000361)
        assert np.isclose(erru, 0.000017)
    
    def test_add_isp2(self):
        dat = misc.PolData('dc_11hs_ep1.pol' , wlmin=4400)
        dat.add_isp(constisp_params=[0.2, 0.01, 0.2, 0.01])
        assert dat.qisp == 0.2
        assert dat.qispr == 0.01
        assert dat.uisp == 0.2
        assert dat.uispr == 0.01
        
    def test_rmv_isp1(self):
        dat = misc.PolData('dc_11hs_ep1.pol' , wlmin=4400)
        dat.add_isp(linearisp_params=[[0.000263, 0.000019], [-2.68, 0.14], 
                                      [0.000361, 0.000017], [-3.94, 0.12],
                                      [-0.00000261,-0.000001934]])
        dat.rmv_isp()
        assert np.isclose(dat.p[0], 0.28586266553)
        assert dat.p0 is not None
    
    def test_rmv_isp2(self):
        dat = misc.PolData('dc_11hs_ep1.pol' , wlmin=4400)
        dat.add_isp(constisp_params=[0.2, 0.01, 0.2, 0.01])
        dat.rmv_isp()
        assert np.isclose(dat.p[0], 2.9419019437)
        assert dat.p0 is not None

