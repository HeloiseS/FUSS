import numpy as np
import FUSS.isp as isp
import FUSS.polmisc as F

FILENAME_POL = 'dc_11hs_ep1.pol'
FILENAME_SPCTR = '11hs_ep1_flux.txt'
poldat = F.PolData(FILENAME_POL, wlmin =4400)

def test_from_range():
    isp_values = isp.from_range(FILENAME_POL, wlmin = 5500, wlmax = 6000)
    assert int(isp_values[0]) == 2 
    # work for FILENAME_POL = 'dc_11hs_ep1.pol' and wlmin = 5500, wlmax = 6000
    
# This technique didn't even work as it was too reliant on user input.     
#def test_from_emline():
#    isp_values = isp.from_emline(FILENAME_POL, FILENAME_SPCTR, wlmin = 5500)

def test_linear_isp1():
    isp_values = isp.linear_isp(poldat.wlp, [0.000263, 0.000019], [-2.68, 0.14], [0.000361, 0.000017], [-3.94, 0.12], -0.00000261,-0.000001934)
    qisp, qisp_r, uisp, uisp_r = isp_values
    assert np.isclose(qisp[0], -1.518948848920000)
    # works for FILENAME_POL = 'dc_11hs_ep1.pol' provided in tests directory

def test_linear_isp2():
    newvals, ispval  = isp.linear_isp(poldat.wlp, [0.000263, 0.000019], [-2.68, 0.14],
     [0.000361, 0.000017], [-3.94, 0.12], -0.00000261, -0.000001934, 
     q = poldat.q, qr = poldat.qr, u = poldat.u, ur = poldat.ur )
    wl, p, pr, q, qr, u, ur, a, ar = newvals
    assert np.isclose(p[0], 0.28586266553)
    
def test_const_isp():
    newvals = isp.const_isp(poldat.wlp, 0.2, 0.01, 0.2, 0.01,  # arbitrary values of isp for which I know wht p[0] is
    poldat.q, poldat.qr, poldat.u, poldat.ur)
    wl, p, pr, q, qr, u, ur, a, ar = newvals
    assert np.isclose(p[0], 2.941901943733591)
    # only works for FILENAME_POL = 'dc_11hs_ep1.pol' provided in tests directory
    # and the arbitrary values of ISP q = u = 0.2 +/- 0.01 that I've picked.




