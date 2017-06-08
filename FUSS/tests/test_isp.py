import FUSS.isp as isp

FILENAME_POL = 'dc_11hs_ep1.pol'

def test_from_range():
    isp_values = isp.from_range(FILENAME_POL, wlmin = 5500, wlmax = 6000)
    assert int(isp_values[0]) == 2 
    # work for FILENAME_POL = 'dc_11hs_ep1.pol' and wlmin = 5500, wlmax = 6000

