import os
from astropy.io import fits

try:
    os.remove('list_calibrate')
except:
    print 'kitten'
    
try:
    os.system('rm c_dSC* -f')
except:
    print 'kitten'

for filename in os.listdir("."):
    if "fits" in filename:
        if "dSCIENCE" in filename:
            with fits.open(filename, ignore_missing_end=True) as hdus:
                head = hdus[0].header['HIERARCH ESO DPR TYPE']
                if head == 'OBJECT' or head == 'STD':
                    
                    airm_i = float(hdus[0].header['HIERARCH ESO TEL AIRM START'])
                    airm_f = float(hdus[0].header['HIERARCH ESO TEL AIRM END'])
                    airm = (airm_i + airm_f)/2
                    exptime = hdus[0].header['EXPTIME']
                    with open('list_airmass', 'a') as f:
                        f.write(filename+" "+str(airm)+" "+str(exptime)+"\n")

os.system('cat list_airmass | sort > list_calibrate')
os.system('rm list_airmass')
               
