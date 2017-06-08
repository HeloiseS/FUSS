"""
08 - June - 2017 / H. F. Stevance / fstevance1@sheffield.ac.uk

datred.py is a module created as part of the FUSS package to help with the data reduction of spectropolarimetric
data (at the present time only used with FORS2 data)

Pre-requisites:
-----------------
os, astropy.io, numpy, math, matplotlib.pyplot, pysynphot, scipy.special

Attributes:
----------
zero_angles : string
    Path to the text file containing the chromatic zero angles for FORS2 (needs updating for you own system).
    Can be found at: http://www.eso.org/sci/facilities/paranal/instruments/fors/inst/pola.html

Functions:
---------
sort_red()
    Creates back-up of the compressed files, uncompresses them, sorts them and re-names them according to the naming
    convention required to use with my .cl files during data reduction in IRAF. Returns: None

info()
    Creates a text file containing useful information on the calibration and data files. Must use in folder containing
    the uncompressed FITS files. Returns: None
-> Output File Format: Filename, ESO label, Retarder Plate Angle, Exposure time,Airmass, Grism, Bin, number of Pixels,
1/Gain, Read Out Noise, Date.

hwrpangles()
    Creates the file used by lin_specpol to know which images correspond to which HWRP angle. It creates separate file
    for the CCSN, zero pol std, and polarised std. The output files are made of 4 columns containing the numbers of
    images corresponding to the 0, 22.5, 45 and 67.5 degree retarder plate angles. 1 set per line.

lin_specpol():
    Calculates the Stokes parameters and P.A of a data set and writes them out in a text file. Also produces
    a plot showing p, q, u, P.A and Delta epsilon_q and Delta epsilon_u. The plots is not automatically saved.

circ_specpol():
    Calculates the circular polarisation v and the delta epsilon. Plot not automatically saved.

flux_spectrum():
    Combines all the flux calibrated apertures to create the flux spectrum.
"""
# TODO: internal functions name _internal_function()
import os
from astropy.io import fits
import numpy as np
import math as m
import matplotlib.pyplot as plt
import pysynphot as S
from scipy import special as special
from astropy.utils.data import get_pkg_data_filename

# ###### LOCATION OF FILE CONTAINING CHROMATIC ZERO ANGLES ######### #
zero_angles = get_pkg_data_filename('data/theta_fors2.txt')
# ################################################################## #


def sort_red():
    """
    Creates back-up of the compressed files, uncompresses them, sorts them and re-names them.

    Notes
    -----
    For the .cl files to work properly, the naming convention used by sort_red is essential.

    """
    # Creating backups of the original uncompressed files
    os.system('mkdir backup')
    os.system('cp * backup')

    os.system('mkdir txt')
    os.system('mv *.txt txt')

    os.system('mkdir FITS')
    os.system('mv *.fits* FITS')

    os.chdir('FITS')
    os.system('uncompress *.Z')
    os.system('rm -f *.Z')

    # Instead of creating an actual file I could I have put the filenames in a list and sorted the list
    # It works and I won't change it cuz I am lazy.
    for filename in os.listdir("."):

        with open("files.txt","a") as f:
            f.write(filename+'\n')
            f.close()

    os.system('sort files.txt > ordered_files.txt')

    with open('ordered_files.txt') as f:
        content = f.readlines()

    filename_list = []

    for i in xrange(len(content)):
        filename_list.append(content[i][:-1])

    os.system('mkdir Data_reduction')
    os.system('mkdir Other')

    bias=0
    arc=0
    flat=0
    sky=0
    sci=0
    slit=0
    os.system('mkdir CHIP2')

    # We don't use CHIP 2 with FORS2 so I put all those in a separate folder
    for filename in filename_list:
        if "fits" in filename:
            try:
                hdulist = fits.open(filename)
                chip = hdulist[0].header['EXTNAME']
                if chip == 'CHIP2':
                    os.system('mv '+filename+' CHIP2')
                hdulist.close()
            except:
                print "Not CHIP2 "+filename

    for filename in filename_list:
        # Here we are renaming the files. The naming convention used here is assumed throughout the rest of the datred sub-module.
        if "fits" in filename:
            try:
                hdulist = fits.open(filename)
                head = hdulist[0].header['HIERARCH ESO DPR TYPE']
                if head == 'BIAS':
                    bias=bias+1
                    new_name='BIAS_'+str(bias).zfill(2)+".fits"
                    os.rename(filename, new_name)
                    os.system('mv '+new_name+' Data_reduction')
                if head == 'FLAT,LAMP':
                    flat=flat+1
                    new_name='FLAT_'+str(flat).zfill(2)+".fits"
                    os.rename(filename, new_name)
                    os.system('mv '+new_name+' Data_reduction')
                if head == 'WAVE,LAMP':
                    arc=arc+1
                    new_name='ARC_'+str(arc).zfill(2)+".fits"
                    os.rename(filename, new_name)
                    os.system('mv '+new_name+' Data_reduction')
                if head == 'OBJECT':
                    sci=sci+1
                    new_name='SCIENCE_'+str(sci).zfill(2)+".fits"
                    os.rename(filename, new_name)
                    os.system('mv '+new_name+' Data_reduction')
                if head == 'STD':
                    sci=sci+1
                    new_name='SCIENCE_'+str(sci).zfill(2)+".fits"
                    os.rename(filename, new_name)
                    os.system('mv '+new_name+' Data_reduction')
                if head == 'SKY':
                    sky=sky+1
                    new_name='SKY_'+str(sky).zfill(2)+".fits"
                    os.rename(filename, new_name)
                    os.system('mv '+new_name+' Other')
                if head == 'SLIT':
                    slit=slit+1
                    new_name='SLIT_'+str(slit).zfill(2)+".fits"
                    os.rename(filename, new_name)
                    os.system('mv '+new_name+' Other')

                hdulist.close()
            except:
                print "Could not sort this file (type?) "+filename

    os.chdir('Data_reduction')

    # If the images have 4096 pixels in x direction they're defo not our science images and if they are calibration
    # images they can't calibrate our science because it's the wrong size.
    os.system('mkdir wrong_size')
    for filename in os.listdir("."):
        if "fits" in filename:
            try:
                hdulist = fits.open(filename)
                size_x = hdulist[0].header['HIERARCH ESO DET OUT1 NX']

                if size_x == 4096:
                    os.system('mv '+filename+' wrong_size')
            except:
                print "Could not sort this file (size?) "+filename
            hdulist.close()

    print '\nAll Done! :D \n'


def info():
    """
    Creates table containing useful information on the images (taken from the headers).

    Notes
    -----
    Use in folder containing the uncompressed FITS files.

    Output File Format: Filename, ESO label, Retarder Plate Angle, Exposure time,Airmass, Grism, Bin, umber
    of Pixels, 1/Gain, Read Out Noise, Date.

    """
    try:
        os.remove('image_info.txt')
    except:
        print 'kitten'

    for filename in os.listdir("."):
        # Here we go through the headers and take out the information that we need to create the fits files.
        # I think the header names and variable names are explicit enough
        if "fits" in filename:
            hdulist = fits.open(filename)

            exptime = hdulist[0].header['EXPTIME']

            binx = hdulist[0].header['HIERARCH ESO DET WIN1 BINX']
            biny = hdulist[0].header['HIERARCH ESO DET WIN1 BINY']
            binning = str(binx)+"x"+str(biny)
            size_x = hdulist[0].header['HIERARCH ESO DET OUT1 NX']
            one_over_gain = hdulist[0].header['HIERARCH ESO DET OUT1 CONAD']
            ron = hdulist[0].header['HIERARCH ESO DET OUT1 RON']  # Read Out Noise
            try:
                grism = hdulist[0].header['HIERARCH ESO INS GRIS1 NAME']
            except:
                grism = 'None'
            '''
            try:
                OSF = hdulist[0].header['HIERARCH ESO INS FILT1 NAME']
            except:
                OSF = 'None'
            '''
            try:
                angle = hdulist[0].header['HIERARCH ESO INS RETA2 ROT']
            except:
                try:
                    angle = str(hdulist[0].header['HIERARCH ESO INS RETA4 ROT'])+'*'
                except:
                    angle = 'None'

            date = hdulist[0].header['DATE-OBS']

            try:
                name = hdulist[0].header['HIERARCH ESO OBS NAME']

            except:
                name = 'None'

            if "fits" and "SCIENCE" in filename:
                head = hdulist[0].header['HIERARCH ESO DPR TYPE']
                if head == 'OBJECT' or head == 'STD':

                    airm_i = float(hdulist[0].header['HIERARCH ESO TEL AIRM START'])
                    airm_f = float(hdulist[0].header['HIERARCH ESO TEL AIRM END'])
                    airm = (airm_i + airm_f)/2
                else:
                    airm = 'None'
            else:
                airm = 'None'

            with open('image_info.txt', 'a') as f:
                f.write(filename[:-5].ljust(15)+str(name).ljust(23)+str(angle).ljust(7)+str(exptime).ljust(10)+
                        str(airm).ljust(8)+str(grism).ljust(12)+binning.ljust(5)+str(size_x).ljust(5)+
                        str(one_over_gain).ljust(7)+str(ron).ljust(5)+str(date)+"\n")

    os.system("sort image_info.txt -o image_info.txt")  # To have filenames written out  alphabetically

    with open('image_info.txt', 'a') as f:  # Just putting labels on columns at end of file
        f.write('=========================================================================='
                '=============================================== \n')
        f.write('Filename'.ljust(15)+'ESO label'.ljust(23)+'Angle'.ljust(7)+'EXP.TIME'.ljust(10)+ 'Airmass'.ljust(8)+
                'Grism'.ljust(12)+'bin'.ljust(5)+'#Pix'.ljust(5)+'1/gain'.ljust(7)+'RON'.ljust(5)+'Date'+"\n")


def hwrpangles(sn_name = 'CCSN', zeropol_name = 'Zero_', polstd_name='NGC2024'):
    """
    Creates the file used by lin_specpol to know which images correspond to which HWRP angle.

    Notes
    -----
    Separate files are created for the CCSN, zero pol std, and polarised std. The output files are made of 4 columns
    containing the numbers of images corresponding to the 0, 22.5, 45 and 67.5 degree retarder plate angles.
    1 set of retarder plate angles per line.

    Parameters
    ----------
    sn_name : string, optional
        A string that is unique to the ESO name of the SN, e.g  Default: 'CCSN'
    zeropol_name : string, optional
        A string that is unique to the ESO name of the zero pol std, e.g Default:  'Zero_'
    polstd_name : string, optional
        Same for polarised std, e.g Default: 'NGC2024'. This one is the most likely to change from one observation set
        to an other. See the info file created using info() tu know what string to give hwrpangles() to be able to
        differentiate between the SN and the standards.

    """

    # We want to make sure that python reads the images in the right order, we also only need to look at the headers
    # of the original images with name format SCIENCE_##.fits where ## is a 2 dgit number.

    ls_filenames = []
    for filename in os.listdir("."):
        if 'SCIENCE' in filename and 'ms' not in filename and 'dSCIENCE' not in filename and 'cal' not in filename:
            ls_filenames.append(filename)
    sorted_files = sorted(ls_filenames)

    # We are sorting the image names in 3 list to distinguish the CCSN from the polarised std from the zero pol std.
    zeropol = []
    sn = []
    polstd= []

    for filename in sorted_files:
        # Here we go through the headers and take out the information that we need to create the fits files.
        # I think the header names and variable names are explicit enough
        if "fits" in filename:
            hdulist = fits.open(filename)
            try:
                name = hdulist[0].header['HIERARCH ESO OBS NAME']
            except:
                name = 'None'

            if sn_name in name:
                sn.append(filename)
            elif zeropol_name in name:
                zeropol.append(filename)
            elif polstd_name in name:
                polstd.append(filename)

    # Now we go through the list of CCSN, zero pol and polarised std images separately and look at the HWRP angle
    # of each frame. We then create the files recording which image number corresponds to which HWRP angle:
    # 1 file for the SN, 1 file for the pol std, 1 file for the zero pol std.

    for list_names in [sn, zeropol, polstd]:  # This loop goes through list of filenames for CCSN, zero pol and pol std
        if list_names == sn:
            output_name = 'hwrpangles'
        elif list_names == zeropol:
            output_name = 'hwrpa_zeropol'
        elif list_names == polstd:
            output_name = 'hwrpa_polstd'
        ls_0=[]
        ls_1=[]
        ls_2=[]
        ls_3=[]
        ls_v_0=[]
        ls_v_1=[]

        for name in list_names:  # This loop goes through each image name within the ccsn, pol or zero pol std list
            hdulist = fits.open(name)
            # So it doesn't crash if there's no info on the HWRP in headers.
            try:
                angle = hdulist[0].header['HIERARCH ESO INS RETA2 ROT']
            except:
                try:
                    angle = str(hdulist[0].header['HIERARCH ESO INS RETA4 ROT'])+'*'
                    print angle
                except:
                    angle = 'None'

            if angle == 0.0:
                ls_0.append(name[8:-5])
            if angle == 22.5:
                ls_1.append(name[8:-5])
            if angle == 45.0:
                ls_2.append(name[8:-5])
            if angle == 67.5:
                ls_3.append(name[8:-5])
            if angle == '45.0*':
                ls_v_0.append(name[8:-5])
            if angle == '315.0*':
                ls_v_1.append(name[8:-5])
        print ls_v_0, ls_v_1

        try:
            os.remove(output_name+'.txt')
        except:
            'kitten'
        try:
            os.remove(output_name+'_v.txt')
        except:
            'kitten'

        for i in xrange(len(ls_0)):
            with open(output_name+'.txt', 'a') as f:
                f.write(str(ls_0[i])+' '+str(ls_1[i])+' '+str(ls_2[i])+' '+str(ls_3[i])+'\n')
            f.close()
        for i in xrange(len(ls_v_0)):
            with open(output_name+'_v.txt', 'a') as fv:
                fv.write(str(ls_v_0[i])+' '+str(ls_v_1[i])+'\n')
            fv.close()

    return 


# ################# SPECPOL. IT'S BIG WITH NESTED FUNCTIONS - MAYBE BAD CODE, BUT IT WORKS ####################### #


def lin_specpol(oray='ap2', hwrpafile = 'hwrpangles.txt', e_min_wl = 3775, bayesian_pcorr=True, p0_step = 0.01):
    """
    Calculates the Stokes parameters and P.A of a data set and writes them out in a text file.

    Notes
    -----
    A plot showing p, q, u, P.A and Delta epsilon_q and Delta epsilon_u is produced. The plots are not automatically
    saved.

    Parameters
    ----------
    oray : string, optional
        Which aperture corresponds to the ordinary ray: 'ap1' or 'ap2'. Default is 'ap2'.
    hwrpafile : string, optional
        The file telling lin_specpol() which image corresponds to which HWRP angle. Created by hwrpangles().
        Default is 'hwrpangles.txt'
    e_min_wl : string, optional
        The first wavelength of the range within which Delta epsilons will be calculated. Default is 3775 (ang).
    bayesian_pcorr : bool, optional
        Turns on or off the bayesian p debiasing method (J. L. Quinn 2012). If False then the step function method will
        be used (wang et al 1997). Default is True.
    p0_step : float, optional
        Step size (and starting point) of the p0 distribution. If the step is larger that an observed value of p then
        the code will fail, and you should decrease the step size. Also increases the run time significantly.
        Default is 0.01
    """
    if oray=='ap2':
        eray='ap1'
    elif oray=='ap1':
        eray='ap2'

    #########################################
    #                LIN_SPECPOL            #
    #########################################
        
    def get_data(ls_0, ls_1, ls_2, ls_3):
        """
        This takes the flux data from the text files given by IRAF and sorts them in lists for later use.

        Notes
        -----
        For lin_specpol() use only.

        /!\ start wavelength and dispersion for each data file should be the same /!\

    `   Parameters
        ----------
        ls_0 : list of ints
            list of file number for files containing dat at 0 deg
        ls_1 : list of ints
            list of file number for files containing dat at 22.5 deg
        ls_2 : list of ints
            list of file number for files containing dat at 45 deg
        ls_3: list of ints
            list of file number for files containing dat at 67.5 deg

        Returns
        -------
        Lists for wl, o and r ray for each angle and errors for o and ray for each angle:
        wl, ls_fo0, ls_fe0, ls_fo0_err, ls_fe0_err, ls_fo1, ls_fe1, ls_fo1_err, ls_fe1_err, ls_fo2,
        ls_fe2, ls_fo2_err, ls_fe2_err, ls_fo3, ls_fe3, ls_fo3_err, ls_fe3_err
        """

        # Need to do this because python doesn't read files in alphabetical order but in order they
        # are written on the disc
        list_file=[]
        for name in os.listdir('.'):
            list_file.append(name)
        sorted_files = sorted(list_file)

        ls_fo0=[]
        ls_fe0=[]
        ls_fo0_err=[]
        ls_fe0_err=[]

        ls_fo1=[]
        ls_fe1=[]
        ls_fo1_err=[]
        ls_fe1_err=[]

        ls_fo2=[]
        ls_fe2=[]
        ls_fo2_err=[]
        ls_fe2_err=[]

        ls_fo3=[]
        ls_fe3=[]
        ls_fo3_err=[]
        ls_fe3_err=[]

        for filename in sorted_files:
            # This condition is related to the naming convention I have adopted.
            if 'c_' not in filename:
                # The following compares the number in the filename to the number in ls_0 to see if the image
                # correspond to a 0 deg HWRP angle set up. The naming convention is crucial for this line to work
                # as it keeps the number in the filename in the location: filename[-10:-8] or filename[-14:-12] for
                # flux and flux_error files, respectively.

                # Francesco: Use regex to look for \w*?_\d{2}*? check the regex codes
                if filename[-10:-8] in ls_0 or filename[-14:-12] in ls_0:
                    # Now we put the filename in the right list, oray or eray and flux or error on flux.
                    if oray in filename:
                        if 'err' not in filename:
                            wl, fo = np.loadtxt(filename, unpack=True, usecols=(0,1))
                            ls_fo0.append(fo)
                        else:
                            wl, fo = np.loadtxt(filename, unpack=True, usecols=(0,1))
                            ls_fo0_err.append(fo)

                    if eray in filename:
                        if 'err' not in filename:
                            wl, fe = np.loadtxt(filename, unpack=True, usecols=(0,1))
                            ls_fe0.append(fe)
                        else:
                            wl, fe = np.loadtxt(filename, unpack=True, usecols=(0,1))
                            ls_fe0_err.append(fe)

                # Same thing as the first loop but for 22.5 HWRP
                if filename[-10:-8] in ls_1 or filename[-14:-12] in ls_1:
                    if oray in filename:
                        if 'err' not in filename:
                            wl, fo = np.loadtxt(filename, unpack=True, usecols=(0,1))
                            ls_fo1.append(fo)
                        else:
                            wl, fo = np.loadtxt(filename, unpack=True, usecols=(0,1))
                            ls_fo1_err.append(fo)

                    if eray in filename:
                        if 'err' not in filename:
                            wl, fe = np.loadtxt(filename, unpack=True, usecols=(0,1))
                            ls_fe1.append(fe)
                        else:
                            wl, fe = np.loadtxt(filename, unpack=True, usecols=(0,1))
                            ls_fe1_err.append(fe)

                # Same thing as the first loop but for 45 HWRP
                if filename[-10:-8] in ls_2 or filename[-14:-12] in ls_2:
                    if oray in filename:
                        if 'err' not in filename:
                            wl, fo = np.loadtxt(filename, unpack=True, usecols=(0,1))
                            ls_fo2.append(fo)
                        else:
                            wl, fo = np.loadtxt(filename, unpack=True, usecols=(0,1))
                            ls_fo2_err.append(fo)

                    if eray in filename:
                        if 'err' not in filename:
                            wl, fe = np.loadtxt(filename, unpack=True, usecols=(0,1))
                            ls_fe2.append(fe)
                        else:
                            wl, fe = np.loadtxt(filename, unpack=True, usecols=(0,1))
                            ls_fe2_err.append(fe)

                # Same thing as the first loop but for 67.5 HWRP
                if filename[-10:-8] in ls_3 or filename[-14:-12] in ls_3:
                    if oray in filename:
                        if 'err' not in filename:
                            wl, fo = np.loadtxt(filename, unpack=True, usecols=(0,1))
                            ls_fo3.append(fo)
                        else:
                            wl, fo = np.loadtxt(filename, unpack=True, usecols=(0,1))
                            ls_fo3_err.append(fo)

                    if eray in filename:
                        if 'err' not in filename:
                            wl, fe = np.loadtxt(filename, unpack=True, usecols=(0,1))
                            ls_fe3.append(fe)
                        elif 'err' in filename:
                            wl, fe = np.loadtxt(filename, unpack=True, usecols=(0,1))
                            ls_fe3_err.append(fe)

        return wl, ls_fo0, ls_fe0, ls_fo0_err, ls_fe0_err, ls_fo1, ls_fe1, ls_fo1_err, ls_fe1_err, ls_fo2, ls_fe2, ls_fo2_err, ls_fe2_err, ls_fo3, ls_fe3, ls_fo3_err, ls_fe3_err

    def norm_flux(fo, fe, fo_r, fe_r):
        """
        For lin_specpol() use only. Finds normalised flux and error.

        Parameters
        ----------
        fo : array
            Array containing the ordinary spectrum
        fe : array
            Array containing the extra-ordinary spectrum
        fo_r : array
            Array containing the ordinary error spectrum
        fe_r : array
            Array containing the extra-ordinary error spectrum

        Returns
        -------
        arrays
            The F and F_r arrays (the normalised flux and error on normalised flux)

        """
        F = (fo - fe)/(fo + fe)
        F_r = np.array([])
        for i in range(len(fo)):
            F_ri = m.fabs(F[i])*np.sqrt(((fo_r[i]**2)+(fe_r[i]**2))*((1/(fo[i]-fe[i])**2)+(1/(fo[i]+fe[i])**2)))
            F_r=np.append(F_r,F_ri)
        return F, F_r

    def specpol(wl, fo0, fe0, fo0_r, fe0_r, fo1, fe1, fo1_r, fe1_r, fo2, fe2, fo2_r, fe2_r, fo3, fe3, fo3_r, fe3_r):
        """
        Finds the p, q, u, theta and errors on these quantities for a set of spectropolarimetric data.

        Notes
        -----
        For lin_specpol() use only

        Parameters
        ----------
        wl : array
            Wavelengths
        fo0 : array
            o ray at 0 deg
        fe0 : array
            e ray at 0 deg
        fo0_r : array
            error on o ray at 0 deg
        fe0_r : array
            error on e ray at 0 deg
        fo1 : array
        fe1 : array
        fo1_r : array
        fe1_r : array
        fo2 : array
        fe2 : array
        fo2_r : array
        fe2_r : array
        fo3 : array
        fe3 : array
        fo3_r : array
        fe3_r : array
        delta_e : array
            Epsilon at each wl bin.
        avg_e : float
            Average epsilon
        stdv_e : float
            Standard dev of epsilon

        Returns
        -------
        arrays
            p, q, and u in percent, with associated errors, as well as theta in degrees ( 0 < theta < 180) and its
            errors.
        """

        # Calculating the normalised fluxes for all 4 HWRP angles.
        F0,F0_r = norm_flux(fo0, fe0, fo0_r, fe0_r)
        F1,F1_r = norm_flux(fo1, fe1, fo1_r, fe1_r)
        F2,F2_r = norm_flux(fo2, fe2, fo2_r, fe2_r)
        F3,F3_r = norm_flux(fo3, fe3, fo3_r, fe3_r)

        # Now Stokes parameters and degree of pol.
        q = 0.5*(F0-F2)
        u = 0.5*(F1-F3)
        q_r = 0.5*np.sqrt(F0_r**2 + F2_r**2)
        u_r = 0.5*np.sqrt(F1_r**2 + F3_r**2)
        p = np.sqrt(q*q + u*u)
        p_r = (1/p) * np.sqrt( (q*q_r)**2 + (u*u_r)**2 )

        # Arrays where we're going to store the values of p and Stokes parameters and P.A
        # after we've applied corrections.
        pf = np.array([])
        qf = np.array([])
        uf = np.array([])
        theta = np.array([])

        # We take our chromatic zero-angles and interpolate them to match the wavlength bins of our data.
        wl2, thetaz = np.loadtxt(zero_angles, unpack = True, usecols =(0,1))
        theta0 = np.interp(wl, wl2, thetaz)

        # Now we apply corrections to the P.A
        for t in range(len(q)):
            theta_t = 0.5*m.atan2(u[t],q[t])
            theta_r = 0.5* np.sqrt( ( (u_r[t]/u[t])**2 + (q_r[t]/q[t])**2) * ( 1/(1+(u[t]/q[t])**2) )**2 )
            theta_t = (theta_t*180.0) /m.pi
            theta_r = (theta_r*180.0) /m.pi
            if theta_t < 0:
                theta_t = 180 + theta_t  # Making sure P.A is within limit 0<theta<180 deg

            theta_cor = theta_t - theta0[t]
            theta_cor_rad = (theta_cor/180.0)*m.pi
            theta = np.append(theta, theta_cor)
            q_t = p[t]*m.cos(2*theta_cor_rad) # Re-calculating Stokes parameters
            u_t = p[t]*m.sin(2*theta_cor_rad)
            qf = np.append(qf, q_t*100) # Filling our arrays of final Stokes parameters and p.
            uf = np.append(uf, u_t*100)
            pf = np.append(pf, np.sqrt(q_t**2 + u_t**2)*100)


        # Now calculating epsilon q and epsilon u and Delta epsilon.
        eq = 0.5*F0 + 0.5*F2
        eu = 0.5*F1 + 0.5*F3
        delta_e = eq-eu
        stdv_dequ=[]
        for i in xrange(len(wl)):
            if wl[i] > e_min_wl:
                dequ=eq[i]-eu[i]
                stdv_dequ.append(dequ)

        stdv_e = np.std(stdv_dequ)
        avg_e = np.average(stdv_dequ)

        return pf, p_r*100, qf, q_r*100, uf, u_r*100, theta, theta_r, delta_e, avg_e, stdv_e

    #########################################
    #            LIN SPECPOL MAIN           #
    #########################################

    # list of files corresponding to each angle (0, 22.5, 45, 67.5)
    ls_0, ls_1, ls_2, ls_3 = np.genfromtxt(hwrpafile, dtype='str', unpack = True, usecols = (0, 1, 2, 3))

    # Now getting the data from the files in lists that will be used by the specpol() function.
    wl, ls_fo0, ls_fe0, ls_fo0_err, ls_fe0_err, ls_fo1, ls_fe1, ls_fo1_err, ls_fe1_err, ls_fo2, ls_fe2, ls_fo2_err, ls_fe2_err, ls_fo3, ls_fe3, ls_fo3_err, ls_fe3_err = get_data(ls_0, ls_1, ls_2, ls_3)

    qls=[]
    qrls=[]
    uls=[]
    urls=[]
    delta_es=[]
    avg_es=[]
    stdv_es=[]

    # Each index in those lists refers to 1 set of values for all wavelength bins wl. Each set of Stokes parameters
    # is stored in lists defined above so that we can then take the average for each bin, and the standard deviation
    # which will be used as the error.

    for i in range(len(ls_fo0)):
        p, pr, q, qr, u, ur, theta, thetar, delta_e, avg_e, stdv_e = specpol(wl, ls_fo0[i], ls_fe0[i], ls_fo0_err[i], ls_fe0_err[i], ls_fo1[i], ls_fe1[i], ls_fo1_err[i], ls_fe1_err[i], ls_fo2[i], ls_fe2[i], ls_fo2_err[i], ls_fe2_err[i], ls_fo3[i], ls_fe3[i], ls_fo3_err[i], ls_fe3_err[i])
        qls.append(q)
        qrls.append(qr)
        uls.append(u)
        urls.append(ur)
        delta_es.append(delta_e)
        avg_es.append(avg_e)
        stdv_es.append(stdv_e)

    # Where we'll put the final values of the Stokes parameters and their errors.
    qf=np.array([])
    uf=np.array([])
    qrf=np.array([])
    urf=np.array([])

    for num in range(len (qls[0])):
        # num indexes the bins each list of Stokes parameters values
        q_to_avg=[]
        u_to_avg=[]
        qr_to_sum=np.array([])
        ur_to_sum=np.array([])
        for s in range(len(qls)):
            # s indexes the data set from which we are taking a particular Stoke parameter
            # We want to average values fo all data sets at each wavelength bins. For example say I have
            # 3 data sets, at 5000 A say, I am gonna take the 3 values of q in each data set at 5000 A and
            # average them. Do the same accross the whole spectrum and with each Stoke parameter to get final results.
            q_to_avg.append(qls[s][num])
            u_to_avg.append(uls[s][num])

            # For the next 4 lines I am calculating the error on the mean and putting it in final list of errors on
            # Stokes parameters
            qr_to_sum=np.append(qr_to_sum, 1/((qrls[s][num])**2))
            ur_to_sum=np.append(ur_to_sum, 1/((urls[s][num])**2))
        qrf= np.append(qrf, np.sqrt(1/np.sum(qr_to_sum)))
        urf= np.append(urf, np.sqrt(1/np.sum(ur_to_sum)))

        # My final Stokes parameters, doing a weighted average. The weights are 1/stdv**2
        qf= np.append(qf, np.average(q_to_avg, weights=1/(qr_to_sum**2)))
        uf= np.append(uf, np.average(u_to_avg, weights=1/(ur_to_sum**2)))

    # Once I have my final Stokes parameters I can calculate the final degree of polarisation (and error).
    pf = np.sqrt(qf**2 + uf**2)
    prf = (1/pf) * np.sqrt( (qf*qrf)**2 + (uf*urf)**2 )

    # And finally the P.A !
    thetaf=np.array([])
    thetarf=np.array([])
    for t in range(len(qrf)):
        thetaf_t = 0.5*m.atan2(uf[t],qf[t])
        thetarf_t = 0.5* np.sqrt( ( (urf[t]/uf[t])**2 + (qrf[t]/qf[t])**2) * ( 1/(1+(uf[t]/qf[t])**2) )**2 )
        thetaf_t = (thetaf_t*180.0) /m.pi
        thetarf_t = (thetarf_t*180.0) /m.pi
        if thetaf_t < 0:
            thetaf_t = 180 + thetaf_t # Again need to make sure the range is within 0-180 deg
        thetaf=np.append(thetaf, thetaf_t)
        thetarf=np.append(thetarf, thetarf_t)

    # #### P Bias Correction #### #

    #  If bayesian_pcorr is False, P will be debiased as in Wang et al. 1997 using a step function
    if bayesian_pcorr is False:
        print "Step Func - p correction"
        pfinal = np.array([])
        for ind in range(len(pf)):
            condition = pf[ind] - prf[ind]
            if condition > 0:
                p_0i = pf[ind]-((float(prf[ind]**2))/float(pf[ind]))
            elif condition < 0:
                p_0i = pf[ind]

            pfinal = np.append(pfinal, p_0i)

    #  If bayesian_pcorr is True, P will be debiased using the Bayesian method described by J. L. Quinn 2012
    #  the correceted p is pbar_{0,mean} * sigma. pbar_{0,mean} is given by equation 47 of J. L. Quinn 2012

    if bayesian_pcorr is True:
        print "Bayesian - p correction"
        sigma = (qrf + urf)/2
        pbar = pf/sigma
        pfinal = np.array([])
        for j in range(len(pbar)):
            p0 = np.arange(p0_step, pbar[j], p0_step)
            rho = np.array([])
            for i in range(len(p0)):
                tau = (sigma[j]**2)*2*p0[i]
                pp0 = pbar[j]*p0[i]
                RiceDistribution = pbar[j]*np.exp(-((pbar[j]**2 + p0[i]**2)/2)) * special.iv(0, pp0)
                rhoi = RiceDistribution * tau
                rho = np.append(rho, rhoi)

            p0mean = np.average(p0, weights=rho)
            pfinal = np.append(pfinal, p0mean*sigma[j])  # !!!! need to multiply by sigma to get p0 and not p0/bar.

    # ###### CREATING THE TEXT FILE ###### #
    pol_file = raw_input('What do you want to name the polarisation file? ')
    try:
        os.remove(pol_file)
    except:
         print 'kittens'
    for l in xrange(len(wl)):
        with open(pol_file, 'a') as pol_f:
            pol_f.write(str(wl[l])+'    '+str(pfinal[l])+'    '+str(prf[l])+'    '+str(qf[l])+'    '+str(qrf[l])+'    '+str(uf[l])+'    '+str(urf[l])+'    '+str(thetaf[l])+'    '+str(thetarf[l])+'\n')

    # ###### MAKING PLOTS ########
    # Just to check that everything looks right.

    f, axarr = plt.subplots(5, 1, figsize=(10, 20), sharex=True)
    plt.subplots_adjust(hspace=0)

    # First axis is p
    axarr[0].errorbar(wl, pf, yerr=prf, c='#D92F2F')
    axarr[0].axhline(0,0, ls='--', c='k')
    pmax=-1000
    for i in range(len(wl)):
        if wl[i]>4500 and pf[i]>pmax:
            pmax=pf[i]
    axarr[0].set_ylim([-0.1,pmax+0.4])
    axarr[0].set_ylabel('p(%)', fontsize=14)

    # Then q
    axarr[1].errorbar(wl, qf, yerr=qrf, c='#D92F2F')
    axarr[1].axhline(0,0, ls='--', c='k')
    qmax=-1000
    qmin=1000
    for i in range(len(wl)):
        if wl[i]>4500 and qf[i]>qmax:
            qmax=qf[i]
        if wl[i]>4500 and qf[i]<qmin:
            qmin=qf[i]
    axarr[1].set_ylim([qmin-0.3,qmax+0.3])
    axarr[1].set_ylabel('q(%)', fontsize=14)

    # And u
    axarr[2].errorbar(wl, uf, yerr=urf, c='#D92F2F')
    axarr[2].axhline(0,0, ls='--', c='k')
    umax=-1000
    umin=1000
    for i in range(len(wl)):
        if wl[i]>4500 and uf[i]>umax:
            umax=uf[i]
        if wl[i]>4500 and uf[i]<umin:
            umin=uf[i]
    axarr[2].set_ylim([umin-0.3,umax+0.3])
    axarr[2].set_ylabel('u(%)', fontsize=14)

    # P.A
    axarr[3].errorbar(wl, thetaf, yerr=thetarf, c='#D92F2F')
    axarr[3].axhline(0,0, ls='--', c='k')
    axarr[3].set_ylim([-0,180])
    axarr[3].set_ylabel('theta', fontsize=14)

    # And finally the Delta epsilons of each data set.
    for i in range(len(delta_es)):
        axarr[4].plot(wl, delta_es[i], alpha= 0.8)
        print "Average Delta epsilon =",avg_es[i],"STDV =",stdv_es[i]

    axarr[4].set_ylabel(r"$\Delta \epsilon", fontsize = 16)
    axarr[4].set_ylim([-0.4,0.4])
    plt.xlim([3500,10000])

    plt.show()

    return pf, prf, qf, qrf, uf, urf, thetaf, thetarf, delta_es


# #####################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########################### #


def circ_specpol(oray='ap2', hwrpafile = 'hwrpangles_v.txt', e_min_wl = 3775):
    """
    Calculates the circular polarisation v and epsilon_v for all data sets. The plot is not automatically saved.

    Parameters
    ----------
    oray : string
        Which aperture corresponds to the ordinary ray: 'ap1' or 'ap2'. Default is 'ap2'.
    hwrpafile : string
        The file telling circ_specpol() which image corresponds to which HWRP angle. Created by hwrpangles().
        Default is 'hwrpangles_v.txt',
    e_min_wl : int
        The first wavelength of the range within which Delta epsilons will be calculated. Default is 3775 (ang).

    """
    if oray=='ap2':
        eray='ap1'
    elif oray=='ap1':
        eray='ap2'

    #########################################
    #                CIRC_SPECPOL            #
    #########################################
        
    def get_data(ls_0, ls_1):
        """
        This takes the flux data from the text files given by IRAF and sorts them in lists for later use.

        Notes
        -----
        For lin_specpol() use only.
        /!\ start wavelength and dispersion for each data file should be the same /!\

        Parameters
        ----------
        ls_0 : list of ints
            list of file number for files containing dat at 45 deg
        ls_1 : list of ints
            list of file number for files containing dat at 315 deg

        Returns
        -------
        lists of lists
            lists for wl, o and e ray for each angle and errors for o and ray for each angle:
        wl, ls_fo0, ls_fe0, ls_fo0_err, ls_fe0_err, ls_fo1, ls_fe1, ls_fo1_err, ls_fe1_err, ls_fo2,
        ls_fe2, ls_fo2_err, ls_fe2_err, ls_fo3, ls_fe3, ls_fo3_err, ls_fe3_err
        """

        # Need to do this because python doesn't read files in alphabetical order but in order they
        # are written on the disc
        list_file=[]
        for name in os.listdir('.'):
            list_file.append(name)
        sorted_files = sorted(list_file)

        ls_fo0=[]
        ls_fe0=[]
        ls_fo0_err=[]
        ls_fe0_err=[]

        ls_fo1=[]
        ls_fe1=[]
        ls_fo1_err=[]
        ls_fe1_err=[]

        for filename in sorted_files:
            # This condition is related to the naming convention I have adopted.
            if 'c_' not in filename:
                # The following compares the number in the filename to the number in ls_0 to see if the image
                # correspond to a 0 deg HWRP angle set up. The naming convention is crucial for this line to work
                # as it keeps the number in the filename in the location: filename[-10:-8] or filename[-14:-12] for
                # flux and flux_error files, respectively.
                if filename[-10:-8] in ls_0 or filename[-14:-12] in ls_0:
                    # Now we put the filename in the right list, oray or eray and flux or error on flux.
                    if oray in filename:
                        if 'err' not in filename:
                            wl, fo = np.loadtxt(filename, unpack=True, usecols=(0,1))
                            ls_fo0.append(fo)
                        else:
                            wl, fo = np.loadtxt(filename, unpack=True, usecols=(0,1))
                            ls_fo0_err.append(fo)

                    if eray in filename:
                        if 'err' not in filename:
                            wl, fe = np.loadtxt(filename, unpack=True, usecols=(0,1))
                            ls_fe0.append(fe)
                        else:
                            wl, fe = np.loadtxt(filename, unpack=True, usecols=(0,1))
                            ls_fe0_err.append(fe)

                # Same thing as the first loop but for 22.5 HWRP
                if filename[-10:-8] in ls_1 or filename[-14:-12] in ls_1:
                    if oray in filename:
                        if 'err' not in filename:
                            wl, fo = np.loadtxt(filename, unpack=True, usecols=(0,1))
                            ls_fo1.append(fo)
                        else:
                            wl, fo = np.loadtxt(filename, unpack=True, usecols=(0,1))
                            ls_fo1_err.append(fo)

                    if eray in filename:
                        if 'err' not in filename:
                            wl, fe = np.loadtxt(filename, unpack=True, usecols=(0,1))
                            ls_fe1.append(fe)
                        else:
                            wl, fe = np.loadtxt(filename, unpack=True, usecols=(0,1))
                            ls_fe1_err.append(fe)

        return wl, ls_fo0, ls_fe0, ls_fo0_err, ls_fe0_err, ls_fo1, ls_fe1, ls_fo1_err, ls_fe1_err

    def normalized_V(wl, fo, fo_err, fe, fe_err):
        """
        Finds the normalized v

        Parameters
        ----------
        wl : array
            Array containing the wavelengths
        fo : array
            Array containing the values of the ordinary flux
        fo_err : array
            Array containing the values of the uncertainties on the ordinary flux
        fe : array
            Array containing the values of the extra-ordinary flux
        fe_err : array
            Array containing the values of the uncertainties on the extra-ordinary flux

        Returns
        -------
        arrays
            2 Arrays containing the values of the normalized V Stokes parameter and its errors at each wavelength
        """

        v = np.array([])
        v_err = np.array([])

        for i in xrange(len(wl)):
            F = (fo[i]-fe[i])/(fo[i]+fe[i])
            F_err = m.fabs(F)*np.sqrt(((fo_err[i]**2)+(fe_err[i]**2))*((1/(fo[i]-fe[i])**2)+(1/(fo[i]+fe[i])**2)))
            v = np.append(v, F)
            v_err = np.append(v_err, F_err)
            
        return v, v_err

    def v_1set(wl, ls_fo0, ls_fe0, ls_fo0_err, ls_fe0_err, ls_fo1, ls_fe1, ls_fo1_err, ls_fe1_err):
        """
        Calculate v and the error on v for a set of observations.

        Notes
        -----
        /!\ Requires the functions GET_DATA() and NORMALIZED_V() /!\

        Parameters
        ----------
        num2 : int
            Number in file name corresponding to observation at +45 deg
        num1 : int
            Number in file name corresponding to observation at -45 deg

        Returns
        -------
        arrays
            wl, v, v_err
        """
        v0, v0_err = normalized_V(wl, ls_fo0, ls_fo0_err, ls_fe0, ls_fe0_err)
        v1, v1_err = normalized_V(wl, ls_fo1, ls_fo1_err, ls_fe1, ls_fe1_err)

        v = []
        v_err = []
        eps=np.array([])

        for i in xrange(len(v0)):

            v_i = 0.5*(v1[i]-v0[i])  # v1 is 45 and v0 is -45 here
            eps_el = 0.5*(v1[i]+v0[i])
            eps=np.append(eps, eps_el*100)
            v_i_err = 0.5*np.sqrt(v1_err[i]**2 + v0_err[i]**2)
            v_i = v_i*100
            v_i_err = v_i_err*100
            v.append(v_i)
            v_err.append(v_i_err)
            #eps= eps*10
            
        cond = (wl>e_min_wl)
        eps_crop = eps[cond]
        eps_avg = np.average(eps_crop)
        eps_std = np.std(eps_crop)

        return v, v_err, eps, eps_avg, eps_std

    def wghtd_mean(values, err):
        """
        Weighted mean and error on the mean

        Parameters
        ----------
        values : array
            Array of the values
        err : array
            Array of the errors on the values, must have same dimension as 'values'.

        Returns
        -------
        floats
            mean, err_mean
        """
        num_to_sum = []
        den_to_sum = []
        for i in xrange(len(values)):
            numerator = values[i] / (err[i]*err[i])
            denominator = 1 / (err[i]*err[i])
            num_to_sum.append(numerator)
            den_to_sum.append(denominator)
        mean = np.sum(num_to_sum) / np.sum(den_to_sum)
        err_mean = np.sqrt( 1/(np.sum(den_to_sum)) )

        return mean, err_mean

    #########################################
    #                SPECPOL MAIN           #
    #########################################

    # list of files corresponding to each angle (0)
    ls_0, ls_1 = np.genfromtxt(hwrpafile, dtype='str', unpack = True, usecols = (0, 1))

    # Now getting the data from the files in lists that will be used by the specpol() function.
    wl, ls_fo0, ls_fe0, ls_fo0_err, ls_fe0_err, ls_fo1, ls_fe1, ls_fo1_err, ls_fe1_err = get_data(ls_0, ls_1)

    v_ls = []
    verr_ls = []
    ls_eps = []
    avg_eps = []
    stdv_eps = []
    
    for i in range(len(ls_fo0)):
        v, verr, eps, eps_avg, eps_std = v_1set(wl, ls_fo0[i], ls_fe0[i], ls_fo0_err[i], ls_fe0_err[i], ls_fo1[i], ls_fe1[i], ls_fo1_err[i], ls_fe1_err[i])
        v_ls.append(v)
        verr_ls.append(verr)
        ls_eps.append(eps)
        avg_eps.append(eps_avg)
        stdv_eps.append(eps_std)
    
    vf = np.array([])
    vf_err = np.array([])

    for num in range(len (v_ls[0])):
        # num indexes the bins each list of Stokes parameters values
        v_to_avg=[]
        verr_to_avg=[]
        for s in range(len(v_ls)):
            # s indexes the data set from which we are taking a particular Stoke parameter
            # We want to average values fo all data sets at each wavelength bins. For example say I have
            # 3 data sets, at 5000 A say, I am gonna take the 3 values of q in each data set at 5000 A and
            # average them. Do the same accross the whole spectrum and with each Stoke parameter to get final results.
            v_to_avg.append(v_ls[s][num])
            verr_to_avg.append(verr_ls[s][num])
        vi, vi_err = wghtd_mean(v_to_avg, verr_to_avg)
        vf=np.append(vf, vi)
        vf_err=np.append(vf_err,vi_err)

    # ###### CREATING THE TEXT FILE ###### #
    pol_file = raw_input('What do you want to name the polarisation file? ')

    try:
        os.remove(pol_file)
    except:
         print 'kittens'
    for l in xrange(len(wl)):
        with open(pol_file, 'a') as pol_f:
            pol_f.write(str(wl[l])+'    '+str(vf[l])+'    '+str(vf_err[l])+'\n')

    # ###### MAKING PLOT ########
    # Just to check that everything looks right.

    f, axarr = plt.subplots(2, 1, figsize=(10, 10), sharex=True)
    plt.subplots_adjust(hspace=0)

    # First axis is v
    axarr[0].errorbar(wl, vf, yerr=vf_err, c='#D92F2F')
    axarr[0].axhline(0,0, ls='--', c='k')
    vmax=-1000
    vmin=10000
    for i in range(len(wl)):
        if wl[i]>4500 and vf[i]>vmax:
            vmax=vf[i]
        if wl[i]>4500 and vf[i]<vmin:
            vmin=vf[i]

    axarr[0].set_ylim([vmin-0.4,vmax+0.4])
    axarr[0].set_ylabel('v(%)', fontsize=14)
    
    # And then the Delta epsilons of each data set.
    for i in range(len(ls_eps)):
        axarr[1].plot(wl, ls_eps[i], alpha= 0.8)
        print "Average Delta epsilon =",avg_eps[i],"STDV =",stdv_eps[i]

    axarr[1].set_ylabel(r"$\Delta \epsilon", fontsize = 16)
    axarr[1].set_ylim([-6,0])
    plt.xlim([3500,10000])

    plt.show()

    return wl, vf, vf_err, ls_eps

# ########################## V BAND POL ###############################


def lin_vband(oray = 'ap2', hwrpafile = 'hwrpangles.txt'):
    """
    This creates synthetic V band linear polarimetry data from the spectropolarimetric data.

    Parameters
    ----------
    hwrpafile : string, optional
        The file telling lin_specpol() which image corresponds to which HWRP angle. Created by hwrpangles().
        Default is 'hwrpangles.txt'.
    oray : string, optional
        Which aperture is the oridnary ray. Shoudl be either 'ap1' or 'ap2'. Default is 'ap2'

    """
    S.setref(area=10053097)  # area in cm^2 - has to be set cuz different to HST

    if oray=='ap2':
        eray='ap1'
    elif oray=='ap1':
        eray='ap2'

    def v_counts(filename, bp_v):
        # Only used by lin_vband(). Finds v counts from spectrum in filename. bp_v is tthe bandpass of V band
        sp = S.FileSpectrum(filename)
        obs = S.Observation(sp, bp_v)

        return obs.wave, obs.flux

    def norm_flux(fo, fe,fo_r, fe_r):
        # normalised flux. Repeat of the one in lin_specpol()
        F = (fo - fe)/(fo + fe)
        F_r = m.fabs(F)*np.sqrt(((fo_r**2)+(fe_r**2))*((1/(fo-fe)**2)+(1/(fo+fe)**2)))

        return F, F_r

    def vpol(fo0, fo1, fo2, fo3, fe0, fe1, fe2, fe3, fo0_r, fe0_r,fo1_r,fe1_r,fo2_r,fe2_r,fo3_r,fe3_r, wl):
        """
        Similar to specpol() in lin_specpol(), but for V-band polarisation.
        """

        # Normalised fluxes
        F0,F0_r = norm_flux(fo0, fe0,fo0_r, fe0_r)
        F1,F1_r = norm_flux(fo1, fe1,fo1_r, fe1_r)
        F2,F2_r = norm_flux(fo2, fe2,fo2_r, fe2_r)
        F3,F3_r = norm_flux(fo3, fe3,fo3_r, fe3_r)

        # Stokes parameters
        q = 0.5*(F0-F2)
        u = 0.5*(F1-F3)
        q_r = 0.5*np.sqrt(F0_r**2 + F2_r**2)
        u_r = 0.5*np.sqrt(F1_r**2 + F3_r**2)

        # degree of polarisation
        p = np.sqrt(q*q + u*u)
        p_r = (1/p) * np.sqrt( (q*q_r)**2 + (u*u_r)**2 )

        # Interpolation of chromatic zero angle values
        wl2, thetaz = np.loadtxt(zero_angles, unpack = True, usecols =(0,1))
        theta0 = np.interp(wl, wl2, thetaz)

        # Finding Polarisation angle
        theta = 0.5*m.atan2(u,q)
        theta_r = 0.5* np.sqrt( ( (u_r/u)**2 + (q_r/q)**2) * ( 1/(1+(u/q)**2) )**2 )
        theta = (theta*180.0) /m.pi
        theta_r = (theta_r*180.0) /m.pi
        if theta < 0:
            theta = 180 + theta
        theta_cor = theta - theta0  # Correction of chromatic zero angle
        theta_cor_rad = (theta_cor/180.0)*m.pi

        # Correction of Stokes parameters and p from P.A correction
        q = p*m.cos(2*theta_cor_rad)
        u = p*m.sin(2*theta_cor_rad)
        p = np.sqrt(q*q + u*u)

        return p*100, p_r*100, q*100, q_r*100, u*100, u_r*100, theta, theta_r

    # ############################### #
    bp_v=S.ObsBandpass("johnson,v")
    wl_v= bp_v.avgwave()

    # list of files corresponding to each angle (0, 22.5, 45, 67.5)
    list_0, list_1, list_2, list_3 = np.genfromtxt(hwrpafile, dtype='str', unpack = True, usecols = (0, 1, 2, 3))
    
    o0 = np.array([])  # 0 deg
    e0 = np.array([])
    o1 = np.array([])  # 22.5 deg
    e1 = np.array([])
    o2 = np.array([])  # 45 deg
    e2 = np.array([])
    o3 = np.array([])  # 67.5deg
    e3 = np.array([])
    o0_r = np.array([])  # 0 deg
    e0_r = np.array([])
    o1_r = np.array([])  # 22.5 deg
    e1_r = np.array([])
    o2_r = np.array([])  # 45 deg
    e2_r = np.array([])
    o3_r = np.array([])  # 67.5deg
    e3_r = np.array([])

    list_file=[]
    for name in os.listdir('.'):
        list_file.append(name)

    sorted_files = sorted(list_file)

    for filename in sorted_files:

        if 'dSC'in filename and 'c_' not in filename and '.fits' not in filename:
            vflux=v_counts(filename, bp_v) # finds counts across Vband given spectrum
            counts = np.sum(vflux[1]) # sum over all bins to get total number of counts in Vband
            counts_r = np.sqrt(np.sum(vflux[1]**2))  # Poisson noise

            # very similar to what's done in get_data() in lin_specpol() so refer to that
            if filename[-10:-8] in list_0 or filename[-14:-12] in list_0:
                if oray in filename:      
                    if 'err' not in filename:
                        o0=np.append(o0, counts)
                    else:
                        o0_r=np.append(o0_r, counts_r)
                
                if eray in filename:
                   # print filename
                    if 'err' not in filename:
                        e0=np.append(e0, counts)
                    else:
                        e0_r=np.append(e0_r, counts_r)

            if filename[-10:-8] in list_1 or filename[-14:-12] in list_1:
                if oray in filename:
                   # print filename
                    if 'err' not in filename:
                        o1=np.append(o1, counts)
                    else:
                        o1_r=np.append(o1_r, counts_r)

                if eray in filename:
                   # print filename
                    if 'err' not in filename:
                        e1=np.append(e1, counts)
                    else:
                        e1_r=np.append(e1_r, counts_r)

            if filename[-10:-8] in list_2 or filename[-14:-12] in list_2:
                if oray in filename:
                   # print filename
                    if 'err' not in filename:
                        o2=np.append(o2, counts)
                    else:
                        o2_r=np.append(o2_r, counts_r)

                if eray in filename:
                   # print filename
                    if 'err' not in filename:
                        e2=np.append(e2, counts)
                    else:
                        e2_r=np.append(e2_r, counts_r)

            if filename[-10:-8] in list_3 or filename[-14:-12] in list_3:
                if oray in filename:
                   # print filename
                    if 'err' not in filename:
                        o3=np.append(o3, counts)
                    else:
                        o3_r=np.append(o3_r, counts_r)

                if eray in filename:
                   # print filename
                    if 'err' not in filename:
                        e3=np.append(e3, counts)
                    else:
                        e3_r=np.append(e3_r, counts_r)

    p_ls=[]
    pr_ls=[]
    q_ls=[]
    qr_ls=np.array([])
    u_ls=[]
    ur_ls=np.array([])
    if len(o0) > 1:
        for i in range(len(list_0)):
            qr_to_sum=np.array([])
            ur_to_sum=np.array([])
            p, p_r, q, q_r, u, u_r, theta, theta_r = vpol(o0[i], o1[i],o2[i],o3[i],e0[i], e1[i],e2[i],e3[i],o0_r[i],e0_r[i],o1_r[i],e1_r[i],o2_r[i],e2_r[i],o3_r[i],e3_r[i], wl_v)
            p_ls.append(p)
            q_ls.append(q)
            u_ls.append(u)
            pr_ls=np.append(pr_ls,p)
            qr_ls=np.append(qr_ls,q)
            ur_ls=np.append(ur_ls,u)

            qr_to_sum=np.append(qr_to_sum, 1/((q_r)**2))
            ur_to_sum=np.append(ur_to_sum, 1/((u_r)**2))

        qavg=np.average(q_ls, weights=1/(qr_ls**2))
        uavg=np.average(u_ls, weights=1/(ur_ls**2))
        qavg_r=np.sqrt(1/np.sum(qr_to_sum))
        uavg_r=np.sqrt(1/np.sum(ur_to_sum))
        pavg= np.sqrt(qavg**2 + uavg**2)
        pavg_r = (1/pavg) * np.sqrt( (qavg*qavg_r)**2 + (uavg*uavg_r)**2 )
        theta_v = (0.5*m.atan2(uavg,qavg))*180/m.pi
        if theta_v <0:
            theta_v = 180+theta_v
        theta_vr = (0.5* np.sqrt( ( (uavg_r/uavg)**2 + (qavg_r/qavg)**2) * ( 1/(1+(uavg/qavg)**2) )**2 ))*180/m.pi
        
    elif len(o0) == 1:
        pavg, pavg_r, qavg, qavg_r, uavg, uavg_r, theta_v, theta_vr = vpol(o0, o1,o2,o3,e0, e1,e2,e3,o0_r,e0_r,o1_r,e1_r,o2_r,e2_r,o3_r,e3_r, wl_v)
        if theta_v <0:
            theta_v = 180+theta_v
        
    print pavg, pavg_r, qavg, qavg_r, uavg, uavg_r, theta_v, theta_vr
    return pavg, pavg_r, qavg, qavg_r, uavg, uavg_r, theta_v, theta_vr


#############################################################################

def flux_spectrum():
    """
    Combines all the flux calibrated apertures to create the flux spectrum.

    Notes
    -----
    Creates a text file with 3 columns columns: wavelength flux errors
    """
    flux = []
    flux_err = []
    i = 0
    output = raw_input('What do you want to call the output file? ')

    for filename in os.listdir("."):
        # Putting flux in from each file in list.
        if "1D_c_" in filename and "err" not in filename:
            wl, f = np.loadtxt(filename, unpack=True, usecols=(0,1))
            if i == 0:
                flux.append(wl)
                flux_err.append(wl)
            flux.append(f)
            i = i+1
        # Putting ERROR on flux from each error file in error list.
        elif "1D_c_" in filename and "err" in filename:
            wl, f = np.loadtxt(filename, unpack=True, usecols=(0,1))
            if i == 0:
                flux.append(wl)
                flux_err.append(wl)

            flux_err.append(f)
            i = i+1

    for x in xrange(len(flux[:][0])):
        if x == 0:
            try:
                os.remove(output)
            except:
                print "kittens"
        sum_flux = 0
        error_sqrd = 0
        for i in xrange(len(flux)):
            if i == 0:
                wl=flux[i][x]
            else:
                sum_flux = sum_flux + flux[i][x]  # Summing flux in given wavelength bin
                error_sqrd = error_sqrd + flux_err[i][x]*flux_err[i][x]  # Adding the square of the errors...
        error = np.sqrt(error_sqrd)  # ... and taking the square-root => error of the bin

        # Writing out the wl, flux and error in each bin out in a text file.
        with open(output,"a") as f:
            f.write(str(wl)+' '+str(sum_flux)+' '+str(error)+'\n')

    return
