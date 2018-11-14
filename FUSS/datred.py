"""
2 - Jan - 2018 / H. F. Stevance / fstevance1@sheffield.ac.uk

datred.py is a module created as part of the FUSS package to help with the data reduction of spectropolarimetric
data (at the present time only used with FORS2 data)

Pre-requisites
--------------
os, astropy.io, numpy, math, matplotlib.pyplot, pysynphot, scipy.special, pandas

Variable
--------
zero_angles : string
    Path to the text file containing the chromatic zero angles for FORS2 (needs updating for you own system).
    Can be found at: http://www.eso.org/sci/facilities/paranal/instruments/fors/inst/pola.html

Functions
---------
sort_red
rebin
pol_deg
pol_ang
F_from_oeray

Classes
-------
Meta
SpecPol
LinearSpecPol
CircularSpecPol

"""

from __future__ import division
from __future__ import print_function
import os
from astropy.io import fits
import numpy as np
import math as m
import matplotlib.pyplot as plt
import pysynphot as S
from scipy import special as special
from astropy.utils.data import get_pkg_data_filename
import re
import pandas as pd
import sys


if sys.version_info.major < 3:
    range = xrange
    input = raw_input

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

    filename_list = sorted([filename for filename in os.listdir(".")])

    os.system('mkdir Data_reduction')
    os.system('mkdir Other')

    bias = 0
    arc = 0
    flat = 0
    sky = 0
    sci = 0
    slit = 0
    os.system('mkdir CHIP2')

    # We don't use CHIP 2 with FORS2 so I put all those in a separate folder
    for filename in filename_list:
        if "fits" in filename:
            try:
                hdulist = fits.open(filename)
                chip = hdulist[0].header['EXTNAME']
                if chip == 'CHIP2':
                    os.system('mv ' + filename + ' CHIP2')
                hdulist.close()
            except:
                print("Not CHIP2 " + filename)

    for filename in filename_list:
        # Here we are renaming the files. The naming convention used here is assumed
        # throughout the rest of the datred sub-module.
        if "fits" in filename:
            try:
                hdulist = fits.open(filename)
                head = hdulist[0].header['HIERARCH ESO DPR TYPE']
                if head == 'BIAS':
                    bias = bias + 1
                    new_name = 'BIAS_' + str(bias).zfill(3) + ".fits"
                    os.rename(filename, new_name)
                    os.system('mv ' + new_name + ' Data_reduction')
                if head == 'FLAT,LAMP':
                    flat = flat + 1
                    new_name = 'FLAT_' + str(flat).zfill(3) + ".fits"
                    os.rename(filename, new_name)
                    os.system('mv ' + new_name + ' Data_reduction')
                if head == 'WAVE,LAMP':
                    arc = arc + 1
                    new_name = 'ARC_' + str(arc).zfill(3) + ".fits"
                    os.rename(filename, new_name)
                    os.system('mv ' + new_name + ' Data_reduction')
                if head == 'OBJECT':
                    sci = sci + 1
                    new_name = 'SCIENCE_' + str(sci).zfill(3) + ".fits"
                    os.rename(filename, new_name)
                    os.system('mv ' + new_name + ' Data_reduction')
                if head == 'STD':
                    sci = sci + 1
                    new_name = 'SCIENCE_' + str(sci).zfill(3) + ".fits"
                    os.rename(filename, new_name)
                    os.system('mv ' + new_name + ' Data_reduction')
                if head == 'SKY':
                    sky = sky + 1
                    new_name = 'SKY_' + str(sky).zfill(3) + ".fits"
                    os.rename(filename, new_name)
                    os.system('mv ' + new_name + ' Other')
                if head == 'SLIT':
                    slit = slit + 1
                    new_name = 'SLIT_' + str(slit).zfill(3) + ".fits"
                    os.rename(filename, new_name)
                    os.system('mv ' + new_name + ' Other')

                hdulist.close()
            except:
                print("Could not sort this file (type?) " + filename)

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
                    os.system('mv ' + filename + ' wrong_size')
            except:
                print("Could not sort this file (size?) " + filename)
            hdulist.close()

    print('\nAll Done! :D \n')


def rebin(wl, f, r, bin_siz=15):
    """
    To rebin my flux spectra

    Parameters
    ----------
    wl : array
        1D array containing the wavelengths to be rebinned
    f : array
        1D array containing the fluxes to be rebinned
    r : array
        1D array containing the errors on the fluxes to be rebinned
    bin_siz : int
        Size of the new bins in Angstrom.

    Returns
    -------
        tuple of 1D arrays: wl, f, err all rebinned to the new bin size
    if bin_siz is None:
        print "No binning"
        return wl, f, r
    """
    wl = np.array(wl)
    f = np.array(f)
    r = np.array(r)
    small_bin_sizes = []

    bins_f = np.zeros(int((max(wl) - min(wl)) / bin_siz) + 1)  # new flux bins, empty for now
    bins_w = np.zeros(int((max(wl) - min(wl)) / bin_siz) + 1)  # new error bins, empty for now

    weights = 1 / (r ** 2)

    for i in range(len(wl) - 1):
        # n = int((wl[i]-min(wl))/bin_siz) # n is the number of the new bin
        small_bin_sizes.append((wl[i + 1] - wl[i]))  # filling list of small bin sizes

    bin_centers = [(min(wl) + bin_siz / 2) + bin_siz * n for n in range(len(bins_f))]  # finding the new bin centers

    bin_edges = [bin_centers[0] - bin_siz / 2] + [bin1 + bin_siz / 2 for bin1 in
                                                  bin_centers]  # finding the new bin edges

    ind_edge = []  # in this list I'll put the index of the array wl corresponding to the wavelength values
    # that are close to the bin edges.

    for edge in bin_edges:
        i_wl_at_edge = min(range(len(wl[:-1])), key=lambda i: abs(edge - wl[i]))
        # this is to find the small bin that is closest to the edge of the new bin
        # print wl[i_wl_at_edge], small_bin_sizes[i_wl_at_edge]
        ind_edge.append(i_wl_at_edge)

    for i in range(len(wl)):
        n = int((wl[i] - min(wl)) / bin_siz)
        if i in ind_edge:
            j = ind_edge.index(i)  # finding index j of the wavelength index i I am interested in
            edge = bin_edges[j]  # the edge to compare to wl[i] will then be at bin_edges[j]

            if wl[i] < edge:
                frac_overlap = (wl[i] + small_bin_sizes[i] / 2 - edge) / (small_bin_sizes[i])
                try:
                    bins_f[n] += f[i] * weights[i] * (1 - frac_overlap)
                    bins_w[n] += weights[i] * (1 - frac_overlap)
                    bins_f[n + 1] += f[i] * weights[i] * frac_overlap
                    bins_w[n + 1] += weights[i] * frac_overlap

                except IndexError:
                    print("Index Error at ", wl[i])
                    pass

            elif wl[i] > edge:
                frac_overlap = (wl[i] + small_bin_sizes[i] / 2 - edge) / (small_bin_sizes[i])
                try:
                    bins_f[n] += f[i] * weights[i] * frac_overlap
                    bins_w[n] += weights[i] * frac_overlap
                    bins_f[n + 1] += f[i] * weights[i] * (1 - frac_overlap)
                    bins_w[n + 1] += weights[i] * (1 - frac_overlap)

                except IndexError:
                    print("Index Error at ", wl[i])
                    pass

        else:
            try:
                bins_f[n] += f[i] * weights[i]
                bins_w[n] += weights[i]
            except IndexError:
                print("Index Error at ", wl[i])
                pass

    for i in range(len(bin_centers)):
        if bins_w[i] == 0.0:
            print(bin_centers[i], bins_w[i], "\nCAREFUL! BIN WLGTH == 0!")

    bins_f[:-1] /= bins_w[:-1]  # normalise weighted values by sum of weights to get weighted average

    bins_err = np.sqrt(1 / bins_w[:-1])

    return bin_centers[:-1], bins_f[:-1], bins_err


def pol_deg(q, u, q_r=None, u_r=None):
    """
    Finds the degree of polarisation p from Stokes parameters q and u. Does debiasing using step function.

    Notes
    -----
    q, u, q_r and u_r must have the same dimension

    Parameters
    ----------
    q : float or 1D numpy.array
        Stokes q.
    u : float or 1D numpy.array
        Stokes u.
    q_r : float or 1D numpy.array
        Errors on Stokes q.
    u_r : float or 1D numpy.array
        Errors on Stokes u.

    Returns
    -------
    tuple (p, error on p) if errors on q and u are given.
    Only p if errors are not given.

    """
    p = np.sqrt(q * q + u * u)
    if q_r is not None and u_r is not None:
        p_r = (1 / p) * np.sqrt((q * q_r) ** 2 + (u * u_r) ** 2)

        # Here we are debiasing using a step function. See Wang et al. (1997)
        # (Polarimetry of the Type IA Supernova SN 1996X -- their equation 3)
        try:
            for i in range(len(p)):
                if p[i] - p_r[i] > 0:
                    p[i] -= (p_r[i]**2)/p[i]
        except TypeError:
            if p - p_r > 0:
                    p -= (p_r**2)/p   
            
        return p, p_r
    else:
        return p


def pol_ang(q, u, q_r=None, u_r=None):
    """
    Calculates the polarisation angle

    Parameters
    ----------
    q : float or 1D numpy.array
        Stokes q.
    u : float or 1D numpy.array
        Stokes u.
    q_r : float or 1D numpy.array
        Errors on Stokes q.
    u_r : float or 1D numpy.array
        Errors on Stokes u.

    Returns
    -------
    tuple (theta, error on theta) if errors on q and u are given.
    Only theta if errors are not given.
4
    """

    if isinstance(q, float):
        theta = 0.5 * m.atan2(u, q)
        theta = (theta * 180.0) / m.pi
        if theta < 0:
            theta += 180  # Making sure P.A is within limit 0<theta<180 deg

        if q_r is not None and u_r is not None:
            theta_r = 0.5 * np.sqrt(((u_r / u) ** 2 + (q_r / q) ** 2) * (1 / (1 + (u / q) ** 2)) ** 2)
            theta_r = (theta_r * 180.0) / m.pi
            if theta_r > 180:
                theta_r = 180
            return theta, theta_r
        else:
            return theta

    else:
        theta = np.array([])
        theta_r = np.array([])
        for t in range(len(q)):
            theta_t = 0.5 * m.atan2(u[t], q[t])
            theta_t = (theta_t * 180.0) / m.pi
            if theta_t < 0:
                theta_t += 180  # Making sure P.A is within limit 0<theta<180 deg
            theta = np.append(theta, theta_t)
            if q_r is not None and u_r is not None:
                theta_tr = 0.5 * np.sqrt(
                        ((u_r[t] / u[t]) ** 2 + (q_r[t] / q[t]) ** 2) * (1 / (1 + (u[t] / q[t]) ** 2)) ** 2)
                theta_tr = (theta_tr * 180.0) / m.pi
                if theta_tr > 180:
                    theta_tr = 180
                theta_r = np.append(theta_r, theta_tr)
                

        if q_r is not None and u_r is not None:
            return theta, theta_r
        else:
            return theta


def F_from_oeray(fo, fe, fo_r, fe_r):
    """
    Normalised flux (F) from ordinary and extra-ordinary rays

    Notes
    -----
    All the arrays should be the same length.

    Parameters
    ----------
    fo : 1D numpy.array
        flux of the ordinary ray
    fe : 1D numpy.array
        flux of the extra-ordinary ray
    fo_r : 1D numpy.array
        error on the flux of the ordinary ray
    fe_r : 1D numpy.array
        error on the flux of the extra-ordinary ray

    Returns
    -------
    Normalised flux : np.array
    """

    F = (fo - fe) / (fo + fe)
    F_r = abs(F) * np.sqrt( ((fo_r ** 2) + (fe_r ** 2)) * ( (1 / (fo - fe) ** 2) + (1 / (fo + fe) ** 2)) )

    return F, F_r



class Meta(object):
    """
    To find the meta data of the images we are working on, create a data frame and save it to a tab separated file.

    Examples
    --------
    >> import FUSS.datred as r
    >> path = "/home/heloise/Data/11hs/88.D-0761/2011-12-23/FITS/Data_reduction/bad_binning/"
    >> metadataframe = r.Meta(path=path).dataThe Meta Data file metadata already exists. Would you like to replace it? (Y/n)N
    File not replaced
    >> metadataframe.iloc[-5:-1]
          Flag         Filename     ESO label Angle Pol. Type Exp. Time Airmass  \
    24  polstd  SCIENCE_84.fits  PMOS_NGC2024     0       lin   19.9985   1.095
    25  polstd  SCIENCE_85.fits  PMOS_NGC2024  22.5       lin   19.9922  1.0945
    26  polstd  SCIENCE_86.fits  PMOS_NGC2024    45       lin   19.9986   1.094
    27  polstd  SCIENCE_87.fits  PMOS_NGC2024  67.5       lin   20.0026   1.093

       1/gain  RON      Grism  Bin                     Date
    24    0.7  2.9  GRIS_300V  2x2  2011-12-23T03:44:59.590
    25    0.7  2.9  GRIS_300V  2x2  2011-12-23T03:46:04.915
    26    0.7  2.9  GRIS_300V  2x2  2011-12-23T03:47:09.870
    27    0.7  2.9  GRIS_300V  2x2  2011-12-23T03:48:15.905
    >> metadataframe.info()
    <class 'pandas.core.frame.DataFrame'>
    Int64Index: 29 entries, 0 to 28
    Data columns (total 12 columns):
    Flag         8 non-null object
    Filename     29 non-null object
    ESO label    29 non-null object
    Angle        29 non-null object
    Pol. Type    23 non-null object
    Exp. Time    29 non-null object
    Airmass      29 non-null object
    1/gain       29 non-null object
    RON          29 non-null object
    Grism        29 non-null object
    Bin          29 non-null object
    Date         29 non-null object
    dtypes: object(12)
    memory usage: 2.9+ KB


    Notes
    -----
    The files whose headers we will look at must have ".fits" in the name, cannot be "dSCIENCE" or "cal_" files.
    (This is to do with the notation used throughout the FUSS package both in the IRAF and python scripts)

    Attributes
    ----------
    clobber_flag : bool or None
        Whether or not to clobber the output_file if it exists. Initiated as "None".
    output : str
        Name of the output file
    path : str
        Path to the location of the images the metadata should be read from. Default: "./".
    target_flag : str
        String present in ESO labels that uniquely identifies the target.
    zeropol_flag : str
        String present in ESO labels that uniquely identifies the zero polarisation standard.
    polstd_flag : str
        String present in ESO labels that uniquely identifies the polarised standard.
    sorted_files : list of str
        Sorted list containing the names of the files whose headers we want to look at.
    data : pandas.core.frame.DataFrame
        Data frame containing the meta data.
        Columns: ['Flag', 'Filename', 'ESO label', 'Angle', 'Pol. Type', 'Exp. Time', 'Airmass', '1/gain', 'RON', 'Grism', 'Bin', 'Date']

    Methods
    -------
    _write_out()
        Writes to file
    _clobber_flag()
        Sets the clobber file to True or False if it has not been set on initiation.
    _flag()
        Writes the corresponding flags to each rows. Either 'tar', 'zpol' or 'polstd'
    _read_headers()
        Reads the fits headers and fills the info attribute data frame.


    """
    def __init__(self, path = "./", output_file = 'metadata',
                 target_flag='CCSN', zpol_flag='Zero_', polstd_flag='NGC2024',
                 clobber = None, make_file = True):

        """
        Initiates MetaData object

        Parameters
        ----------
        path : str, optional
            path to the location of the images the metadata should be read from. Default: "./".
        output_file : str, optional
            path to the output file where the data frame containing metadata will be written. Default: 'metadata'.
        target_flag : str, optional
            String present in ESO labels that uniquely identifies the target.
        zpol_flag : str, optional
            String present in ESO labels that uniquely identifies the zero polarisation standard.
        polstd_flag : str, optional
            String present in ESO labels that uniquely identifies the polarised standard.
        clobber : bool, optional
            If a file "output_file" already exists, whether to replace it or not. Default: False.

        Returns
        -------
        None
        """

        self.clobber_flag = clobber
        self.output = output_file
        self.path = path
        self.target_flag = target_flag
        self.zeropol_flag = zpol_flag
        self.polstd_flag = polstd_flag

        # We want to make sure that python reads the images in the right order, we also only need to look at the headers
        # of the original images with name format SCIENCE_##.fits where ## is a 2 or 3 digit number.
        self.sorted_files = sorted([filename for filename in os.listdir(self.path)
                               if 'fits' in filename and 'ms' not in filename
                               and 'dSCIENCE' not in filename and 'cal' not in filename])

        self.data = pd.DataFrame(columns = ['Flag', 'Filename', 'ESO label', 'Angle', 'Pol. Type',
                                            'Exp. Time', 'Airmass', '1/gain',
                                            'RON', 'Grism', 'Bin', 'Date'])

        self._read_headers()  # read headers and fills the dataframce
        self._flag()  # to flag targets and pol standards
        if make_file is False:
            return

        self._write_out()  # writes it out to a file
        return

    def _write_out(self):
        """
        Writes to file
        """
        if not os.path.isfile(self.output):
            self.data.to_csv(self.output, sep='\t', index=False)
            print("File Created")
            return
        else:
            self._clobber_flag()
            if self.clobber_flag is True:
                os.remove(self.output)
                self.data.to_csv(self.output, sep='\t', index=False)
                print("File replaced")
            if self.clobber_flag is False:
                print("File not replaced")

    def _clobber_flag(self):
        """
        Sets the clobber file to True or False if it has not been set on initiation.
        """
        while self.clobber_flag is None:
            clobber = input("The Meta Data file "+self.output+" already exists. Would you like to replace it? (Y/n)")
            if clobber in ['n', 'N', 'no', 'no']:
                self.clobber_flag = False
            elif clobber in ['', ' ', 'Y', 'y', 'yes', 'Yes', 'ye', 'Ye']:
                self.clobber_flag = True

    def _flag(self):
        """
        Writes the corresponding flags to each rows. Either 'tar', 'zpol' or 'polstd'
        """
        for i in range(len(self.data['ESO label'])):
            if self.target_flag in self.data.loc[i,'ESO label']:
                self.data.loc[i, 'Flag'] = 'tar'
            elif self.zeropol_flag in self.data.loc[i,'ESO label']:
                self.data.loc[i, 'Flag'] = 'zpol'
            elif self.polstd_flag in self.data.loc[i,'ESO label']:
                self.data.loc[i, 'Flag'] = 'polstd'

    def _read_headers(self):
        """
        Reads the fits headers and fills the info attribute data frame.
        """
        for i in range(len(self.sorted_files)):
            filename = self.sorted_files[i]
            self.data.loc[i, 'Filename'] = filename
            # Here we go through the headers and take out the information that we need to create the fits files.
            # I think the header names and variable names are explicit enough
            hdulist = fits.open(self.path+filename)

            exptime = hdulist[0].header['EXPTIME']
            self.data.loc[i,'Exp. Time'] = exptime

            binx = hdulist[0].header['HIERARCH ESO DET WIN1 BINX']
            biny = hdulist[0].header['HIERARCH ESO DET WIN1 BINY']
            binning = str(binx) + "x" + str(biny)
            self.data.loc[i,'Bin'] = binning

            #size_x = hdulist[0].header['HIERARCH ESO DET OUT1 NX']

            one_over_gain = hdulist[0].header['HIERARCH ESO DET OUT1 CONAD']
            self.data.loc[i,'1/gain'] = one_over_gain

            ron = hdulist[0].header['HIERARCH ESO DET OUT1 RON']  # Read Out Noise
            self.data.loc[i,'RON'] = ron

            try:
                grism = hdulist[0].header['HIERARCH ESO INS GRIS1 NAME']
            except KeyError:
                grism = 'None'
            self.data.loc[i,'Grism'] = grism
            poltype = None
            try:
                angle = hdulist[0].header['HIERARCH ESO INS RETA2 ROT']
                poltype = 'lin'
            except KeyError:
                try:
                    angle = hdulist[0].header['HIERARCH ESO INS RETA2 POSANG']
                    poltype = 'lin'
                except KeyError:
                    try:
                        angle = str(hdulist[0].header['HIERARCH ESO INS RETA4 ROT'])
                        poltype = 'circ'
                    except KeyError:
                        angle = 'None'

            self.data.loc[i,'Angle'] = angle
            self.data.loc[i,'Pol. Type'] = poltype

            date = hdulist[0].header['DATE-OBS']
            self.data.loc[i,'Date'] = date

            try:
                esoname = hdulist[0].header['HIERARCH ESO OBS NAME']
            except KeyError:
                esoname = 'None'
            self.data.loc[i,'ESO label'] = esoname

            if "fits" and "SCIENCE" in filename:
                head = hdulist[0].header['HIERARCH ESO DPR TYPE']
                if head == 'OBJECT' or head == 'STD':

                    airm_i = float(hdulist[0].header['HIERARCH ESO TEL AIRM START'])
                    airm_f = float(hdulist[0].header['HIERARCH ESO TEL AIRM END'])
                    airm = (airm_i + airm_f) / 2
                else:
                    airm = 'None'
            else:
                airm = 'None'

            self.data.loc[i,'Airmass'] = airm


# ################# LINEAR SPECPOL ####################### #
class SpecPol(object):
    """
    Base class for LinearSpecPol and CircularSpecPol.

    Notes
    -----
    SpecPol is not particularly useful on its own, it is simply a basic framework of methods and attributes we
    need for both linear and circular polarimetry calculations.

    Attributes
    ---------
    metadata : pandas.core.frame.DataFrame
        Data frame containing the metadata as built by MetaData. Can be initiated from tab separated file containing
        the data frame or from the data frame directly.
    oray : str
        Which aperture flag that corresponds to the ordinary ray files. In our naming convention the flag can be either
        'ap1' or 'ap2'.
    eray : str
        Which aperture flag that corresponds to the extra-ordinary ray files. In our naming convention the flag can be
        either 'ap1' or 'ap2'.
    bin_size : int or None
        Size of the bins in Angstrom if rebinning. Otherwise None.
    snrplot : bool
        Whether or not to plot the signal to noise ratio plots (expectd and calculated -- to check how well a job the
        binning did). Default is False.
    pol_file : Initiated as None
        Name of the output polarisation file. Will be defined through user input.
    flag : Initiated as None
        Flag to focus on in the metadata: 'tar', 'zpol' or 'polstd'. Defined in the calculate() method of LinearSpecPol
        and CircularSpecPol
    """

    def __init__(self, oray='ap2', metadata = 'metadata', bin_size=None,
                 snrplot=False):

        if oray == 'ap2':
            self.oray, self.eray = 'ap2', 'ap1'
        elif oray == 'ap1':
            self.oray, self.eray = 'ap1', 'ap2'

        if isinstance(metadata, str):
            assert os.path.isfile(metadata), metadata+" is not a valid file path."
            self.metadata = pd.read_csv(metadata, sep='\t')
        elif isinstance(metadata, pd.core.frame.DataFrame):
            self.metadata = metadata

        self.bin_size=bin_size
        self.snrplot = snrplot
        self.pol_file = None
        self.flag = None

    def _flux_diff_from_file(self, files, check_bin=False):
        """
        Calculates the normalised flux differences from files.

        Parameters
        ----------
        files : list of str
            list of files
        check_bin : bool
            Whether to check the binning with snr plots.

        """

        # Extracting polarised fluxes of o and e ray and their errors according to filenames
        for filename in files:
            if self.oray in filename:
                if 'err' not in filename:
                    self.wl, fo = np.loadtxt(filename, unpack=True, usecols=(0, 1))
                else:
                    self.wl, fo_r = np.loadtxt(filename, unpack=True, usecols=(0, 1))

            if self.eray in filename:
                if 'err' not in filename:
                    self.wl, fe = np.loadtxt(filename, unpack=True, usecols=(0, 1))
                else:
                    self.wl, fe_r = np.loadtxt(filename, unpack=True, usecols=(0, 1))

        # BINNING
        if self.bin_size is None:
            self.wl_bin, fo_bin, fo_bin_r, fe_bin, fe_bin_r = self.wl, fo, fo_r, fe, fe_r

        else:
            print("Binning to ", self.bin_size, "Angstrom")
            self.wl_bin, fo_bin, fo_bin_r = rebin(self.wl, fo, fo_r, bin_siz=self.bin_size)
            self.wl_bin, fe_bin, fe_bin_r = rebin(self.wl, fe, fe_r, bin_siz=self.bin_size)

            # To perform a few checks on the binning (in particular the SNR yielded by binning)
            if check_bin is True:
                self._check_binning(fo, fe, fo_r, fe_r, fo_bin, fe_bin, fo_bin_r, fe_bin_r)

        # Finding flux difference
        F, F_r = F_from_oeray(fo_bin, fe_bin, fo_bin_r, fe_bin_r)

        return self.wl_bin, F, F_r

    def _list_files(self, angle, linpol = True, fileloc='.'):
        """
        Creates list of files of interest.

        Parameters
        ----------
        angle : float
            Which HWRP angle we are interested in.
        linpol : bool, optional
            Whether we are doing linear (True) or circular (False) polarisation. Default is True.
        fileloc : str
            Path to location of files. Usually will be '.' so this is the default.
        """
        if linpol is True:
            poltype = 'lin'
        else:
            poltype = 'circ'

        angle = float(angle)
        root_list = [str(self.metadata.loc[i,"Filename"])[:-5]\
                  for i in xrange(len(self.metadata["Filename"])) \
                  if self.metadata.loc[i, 'Flag'] == self.flag \
                  if float(self.metadata.loc[i, 'Angle']) == angle \
                    if self.metadata.loc[i, 'Pol. Type'] == poltype]

        mylist = []
        sorted_files = sorted([filename for filename in os.listdir(fileloc) if '.txt' in filename and '1D' in filename])
        for filename in sorted_files:
            for root in root_list:
                if root in filename:
                    mylist.append(filename)
        return mylist

    def _check_binning(self, fo, fe, fo_r, fe_r, bin_fo, bin_fe, bin_fo_r, bin_fe_r):
        """
        This performs a few checks on the binning: Calculates expected values of SNR and compares to the values of
        SNR after binning: does the median and the values at central wavelength. Also produces a plot of whole spectrum
        """
        snr_not_binned = np.array((fo + fe) / np.sqrt(fo_r ** 2 + fe_r ** 2))
        snr_expected = snr_not_binned * np.sqrt(self.bin_size / (self.wl[1] - self.wl[0]))
        ind_not_binned_central_wl = int(np.argwhere(self.wl == min(self.wl, key=lambda x: abs(x - 6204)))[0])
        snr_not_binned_central_wl = snr_not_binned[ind_not_binned_central_wl]
        snr_central_expected = snr_not_binned_central_wl * np.sqrt(self.bin_size / (self.wl[1] - self.wl[0]))

        snr_binned = np.array((bin_fo + bin_fe) / np.sqrt(bin_fo_r ** 2 + bin_fe_r ** 2))
        ind_central_wl = int(np.argwhere(self.wl_bin == min(self.wl_bin, key=lambda x: abs(x - 6204)))[0])
        snr_central_wl = (bin_fo[ind_central_wl] + bin_fe[ind_central_wl]) / \
                         np.sqrt(bin_fo_r[ind_central_wl] ** 2 + bin_fe_r[ind_central_wl] ** 2)

        print("\n======== BEFORE BINNING ======")
        print("MEDIAN SNR ")
        print(np.median(snr_not_binned))

        print("CENTRAL SNR at (", self.wl[ind_not_binned_central_wl], " A)")
        print(snr_not_binned_central_wl)

        print("======== AFTER BINNING ======")

        print("MEDIAN SNR / EXPECTED ")
        print(np.median(snr_binned), np.median(snr_expected))

        print("CENTRAL SNR / EXPECTED (at ", self.wl[ind_not_binned_central_wl], " A)")
        print(snr_central_wl, snr_central_expected)
        print("\n")

        if self.snrplot is True:
            plt.plot(self.wl, snr_expected, marker='o', label='Expected')
            plt.plot(self.wl_bin, snr_binned, marker='x', label='Calculated after binning')
            plt.legend()
            plt.title("SNR")
            plt.show()

        return


class LinearSpecPol(SpecPol):
    """
    class for Linear sepctropolarimetry reduction. Inherits from SpecPol

    Examples
    --------
    >> import FUSS.datred as r
    >> import pandas as pd
    >> metadataframe = pd.read_csv('metadata', sep='\t') # or can just feed a freshly created Meta.data data frame
    >> pol = r.LinearSpecPol(metadata = metadataframe, bin_size = 15)
    >> poldataframe = pol.calculate() # One full set of specpol data (4 images and 4 files to rebin per image)
    Binning to  15 Angstrom
    Index Error at  9316.2 # This can happen on the last bin and isn't a problem. If it is in the middle of the spcetrum
    Index Error at  9326.1 # then you should look into it.
    Index Error at  9316.2
    Index Error at  9326.1

    ======== BEFORE BINNING ======
    MEDIAN SNR
    86.6505323611
    CENTRAL SNR at ( 6204.3  A)
    117.809391701
    ======== AFTER BINNING ======
    MEDIAN SNR / EXPECTED
    185.294052885 184.73955572
    CENTRAL SNR / EXPECTED (at  6204.3  A)
    250.397252588 251.17046704


    Binning to  15 Angstrom
    Index Error at  9316.2
    Index Error at  9326.1
    Index Error at  9316.2
    Index Error at  9326.1
    Binning to  15 Angstrom
    Index Error at  9316.2
    Index Error at  9326.1
    Index Error at  9316.2
    Index Error at  9326.1
    Binning to  15 Angstrom
    Index Error at  9316.2
    Index Error at  9326.1
    Index Error at  9316.2
    Index Error at  9326.1

    What do you want to name the polarisation file? [filename] # the file created will be filename.pol
    >> poldataframe.info()
    <class 'pandas.core.frame.DataFrame'>
    RangeIndex: 421 entries, 0 to 420
    Data columns (total 10 columns):
    wl            421 non-null float64
    p             421 non-null float64
    p_r           421 non-null float64
    q             421 non-null float64
    q_r           421 non-null float64
    u             421 non-null float64
    u_r           421 non-null float64
    theta         421 non-null float64
    theta_r       421 non-null float64
    delta_eps0    421 non-null float64
    dtypes: float64(10)
    memory usage: 33.0 KB

    Attributes
    ----------
    wl : 1D np.array
        Initiated as None. Will contains original wavelength bins
    wl_bin : 1D np.array
        Initiated as None. Will contain the binned wavelength. If bin_size is None then wl_bin will == wl.
    poldata : pandas.core.frame.DataFrame
       Initiated with columns: 'wl', 'p', 'p_r', 'q', 'q_r', 'u', 'u_r', 'theta', 'theta_r' and dtype='float64'
       delta_epsilon# (# being a number) will be added for delta epsilon spectra corresponding to each set of 4
       Half Wave Retarder Plate Angles.
    flag : str
        Flag corresponding to object of interest  in the metadata data frame column "Flag".
        Initiated as None and defined in the method calculate() to specify which object to do the specpol on
        (the target, the zero polarisation standard or the polarisation standard -- usually 'tar', 'zpol', 'polstd',
        respectively)

    Methods
    -------
    calculate()
        Calculates the Stokes parameters, degree and angle polarisation from ordinary and extra-ordinary ray fluxes.
        Creates an output file containing the poldata data frame (filled by this method)
    plot()
        Plots the Stokes parameters, degree and angle of polarisation calculated.
    _get_data()
    _linspecpol()
    + inherited from SpecPol
    """
    def __init__(self, oray='ap2', metadata = 'metadata',
                 bin_size=None, snrplot=False):
        SpecPol.__init__(self, oray, metadata,bin_size, snrplot)
        self.wl= None
        self.wl_bin = None
        self.poldata = pd.DataFrame(columns = ['wl', 'p', 'p_r', 'q', 'q_r',
                                               'u', 'u_r', 'theta', 'theta_r'], dtype='float64')
        self.flag = None

    def calculate(self, flag='tar'):
        """
        Calculates the Stokes parameters, degree and angle polarisation, their errors and delta_epsilon for each set
        from the ordinary and extra-ordinary ray fluxes.

        Notes
        -----
        Creates an output file containing the poldata data frame (filled by this method)

        Parameters
        ----------
        flag : str
            Which flag to look for in the "Flag" column of metadata. Default is 'tar' which identifies target data.
            'zpol' can be used for zero pol. standard and 'polstd' for the polarisation standard.
            Custom flags can be input if you have written your own in the "Flag" column of the metadata Data Frame

        Returns
        -------
        poldata : pandas.core.frame.DataFrame
            Contains the polarisation data frame.
            columns = ['wl', 'p', 'p_r', 'q', 'q_r', 'u', 'u_r', 'theta', 'theta_r'], dtype='float64'
        """
        self.flag = flag
        # Now getting the data from the files in lists that will be used by the specpol() function.
        ls_F0, ls_F0_r, ls_F1, ls_F1_r, ls_F2, ls_F2_r, ls_F3, ls_F3_r = self._get_data()

        self.poldata["wl"] = self.wl_bin # even if have not rebinned since in that case we do self.wl_bin = self.wl

        q_ls = []
        q_r_ls = []
        u_ls = []
        u_r_ls = []

        for i in range(len(ls_F0)):
            p, pr, q, q_r, u, u_r, theta, theta_r, delta_e= self._linspecpol(ls_F0[i], ls_F0_r[i], ls_F1[i],
                                                                                       ls_F1_r[i], ls_F2[i], ls_F2_r[i],
                                                                                       ls_F3[i], ls_F3_r[i])
            q_ls.append(q)
            q_r_ls.append(q_r)
            u_ls.append(u)
            u_r_ls.append(u_r)

            self.poldata["delta_eps"+str(i)] = delta_e

        for num in range(len(q_ls[0])):
            # num indexes the bins each list of Stokes parameters values
            q_to_avg = []
            u_to_avg = []
            q_r_to_sum = np.array([])
            u_r_to_sum = np.array([])
            for s in range(len(q_ls)):
                # s indexes the data set from which we are taking a particular Stoke parameter
                # We want to average values fo all data sets at each wavelength bins. For example say I have
                # 3 data sets, at 5000 A say, I am gonna take the 3 values of q in each data set at 5000 A and
                # average them. Do the same accross the whole spectrum and with each Stoke parameter to get final results.
                q_to_avg.append(q_ls[s][num])
                u_to_avg.append(u_ls[s][num])

                # For the next 4 lines I am calculating the error on the mean and putting it in final list of errors on
                # Stokes parameters
                q_r_to_sum = np.append(q_r_to_sum, 1 / ((q_r_ls[s][num]) ** 2))
                u_r_to_sum = np.append(u_r_to_sum, 1 / ((u_r_ls[s][num]) ** 2))


            self.poldata.loc[num, 'q'] = np.average(q_to_avg, weights=q_r_to_sum)
            self.poldata.loc[num, 'u'] = np.average(u_to_avg, weights=u_r_to_sum)
            self.poldata.loc[num, 'q_r'] = np.sqrt(1 / np.sum(q_r_to_sum))
            self.poldata.loc[num, 'u_r'] = np.sqrt(1 / np.sum(u_r_to_sum))

        # Once I have my final Stokes parameters I can calculate the final debiases degree of polarisation (and error).
        self.poldata['p'], self.poldata['p_r'] = pol_deg(self.poldata['q'].values, self.poldata['u'].values,
                                                         self.poldata['q_r'].values, self.poldata['u_r'].values)

        # And finally the P.A !
        self.poldata['theta'], self.poldata['theta_r'] = pol_ang(self.poldata['q'].values, self.poldata['u'].values,
                                                        self.poldata['q_r'].values, self.poldata['u_r'].values)

        # ###### CREATING THE TEXT FILE ###### #
        self.pol_file = input('What do you want to name the polarisation file? ')

        try:
            os.remove(self.pol_file + ".pol")
        except OSError:
            pass

        self.poldata.to_csv(self.pol_file+".pol", index=False, sep="\t")

        return self.poldata

    def plot(self):
        if self.poldata is None:
            return "You should run the `calculate`method before trying to make plots "

        wavelength = self.poldata['wl'].values
        p = self.poldata['p'].values
        p_r = self.poldata['p_r'].values
        q = self.poldata['q'].values
        q_r = self.poldata['q_r'].values
        u = self.poldata['u'].values
        u_r = self.poldata['u_r'].values
        theta = self.poldata['theta'].values
        theta_r = self.poldata['theta_r'].values

        f, axarr = plt.subplots(5, 1, figsize=(8, 8), sharex=True)
        plt.subplots_adjust(hspace=0)

        # First axis is p
        axarr[0].errorbar(wavelength, p,yerr=p_r, c='#D92F2F')
        axarr[0].axhline(0, 0, ls='--', c='k')
        pmax = -1000
        for i in range(len(wavelength)):
            if wavelength[i] > 4500 and p[i] > pmax:
                pmax = p[i]

        axarr[0].set_ylim([-0.1, pmax + 0.4])
        axarr[0].set_ylabel('p(%)', fontsize=14)

        # Then q
        axarr[1].errorbar(wavelength, q, yerr=q_r, c='#D92F2F')
        axarr[1].axhline(0, 0, ls='--', c='k')
        qmax = -1000
        qmin = 1000
        for i in range(len(wavelength)):
            if wavelength[i] > 4500 and q[i] > qmax:
                qmax = q[i]
            if wavelength[i] > 4500 and q[i] < qmin:
                qmin = q[i]
        axarr[1].set_ylim([qmin - 0.3, qmax + 0.3])
        axarr[1].set_ylabel('q(%)', fontsize=14)

        # And u
        axarr[2].errorbar(wavelength, u, yerr=u_r, c='#D92F2F')
        axarr[2].axhline(0, 0, ls='--', c='k')
        umax = -1000
        umin = 1000
        for i in range(len(wavelength)):
            if wavelength[i] > 4500 and u[i] > umax:
                umax = u[i]
            if wavelength[i] > 4500 and u[i] < umin:
                umin = u[i]
        axarr[2].set_ylim([umin - 0.3, umax + 0.3])
        axarr[2].set_ylabel('u(%)', fontsize=14)

        # P.A
        axarr[3].errorbar(wavelength, theta, yerr=theta_r, c='#D92F2F')
        axarr[3].axhline(0, 0, ls='--', c='k')
        axarr[3].set_ylim([-0, 180])
        axarr[3].set_ylabel('theta', fontsize=14)

        # And finally the Delta epsilons of each data set.
        delta_cols = [col for col in self.poldata.columns if 'delta' in col]
        for column in delta_cols:
            axarr[4].plot(wavelength, self.poldata[column], alpha=0.8)
            print("MEAN Delta epsilon =", self.poldata[column].values.mean(),
                  "STDV =", self.poldata[column].values.std() )

        axarr[4].set_ylabel(r"$\Delta \epsilon$", fontsize=16)
        axarr[4].set_ylim([-4.0, 4.0])
        plt.xlim([3500, 10000])

        save_cond = input("do you want to save the plot?(Y/n): ")
        if save_cond == "y" or save_cond == "Y" or save_cond == "":
            plt.savefig(self.pol_file + ".png")
            print("Plot saved")
        else:
            print("Plot not saved")

        plt.show()

    def _get_data(self):
        """
        This takes the flux data from the text files given by IRAF and sorts them in lists for later use.

        Returns
        -------
        tuple of 8 lists
            Each list corresponds to the normalised flux difference for a given HWRP angle or its error.
            Each lsit contains N 1D numpy arrays containing N noral;ised flux differences spectra (or error), one
            for each set of spectropolarimetric data.

         """

        # Need to do this because python doesn't read files in alphabetical order but in order they
        # are written on the disc
        #sorted_files = sorted([filename for filename in os.listdir(".")
        #                       if 'dSCIENCE' in filename and 'fits' not in filename and 'c_' not in filename])
        files_0_deg = SpecPol._list_files(self, 0.0)
        files_22_deg = SpecPol._list_files(self, 22.5)
        files_45_deg = SpecPol._list_files(self, 45.0)
        files_67_deg = SpecPol._list_files(self, 67.5)
        errormessage = "It seems you don't have the same number of images for each retarder plate angle. " \
                       "This may not be code breaking but only complete sets of retarder plate angles can be used."

        assert len(files_0_deg) == len(files_22_deg) == len(files_45_deg) == len(files_67_deg), errormessage

        ls_F0 = []
        ls_F0_r = []

        ls_F1 = []
        ls_F1_r = []

        ls_F2 = []
        ls_F2_r = []

        ls_F3 = []
        ls_F3_r = []

        for file_list in [files_0_deg, files_22_deg, files_45_deg, files_67_deg]:
            nbre_sets = len(file_list) / 4
            nbre_sets_remainder = len(file_list) % 4
            assert nbre_sets_remainder == 0, "There should be 4 data files for each image (2 apertures x (flux + err))"
            print("4 Files per images... All good here")

        for i in range(int(nbre_sets)):
            step = i * 4
            files_0_deg_subset = files_0_deg[0 + step:4 + step]
            files_22_deg_subset = files_22_deg[0 + step:4 + step]
            files_45_deg_subset = files_45_deg[0 + step:4 + step]
            files_67_deg_subset = files_67_deg[0 + step:4 + step]
            if i == 0:
                check = True
            wl0, F0, F0_r = SpecPol._flux_diff_from_file(self, files_0_deg_subset, check_bin=check)
            ls_F0.append(F0)
            ls_F0_r.append(F0_r)
            wl1, F1, F1_r = SpecPol._flux_diff_from_file(self, files_22_deg_subset)
            ls_F1.append(F1)
            ls_F1_r.append(F1_r)
            wl2, F2, F2_r = SpecPol._flux_diff_from_file(self, files_45_deg_subset)
            ls_F2.append(F2)
            ls_F2_r.append(F2_r)
            wl3, F3, F3_r = SpecPol._flux_diff_from_file(self, files_67_deg_subset)
            ls_F3.append(F3)
            ls_F3_r.append(F3_r)

        assert len(wl0) == len(wl1) == len(wl2) == len(wl3), "Wavelength bins not homogeneous. This will be an issue."

        return ls_F0, ls_F0_r, ls_F1, ls_F1_r, ls_F2, ls_F2_r, ls_F3, ls_F3_r

    def _linspecpol(self, F0, F0_r, F1, F1_r, F2, F2_r, F3, F3_r):
        """
        Calculates the spectropolarimetric data from the normalised flux differences.

        Parameters
        ----------
        Eight 1D numpy arrays
            The normalised flux differences (and errors) for each HWRP angle.

        Returns
        -------
        arrays
            p, q, and u in percent, with associated errors, as well as theta in degrees ( 0 < theta < 180) and its
            errors, and delta epsilon (see Maund 2008)
        """

        # Now Stokes parameters and degree of pol.
        q = 0.5 * (F0 - F2)
        u = 0.5 * (F1 - F3)
        q_r = 0.5 * np.sqrt(F0_r ** 2 + F2_r ** 2)
        u_r = 0.5 * np.sqrt(F1_r ** 2 + F3_r ** 2)
        p, p_r = pol_deg(q, u, q_r, u_r)

        # Arrays where we're going to store the values of p and Stokes parameters and P.A
        # after we've applied corrections.
        pf = np.array([])
        qf = np.array([])
        uf = np.array([])
        theta = np.array([])

        # We take our chromatic zero-angles and interpolate them to match the wavlength bins of our data.
        wl2, thetaz = np.loadtxt(zero_angles, unpack=True, usecols=(0, 1))
        theta0 = np.interp(self.wl_bin, wl2, thetaz)
        # Now we apply corrections to the P.A
        for t in range(len(q)):
            theta_t = 0.5 * m.atan2(u[t], q[t])
            theta_r = 0.5 * np.sqrt(((u_r[t] / u[t]) ** 2 + (q_r[t] / q[t]) ** 2) * (1 / (1 + (u[t] / q[t]) ** 2)) ** 2)
            theta_t = (theta_t * 180.0) / m.pi
            theta_r = (theta_r * 180.0) / m.pi
            if theta_t < 0:
                theta_t += 180  # Making sure P.A is within limit 0<theta<180 deg

            theta_cor = theta_t - theta0[t]
            theta_cor_rad = (theta_cor / 180.0) * m.pi
            theta = np.append(theta, theta_cor)
            q_t = p[t] * m.cos(2 * theta_cor_rad)  # Re-calculating Stokes parameters
            u_t = p[t] * m.sin(2 * theta_cor_rad)
            qf = np.append(qf, q_t * 100)  # Filling our arrays of final Stokes parameters and p.
            uf = np.append(uf, u_t * 100)
            pf = np.append(pf, np.sqrt(q_t ** 2 + u_t ** 2) * 100)

        # Now calculating epsilon q and epsilon u and Delta epsilon.
        eq = (0.5 * F0 + 0.5 * F2) * 100  # in percent
        eu = (0.5 * F1 + 0.5 * F3) * 100
        delta_e = eq - eu

        return pf, p_r * 100, qf, q_r * 100, uf, u_r * 100, theta, theta_r, delta_e


class CircularSpecPol(SpecPol):
    """
    class for Circular sepctropolarimetry reduction. Inherits from SpecPol

    Examples
    --------
    Similar to LinearSpecpol examples.

    Attributes
    ----------
    wl : 1D np.array
        Initiated as None. Will contains original wavelength bins
    wl_bin : 1D np.array
        Initiated as None. Will contain the binned wavelength. If bin_size is None then wl_bin will == wl.
    poldata : pandas.core.frame.DataFrame
       Initiated with columns: 'wl', 'p', 'p_r', 'q', 'q_r', 'u', 'u_r', 'theta', 'theta_r' and dtype='float64'
       delta_epsilon# (# being a number) will be added for delta epsilon spectra corresponding to each set of 4
       Half Wave Retarder Plate Angles.
    flag : str
        Flag corresponding to object of interest  in the metadata data frame column "Flag".
        Initiated as None and defined in the method calculate() to specify which object to do the specpol on
        (the target, the zero polarisation standard or the polarisation standard -- usually 'tar', 'zpol', 'polstd',
        respectively)

    Methods
    -------
    calculate()
        Calculates the Stokes parameters, degree and angle polarisation from ordinary and extra-ordinary ray fluxes.
        Creates an output file containing the poldata data frame (filled by this method)
    plot()
        Plots the Stokes parameters, degree and angle of polarisation calculated.
    _get_data()
    _circspecpol()
    + inherited from SpecPol
    """
    def __init__(self, oray='ap2', metadata = 'metadata',
                 bin_size=None, snrplot=False):
        SpecPol.__init__(self, oray, metadata, bin_size, snrplot)
        self.wl= None
        self.wl_bin = None
        self.poldata = pd.DataFrame(columns = ['wl', 'v', 'v_r'], dtype='float64')

    def calculate(self, flag='tar'):
        """
        Calculates the Stokes v, its error and epsilon from the ordinary and extra-ordinary ray fluxes.

        Notes
        -----
        Creates an output file containing the poldata data frame (filled by this method)

        Parameters
        ----------
        flag : str
            Which flag to look for in the "Flag" column of metadata. Default is 'tar' which identifies target data.
            'zpol' can be used for zero pol. standard and 'polstd' for the polarisation standard.
            Custom flags can be input if you have written your own in the "Flag" column of the metadata Data Frame

        Returns
        -------
        poldata : pandas.core.frame.DataFrame
            Contains the polarisation data frame. columns = ['wl', 'v', 'v_r'], dtype='float64'
        """

        self.flag = flag
        # Now getting the data from the files in lists that will be used by the specpol() function.
        ls_F0, ls_F0_r, ls_F1, ls_F1_r = self._get_data()

        self.poldata["wl"] = self.wl_bin # even if have not rebinned since in that case we do self.wl_bin = self.wl

        v_ls = []
        v_r_ls = []

        for i in range(len(ls_F0)):
            v, v_r, epsilon = self._circspecpol(ls_F0[i], ls_F0_r[i], ls_F1[i], ls_F1_r[i])
            v_ls.append(v)
            v_r_ls.append(v_r)

            self.poldata["epsilon"+str(i)] = epsilon

        for num in range(len(v_ls[0])):
            # num indexes the bins each list of Stokes parameters values
            v_to_avg = []
            v_r_to_sum = np.array([])
            for s in range(len(v_ls)):
                v_to_avg.append(v_ls[s][num])
                v_r_to_sum = np.append(v_r_to_sum, (1 / ((v_r_ls[s][num]) ** 2)))

            self.poldata.loc[num, 'v'] = np.average(v_to_avg, weights = v_r_to_sum)
            self.poldata.loc[num, 'v_r'] = np.sqrt(1 / np.sum(v_r_to_sum))


        # ###### CREATING THE TEXT FILE ###### #
        self.pol_file = input('What do you want to name the polarisation file? ')

        try:
            os.remove(self.pol_file + ".pol")
        except OSError:
            pass

        self.poldata.to_csv(self.pol_file+".pol", index=False, sep="\t")

        return self.poldata

    def _get_data(self):
        """
        This takes the flux data from the text files given by IRAF and sorts them in lists for later use.

        Returns
        -------
        tuple of 4 lists
            Each list corresponds to the normalised flux difference for a given HWRP angle or its error.
            Each lsit contains N 1D numpy arrays containing N noral;ised flux differences spectra (or error), one
            for each set of spectropolarimetric data.

         """
        files_45_deg = SpecPol._list_files(self, 45.0, linpol=False)
        files_315_deg = SpecPol._list_files(self, 315, linpol=False)

        errormessage = "It seems you don't have the same number of images for each retarder plate angle. " \
                       "This may not be code breaking but only complete sets of retarder plate angles can be used."

        assert len(files_45_deg) == len(files_315_deg), errormessage

        ls_F0 = []
        ls_F0_r = []

        ls_F1 = []
        ls_F1_r = []

        for file_list in [files_45_deg, files_315_deg]:
            nbre_sets = len(file_list) / 4
            nbre_sets_remainder = len(file_list) % 4
            assert nbre_sets_remainder == 0, "There should be 4 data files for each image (2 apertures x (flux + err))"
            print("4 Files per images... All good here")

        for i in range(int(nbre_sets)):
            step = i * 4
            files_45_deg_subset = files_45_deg[0 + step:4 + step]
            files_315_deg_subset = files_315_deg[0 + step:4 + step]

            if i == 0:
                check = True
            wl0, F0, F0_r = SpecPol._flux_diff_from_file(self, files_45_deg_subset, check_bin=check)
            ls_F0.append(F0)
            ls_F0_r.append(F0_r)
            wl1, F1, F1_r = SpecPol._flux_diff_from_file(self, files_315_deg_subset)
            ls_F1.append(F1)
            ls_F1_r.append(F1_r)

        assert len(wl0) == len(wl1), "Wavelength bins not homogeneous. This will be an issue."

        return ls_F0, ls_F0_r, ls_F1, ls_F1_r

    def _circspecpol(self, F0, F0_r, F1, F1_r):
        """
        Calulates the circular polarisation from the normalised flux differences

        Parameters
        ----------
        Four 1D numpy arrays
            The normalised flux differences (and errors) for each HWRP angle.

        Returns
        -------
        v, v_r and epsilon
        """
        # Now Stokes parameters and degree of pol.
        v = 0.5 * (F0 - F1)
        v_r = 0.5*np.sqrt(F0_r**2 + F1_r**2)

        # Now calculating epsilon q and epsilon u and Delta epsilon.
        epsilon = 0.5 * (F0 + F1)
        return v*100, v_r * 100, epsilon*100

    def plot(self):

        f, axarr = plt.subplots(2, 1, figsize=(10, 10), sharex=True)
        plt.subplots_adjust(hspace=0)

        # First axis is v
        wl = self.poldata["wl"]
        v = self.poldata["v"]
        v_r = self.poldata['v_r']
        axarr[0].errorbar(wl, v, yerr=v_r, c='#D92F2F')
        axarr[0].axhline(0, 0, ls='--', c='k')
        vmax = -1000
        vmin = 10000
        for i in range(len(wl)):
            if wl[i] > 4500 and v[i] > vmax:
                vmax = v[i]
            if wl[i] > 4500 and v[i] < vmin:
                vmin = v[i]

        axarr[0].set_ylim([vmin - 0.4, vmax + 0.4])
        axarr[0].set_ylabel('v(%)', fontsize=14)

        # And then the Delta epsilons of each data set.

        delta_cols = [col for col in self.poldata.columns if 'epsilon' in col]
        for column in delta_cols:
            axarr[1].plot(wl, self.poldata[column], alpha=0.8)
            print("MEAN Delta epsilon =", self.poldata[column].values.mean(),
                  "STDV =", self.poldata[column].values.std() )

        axarr[1].set_ylabel(r"$\Delta \epsilon$", fontsize=16)
        axarr[1].set_ylim([-6, 6])
        plt.xlim([3500, 10000])

        save_cond = input("do you want to save the plot?(Y/n): ")
        if save_cond == "y" or save_cond == "Y" or save_cond == "":
            plt.savefig(self.pol_file+ ".png")
            print("Plot saved")
        else:
            print("Plot not saved")

        plt.show()


def mk_flx_spctr(metadata = 'metadata', fileloc='.', flag = 'tar', output = None, header = True, front="1D_c_d"):
    """
    Combines all the flux calibrated apertures to create the flux spectrum.

    Notes
    -----
    Creates a text file with 3 columns columns: wavelength flux errors
    """
    if isinstance(metadata, str):
        metadataframe = pd.read_csv(metadata, sep='\t')

    assert "Filename" in list(metadataframe), "Are you sure "+metadata+" is a (or location to a) Meta Data frame? "

    fluxdata = pd.DataFrame(columns=['wl', 'f', 'f_r'], dtype='float64')

    if output is None:
        output = input('What do you want to call the output file? ')
    output += ".flx"

    files = [str(metadataframe.loc[i,"Filename"])[:-5] for i in xrange(len(metadataframe["Filename"])) \
                  if metadataframe.loc[i, 'Flag'] == flag ]

    file_list1 = [front+filename+"_ap1.txt" for filename in files]
    file_list1r = [front+filename+"_ap1_err.txt" for filename in files]
    file_list2 = [front+filename+"_ap2.txt" for filename in files]
    file_list2r = [front+filename+"_ap2_err.txt" for filename in files]
    file_list = file_list1 + file_list1r + file_list2 + file_list2r

    for filename in file_list:

        if 'err' in filename:
            wl, fr = np.loadtxt(filename, unpack=True, usecols=(0, 1))

            if len(fluxdata['wl']) == 0:
                fluxdata['wl'] = wl
                fluxdata.fillna(value=0, inplace=True)

            fluxdata['f_r'] += fr**2

        if 'err' not in filename:
            wl, f = np.loadtxt(filename, unpack=True, usecols=(0, 1))
            if len(fluxdata['wl']) == 0:
                fluxdata['wl'] = wl
                fluxdata.fillna(value=0, inplace=True)

            fluxdata['f'] += f

    fluxdata['f_r'] = np.sqrt(fluxdata['f_r'])

    try:
        os.remove(output)
    except OSError:
        pass

    fluxdata.to_csv(output, sep='\t', header = header, index=False)

    return fluxdata

# ##################  THE FOLLOWING IS FOR BACKWARDS COMPATIBILITY ONLY##########################

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
    except OSError:
        pass

    for filename in os.listdir("."):
        # Here we go through the headers and take out the information that we need to create the fits files.
        # I think the header names and variable names are explicit enough
        if "fits" in filename:
            hdulist = fits.open(filename)

            exptime = hdulist[0].header['EXPTIME']

            binx = hdulist[0].header['HIERARCH ESO DET WIN1 BINX']
            biny = hdulist[0].header['HIERARCH ESO DET WIN1 BINY']
            binning = str(binx) + "x" + str(biny)
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
                    angle = hdulist[0].header['HIERARCH ESO INS RETA2 POSANG']
                except:
                    try:
                        angle = str(hdulist[0].header['HIERARCH ESO INS RETA4 ROT']) + '*'
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
                    airm = (airm_i + airm_f) / 2
                else:
                    airm = 'None'
            else:
                airm = 'None'

            with open('image_info.txt', 'a') as f:
                f.write(filename[:-5].ljust(15) + str(name).ljust(23) + str(angle).ljust(7)
                        + str(exptime).ljust(10) + str(airm).ljust(8) + str(grism).ljust(12)
                        + binning.ljust(5) + str(size_x).ljust(5) + str(one_over_gain).ljust(7)
                        + str(ron).ljust(5) + str(date) + "\n")

    os.system("sort image_info.txt -o image_info.txt")  # To have filenames written out  alphabetically

    with open('image_info.txt', 'a') as f:  # Just putting labels on columns at end of file
        f.write('=========================================================================='
                '=============================================== \n')
        f.write(
                'Filename'.ljust(15) + 'ESO label'.ljust(23) + 'Angle'.ljust(7) + 'EXP.TIME'.ljust(
                        10) + 'Airmass'.ljust(
                        8) +
                'Grism'.ljust(12) + 'bin'.ljust(5) + '#Pix'.ljust(5) + '1/gain'.ljust(7) + 'RON'.ljust(
                        5) + 'Date' + "\n")


def hwrpangles(sn_name='CCSN', zeropol_name='Zero_', polstd_name='NGC2024'):
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
    # of the original images with name format SCIENCE_##.fits where ## is a 2 or 3 digit number.

    sorted_files = sorted([filename for filename in os.listdir(".")
                           if 'SCIENCE' in filename and 'ms' not in filename
                           and 'dSCIENCE' not in filename and 'cal' not in filename])

    # We are sorting the image names in 3 list to distinguish the CCSN from the polarised std from the zero pol std.
    zeropol = []
    sn = []
    polstd = []

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
        ls_0 = []
        ls_1 = []
        ls_2 = []
        ls_3 = []
        ls_v_0 = []
        ls_v_1 = []

        for name in list_names:  # This loop goes through each image name within the ccsn, pol or zero pol std list
            hdulist = fits.open(name)
            # So it doesn't crash if there's no info on the HWRP in headers.
            try:
                angle = hdulist[0].header['HIERARCH ESO INS RETA2 ROT']
            except:
                try:
                    angle = hdulist[0].header['HIERARCH ESO INS RETA2 POSANG']
                except:
                    try:
                        angle = str(hdulist[0].header['HIERARCH ESO INS RETA4 ROT']) + '*'
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

        try:
            os.remove(output_name + '.txt')
        except OSError:
            pass

        try:
            os.remove(output_name + '_v.txt')
        except OSError:
            pass

        for i in xrange(len(ls_0)):
            with open(output_name + '.txt', 'a') as f:
                f.write(str(ls_0[i]) + ' ' + str(ls_1[i]) + ' ' + str(ls_2[i]) + ' ' + str(ls_3[i]) + '\n')
            f.close()
        for i in xrange(len(ls_v_0)):
            with open(output_name + '_v.txt', 'a') as fv:
                fv.write(str(ls_v_0[i]) + ' ' + str(ls_v_1[i]) + '\n')
            fv.close()

    return





def flux_diff_from_file(files, ordinary_ray, extra_ray, bin_size, check_bin=False, snrplot=False):
    # keeping this as separate function because making instance within the function means they can be forgotten
    # when come out of it and not take up too much memory

    # Extracting polarised fluxes of o and e ray and their errors according to filenames
    for filename in files:
        if ordinary_ray in filename:
            if 'err' not in filename:
                wl, fo = np.loadtxt(filename, unpack=True, usecols=(0, 1))
            else:
                wl, fo_r = np.loadtxt(filename, unpack=True, usecols=(0, 1))

        if extra_ray in filename:
            if 'err' not in filename:
                wl, fe = np.loadtxt(filename, unpack=True, usecols=(0, 1))
            else:
                wl, fe_r = np.loadtxt(filename, unpack=True, usecols=(0, 1))

    # BINNING
    if bin_size is None:
        wl_bin, fo_bin, fo_bin_r, fe_bin, fe_bin_r = wl, fo, fo_r, fe, fe_r

    else:
        print("Binning to ", bin_size, "Angstrom")
        wl_bin, fo_bin, fo_bin_r = rebin(wl, fo, fo_r, bin_siz=bin_size)
        wl_bin, fe_bin, fe_bin_r = rebin(wl, fe, fe_r, bin_siz=bin_size)

        # To perform a few checks on the binning (in particular the SNR yielded by binning)
        if check_bin is True:
            check_binning(wl, fo, fe, fo_r, fe_r, wl_bin, fo_bin, fe_bin, fo_bin_r, fe_bin_r, bin_size, snrplot=snrplot)

    # Creating the oddinary ray - extraordinary ray data object
    data_object = OERay(wl_bin, fo_bin, fo_bin_r, fe_bin, fe_bin_r)
    F, F_r = data_object.norm_flux_diff()  # and calculating the flux difference

    return wl_bin, F, F_r


def check_binning(wl, fo, fe, fo_r, fe_r, bin_wl, bin_fo, bin_fe, bin_fo_r, bin_fe_r, bin_size, snrplot=False):
    snr_not_binned = np.array((fo + fe) / np.sqrt(fo_r ** 2 + fe_r ** 2))
    snr_expected = snr_not_binned * np.sqrt(bin_size / (wl[1] - wl[0]))
    ind_not_binned_central_wl = int(np.argwhere(wl == min(wl, key=lambda x: abs(x - 6204)))[0])
    snr_not_binned_central_wl = snr_not_binned[ind_not_binned_central_wl]
    snr_central_expected = snr_not_binned_central_wl * np.sqrt(bin_size / (wl[1] - wl[0]))

    snr_binned = np.array((bin_fo + bin_fe) / np.sqrt(bin_fo_r ** 2 + bin_fe_r ** 2))
    ind_central_wl = int(np.argwhere(bin_wl == min(bin_wl, key=lambda x: abs(x - 6204)))[0])
    snr_central_wl = (bin_fo[ind_central_wl] + bin_fe[ind_central_wl]) / \
                     np.sqrt(bin_fo_r[ind_central_wl] ** 2 + bin_fe_r[ind_central_wl] ** 2)

    print("\n======== BEFORE BINNING ======")
    print("MEDIAN SNR ")
    print(np.median(snr_not_binned))

    print("CENTRAL SNR at (", wl[ind_not_binned_central_wl], " A)")
    print(snr_not_binned_central_wl)

    print("======== AFTER BINNING ======")

    print("MEDIAN SNR / EXPECTED ")
    print(np.median(snr_binned), np.median(snr_expected))

    print("CENTRAL SNR / EXPECTED (at ", wl[ind_not_binned_central_wl], " A)")
    print(snr_central_wl, snr_central_expected)
    print("\n")

    if snrplot is True:
        plt.plot(wl, snr_expected, marker='o', label='Expected')
        plt.plot(bin_wl, snr_binned, marker='x', label='Calculated after binning')
        plt.legend()
        plt.title("SNR")
        plt.show()

    return


def lin_specpol(oray='ap2', hwrpafile='hwrpangles.txt',
                bin_size=None, e_min_wl=3775,
                bayesian_pcorr=False, p0_step=0.01, snrplot=False):
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
    bin_size : int, optional
        Size of the wavelength bins in Angstrom. Default is 15A.
    e_min_wl : string, optional
        The first wavelength of the range within which Delta epsilons will be calculated. Default is 3775 (ang).
    bayesian_pcorr : bool, optional
        Turns on or off the bayesian p debiasing method (J. L. Quinn 2012). If False then the step function method will
        be used (wang et al 1997). Default is False. Does not work if p large and errors small.
    p0_step : float, optional
        Step size (and starting point) of the p0 distribution. If the step is larger that an observed value of p then
        the code will fail, and you should decrease the step size. Also increases the run time significantly.
        Default is 0.01
    snrplot : Boolean
        If True and bin_size is not None, plots of expected Vs calculated SNR after binning will be showed for the
        first ordinary and extra ordinary ray at 0 degree HWRP angle.
    """
    if oray == 'ap2':
        eray = 'ap1'
    elif oray == 'ap1':
        eray = 'ap2'

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
        sorted_files = sorted([filename for filename in os.listdir(".")
                               if 'dSCIENCE' in filename and 'fits' not in filename and 'c_' not in filename])
        ls_F0 = []
        ls_F0_r = []

        ls_F1 = []
        ls_F1_r = []

        ls_F2 = []
        ls_F2_r = []

        ls_F3 = []
        ls_F3_r = []

        valid1 = re.compile('SCI')
        valid2 = re.compile('STD')
        find_nbr = re.compile('\d{1,3}')  # This is what we'll look for in filename: a number 1-3 digits long
        # The first part searches for the
        files_0_deg = []
        files_22_deg = []
        files_45_deg = []
        files_67_deg = []

        for filename in sorted_files:
            nbr_in_file_name = "PasLa"
            # finding the number in the filename. Searched through filename for a 1-3 digit number and returns it.
            try:
                if valid1.search(filename) or valid2.search(filename):
                    nbr_in_file_name = find_nbr.search(filename[1:]).group()
                    # removing first character as files start with a 1 usually and that messes things up
            except AttributeError:
                print("Couldn't find a number in this filename - passing")
                pass
            # The following compares the number in the filename to the number in ls_0 to see if the image
            # correspond to a 0 deg HWRP angle set up.

            if str(nbr_in_file_name) in ls_0:
                files_0_deg.append(filename)

            # Same thing as the first loop but for 22.5 HWRP
            elif str(nbr_in_file_name) in ls_1:
                files_22_deg.append(filename)

            # Same thing as the first loop but for 45 HWRP
            elif str(nbr_in_file_name) in ls_2:
                files_45_deg.append(filename)

            # Same thing as the first loop but for 67.5 HWRP
            elif str(nbr_in_file_name) in ls_3:
                files_67_deg.append(filename)

        for file_list in [files_0_deg, files_22_deg, files_45_deg, files_67_deg]:
            nbre_sets = len(file_list) / 4
            nbre_sets_remainder = len(file_list) % 4
            assert nbre_sets_remainder == 0, "There should be 4 data files for each image "
            print("4 Files per images... All good here")

        for i in range(int(nbre_sets)):
            step = i * 4
            files_0_deg_subset = files_0_deg[0 + step:4 + step]
            files_22_deg_subset = files_22_deg[0 + step:4 + step]
            files_45_deg_subset = files_45_deg[0 + step:4 + step]
            files_67_deg_subset = files_67_deg[0 + step:4 + step]
            if i == 0:
                check = True
            wl0, F0, F0_r = flux_diff_from_file(files_0_deg_subset, oray, eray, bin_size,
                                                check_bin=check, snrplot=snrplot)

            ls_F0.append(F0)
            ls_F0_r.append(F0_r)
            wl1, F1, F1_r = flux_diff_from_file(files_22_deg_subset, oray, eray, bin_size)
            ls_F1.append(F1)
            ls_F1_r.append(F1_r)
            wl2, F2, F2_r = flux_diff_from_file(files_45_deg_subset, oray, eray, bin_size)
            ls_F2.append(F2)
            ls_F2_r.append(F2_r)
            wl3, F3, F3_r = flux_diff_from_file(files_67_deg_subset, oray, eray, bin_size)
            ls_F3.append(F3)
            ls_F3_r.append(F3_r)

        assert len(wl0) == len(wl1) == len(wl2) == len(wl3), "Wavelength bins not homogenous. This will be an issue."

        return wl0, ls_F0, ls_F0_r, ls_F1, ls_F1_r, ls_F2, ls_F2_r, ls_F3, ls_F3_r

    def specpol(wl, F0, F0_r, F1, F1_r, F2, F2_r, F3, F3_r):
        """
        Finds the p, q, u, theta and errors on these quantities for a set of spectropolarimetric data.

        Notes
        -----
        For lin_specpol() use only

        Parameters
        ----------
        wl : array
            Wavelengths


        Returns
        -------
        arrays
            p, q, and u in percent, with associated errors, as well as theta in degrees ( 0 < theta < 180) and its
            errors.
        """

        # Now Stokes parameters and degree of pol.
        q = 0.5 * (F0 - F2)
        u = 0.5 * (F1 - F3)
        q_r = 0.5 * np.sqrt(F0_r ** 2 + F2_r ** 2)
        u_r = 0.5 * np.sqrt(F1_r ** 2 + F3_r ** 2)
        p, p_r = pol_deg(q, u, q_r, u_r)

        # Arrays where we're going to store the values of p and Stokes parameters and P.A
        # after we've applied corrections.
        pf = np.array([])
        qf = np.array([])
        uf = np.array([])
        theta = np.array([])

        # We take our chromatic zero-angles and interpolate them to match the wavlength bins of our data.
        wl2, thetaz = np.loadtxt(zero_angles, unpack=True, usecols=(0, 1))
        theta0 = np.interp(wl, wl2, thetaz)

        # Now we apply corrections to the P.A
        for t in range(len(q)):
            theta_t = 0.5 * m.atan2(u[t], q[t])
            theta_r = 0.5 * np.sqrt(((u_r[t] / u[t]) ** 2 + (q_r[t] / q[t]) ** 2) * (1 / (1 + (u[t] / q[t]) ** 2)) ** 2)
            theta_t = (theta_t * 180.0) / m.pi
            theta_r = (theta_r * 180.0) / m.pi
            if theta_t < 0:
                theta_t = 180 + theta_t  # Making sure P.A is within limit 0<theta<180 deg

            theta_cor = theta_t - theta0[t]
            theta_cor_rad = (theta_cor / 180.0) * m.pi
            theta = np.append(theta, theta_cor)
            q_t = p[t] * m.cos(2 * theta_cor_rad)  # Re-calculating Stokes parameters
            u_t = p[t] * m.sin(2 * theta_cor_rad)
            qf = np.append(qf, q_t * 100)  # Filling our arrays of final Stokes parameters and p.
            uf = np.append(uf, u_t * 100)
            pf = np.append(pf, np.sqrt(q_t ** 2 + u_t ** 2) * 100)

        # Now calculating epsilon q and epsilon u and Delta epsilon.
        eq = (0.5 * F0 + 0.5 * F2) * 100  # in percent
        eu = (0.5 * F1 + 0.5 * F3) * 100
        delta_e = eq - eu
        stdv_dequ = []
        try:
            for i in xrange(len(wl)):
                if wl[i] > e_min_wl:
                    dequ = eq[i] - eu[i]
                    stdv_dequ.append(dequ)
        except IndexError:
            dequ = eq - eu
            stdv_dequ.append(dequ)

        stdv_e = np.std(stdv_dequ)
        avg_e = np.average(stdv_dequ)

        return pf, p_r * 100, qf, q_r * 100, uf, u_r * 100, theta, theta_r, delta_e, avg_e, stdv_e

    #########################################
    #            LIN SPECPOL MAIN           #
    #########################################

    # list of files corresponding to each angle (0, 22.5, 45, 67.5)
    ls_0, ls_1, ls_2, ls_3 = np.genfromtxt(hwrpafile, dtype='str', unpack=True, usecols=(0, 1, 2, 3))

    # Now getting the data from the files in lists that will be used by the specpol() function.
    wl, ls_F0, ls_F0_r, ls_F1, ls_F1_r, ls_F2, ls_F2_r, ls_F3, ls_F3_r = get_data(ls_0, ls_1, ls_2, ls_3)

    qls = []
    qrls = []
    uls = []
    urls = []
    delta_es = []
    avg_es = []
    stdv_es = []

    for i in range(len(ls_F0)):
        p, pr, q, qr, u, ur, theta, thetar, delta_e, avg_e, stdv_e = specpol(wl, ls_F0[i], ls_F0_r[i], ls_F1[i],
                                                                             ls_F1_r[i], ls_F2[i], ls_F2_r[i],
                                                                             ls_F3[i], ls_F3_r[i])
        qls.append(q)
        qrls.append(qr)
        uls.append(u)
        urls.append(ur)
        delta_es.append(delta_e)
        avg_es.append(avg_e)
        stdv_es.append(stdv_e)

    # Where we'll put the final values of the Stokes parameters and their errors.
    qf = np.array([])
    uf = np.array([])
    qrf = np.array([])
    urf = np.array([])

    for num in range(len(qls[0])):
        # num indexes the bins each list of Stokes parameters values
        q_to_avg = []
        u_to_avg = []
        qr_to_sum = np.array([])
        ur_to_sum = np.array([])
        for s in range(len(qls)):
            # s indexes the data set from which we are taking a particular Stoke parameter
            # We want to average values fo all data sets at each wavelength bins. For example say I have
            # 3 data sets, at 5000 A say, I am gonna take the 3 values of q in each data set at 5000 A and
            # average them. Do the same accross the whole spectrum and with each Stoke parameter to get final results.
            q_to_avg.append(qls[s][num])
            u_to_avg.append(uls[s][num])

            # For the next 4 lines I am calculating the error on the mean and putting it in final list of errors on
            # Stokes parameters
            qr_to_sum = np.append(qr_to_sum, 1 / ((qrls[s][num]) ** 2))
            ur_to_sum = np.append(ur_to_sum, 1 / ((urls[s][num]) ** 2))
        qrf = np.append(qrf, np.sqrt(1 / np.sum(qr_to_sum)))
        urf = np.append(urf, np.sqrt(1 / np.sum(ur_to_sum)))

        qf = np.append(qf, np.average(q_to_avg, weights=qr_to_sum))
        uf = np.append(uf, np.average(u_to_avg, weights=ur_to_sum))

    # Once I have my final Stokes parameters I can calculate the final degree of polarisation (and error).
    pf = np.sqrt(qf ** 2 + uf ** 2)
    prf = (1 / pf) * np.sqrt((qf * qrf) ** 2 + (uf * urf) ** 2)

    # And finally the P.A !
    thetaf = np.array([])
    thetarf = np.array([])
    for t in range(len(qrf)):
        thetaf_t = 0.5 * m.atan2(uf[t], qf[t])
        thetarf_t = 0.5 * np.sqrt(
                ((urf[t] / uf[t]) ** 2 + (qrf[t] / qf[t]) ** 2) * (1 / (1 + (uf[t] / qf[t]) ** 2)) ** 2)
        thetaf_t = (thetaf_t * 180.0) / m.pi
        thetarf_t = (thetarf_t * 180.0) / m.pi
        if thetaf_t < 0:
            thetaf_t = 180 + thetaf_t  # Again need to make sure the range is within 0-180 deg
        thetaf = np.append(thetaf, thetaf_t)
        thetarf = np.append(thetarf, thetarf_t)

    # #### P Bias Correction #### #

    #  If bayesian_pcorr is False, P will be debiased as in Wang et al. 1997 using a step function
    if bayesian_pcorr is False:
        print("Step Func - p correction")
        pfinal = np.array([])
        for ind in range(len(pf)):
            condition = pf[ind] - prf[ind]
            if condition > 0:
                p_0i = pf[ind] - ((float(prf[ind] ** 2)) / float(pf[ind]))
            elif condition < 0:
                p_0i = pf[ind]

            pfinal = np.append(pfinal, p_0i)

    # If bayesian_pcorr is True, P will be debiased using the Bayesian method described by J. L. Quinn 2012
    #  the correceted p is pbar_{0,mean} * sigma. pbar_{0,mean} is given by equation 47 of J. L. Quinn 2012

    if bayesian_pcorr is True:
        print("Bayesian - p correction")
        sigma = (qrf + urf) / 2
        pbar = pf / sigma
        pfinal = np.array([])
        for j in range(len(pbar)):
            p0 = np.arange(p0_step, pbar[j], p0_step)
            rho = np.array([])
            for i in range(len(p0)):
                tau = (sigma[j] ** 2) * 2 * p0[i]
                pp0 = pbar[j] * p0[i]

                # print pbar[j]**2, (p0[i]**2)/2, special.iv(0, pp0)
                rice_distribution = pbar[j] * np.exp(-((pbar[j] ** 2 + p0[i] ** 2) / 2)) * special.iv(0, pp0)
                if m.isnan(rice_distribution) is True or m.isinf(rice_distribution) is True:
                    # print pf[j], sigma[j], pbar[j], p0[i-1], pbar[j]*p0[i-1]
                    # print pf[j], sigma[j], pbar[j], p0[i], pp0
                    print("Infinite values encountered. Resulting polarisation may be invalid. Use the Step function method")

                    # return
                rhoi = rice_distribution * tau
                rho = np.append(rho, rhoi)

            p0mean = np.average(p0, weights=rho)
            pfinal = np.append(pfinal, p0mean * sigma[j])  # !!!! need to multiply by sigma to get p0 and not p0/bar.

    # ###### CREATING THE TEXT FILE ###### #
    pol_file = input('What do you want to name the polarisation file? ')
    try:
        os.remove(pol_file + ".pol")
        os.remove(pol_file + ".delta")
    except:
        print('kittens')

    for l in xrange(len(wl)):
        with open(pol_file + ".pol", 'a') as pol_f:
            pol_f.write(
                    str(wl[l]) + '    ' + str(pfinal[l]) + '    ' + str(prf[l]) + '    ' + str(qf[l]) + '    ' + str(
                            qrf[l]) +
                    '    ' + str(uf[l]) + '    ' + str(urf[l]) + '    ' + str(thetaf[l]) + '    ' + str(
                            thetarf[l]) + '\n')
            # writing the file containing delta epsilon
            # because delta_es is a list of lists the file write out is a bit convoluted

    with open(pol_file + ".delta", 'a') as delta_f:
        if len(delta_es) > 1:  # Case were I have more than one list in delta_es if mroe than one set of data
            for deltas in zip(wl, *delta_es):  # first I am zipping my lists. The * to unpack unknwon number of lists
                to_write = str(deltas[0])  # I am creating a variable to store my output string, deltas[0] is the wvlgth
                for delta in deltas[1:]:  # then iterating over the values in each column
                    to_write += '    ' + str(delta)  # to add them onto the string
                to_write += '\n'  # and when that's done I add the end of line character
                delta_f.write(to_write)  # and finally write it into my file
        else:  # if I have only one list in delta_es then zip is going to fail and I can just do the following
            for wlgth, delta in zip(wl, delta_es[0]):  # first I am zipping my lists
                delta_f.write(str(wlgth) + '    ' + str(delta) + '\n')

    # ###### MAKING PLOTS ########
    # Just to check that everything looks right.

    f, axarr = plt.subplots(5, 1, figsize=(8, 8), sharex=True)
    plt.subplots_adjust(hspace=0)

    # First axis is p
    print(wl)
    axarr[0].errorbar(wl, pf, yerr=prf, c='#D92F2F')
    axarr[0].axhline(0, 0, ls='--', c='k')
    pmax = -1000
    for i in range(len(wl)):
        if wl[i] > 4500 and pf[i] > pmax:
            pmax = pfinal[i]

    axarr[0].set_ylim([-0.1, pmax + 0.4])
    axarr[0].set_ylabel('p(%)', fontsize=14)

    # Then q
    axarr[1].errorbar(wl, qf, yerr=qrf, c='#D92F2F')
    axarr[1].axhline(0, 0, ls='--', c='k')
    qmax = -1000
    qmin = 1000
    for i in range(len(wl)):
        if wl[i] > 4500 and qf[i] > qmax:
            qmax = qf[i]
        if wl[i] > 4500 and qf[i] < qmin:
            qmin = qf[i]
    axarr[1].set_ylim([qmin - 0.3, qmax + 0.3])
    axarr[1].set_ylabel('q(%)', fontsize=14)

    # And u
    axarr[2].errorbar(wl, uf, yerr=urf, c='#D92F2F')
    axarr[2].axhline(0, 0, ls='--', c='k')
    umax = -1000
    umin = 1000
    for i in range(len(wl)):
        if wl[i] > 4500 and uf[i] > umax:
            umax = uf[i]
        if wl[i] > 4500 and uf[i] < umin:
            umin = uf[i]
    axarr[2].set_ylim([umin - 0.3, umax + 0.3])
    axarr[2].set_ylabel('u(%)', fontsize=14)

    # P.A
    axarr[3].errorbar(wl, thetaf, yerr=thetarf, c='#D92F2F')
    axarr[3].axhline(0, 0, ls='--', c='k')
    axarr[3].set_ylim([-0, 180])
    axarr[3].set_ylabel('theta', fontsize=14)

    # And finally the Delta epsilons of each data set.
    for i in range(len(delta_es)):
        axarr[4].plot(wl, delta_es[i], alpha=0.8)
        print("Average Delta epsilon =", avg_es[i], "STDV =", stdv_es[i])

    axarr[4].set_ylabel(r"$\Delta \epsilon", fontsize=16)
    axarr[4].set_ylim([-4.0, 4.0])
    plt.xlim([3500, 10000])

    save_cond = input("do you want to save the plot?(Y/n): ")
    if save_cond == "y" or save_cond == "Y" or save_cond == "":
        plt.savefig(pol_file + ".png")
        print("Plot saved")
    else:
        print("Plot not saved")

    plt.show()

    return pfinal, prf, qf, qrf, uf, urf, thetaf, thetarf, delta_es


def circ_specpol(oray='ap2', hwrpafile='hwrpangles_v.txt', bin_size=None, e_min_wl=3775):
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
    if oray == 'ap2':
        eray = 'ap1'
    elif oray == 'ap1':
        eray = 'ap2'

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
        list_file = []
        for name in os.listdir('.'):
            list_file.append(name)
        sorted_files = sorted(list_file)

        ls_fo0 = []
        ls_fe0 = []
        ls_fo0_err = []
        ls_fe0_err = []

        ls_fo1 = []
        ls_fe1 = []
        ls_fo1_err = []
        ls_fe1_err = []

        valid1 = re.compile('SCI')
        valid2 = re.compile('STD')
        find_nbr = re.compile('\d{1,3}')  # This is what we'll look for in filename: a number 1-3 digits long
        # The first part searches for the

        for filename in sorted_files:
            nbr_in_file_name = "PasLa"
            # finding the number in the filename. Searched through filename for a 1-3 digit number and returns it.
            try:
                if valid1.search(filename) or valid2.search(filename):
                    nbr_in_file_name = find_nbr.search(filename[1:]).group()
                    # removing first character as files start with a 1 usually and that messes things up
            except AttributeError:
                print("Couldn't find a number in this filename - passing")
                pass

            # This condition is related to the naming convention I have adopted.
            if 'c_' not in filename:
                # The following compares the number in the filename to the number in ls_0 to see if the image
                # correspond to a 0 deg HWRP angle set up. The naming convention is crucial for this line to work
                # as it keeps the number in the filename in the location: filename[-10:-8] or filename[-14:-12] for
                # flux and flux_error files, respectively.
                if nbr_in_file_name in ls_0:
                    # Now we put the filename in the right list, oray or eray and flux or error on flux.
                    if oray in filename:
                        if 'err' not in filename:
                            wl, fo = np.loadtxt(filename, unpack=True, usecols=(0, 1))
                            ls_fo0.append(fo)
                        else:
                            wl, fo = np.loadtxt(filename, unpack=True, usecols=(0, 1))
                            ls_fo0_err.append(fo)

                    if eray in filename:
                        if 'err' not in filename:
                            wl, fe = np.loadtxt(filename, unpack=True, usecols=(0, 1))
                            ls_fe0.append(fe)
                        else:
                            wl, fe = np.loadtxt(filename, unpack=True, usecols=(0, 1))
                            ls_fe0_err.append(fe)

                # Same thing as the first loop but for 45 HWRP
                if nbr_in_file_name in ls_1:
                    if oray in filename:
                        if 'err' not in filename:
                            wl, fo = np.loadtxt(filename, unpack=True, usecols=(0, 1))
                            ls_fo1.append(fo)
                        else:
                            wl, fo = np.loadtxt(filename, unpack=True, usecols=(0, 1))
                            ls_fo1_err.append(fo)

                    if eray in filename:
                        if 'err' not in filename:
                            wl, fe = np.loadtxt(filename, unpack=True, usecols=(0, 1))
                            ls_fe1.append(fe)
                        else:
                            wl, fe = np.loadtxt(filename, unpack=True, usecols=(0, 1))
                            ls_fe1_err.append(fe)

        return wl, ls_fo0, ls_fe0, ls_fo0_err, ls_fe0_err, ls_fo1, ls_fe1, ls_fo1_err, ls_fe1_err

    def rebin(wl, f, r, bin_siz=bin_size):
        """
        To rebin my flux spectra

        Parameters
        ----------
        wl : array
            1D array containing the wavelengths to be rebinned
        f : array
            1D array containing the fluxes to be rebinned
        r : array
            1D array containing the errors on the fluxes to be rebinned
        bin_siz : int
            Size of the new bins in Angstrom.

        Returns
        -------
            tuple of 1D arrays: wl, f, err all rebinned to the new bin size
        if bin_siz is None:
            print "No binning"
            return wl, f, r
        """
        wl = np.array(wl)
        f = np.array(f)
        r = np.array(r)
        small_bin_sizes = []

        bins_f = np.zeros(int((max(wl) - min(wl)) / bin_siz) + 1)  # new flux bins, empty for now
        bins_w = np.zeros(int((max(wl) - min(wl)) / bin_siz) + 1)  # new error bins, empty for now

        weights = 1 / (r ** 2)

        for i in range(len(wl) - 1):
            n = int((wl[i] - min(wl)) / bin_siz)  # n is the number of the new bin
            small_bin_sizes.append((wl[i + 1] - wl[i]))  # filling list of small bin sizes

        bin_centers = [(min(wl)) + bin_siz * n for n in range(len(bins_f))]  # finding the new bin centers

        bin_edges = [bin_centers[0]] + [bin1 + bin_siz for bin1 in bin_centers]  # finding the new bin edges

        ind_edge = []  # in this list I'll put the index of the array wl corresponding to the wavelength values
        # that are close to the bin edges.

        for edge in bin_edges:
            i_wl_at_edge = min(range(len(wl[:-1])), key=lambda i: abs(edge - wl[i]))
            # this is to find the small bin that is closest to the edge of the new bin
            # print wl[i_wl_at_edge], small_bin_sizes[i_wl_at_edge]
            ind_edge.append(i_wl_at_edge)

        for i in range(len(wl)):
            n = int((wl[i] - min(wl)) / bin_siz)
            if i in ind_edge:
                j = ind_edge.index(i)  # finding index j of the wavelength index i I am interested in
                edge = bin_edges[j]  # the edge to compare to wl[i] will then be at bin_edges[j]

                if wl[i] < edge:
                    frac_overlap = (wl[i] + small_bin_sizes[i] / 2 - edge) / (small_bin_sizes[i])
                    try:
                        bins_f[n] += f[i] * weights[i] * (1 - frac_overlap)
                        bins_w[n] += weights[i] * (1 - frac_overlap)
                        bins_f[n + 1] += f[i] * weights[i] * frac_overlap
                        bins_w[n + 1] += weights[i] * frac_overlap

                    except IndexError:
                        print("Index Error at ", wl[i])
                        pass

                elif wl[i] > edge:
                    frac_overlap = (wl[i] + small_bin_sizes[i] / 2 - edge) / (small_bin_sizes[i])
                    try:
                        bins_f[n] += f[i] * weights[i] * frac_overlap
                        bins_w[n] += weights[i] * frac_overlap
                        bins_f[n + 1] += f[i] * weights[i] * (1 - frac_overlap)
                        bins_w[n + 1] += weights[i] * (1 - frac_overlap)

                    except IndexError:
                        print("Index Error at ", wl[i])
                        pass

            else:
                try:
                    bins_f[n] += f[i] * weights[i]
                    bins_w[n] += weights[i]
                except IndexError:
                    print("Index Error at ", wl[i])
                    pass

        for i in range(len(bin_centers)):
            if bins_w[i] == 0.0:
                print(bin_centers[i], bins_w[i])

        bins_f[:-1] /= bins_w[:-1]  # normalise weighted values by sum of weights to get weighted average
        bins_err = np.sqrt(1 / bins_w[:-1])

        return bin_centers[:-1], bins_f[:-1], bins_err

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
            F = (fo[i] - fe[i]) / (fo[i] + fe[i])
            F_err = m.fabs(F) * np.sqrt(
                    ((fo_err[i] ** 2) + (fe_err[i] ** 2)) * ((1 / (fo[i] - fe[i]) ** 2) + (1 / (fo[i] + fe[i]) ** 2)))
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
        eps = np.array([])

        for i in xrange(len(v0)):
            v_i = 0.5 * (v1[i] - v0[i])  # v1 is 45 and v0 is -45 here
            eps_el = 0.5 * (v1[i] + v0[i])
            eps = np.append(eps, eps_el * 100)
            v_i_err = 0.5 * np.sqrt(v1_err[i] ** 2 + v0_err[i] ** 2)
            v_i = v_i * 100
            v_i_err = v_i_err * 100
            v.append(v_i)
            v_err.append(v_i_err)
            # eps= eps*10

        # cond = (wl>e_min_wl)
        eps_crop = eps[int(np.argwhere(wl > e_min_wl)[0]):]
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
            numerator = values[i] / (err[i] * err[i])
            denominator = 1 / (err[i] * err[i])
            num_to_sum.append(numerator)
            den_to_sum.append(denominator)
        mean = np.sum(num_to_sum) / np.sum(den_to_sum)
        err_mean = np.sqrt(1 / (np.sum(den_to_sum)))

        return mean, err_mean

    #########################################
    #                SPECPOL MAIN           #
    #########################################

    # list of files corresponding to each angle (0)
    ls_0, ls_1 = np.genfromtxt(hwrpafile, dtype='str', unpack=True, usecols=(0, 1))

    # Now getting the data from the files in lists that will be used by the specpol() function.
    wl, ls_fo0, ls_fe0, ls_fo0_err, ls_fe0_err, ls_fo1, ls_fe1, ls_fo1_err, ls_fe1_err = get_data(ls_0, ls_1)

    v_ls = []
    verr_ls = []
    ls_eps = []
    avg_eps = []
    stdv_eps = []

    for i in range(len(ls_fo0)):
        print("Set ", i + 1)
        snr_nb0 = np.array((ls_fo0[i] + ls_fe0[i]) / np.sqrt(ls_fo0_err[i] ** 2 + ls_fe0_err[i] ** 2))  # SNR not binned
        snr_nb1 = np.array((ls_fo1[i] + ls_fe1[i]) / np.sqrt(ls_fo1_err[i] ** 2 + ls_fe1_err[i] ** 2))

        # print (wl[1]-wl[0])
        snr_exp0 = snr_nb0 * np.sqrt(bin_size / (wl[1] - wl[0]))  # expected snr
        # snr_exp1 = snr_nb1*np.sqrt(bin_size/(wl[1]-wl[0]))


        ind_central_wl = int(np.argwhere(wl == min(wl, key=lambda x: abs(x - 6204)))[0])
        snr_nbc0 = snr_nb0[ind_central_wl]
        snr_nbc1 = snr_nb1[ind_central_wl]

        if bin_size is not None:
            print("Rebinning to ", str(bin_size), " Angstrom")
            bin_wl, bin_fo0, bin_fo0_err = rebin(wl, ls_fo0[i], ls_fo0_err[i])
            bin_wl, bin_fe0, bin_fe0_err = rebin(wl, ls_fe0[i], ls_fe0_err[i])
            bin_wl, bin_fo1, bin_fo1_err = rebin(wl, ls_fo1[i], ls_fo1_err[i])
            bin_wl, bin_fe1, bin_fe1_err = rebin(wl, ls_fe1[i], ls_fe1_err[i])

        elif bin_size is None:
            print("Not Rebinning")
            bin_wl, bin_fo0, bin_fo0_err = wl, ls_fo0[i], ls_fo0_err[i]
            bin_wl, bin_fe0, bin_fe0_err = wl, ls_fe0[i], ls_fe0_err[i]
            bin_wl, bin_fo1, bin_fo1_err = wl, ls_fo1[i], ls_fo1_err[i]
            bin_wl, bin_fe1, bin_fe1_err = wl, ls_fe1[i], ls_fe1_err[i]

        # v, verr, eps, eps_avg, eps_std = v_1set(bin_wl, bin_fo0, bin_fe0, bin_fo0_err, bin_fe0_err, bin_fo1, bin_fe1, bin_fo1_err, bin_fe1_err)
        # v_ls.append(v)
        # verr_ls.append(verr)
        # ls_eps.append(eps)
        # avg_eps.append(eps_avg)
        # stdv_eps.append(eps_std)

        ind_central_wl = int(np.argwhere(bin_wl == min(bin_wl, key=lambda x: abs(x - 6204)))[0])
        snr_c0 = (bin_fo0[ind_central_wl] + bin_fe0[ind_central_wl]) / np.sqrt(
                bin_fo0_err[ind_central_wl] ** 2 + bin_fe0_err[ind_central_wl] ** 2)
        snr_c1 = (bin_fo1[ind_central_wl] + bin_fe1[ind_central_wl]) / np.sqrt(
                bin_fo1_err[ind_central_wl] ** 2 + bin_fe1_err[ind_central_wl] ** 2)

        print("\n")

        snr0 = (bin_fo0 + bin_fe0) / np.sqrt(bin_fo0_err ** 2 + bin_fe0_err ** 2)
        # snr1 = (bin_fo1 + bin_fe1)/np.sqrt(bin_fo1_err**2 + bin_fe1_err**2)

        snr_c_exp0 = snr_nbc0 * np.sqrt(bin_size / (wl[1] - wl[0]))
        snr_c_exp1 = snr_nbc1 * np.sqrt(bin_size / (wl[1] - wl[0]))

        print("--------- 0 ------ 22.5 ----- 45 ----- 67.5 -------- (CNTR WL ", int(bin_wl[ind_central_wl]), ")")
        print("CNTR WL THEORETICAL: ", int(snr_c_exp0), int(snr_c_exp1))
        print("CNTR WL CALCULATED:  ", int(snr_c0), int(snr_c1))
        print("CNTR WL RESIDUALS (EXPECTED SNR - OBTAINED)")
        print(snr_c_exp0 - snr_c0, snr_c_exp1 - snr_c1)
        print("\n")

        plt.scatter(wl, snr_exp0, c='k', marker='.', label="expected SN 0 deg")
        # plt.scatter(wl, snr_exp1, c='m', marker='.')
        plt.plot(bin_wl, snr0, c='r', marker='+', ls='--', label="calculated", alpha=0.8)
        plt.legend()
        plt.show()

        plt.plot(bin_wl, bin_fo0 / bin_fo0_err, label="F/err 0 deg")
        plt.show()

    for i in range(len(ls_fo0)):
        print("Rebinning to ", str(bin_size), " Angstrom")
        bin_wl, bin_fo0, bin_fo0_err = rebin(wl, ls_fo0[i], ls_fo0_err[i])
        bin_wl, bin_fe0, bin_fe0_err = rebin(wl, ls_fe0[i], ls_fe0_err[i])
        bin_wl, bin_fo1, bin_fo1_err = rebin(wl, ls_fo1[i], ls_fo1_err[i])
        bin_wl, bin_fe1, bin_fe1_err = rebin(wl, ls_fe1[i], ls_fe1_err[i])

        v, verr, eps, eps_avg, eps_std = v_1set(bin_wl, bin_fo0, bin_fe0, bin_fo0_err, bin_fe0_err, bin_fo1, bin_fe1,
                                                bin_fo1_err, bin_fe1_err)
        v_ls.append(v)
        verr_ls.append(verr)
        ls_eps.append(eps)
        avg_eps.append(eps_avg)
        stdv_eps.append(eps_std)

    vf = np.array([])
    vf_err = np.array([])

    for num in range(len(v_ls[0])):
        # num indexes the bins each list of Stokes parameters values
        v_to_avg = []
        verr_to_avg = []
        for s in range(len(v_ls)):
            # s indexes the data set from which we are taking a particular Stoke parameter
            # We want to average values fo all data sets at each wavelength bins. For example say I have
            # 3 data sets, at 5000 A say, I am gonna take the 3 values of q in each data set at 5000 A and
            # average them. Do the same accross the whole spectrum and with each Stoke parameter to get final results.
            v_to_avg.append(v_ls[s][num])
            verr_to_avg.append(verr_ls[s][num])
        vi, vi_err = wghtd_mean(v_to_avg, verr_to_avg)
        vf = np.append(vf, vi)
        vf_err = np.append(vf_err, vi_err)

    # ###### CREATING THE TEXT FILE ###### #
    pol_file = input('What do you want to name the polarisation file? ')

    try:
        os.remove(pol_file + '.pol')
    except:
        print('kittens')
    for l in xrange(len(bin_wl)):
        with open(pol_file + '.pol', 'a') as pol_f:
            pol_f.write(str(bin_wl[l]) + '    ' + str(vf[l]) + '    ' + str(vf_err[l]) + '\n')

    # ###### MAKING PLOT ########
    # Just to check that everything looks right.

    f, axarr = plt.subplots(2, 1, figsize=(10, 10), sharex=True)
    plt.subplots_adjust(hspace=0)

    # First axis is v
    axarr[0].errorbar(bin_wl, vf, yerr=vf_err, c='#D92F2F')
    axarr[0].axhline(0, 0, ls='--', c='k')
    vmax = -1000
    vmin = 10000
    for i in range(len(bin_wl)):
        if bin_wl[i] > 4500 and vf[i] > vmax:
            vmax = vf[i]
        if bin_wl[i] > 4500 and vf[i] < vmin:
            vmin = vf[i]

    axarr[0].set_ylim([vmin - 0.4, vmax + 0.4])
    axarr[0].set_ylabel('v(%)', fontsize=14)

    # And then the Delta epsilons of each data set.
    for i in range(len(ls_eps)):
        axarr[1].plot(bin_wl, ls_eps[i], alpha=0.8)
        print("Average Delta epsilon =", avg_eps[i], "STDV =", stdv_eps[i])

    axarr[1].set_ylabel(r"$\Delta \epsilon", fontsize=16)
    axarr[1].set_ylim([-6, 0])
    plt.xlim([3500, 10000])

    save_cond = input("do you want to save the plot?(Y/n): ")
    if save_cond == "y" or save_cond == "Y" or save_cond == "":
        plt.savefig(pol_file + ".png")
        print("Plot saved")
    else:
        print("Plot not saved")

    plt.show()

    return bin_wl, vf, vf_err, ls_eps


# ########################## V BAND POL ###############################


def lin_vband(oray='ap2', hwrpafile='hwrpangles.txt'):
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

    if oray == 'ap2':
        eray = 'ap1'
    elif oray == 'ap1':
        eray = 'ap2'

    def v_counts(filename, bp_v):
        # Only used by lin_vband(). Finds v counts from spectrum in filename. bp_v is tthe bandpass of V band
        sp = S.FileSpectrum(filename)
        obs = S.Observation(sp, bp_v)

        return obs.wave, obs.flux

    def norm_flux(fo, fe, fo_r, fe_r):
        # normalised flux. Repeat of the one in lin_specpol()
        F = (fo - fe) / (fo + fe)
        F_r = m.fabs(F) * np.sqrt(((fo_r ** 2) + (fe_r ** 2)) * ((1 / (fo - fe) ** 2) + (1 / (fo + fe) ** 2)))

        return F, F_r

    def vpol(fo0, fo1, fo2, fo3, fe0, fe1, fe2, fe3, fo0_r, fe0_r, fo1_r, fe1_r, fo2_r, fe2_r, fo3_r, fe3_r, wl):
        """
        Similar to specpol() in lin_specpol(), but for V-band polarisation.
        """

        # Normalised fluxes
        F0, F0_r = norm_flux(fo0, fe0, fo0_r, fe0_r)
        F1, F1_r = norm_flux(fo1, fe1, fo1_r, fe1_r)
        F2, F2_r = norm_flux(fo2, fe2, fo2_r, fe2_r)
        F3, F3_r = norm_flux(fo3, fe3, fo3_r, fe3_r)

        # Stokes parameters
        q = 0.5 * (F0 - F2)
        u = 0.5 * (F1 - F3)
        q_r = 0.5 * np.sqrt(F0_r ** 2 + F2_r ** 2)
        u_r = 0.5 * np.sqrt(F1_r ** 2 + F3_r ** 2)

        # degree of polarisation
        p = np.sqrt(q * q + u * u)
        p_r = (1 / p) * np.sqrt((q * q_r) ** 2 + (u * u_r) ** 2)

        # Interpolation of chromatic zero angle values
        wl2, thetaz = np.loadtxt(zero_angles, unpack=True, usecols=(0, 1))
        theta0 = np.interp(wl, wl2, thetaz)

        # Finding Polarisation angle
        theta = 0.5 * m.atan2(u, q)
        theta_r = 0.5 * np.sqrt(((u_r / u) ** 2 + (q_r / q) ** 2) * (1 / (1 + (u / q) ** 2)) ** 2)
        theta = (theta * 180.0) / m.pi
        theta_r = (theta_r * 180.0) / m.pi
        if theta < 0:
            theta = 180 + theta
        theta_cor = theta - theta0  # Correction of chromatic zero angle
        theta_cor_rad = (theta_cor / 180.0) * m.pi

        # Correction of Stokes parameters and p from P.A correction
        q = p * m.cos(2 * theta_cor_rad)
        u = p * m.sin(2 * theta_cor_rad)
        p = np.sqrt(q * q + u * u)

        return p * 100, p_r * 100, q * 100, q_r * 100, u * 100, u_r * 100, theta, theta_r

    # ############################### #
    bp_v = S.ObsBandpass("johnson,v")
    wl_v = bp_v.avgwave()

    # list of files corresponding to each angle (0, 22.5, 45, 67.5)
    list_0, list_1, list_2, list_3 = np.genfromtxt(hwrpafile, dtype='str', unpack=True, usecols=(0, 1, 2, 3))

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

    list_file = []
    for name in os.listdir('.'):
        list_file.append(name)

    sorted_files = sorted(list_file)

    valid1 = re.compile('SCI')
    valid2 = re.compile('STD')
    find_nbr = re.compile('\d{1,3}')  # This is what we'll look for in filename: a number 1-3 digits long
    # The first part searches for the

    for filename in sorted_files:
        nbr_in_file_name = "PasLa"
        # finding the number in the filename. Searched through filename for a 1-3 digit number and returns it.
        try:
            if valid1.search(filename) or valid2.search(filename):
                nbr_in_file_name = find_nbr.search(filename[1:]).group()
                # removing first character as files start with a 1 usually and that messes things up
        except AttributeError:
            print("Couldn't find a number in this filename - passing")
            pass

        if 'dSC' in filename and 'c_' not in filename and '.fits' not in filename:
            vflux = v_counts(filename, bp_v)  # finds counts across Vband given spectrum
            counts = np.sum(vflux[1])  # sum over all bins to get total number of counts in Vband
            counts_r = np.sqrt(np.sum(vflux[1] ** 2))  # Poisson noise

            # very similar to what's done in get_data() in lin_specpol() so refer to that
            if nbr_in_file_name in list_0:
                if oray in filename:
                    if 'err' not in filename:
                        o0 = np.append(o0, counts)
                    else:
                        o0_r = np.append(o0_r, counts_r)

                if eray in filename:
                    # print filename
                    if 'err' not in filename:
                        e0 = np.append(e0, counts)
                    else:
                        e0_r = np.append(e0_r, counts_r)

            if nbr_in_file_name in list_1:
                if oray in filename:
                    # print filename
                    if 'err' not in filename:
                        o1 = np.append(o1, counts)
                    else:
                        o1_r = np.append(o1_r, counts_r)

                if eray in filename:
                    # print filename
                    if 'err' not in filename:
                        e1 = np.append(e1, counts)
                    else:
                        e1_r = np.append(e1_r, counts_r)

            if nbr_in_file_name in list_2:
                if oray in filename:
                    # print filename
                    if 'err' not in filename:
                        o2 = np.append(o2, counts)
                    else:
                        o2_r = np.append(o2_r, counts_r)

                if eray in filename:
                    # print filename
                    if 'err' not in filename:
                        e2 = np.append(e2, counts)
                    else:
                        e2_r = np.append(e2_r, counts_r)

            if nbr_in_file_name in list_3:
                if oray in filename:
                    # print filename
                    if 'err' not in filename:
                        o3 = np.append(o3, counts)
                    else:
                        o3_r = np.append(o3_r, counts_r)

                if eray in filename:
                    # print filename
                    if 'err' not in filename:
                        e3 = np.append(e3, counts)
                    else:
                        e3_r = np.append(e3_r, counts_r)

    p_ls = []
    pr_ls = []
    q_ls = []
    qr_ls = np.array([])
    u_ls = []
    ur_ls = np.array([])
    if len(o0) > 1:
        for i in range(len(list_0)):
            qr_to_sum = np.array([])
            ur_to_sum = np.array([])
            p, p_r, q, q_r, u, u_r, theta, theta_r = vpol(o0[i], o1[i], o2[i], o3[i], e0[i], e1[i], e2[i], e3[i],
                                                          o0_r[i], e0_r[i], o1_r[i], e1_r[i], o2_r[i], e2_r[i], o3_r[i],
                                                          e3_r[i], wl_v)
            p_ls.append(p)
            q_ls.append(q)
            u_ls.append(u)
            pr_ls = np.append(pr_ls, p)
            qr_ls = np.append(qr_ls, q)
            ur_ls = np.append(ur_ls, u)

            qr_to_sum = np.append(qr_to_sum, 1 / ((q_r) ** 2))
            ur_to_sum = np.append(ur_to_sum, 1 / ((u_r) ** 2))

        qavg = np.average(q_ls, weights=1 / (qr_ls ** 2))
        uavg = np.average(u_ls, weights=1 / (ur_ls ** 2))
        qavg_r = np.sqrt(1 / np.sum(qr_to_sum))
        uavg_r = np.sqrt(1 / np.sum(ur_to_sum))
        pavg = np.sqrt(qavg ** 2 + uavg ** 2)
        pavg_r = (1 / pavg) * np.sqrt((qavg * qavg_r) ** 2 + (uavg * uavg_r) ** 2)
        theta_v = (0.5 * m.atan2(uavg, qavg)) * 180 / m.pi
        if theta_v < 0:
            theta_v = 180 + theta_v
        theta_vr = (0.5 * np.sqrt(
                ((uavg_r / uavg) ** 2 + (qavg_r / qavg) ** 2) * (1 / (1 + (uavg / qavg) ** 2)) ** 2)) * 180 / m.pi

    elif len(o0) == 1:
        pavg, pavg_r, qavg, qavg_r, uavg, uavg_r, theta_v, theta_vr = vpol(o0, o1, o2, o3, e0, e1, e2, e3, o0_r, e0_r,
                                                                           o1_r, e1_r, o2_r, e2_r, o3_r, e3_r, wl_v)
        if theta_v < 0:
            theta_v = 180 + theta_v
    print("PAVG         PAVG_err")
    print(pavg, pavg_r)
    print("QAVG         QAVG_err")
    print(qavg, qavg_r)
    print("UAVG         UAVG_err")
    print(uavg, uavg_r)
    print("P.A          P.A_err")
    print(theta_v, theta_vr)
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
    output = input('What do you want to call the output file? ')

    for filename in os.listdir("."):
        # Putting flux in from each file in list.
        if "1D_c_" in filename and "err" not in filename:
            wl, f = np.loadtxt(filename, unpack=True, usecols=(0, 1))
            if i == 0:
                flux.append(wl)
                flux_err.append(wl)
            flux.append(f)
            i = i + 1
        # Putting ERROR on flux from each error file in error list.
        elif "1D_c_" in filename and "err" in filename:
            wl, f = np.loadtxt(filename, unpack=True, usecols=(0, 1))
            if i == 0:
                flux.append(wl)
                flux_err.append(wl)

            flux_err.append(f)
            i = i + 1

    for x in xrange(len(flux[:][0])):
        if x == 0:
            try:
                os.remove(output)
            except:
                print("kittens")
        sum_flux = 0
        error_sqrd = 0
        for i in xrange(len(flux)):
            if i == 0:
                wl = flux[i][x]
            else:
                sum_flux = sum_flux + flux[i][x]  # Summing flux in given wavelength bin
                error_sqrd = error_sqrd + flux_err[i][x] * flux_err[i][x]  # Adding the square of the errors...
        error = np.sqrt(error_sqrd)  # ... and taking the square-root => error of the bin

        # Writing out the wl, flux and error in each bin out in a text file.
        with open(output, "a") as f:
            f.write(str(wl) + ' ' + str(sum_flux) + ' ' + str(error) + '\n')

    return
