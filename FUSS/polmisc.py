# TODO: Do I need to update this dosctring to have something show up in the documentation?
"""
09 - May - 2017 / H. F. Stevance / fstevance1@sheffield.ac.uk

This is the main module of FUSS. It contains general utility functions, a couple of interactive routines and
also defines a new class: PolData, to deal with specpol data.
All this should make dealing with and analysing specpol data easier.

New: It is all (*) loaded by __init__.py

Functions:
----------
get_spctr(): Gets flux data from text file.
get_pol(): Gets pol data from text file.
dopcor(): Doppler Correction.
dopcor_file(): Doppler correction from data from a file output into a new file
ylim_def(): Used to define y limits for plots. Used within FUSS.
rot_data(): To rotate 2D data.
norm_ellipse(): Creates random data where the x and y coordinates are described by 2 different normal distributions.


Interactive routines:
---------------------
ep_date(): Taking a date as reference point, finds epoch from date or date from epoch.
vel(): Finds expansion velocity of element from observed and rest wavelength.

Class PolData():
----------------
Attributes:
    Defined by __init__
    - name: name
    - wlp = wavelength bins of polarisation data
    - p = p
    - pr = Delta p
    - q = q
    - qr = Delta q
    - u = u
    - ur = Delta u
    - a = Polarisation Angle P.A
    - ar = Delta P.A
    - wlf = wavelength bins of flux spectrum
    - f = Flux
    - fr = Delta F

    Defined by find_isp() or add_isp()
    - qisp, qispr, uisp, uispr, aisp, aispr: Stokes parameters and P.A of ISP

    Defined by rmv_isp()
    - p0, p0r, q0, ... , a0r : Original polarisation data before ISP correction
    - Updates p, pr, q, ..., ar with ISP corrected values.

Methods:
    - add_flux_data()
    - flu_n_pol()
    - find_isp()
    - add_isp()
    - rmv_isp()
    - qu_plt()

"""

import numpy as np
import matplotlib.pyplot as plt
import math as m
from scipy.stats import kde
import matplotlib.gridspec as gridspec
from scipy.odr import ODR, Model, Data, RealData, odr, Output
import os
from astropy.io import fits
import datetime as dt
from FUSS import isp as isp


# ################## FUNCTIONS ###################### FUNCTIONS #################### FUNCTIONS ################# #


def get_spctr(filename, wlmin=0, wlmax=100000, err=False, scale=True):
    """
    Imports spectrum.

    Notes
    -----
    Required file format: wl(Angstrom) flux *flux_error* (*optional*)

    Parameters
    ----------
    filename : string
        Name of the ASCII file where the spectrum is.
    wlmin : int, optional
        Lower wavelength cutoff. Default = 0.
    wlmax : int, optional
        Upper wavelength cutoff. Default = 100000.
    err : bool, optional
        If there is an error column, set to True. Default is False.
    scale : bool, optional
        Default is True. Multiplies the spectrum (and error) by the median values of the flux.

    Returns
    -------
    Tuple of 1D Arrays
        => Wavelength, Flux, *flux_error* (optional)

    """

    if err is False:
        flux = np.loadtxt(filename, unpack=True, usecols=(0, 1))
        cond = (flux[0] > wlmin) & (flux[0] < wlmax)
        wl = flux[0][cond]
        f = flux[1][cond]
        if scale is True:
            s = 1 / np.median(f)  # normalising the spectrum
            f = f * s
        return wl, f

    else:
        flux = np.loadtxt(filename, unpack=True, usecols=(0, 1))
        cond = (flux[0] > wlmin) & (flux[0] < wlmax)
        wl = flux[0][cond]
        f = flux[1][cond]
        r = flux[2][cond]
        if scale is True:
            s = 1 / np.median(f)
            f = f * s
            r = r * s
        return wl, f, r


def get_pol(filename, wlmin=0, wlmax=100000):
    """
    Imports values from polarisation files (given by specpol routine).

    Notes
    -----
    Required File format: 9 columns.
    First column must be wavelength in Angstrom.
    The other 8 columns are for stokes parameters, degree of pol and P.A, and associated errors:
        => wl p p_err q q_err u u_err angle angle_err

    Parameters
    ----------
    filename : string
        Name of the ASCII file.
    wlmin : int, optional
        Lower wavelength cutoff. Default = 0.
    wlmax : int, optional
        Upper wavelength cutoff. Default = 100000.

    Returns
    -------
    Tuple of 1D Arrays
        One 1 D array per parameter (so first must be wavelength, order of the rest depends on input file).
        => 9 arrays total.

    """

    pol0 = np.loadtxt(filename, unpack=True, usecols=(0, 1, 2, 3, 4, 5, 6, 7, 8))
    pol = []
    cond = (pol0[0] > wlmin) & (pol0[0] < wlmax)  # pol0[0] should contain the wavelength bins
    for val in pol0:
        # Applies the limits determined by wlmin, wlmax
        valn = val[cond]
        pol.append(valn)

    return pol[0], pol[1], pol[2], pol[3], pol[4], pol[5], pol[6], pol[7], pol[8]


def dopcor(val, z):
    # TODO: update dosctring
    # TODO: Create test
    """
    Doppler Correction.
    :param val: Array containing the data. val[0] MUST BE THE WAVELENGTH. NEED AT LEAST 2 COLUMNS!!
    :param z: Redshift
    :return: Array containing the data with the wavelength column doppler corrected.
    """

    values = np.array(val)  # need this in case val is not an array but a list
    wl0 = values[0]
    wln = np.array([])
    for wl in wl0:
        wl_dopcor = (wl) - (wl * z)
        wln = np.append(wln, wl_dopcor)
    values[0] = wln
    return values


def dopcor_file(filename, z):
    # TODO: update dosctring
    # TODO: Create test
    """
    Doppler Correction of data from a file (filename), into another file (output)
    :param filename: File containing the data to be Doppler corrected, the firsty column MUST contain the wavelength bins
    :param z: Redshift
    :return: 
    """
    output = 'dc_' + filename
    os.system('cp -i ' + filename + ' ' + output)
    f = file(output, 'r+')

    dopcor = []
    for line in f:
        columns = line.split()
        wl = float(columns[0])
        wl_dopcor = (wl) - (wl * z)
        dopcor.append(wl_dopcor)
    f.close()

    f0 = file(filename, 'r')
    f = file(output, 'w')
    i = 0
    for line in f0:
        columns = line.split()
        n_line = line.replace(columns[0], str(dopcor[i]))
        f.write(n_line)
        i = i + 1

    print output + ' created'


def ylim_def(wl, f, wlmin=4500, wlmax=9500):
    # TODO: Do I even use this?? Can I make it less convoluted?
    '''
    YLIM_DEF finds appropriate y limits for a spectrum. Look at values between a given range (Default: 4500-9500A) where
    we don't expect few order of magnitudes discrepancies like we see sometimes at the extremeties of the spectrum, then
    find the max and min value then define ymax and ymin.
    :param wl: Array containing the wavelengths
    :param f: Array containing the flux
    :param wlmin: Start of wavelength range where we search for min and max.
    :param wlmax: End of wavelength range where we search for min and max.
    :return: ymin (defined as flux min - 2*flux min) and ymax (defined as flux max + flux max/10)
    '''

    fmax = -100000
    fmin = 1000
    for i in xrange(len(wl)):
        if wl[i] < wlmax and wl[i] > wlmin:
            if f[i] < fmin:
                fmin = f[i]
                print fmin
            elif f[i] > fmax:
                fmax = f[i]
                print fmax

    # These tweaks to make the y limit okay were determined through testing. May not always
    # be appropriate and might need fixing later.
    if fmin > 0 and fmin < 1:
        ymin = fmin - 1.2 * fmin
    elif fmin > 0 and fmin > 1:
        ymin = fmin - fmin / 5
    elif fmin < 0 and fmin > -1:
        ymin = fmin + 1.2 * fmin
    elif fmin < 0 and fmin < -1:
        ymin = fmin + fmin / 5

    if fmax > 0 and fmax < 1:
        ymax = fmax + 1.2 * fmax
    elif fmax > 0 and fmax > 1:
        ymax = fmax + fmax / 5
    elif fmax < 0 and fmax > -1:
        ymax = fmax - 1.2 * fmax
    elif fmax < 0 and fmin < -1:
        ymax = fmax - fmax / 10

    return ymin, ymax


def rot_data(q, u, theta):
    # TODO: update dosctring
    # TODO: Create test
    """
    Rotates data
    :param q: X coordinate
    :param u: Y coordinate
    :param theta: Angle to rotate data by in rad.
    :return:
    """
    rot_matrix = np.array([[np.cos(theta), -np.sin(theta)],
                           [np.sin(theta), np.cos(theta)]])

    q_rot = np.array([])
    u_rot = np.array([])

    # Applying rotation to all bins and storing result in q_rot and u_rot
    for i in range(len(u)):
        coor = np.array([[q[i]],
                         [u[i]]])
        new_coor_i = np.dot(rot_matrix, coor)
        q_rot = np.append(q_rot, new_coor_i[0])
        u_rot = np.append(u_rot, new_coor_i[1])

    return q_rot, u_rot


def norm_ellipse(xc, yc, a, b, theta, n):
    # TODO: update dosctring
    # TODO: Create test
    """
    Creates ellipsoidal data set normally distributed around (xc,yc)
    :param xc: x coordinate of center of ellipse
    :param yc: y coordinate of center of ellipse
    :param a: standard deviation of Gaussian that describes major axis.
    :param b: standard deviation of Gaussian that describes minor axis
    :param theta: Rotation angle in radians
    :param n: Number of points created.
    :return: Arrays x and y, containing the x and y values of the data created.
    """
    i = 0
    x = np.array([])
    y = np.array([])

    # This creates data within ellipse. The x an y coordinates are defined by normal distribution.
    # That means we get normally distributed points in 2D, also means the ellipse's major and minor axis
    # are aligned with x and y axis or vice versa. So also give possibility to rotate the data set created
    while i <= n:
        x = np.append(x, np.random.normal(xc, a))
        y = np.append(y, np.random.normal(yc, b))
        i += 1
    if theta != 0:
        x, y = rot_data(x, y, theta)  # Applying rotation
    return x, y


def ep_date():
    # TODO: update dosctring
    """
    Interactive Routine. Finds epoch from date or date from epoch given a maximum date.
    :return:
    """

    # ####### Functions used by ep_date ########## #
    def date_input():
        yr = raw_input("Year: ")
        month = raw_input("Month: ")
        day = raw_input("Day: ")
        date = dt.date(int(yr), int(month), int(day))
        return date

    def date_from_epoch():
        ep = dt.timedelta(float(raw_input("\n What epoch (in days) would you like to know the date for: ")))
        print '\nDate at epoch ' + str(ep) + ' days: '
        print vmax + ep
        return vmax + ep

    def ep_from_dates():
        print "\nDate of epoch you want in days"
        date_ep = date_input()
        ep = date_ep - vmax
        print '\nEpoch:'
        print ep
        return ep

    # ############### MAIN ##################### #
    print "\nDate at V-band max"
    vmax = date_input()

    print "\n What do you want to do? \n (1) Get epoch in days. Inputs: Date of epoch" \
          "\n (2) Get date for an epoch in days. Inputs: Epoch in days (can be negative)" \
          "\n (3) Update the V-band max date" \
          "\n (4) Exit"

    to_do = raw_input("#> ")
    while to_do != '4':
        if to_do == '1':
            ep_from_dates()
        if to_do == '2':
            date_from_epoch()
        if to_do == '3':
            print "\nDate at V-band max"
            vmax = date_input()
        if to_do != '1' and to_do != '2' and to_do != '3' and to_do != '4':
            print "Must choose option 1, 2, 3 or 4"

        to_do = raw_input("#> ")

    return "Good Bye"


def vel():
    # TODO: update dosctring
    """
    Interactive routine. Finds the velocity for a given observed wavelength and rest wavelength.
    :return:
    """
    cont = 'y'
    while cont == 'y' or cont == '':
        l_obs = float(raw_input('What is the observed wavelength: '))
        l_emit = float(raw_input('What is the rest wavelength: '))

        c = 299792.458  # Speed of light in km/s

        v = ((l_obs - l_emit) / l_emit) * c
        print v
        cont = raw_input('Continue?(y/n): ')


#  #################################################################################  #
#  ##############  CLASSE ############## POLDATA ########### CLASSE  ###############  #
#  #################################################################################  #

class PolData(object):
    """
    Each instance contains one spectropolarimetric data set.

    Note
    -----
        The attributes wlp, p, pr, q, qr, u, ur, a and ar are 1D arrays and must have the
    same length.
        The attributes wlf, f and fr must have the same length, but it can differ from the
    length of wlp, p, etc...
        When the ISP is removed, the attributes p0, p0r, q0, etc... store the original values
    of p, pr, q, etc..., and the latter are updated to have the ISP corrected values of polarisation.

    Parameters
    ----------
    name : str
        A short handle to make your data object recognisable (e.g. 'ep1', '14ad')
    poldata : str or tuple
        The polarisation data can be imported from a text file containing only the data, where
        the column order is: wavelength p p_err q q_err u u_err a a_err.
        Alternatively a tuple of arrays containing the data can be provided. Make sure the order
        of the arrays in the tuple corresponds to wavelength p p_err q q_err u u_err a a_err.
    wlmin : int, optional
        Minimum wavelength cutoff
    wlmax : int, optional
        Maximum wavelength cutoff

    Attributes
    ----------
    name : str
        A short handle to make your data object recognisable (e.g. 'ep1', '14ad')
    wlp : array
        1D array containing the wavelength bins of the polarisation data.
    p : array
        1D array containing the degree of polarisation in each bin.
    pr : array
        1D array containing the error on p in each bin.
    q : array
        1D array containing Stokes q in each bin.
    qr : array
        1D array containing the error on q in each bin.
    u : array
        1D array containing Stokes u in each bin.
    ur : array
        1D array containing the error on u in each bin.
    a : array
        1D array containing the polarisation angle in each bin.
    ar : array
        1D array containing the error on the polarisation in each bin.
    wlf : array, optional
        1D array containing wavelength bins of the flux spectrum.
    f :  array, optional
        1D array containing the flux in each bin.
    fr : array, optional
        1D array containing the error on the flux in each bin.
    qisp : float, optional
        Stokes q of the ISP.
    qispr : float, optional
        Error on q ISP.
    uisp : float, optional
        Stokes u of the ISP
    uispr : float, optional
        Error on u ISP
    gradq : tuple, optional
        Gradient of Stokes q ISP  and error on the gradient.
    constq : tuple, optional
        Intercept of Stokes q ISP and error on the intercept.
    gradu : tuple, optional
        Gradient of Stokes u ISP  and error on the gradient.
    constu : tuple, optional
        Intercept of Stokes u ISP and error on the intercept.
    p0 : array
        1D array containing the degree of polarisation in each bin BEFORE ISP REMOVAL.
    p0r : array
        1D array containing the error on p in each bin BEFORE ISP REMOVAL.
    q0 : array
        1D array containing Stokes q in each bin BEFORE ISP REMOVAL.
    q0r : array
        1D array containing the error on q in each bin BEFORE ISP REMOVAL.
    u0 : array
        1D array containing Stokes u in each bin BEFORE ISP REMOVAL.
    u0r : array
        1D array containing the error on u in each bin BEFORE ISP REMOVAL.
    a0 : array
        1D array containing the polarisation angle in each bin BEFORE ISP REMOVAL.
    a0r : array
        1D array containing the error on the polarisation in each bin BEFORE ISP REMOVAL.
    """

    def __init__(self, name, poldata, wlmin=None, wlmax=1000000):

        if type(poldata) is str:
            pol0 = get_pol(poldata, wlmin=wlmin, wlmax=wlmax)
        else:
            pol0 = poldata

        self.name = name
        self.wlp = pol0[0]
        self.p = pol0[1]
        self.pr = pol0[2]
        self.q = pol0[3]
        self.qr = pol0[4]
        self.u = pol0[5]
        self.ur = pol0[6]
        self.a = pol0[7]
        self.ar = pol0[8]
        self.wlf = None
        self.f = None
        self.fr = None
        self.qisp = None
        self.qispr = None
        self.uisp = None
        self.uispr = None
        self.pisp = None
        self.pispr = None
        self.aisp = None
        self.aispr = None
        self.gradq = None
        self.constq = None
        self.gradu = None
        self.constu = None
        self.q0 = None
        self.u0 = None
        self.q0r = None
        self.u0r = None
        self.p0 = None
        self.p0r = None
        self.a0 = None
        self.a0r = None

        print " ==== PolData - instance: " + self.name + " ===="
        print "Polarisation data initialised. If you want to add Stokes I use add_flux_data(). " \
              "To find ISP use find_isp(). \n"

    def add_flux_data(self, filename, wlmin=None, wlmax=1000000, err=False):
        """
        Adds flux spectrum data attributes to the PolData.

        Parameters
        ----------
        filename : str
            File containing the flux data. File format: wl, f, fr (no comas)
        wlmin : int
            Minimum wavelength cut off
        wlmax :
            Maximum wavelength cut off
        err : bool
            If false, only imports wavelength and flux, not the error on the flux. Default = False.
        """
        flux = get_spctr(filename, wlmin=wlmin, wlmax=wlmax)
        self.wlf = flux[0]
        self.f = flux[1]
        if err is True:
            self.fr = flux[2]

        print " ==== PolData - instance: " + self.name + " ===="
        print "Flux spectrum added."

    def flu_n_pol(self, save=False):
        """
        Creates plot of p, q, u, theta, and flux.

        Note
        ----
        /!\ The x-axis is SHARED, so limits on polarisation attributes and flux
        attributes should be the same.

        Parameters
        ----------
        save : bool
            Whether to save the plot or not. Saved as [self.name]_fnp.png
        """

        fnp = plt.figure(figsize=(10, 10))
        grid = gridspec.GridSpec(5, 1, hspace=0)
        p_plot = plt.subplot(grid[0])
        q_plot = plt.subplot(grid[1])
        u_plot = plt.subplot(grid[2])
        a_plot = plt.subplot(grid[3])
        f_plot = plt.subplot(grid[4])
        p_plot.errorbar(self.wlp, self.p, yerr=self.pr, color='purple', capsize=0, ecolor='grey')
        q_plot.errorbar(self.wlp, self.q, yerr=self.qr, color='r', alpha=0.8, capsize=0, ecolor='grey')
        u_plot.errorbar(self.wlp, self.u, yerr=self.ur, color='blue', alpha=0.8, capsize=0, ecolor='grey')
        a_plot.errorbar(self.wlp, self.a, yerr=self.ar, color='orange', alpha=0.8, capsize=0, ecolor='grey')

        try:
            f_plot.errorbar(self.wlf, self.f, yerr=self.fr, color='k', alpha=0.5, lw=1.5, capsize=0, ecolor='grey')
        except:
            print 'Flux attributes not defined'

        p_plot.set_ylim(ylim_def(self.wlp, self.p, wlmin=4700))
        p_plot.set_ylabel('p (%)')
        p_plot.set_title(self.name, fontsize=16)

        q_plot.set_ylim(ylim_def(self.wlp, self.q, wlmin=4700))
        q_plot.set_ylabel('q (%)')

        u_plot.set_ylim(ylim_def(self.wlp, self.u, wlmin=4700))
        u_plot.set_ylabel('u (%)')

        a_plot.set_ylim(ylim_def(self.wlp, self.a, wlmin=4700))
        a_plot.set_ylabel('P.A (deg)')

        try:
            f_plot.set_ylim(ylim_def(self.wlf, self.f))
            f_plot.set_ylabel('Flux')
            f_plot.set_xlabel('Wavelength (Ang)', fontsize=14)
        except:
            print 'Flux attributes not defined'

        p_plot.xaxis.set_visible(False)
        q_plot.xaxis.set_visible(False)
        u_plot.xaxis.set_visible(False)
        a_plot.xaxis.set_visible(False)

        if save is True:
            fnp.savefig(self.name + '_fnp.png')
        plt.show()
        return

    def find_isp(self, wlmin, wlmax):
        """
        Estimates ISP

        Notes
        -----
        Simply an average of q and u over a given wavelength range which should correspond to line
        blanketting region.

        Parameters
        ----------
        wlmin : int
            Start of wavelength range.
        wlmax : int
            End of wavelength range.
        """

        ls = [self.q, self.qr, self.u, self.ur]
        cond = (self.wlp > wlmin) & (self.wlp < wlmax)
        crop = []
        for val in ls:
            valn = val[cond]
            crop.append(valn)

        # Values of p, q, u, a and their error for ISP
        self.qisp = np.average(crop[0], weights=1 / (crop[1] ** 2))
        self.qispr = np.std(crop[0])
        self.uisp = np.average(crop[2], weights=1 / (crop[3] ** 2))
        self.uispr = np.std(crop[2])
        self.pisp = np.sqrt(self.qisp ** 2 + self.uisp ** 2)
        self.pispr = (1 / self.pisp) * np.sqrt((self.qisp * self.qispr) ** 2 + (self.uisp * self.uispr) ** 2)
        self.aisp = (0.5 * m.atan2(self.uisp, self.qisp)) * 180.0 / m.pi
        self.aispr = 0.5 * np.sqrt(((self.uispr / self.uisp) ** 2 + (self.qispr / self.qisp) ** 2) * (
        1 / (1 + (self.uisp / self.qisp) ** 2)) ** 2)

        if self.aisp < 0:
            self.aisp = 180 + self.aisp  # Making sure P.A range is 0-180 deg

        print " ==== PolData - instance: " + self.name + " ===="
        print "ISP found: \n qisp = " + str(self.qisp) + " +/- " + str(self.qispr) \
              + "\n usip = " + str(self.uisp) + " +/- " + str(self.uispr) \
              + "\n pisp = " + str(self.pisp) + " +/- " + str(self.pispr) \
              + "\n P.A isp = " + str(self.aisp) + " +/- " + str(self.aispr)
        return self.qisp, self.qispr, self.uisp, self.uispr

    def add_isp(self, constisp_params = None, linearisp_params = None):
        # TODO: I need 2 tests for this function
        """
        Adds parameters of isp to the data.

        Parameters
        ----------
        constisp_params : tuple
            If the isp is constant give the stokes parameters of the isp here in a list:
        [qisp, qisp error, uisp , uisp error]
        linearisp_params : tuple
            Tuple of tuples: [[grad_q, grad_q error],[intercept_q, intercept_q error],
            [grad_u, grad_u error],[intercept_u, intercept_u error]].
            For qisp = grad_q * lambda + intercept_q (and similar equation for u), where lambda is in Angstrom.

        Examples
        --------
            If the ISP is constant across your wavelength range, put its values an associated errors in constisp_params:
            >> PolDataObj.add_isp(constisp_params=[0.14, 0.04, 0.08, 0.03])

            If the isp changes linearly with wavelength, give the parameters for the lines of q and u ISP here.
            >> PolDataObj.add_isp(linearisp_params=[[0.00035, 0.00003],[2.45, 0.19]])
        """

        if linearisp_params is None:
            self.qisp, self.qispr, self.uisp, self.uispr = constisp_params
            # Values of p, q, u, a and their error for ISP
            self.pisp = np.sqrt(self.qisp ** 2 + self.uisp ** 2)
            self.pispr = (1 / self.pisp) * np.sqrt((self.qisp * self.qispr) ** 2 + (self.uisp * self.uispr) ** 2)
            self.aisp = (0.5 * m.atan2(self.uisp, self.qisp)) * 180.0 / m.pi
            self.aispr = 0.5 * np.sqrt(((self.uispr / self.uisp) ** 2 + (self.qispr / self.qisp) ** 2) * (
            1 / (1 + (self.uisp / self.qisp) ** 2)) ** 2)
            self.aispr = (self.aispr * 180.0) / m.pi
            if self.aisp < 0:
                self.aisp = 180 + self.aisp  # Making sure P.A range is 0-180 deg

            print " ==== PolData - instance: " + self.name + " ===="
            print "ISP Added: \n qisp = " + str(self.qisp) + " +/- " + str(self.qispr) \
                  + "\n usip = " + str(self.uisp) + " +/- " + str(self.uispr) \
                  + "\n pisp = " + str(self.pisp) + " +/- " + str(self.pispr) \
                  + "\n P.A isp = " + str(self.aisp) + " +/- " + str(self.aispr) + "\n"
            self.gradq = None # this will be used as a condition for the method of isp removal in rmv_isp
        elif constisp_params is None:
            self.gradq, self.constq, self.gradu, self.constu = linearisp_params
            self.qisp = None # this will be used as a condition for the method of isp removal in rmv_isp
        return

    def rmv_isp(self, bayesian_pcorr=True, p0_step=0.01):
        # TODO: I need 2 tests for this. Maybe will need 14ad data for the constant case and 11hs for the linear case
        """
        Removes ISP and updates q, qr, u, ur, p, pr, a and ar.

        Note
        -----
        Stores the original non ISP corrected degree of polarisation, Stokes parameters, polarisation angle,
        and associated errors in p0, p0r, q0, q0r, u0, u0r, a0, and a0r, and updates p, pr, q, qr, u, ur, a and ar.
        """

        # Storing original values  of Stokes parameters and their errors in newly defined
        # attributes.
        self.q0 = self.q
        self.u0 = self.u
        self.q0r = self.qr
        self.u0r = self.ur
        # Storing original degree of polarisation and it's error in new variable and updating p and pr
        self.p0 = self.p
        self.p0r = self.pr
        # Same as before but for the P.A
        self.a0 = self.a
        self.a0r = self.ar

        if self.qisp is None:
            new_stokes, __ = isp.linear_isp(self.wlp, self.gradq, self.constq,
                                            self.gradu, self.constu,
                                            self.q, self.qr,
                                            self.u, self.ur,
                                            bayesian_pcorr=bayesian_pcorr, p0_step=p0_step)

        elif self.gradq is None:
            new_stokes = isp.const_isp(self.wlp, self.qisp, self.qispr,
                                       self.uisp, self.uispr,
                                       self.q, self.qr,
                                       self.u, self.ur,
                                       bayesian_pcorr=bayesian_pcorr, p0_step=p0_step)

        self.p = new_stokes[1]
        self.pr =new_stokes[2]
        self.q = new_stokes[3] # new_stokes[0] is wavelength bins
        self.qr = new_stokes[4]
        self.u = new_stokes[5]
        self.ur = new_stokes[6]
        self.a = new_stokes[7]
        self.ar = new_stokes[8]

    def qu_plt(self, subplot_loc=111, wlmin=None, wlmax=100000,
               qlim=[-3.0, 3.0], ulim=[-3.0, 3.0], textloc=[-2.7, -2.7], cisp='k', fs=16,
               ls=14, isp=False, wlrest=None, colorbar=True, colorbar_labelsize=14, size_clbar=0.05, line_color=None,
               marker='.', lambda_xshift=1.7, fit=True,
               qlab_vis=True, ulab_vis=True,
               qticks_vis=True, uticks_vis=True):
        # TODO: anyway to use *args here? how does that even work?
        """
        Plots the QU plane corresponding to the imported data.

        Parameters
        ----------
        subplot_loc : int or matplotlib.gridspec.GridSpec, optional
            Location of the subplot.  Can be a 3 digit integer or a gridspec location ifcreated a grid using gridspec.
            Default = 111.
        wlmin : int, optional
            Min wavelength cut off. Default None.
        wlmax : int, optional
            Max wavelength cut off. Default 100000.
        qlim : tuple, optional
            [min q, max q]. Default = [-3.0, 3.0]
        ulim : tuple, optional
            [min u, max u]. Default = [-3.0, 3.0]
        textloc : tuple, optional
            Location of name of qu-plot. Default = [-2.7, -2.7]
        cisp : string, optional
            Color of ISP marker. Default = 'k'
        fs : int, optional
            Font size. Applies to text on plot and axis labels, not graduations on the axes. Default = 16
        ls : int, optional
            Label size. Size of the tick numbers on axes. Default = 14.
        isp : bool, optional
            Whether to plot ISP. Default False.
        wlrest :int, optional
            If plotting qu plot of a line, rest wavelength of that line. Otherwise leave default value: None.
        colorbar : bool, optional
            Default is True. If False the colorbar is not plotted.
        colorbar_labelsize : int, optional
            Label size of the color bar ticks. Default 15.
        size_clbar : float, optional
            Modifies the size of the colour bar. Also screws with the plot somehow. Default = 0.05.
        line_color : string, optional
            If want a solid colour for the lines between the markers. Default is None and gives lines cycling through
            rainbow colors to match the color of the point they are associated with.
        marker : string, optional
            Type of marker to be used. Default is '.'
        lambda_xshift : float, optional
            Position of the colourbar label define as qmax + shift. This is the shift value. Default is 1.7.
        fit : bool, optional
            If False the dominant axis will not be plotted. Its parameters will still be calculated and returned.
            Default is True.
        qlab_vis : bool, optional
            If False, the q label is not plotted. Default is True.
        ulab_vis : bool, optional
            If False, the u label is not plotted. Default is True.
        qticks_vis : bool, optional
            If False, all q tick labels are invisible. Default is True.
        uticks_vis : bool, optional
            If False, all u tick labels are invisible. Default is True.

        Returns
        ------
            matplotlib.axes._subplots.AxesSubplot
            The axis the qu plane is plotted on. That way can plot other things on top, e.g line or ellipse or else.
        """

        # ###################       FITTING THE DATA WITH DOM AXIS         ########################### #
        func = lambda beta,x: beta[0] + beta[1] * x # Expression of the line that we want to fit to the data

        data = RealData(self.q, self.u, self.qr, self.ur)
        model = Model(func)
        odr = ODR(data, model, [0, 0])

        # Given the levels of pol in SNE, I don't expect to ever have to plot a q-u plot with limits [-10,10]
        # The following are just q values from -10 to 10 that will be used to plot the line fit
        q_n = np.arange(-10, 10, 0.1)

        qu = plt.subplot(subplot_loc, aspect='equal')

        odr.set_job(fit_type=0)  # fit_type = 0 => explicit ODR.
        output = odr.run()

        print " ==== QUplot - instance: " + self.name + " ===="
        print "Dom. Axis = a*x + b"
        print "a = " + str(output.beta[1]) + " +/- " + str(output.sd_beta[1])
        print "b = " + str(output.beta[0]) + " +/- " + str(output.sd_beta[0]) + "\n"

        u_n = func(output.beta, q_n)  # Based on fit, get the u values for each q

        if fit is True:
            qu.plot(q_n, u_n, 'k--', linewidth=2, zorder=1000)
            # the zorder is high to sit on top of the scatter created belox

        cond = (self.wlp > wlmin) & (self.wlp < wlmax)
        wl_crop = self.wlp[cond]
        q_crop = self.q[cond]
        qr_crop = self.qr[cond]
        u_crop = self.u[cond]
        ur_crop = self.ur[cond]

        # #################### CREATING THE PLOT ########################
        plt.set_cmap('plasma') # Plasma is perceptually uniform. Otherwise use 'jet'
        if wlrest is None:
            # Defining the min and max wavelength, which are going to be the beginning and end of the colour map
            wlmin = min(wl_crop)
            wlmax = max(wl_crop)
            sc = qu.scatter(q_crop, u_crop, s=100,
                            vmin=wlmin, vmax=wlmax,
                            c=wl_crop, marker=marker,
                            zorder=600, lw=0)

        else:
            vel = np.array([])
            c = 299792.0

            for i in xrange(len(wl_crop)):
                v = c * ((wl_crop[i] - wlrest) / wlrest)
                vel = np.append(vel, v)

            # Defining the min and max VELOCITIES, which are going to be the beginning and end of the colour map
            velmin = min(vel)
            velmax = max(vel)
            print velmin, velmax
            sc = qu.scatter(q_crop, u_crop, s=100,
                            vmin=velmin, vmax=velmax,
                            c=vel, marker=marker,
                            zorder=600, lw=0)


        # ################## Plotting Points ###############################
        # vmin and vmax are the start and end of the colour map. c = wl because we're defining the colourmap using the
        # wavelengths wl. zorder doesn't have to be 600, it just needs to be below that of the fitting line we did above
        # and greater than the zorder of the error bars, because otherwise it doesn't look nice.

        clbar = plt.colorbar(sc, fraction=size_clbar)  # Plotting to colour map. Need to do that to get a rainbow.
        clbar.ax.tick_params(labelsize=colorbar_labelsize)

        if colorbar is False:
            clbar.remove()  # Removing Colormap from plot (but still exists so we can plot rainbows)
        elif colorbar is True:
            if wlrest is None:
                qu.text(qlim[1] + lambda_xshift, (ulim[1] + ulim[0]) / 2, r'$\lambda (\AA)$', fontsize=fs)
            else:
                qu.text(qlim[1] + lambda_xshift, (ulim[1] + ulim[0]) / 2, 'Velocity (km/s)', rotation='vertical',
                        fontsize=fs)

        a, b, c = qu.errorbar(q_crop, u_crop, xerr=qr_crop, yerr=ur_crop, marker='.', capsize=0,
                              zorder=500, linestyle='None', alpha=0.4)  # Plotting error bars

        # Convert my wavelengths into the colour map plotted earlier applying the colourbar to "c",
        #  that is, the errorbars, there are 2 components (c[0] and c[1]) because I have error bars in both x and y.
        if wlrest is None:
            clmap = clbar.to_rgba(wl_crop)
        else:
            clmap = clbar.to_rgba(vel)
        c[0].set_color(clmap)
        c[1].set_color(clmap)

        # The following loop cycles through our colormap. Without this the lines we are about to create to connect
        # the points of the scatter plot will not have colours corresponding to the points they are linking.
        qu.set_prop_cycle(plt.cycler('color', clmap))
        for i in range(len(wl_crop) - 1):
            qu.plot(q_crop[i:i + 2], u_crop[i:i + 2], c=line_color,
                    alpha=1)  # Here we create line for each pair of points
            # Note that it's "i+2" in order for the last point to be i+1 -because it's up to point i+2, excluding i+2.

        # To mark ISP with errorbars
        if isp is True:
            plt.errorbar(self.qisp, self.uisp, xerr=self.qispr, yerr=self.uispr, fmt='o', color=cisp, elinewidth=2.5,
                         capthick=2.5, zorder=5000)

        plt.axvline(0, color='k', linestyle='-.')
        plt.axhline(0, color='k', linestyle='-.')

        qu.tick_params(axis='both', which='major', labelsize=ls)

        # Now fiddling with the ticks: If ticks are made to be visible then sent every other tick to be invisible
        # so bring so space to the axes. If ticks are set to be invisible... well make them invisible.
        xticks = qu.xaxis.get_major_ticks()
        yticks = qu.yaxis.get_major_ticks()

        ''' Didn't work to resize my tick labels :(
        for xtick in xticks:
            xtick.label1.set_fontsize(ticklabelsize)
            
        for ytick in yticks:
            ytick.label1.set_fontsize(ticklabelsize)
        '''
        if qticks_vis is False:
            for i in range(0, len(xticks)):
                xticks[i].label1.set_visible(False)
        else:
            for i in range(0, len(xticks), 2):
                xticks[i].label1.set_visible(False)

        if uticks_vis is False:
            for i in range(0, len(yticks)):
                yticks[i].label1.set_visible(False)
        else:
            for i in range(0, len(yticks), 2):
                yticks[i].label1.set_visible(False)

        if qlab_vis is True:
            qu.set_xlabel('q (%)', fontsize=fs)

        if ulab_vis is True:
            qu.set_ylabel('u (%)', labelpad=-1, fontsize=fs)

        qu.text(textloc[0], textloc[1], self.name, fontsize=fs)

        qu.set_xlim(qlim)  # Setting some limits.
        qu.set_ylim(ulim)
        return qu
