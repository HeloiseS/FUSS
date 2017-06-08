from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
import FUSS as F
import interactive_graph as ig
from FUSS import statistics as Fstat
import math as m
from scipy import special as special


def from_emline(filename_pol, filename_spctr, wlmin=4400, cont2ranges = False):
    """
    This function finds isp from one emission line. Requires interactive_range.def_ranges()
    :param filename_pol: name of the file containing the polarisation data
    :param filename_spctr: name of the file containing the flux data
    :param wlmin: minimum wavelength cutoff for the data imported from the polarisation and flux spectrum files.
    Default is 4400.
    :param cont2ranges: Boolean. If the continuum is the be defined by 2 ranges of values on either side of the line,
    set to True. If False, then the user should indicate the continuum by just two points on either side of the line.
    Default is False.
    :return: pol_isp, pol_cont
    """
    # importing the data
    flux = F.get_spctr(filename_spctr, wlmin=wlmin, scale = False, err = True)
    pol = F.PolData('pol', filename_pol , wlmin=wlmin )

    scale = np.median(flux[1])  # scale factor used for plotting later

    # Need to define figure and plot the spectrum before calling ig.def_ranges()
    # not calling plot.show() though because that is done in igr.def_ranges()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(flux[0], flux[1])

    cont_ranges = ig.def_ranges(fig, flux, err=True)

    if cont2ranges is True:
        ###################################################################
        # Defining continuum (should need only 2 ranges so only considers #
        # the first 2 ranges defined with def_range)                      #
        ###################################################################
        cont_ranges[0].average()
        cont_ranges[1].average()

    ################################################################
    # Defining emission line region. Only 1 range defined by user  #
    ################################################################

    # need to plot again otherwise the mouse click function does not work.
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(flux[0], flux[1])  # plotting flux spectrum

    if cont2ranges is True:
        ax.plot(cont_ranges[0].x, cont_ranges[0].y, lw=2)  # plotting the first range
        ax.plot(cont_ranges[1].x, cont_ranges[1].y, lw=2)  # plotting the second range
        # plotting the line defining continuum (according to the 2 ranges picked)
        ax.plot([cont_ranges[0].middle,cont_ranges[1].middle], [cont_ranges[0].avg,cont_ranges[1].avg], lw=2, ls='--')
    else:
        ax.scatter([cont_ranges[0].x[0],cont_ranges[0].x[-1]], [cont_ranges[0].y[0], cont_ranges[0].y[-1] ], marker='o', c='r')  # plotting the first range
        ax.plot([cont_ranges[0].x[0],cont_ranges[0].x[-1]], [cont_ranges[0].y[0], cont_ranges[0].y[-1] ], marker='o', c='r')  # plotting the first range

    # plotting q and u just to see depolarisation regions. Scaled to fit on graph
    ax.plot(pol.wlp, pol.q*scale)
    ax.plot(pol.wlp, pol.u*scale)

    emission_range = ig.def_ranges(fig, flux, err = True)

    start=emission_range[0].start
    end=emission_range[0].end

    if cont2ranges is True:
        # To find the continuum flux we just interpolate between the averages of the first and second continuum ranges
        Fcont = np.interp(emission_range[0].x, [cont_ranges[0].middle, cont_ranges[1].middle], [cont_ranges[0].avg,cont_ranges[1].avg])
    else:
        Fcont = np.interp(emission_range[0].x, [cont_ranges[0].x[0],cont_ranges[0].x[-1]], [cont_ranges[0].y[0], cont_ranges[0].y[-1]])

    # Total flux of emission line is just array of all values of flux at each wavelength bin
    Ftot = emission_range[0].y
    Ftot_r = emission_range[0].yr
    # Line flux is total flux - continuum flux at each wavelength bin
    Fline = Ftot-Fcont
    Fline_r = np.array(emission_range[0].yr)

    # interpolating values of stokes parameters to match the wavelength bins of the flux so we can do
    # operations with all of these quantities.
    qtot = np.interp(emission_range[0].x, pol.wlp, pol.q)
    qtot_r = np.interp(emission_range[0].x, pol.wlp, pol.qr)
    utot = np.interp(emission_range[0].x, pol.wlp, pol.u)
    utot_r = np.interp(emission_range[0].x, pol.wlp, pol.ur)

    # qtot*Ftot/Fcont = yq (similar equation for u)
    # Fline/Fcont = x
    yq = (qtot * Ftot)/Fcont
    yqr = yq*np.sqrt( (qtot_r/qtot)**2 + (Ftot_r/Ftot)**2)
    yu = (utot * Ftot)/Fcont
    yur = yu*np.sqrt( (utot_r/utot)**2 + (Ftot_r/Ftot)**2)
    x = Fline/Fcont
    xr = Fline_r/Fcont

    qisp, qisp_r, qcont, qcont_r = Fstat.odr_fit(x, xr, yq, yqr)
    uisp, uisp_r, ucont, ucont_r = Fstat.odr_fit(x, xr, yu, yur)

    qfit = x*qisp + qcont
    ufit = x*uisp + ucont

    plt.errorbar(x, yq, xerr=xr, yerr=yqr)
    plt.errorbar(x, yu, xerr=xr, yerr=yur)
    plt.plot(x, qfit)
    plt.plot(x, ufit)
    plt.show()

    pisp = np.sqrt(qisp**2 + uisp**2)
    pisp_r = (1/pisp)*np.sqrt((qisp*qisp_r)**2 + (uisp*uisp_r)**2 )
    pol_isp = [pisp, pisp_r, qisp, qisp_r, uisp, uisp_r]

    pcont = np.sqrt(qcont**2 + ucont**2)
    pcont_r = (1/pcont)*np.sqrt((qcont*qcont_r)**2 + (ucont*ucont_r)**2 )
    pol_cont = [pcont, pcont_r, qcont, qcont_r, ucont, ucont_r]
    
    
    emline_wl = [(start+end)/2, end-((start+end)/2)]
    if cont2ranges is True:
        print "-------------------------- ISP from emission line ----------------------"
        print "For the emission line in range {0:.0f} - {1:.0f} Ang".format(start, end)
        print "With continuum defined by the ranges:"
        print "{0:.0f} - {1:.0f} | center: {2:.1f}".format(min(cont_ranges[0].x), max(cont_ranges[0].x), cont_ranges[0].middle)
        print "{0:.0f} - {1:.0f} | center: {2:.1f}".format(min(cont_ranges[1].x), max(cont_ranges[1].x), cont_ranges[1].middle)
        print "\nWe find:"
        print "ISP: p = {0:.3f} +/- {1:.3f} | q = {2:.3f} +/- {3:.3f} | u = {4:.3f} +/- {5:.3f}" .format(pisp, pisp_r,qisp, qisp_r,uisp,uisp_r)
        print "Continuum: p = {0:.3f} +/- {1:.3f} | q = {2:.3f} +/- {3:.3f} | u = {4:.3f} +/- {5:.3f}" .format(pcont, pcont_r,qcont, qcont_r,ucont,ucont_r)
    else:
        print "-------------------------- ISP from emission line ----------------------"
        print "For the emission line in range {0:.0f} - {1:.0f} Ang".format(start, end)
        print "With continuum defined by the points at:"
        print "{0:.0f} and {1:.0f}".format(cont_ranges[0].x[0], cont_ranges[0].x[-1])
        print "\nWe find:"
        print "ISP: p = {0:.3f} +/- {1:.3f} | q = {2:.3f} +/- {3:.3f} | u = {4:.3f} +/- {5:.3f}" .format(pisp, pisp_r,qisp, qisp_r,uisp,uisp_r)
        print "Continuum: p = {0:.3f} +/- {1:.3f} | q = {2:.3f} +/- {3:.3f} | u = {4:.3f} +/- {5:.3f}" .format(pcont, pcont_r,qcont, qcont_r,ucont,ucont_r)

    return emline_wl, pol_isp, pol_cont


def from_range(filename_pol, wlmin=None, wlmax=None):
    """
    Estimates ISP from polarisation within a range either defined from parameters or interactively.

    Notes
    -----
    If wlmin and wlmax are not given a plot will be displayed for the user to indicate the location of the range.

    Parameters
    ----------
    filename_pol : string
        Name of the text file were the polarisation data is located.
    wlmin : int
        Start of wavelength range. Default is None.
    wlmax : int
        End of wavelength range. Default is None.

    Returns
    -------
    tuple of floats
        pisp, pispr, qisp, qispr, uisp, uispr
    """
    pol = F.PolData('pol', filename_pol , wlmin=3500 )
    ls = [pol.q, pol.qr, pol.u, pol.ur]
    crop = []
    cond = (pol.wlp > wlmin) & (pol.wlp < wlmax)
    if wlmin is not None:
        for val in ls:
            valn = val[cond]
            crop.append(valn)
    else:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(pol.wlp, pol.q)  # plotting flux spectrum
        ax.plot(pol.wlp, pol.u)  # plotting flux spectrum
        isp_range = ig.def_ranges(fig, [pol.wlp,pol.q], err = False)
        for val in ls:
            valn = val[cond]
            crop.append(valn)

    # Values of p, q, u, a and their error for ISP
    qisp = np.average(crop[0], weights=1 / (crop[1] ** 2))
    qispr = np.std(crop[0])
    uisp = np.average(crop[2], weights=1 / (crop[3] ** 2))
    uispr = np.std(crop[2])
    pisp = np.sqrt(qisp ** 2 + uisp ** 2)
    pispr = (1 / pisp) * np.sqrt((qisp * qispr) ** 2 + (uisp * uispr) ** 2)
    aisp = (0.5 * m.atan2(uisp, qisp)) * 180.0 / m.pi
    aispr = 0.5 * np.sqrt(((uispr / uisp) ** 2 + (qispr / qisp) ** 2) * (
    1 / (1 + (uisp / qisp) ** 2)) ** 2)

    if aisp < 0:
        aisp = 180 + aisp  # Making sure P.A range is 0-180 deg

    if wlmin is None:
        print "Range: {0:.0f} - {1:.0f}".format(isp_range[0].start, isp_range[0].end)
    else:
        print "Range: {0:.0f} - {1:.0f}".format(wlmin, wlmax)

    print "ISP found: \n qisp = " + str(qisp) + " +/- " + str(qispr) \
          + "\n usip = " + str(uisp) + " +/- " + str(uispr) \
          + "\n pisp = " + str(pisp) + " +/- " + str(pispr) \
          + "\n P.A isp = " + str(aisp) + " +/- " + str(aispr)

    return pisp, pispr, qisp, qispr, uisp, uispr


def debias_p(p, pr, q=None, qr=None, u=None, ur=None, bayesian_pcorr = True, p0_step = 0.01):
    #  If bayesian_pcorr is False, P will be debiased as in Wang et al. 1997 using a step function
    if bayesian_pcorr is False:
        print "Step Func - p correction"
        pfinal = np.array([])
        for ind in range(len(p)):
            condition = p[ind] - pr[ind]
            if condition > 0:
                p_0i = p[ind]-((float(pr[ind]**2))/float(p[ind]))
            elif condition < 0:
                p_0i = p[ind]

            pfinal = np.append(pfinal, p_0i)
        return pfinal

    #  If bayesian_pcorr is True, P will be debiased using the Bayesian method described by J. L. Quinn 2012
    #  the correceted p is pbar_{0,mean} * sigma. pbar_{0,mean} is given by equation 47 of J. L. Quinn 2012
    if bayesian_pcorr is True:
        print "Bayesian - p correction"
        sigma = (qr + ur)/2
        pbar = p/sigma
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
        return pfinal


def linear_isp(wlp, gradq, constq, gradu, constu, q=None, qr=None, u=None, ur=None, bayesian_pcorr=False, p0_step = 0.01):
    """
    Calculates a linear and can also remove it from polarisation data if provided
    :param wlp: 1D Array of wavelength bins of final desired isp (often the wavelength bins fo your pol data)
    :param gradq: [gradient of q isp, error on gradient]
    :param constq: [intercept of q isp, error on intercept]
    :param gradu: [gradient of u isp, error on gradient]
    :param constu: [intercept of u isp, error on intercept]
    :param q: Stokes q of pol data you want to correct for isp. Default is None.
    :param qr: Err on Stokes q of pol data you want to correct for isp. Default is None.
    :param u: Stokes u of pol data you want to correct for isp. Default is None.
    :param ur: Err on Stokes u of pol data you want to correct for isp. Default is None.
    :param bayesian_pcorr: Boolean. If the data is being corrected for ISP and p debiasing should be done using
    the bayesian aproach set to True. Default is False.
    :return: If polarisation data not provided, only returns the array isp = [qisp, qisp_r, uisp, uisp_r],
    if data provided returns new_stokes = [wlp, pisp, pispr error, qisp, qisp error, uisp, uisp error,
    pol angle, pol angle error], and then isp = [qisp, qisp_r, uisp, uisp_r].
    """

    qisp = np.array([])
    qisp_r = np.array([])
    uisp = np.array([])
    uisp_r = np.array([])

    newq = np.array([])
    newqr = np.array([])
    newu = np.array([])
    newur = np.array([])

    for wl in wlp:
        qisp = np.append(qisp, gradq[0]*wl+constq[0])
        uisp = np.append(uisp, gradu[0]*wl+constu[0])
        qisp_r = np.append(qisp_r, np.sqrt((gradq[1]*wl)**2 + constq[1]**2))
        uisp_r = np.append(uisp_r, np.sqrt((gradu[1]*wl)**2 + constu[1]**2))

    isp = [qisp, qisp_r, uisp, uisp_r]

    if q is None:
        return isp

    for i in range(len(wlp)):
        newq = np.append(newq, q[i]-qisp[i])
        newqr = np.append(newqr, np.sqrt(qr[i]**2+qisp_r[i]**2))
        newu = np.append(newu, u[i]-uisp[i])
        newur = np.append(newur, np.sqrt(ur[i]**2+uisp_r[i]**2))

    newp = np.sqrt(newq**2 + newu**2)
    newpr = (1 / newp) * np.sqrt((newq * newqr) ** 2 + (newu * newur) ** 2)

    newa = np.array([])
    newar = np.array([])
    for i in range(len(wlp)):
        thetai = 0.5 * m.atan2(newu[i], newq[i])
        thetai_r = 0.5 * np.sqrt(((newur[i] / newu[i]) ** 2 + (newqr[i] / newq[i]) ** 2) * (
        1 / (1 + (newu[i] / newq[i]) ** 2)) ** 2)
        thetai = (thetai * 180.0) / m.pi
        thetai_r = (thetai_r * 180.0) / m.pi
        if thetai < 0:
            thetai = 180 + thetai

        newa = np.append(newa, thetai)
        newar = np.append(newar, thetai_r)

    if bayesian_pcorr is False:
        newp_debias = debias_p(newp, newpr, bayesian_pcorr=False)
    elif bayesian_pcorr is True:
        newp_debias = debias_p(newp, newpr, newq, newqr, newu, newur, bayesian_pcorr=True, p0_step=p0_step)

    new_stokes = [wlp, newp_debias, newpr, newq, newqr, newu, newur, newa, newar]

    return new_stokes, isp

def const_isp(wlp, qisp, qispr, uisp, uispr, q, qr, u, ur, bayesian_pcorr=True, p0_step=0.01):
    """
    Removes single valued (constant with wavelength) isp from data
    :param wlp: wavelength bins of the data
    :param qisp: stokes params of the isp
    :param qispr:
    :param uisp:
    :param uispr:
    :param q: stokes params of the data
    :param qr:
    :param u:
    :param ur:
    :return: new_stokes = [wlp, new p, new p error, new q, new q error, new u, new u error,
    new pol angle, new pol angle error]
    """
    newq = q - qisp
    newu = u - uisp
    newqr = np.sqrt(qr ** 2 + qispr ** 2)
    newur = np.sqrt(ur ** 2 + uispr ** 2)

    newp = np.sqrt(newq ** 2 + newu ** 2)
    newpr = (1 / newp) * np.sqrt((newq * newqr) ** 2 + (newu * newur) ** 2)

    newa = np.array([])
    newar = np.array([])
    for i in range(len(wlp)):
        thetai = 0.5 * m.atan2(newu[i], newq[i])
        thetai_r = 0.5 * np.sqrt(((newur[i] / newu[i]) ** 2 + (newqr[i] / newq[i]) ** 2) * (
        1 / (1 + (newu[i] / newq[i]) ** 2)) ** 2)
        thetai = (thetai * 180.0) / m.pi
        thetai_r = (thetai_r * 180.0) / m.pi
        if thetai < 0:
            thetai = 180 + thetai

        newa = np.append(newa, thetai)
        newar = np.append(newar, thetai_r)

    if bayesian_pcorr is False:
        newp_debias = debias_p(newp, newpr, bayesian_pcorr=False)
    elif bayesian_pcorr is True:
        newp_debias = debias_p(newp, newpr, newq, newqr, newu, newur, bayesian_pcorr=True, p0_step=p0_step)

    new_stokes =[wlp, newp_debias, newpr, newq, newqr, newu, newur, newa, newar]

    return new_stokes
