import matplotlib.pyplot as plt
import numpy as np
import FUSS as F
import interactive_graph as ig
from FUSS import stat as Fstat

def isp_emline(filename_pol, filename_spctr, wlmin=4400, cont2ranges = False):
    """
    This function finds isp from one emission line. Requires interactive_range.def_ranges()
    :param filename_pol: name of the file containing the polarisation data
    :param filename_spctr: name of the file containing the flux data
    :param wlmin: minimum wavelength cutoff for the data imported from the polarisation and flux spectrum files.
    Default is 4400.
    :param cont2ranges: Boolean. If the continuum is the be defined by 2 ranges of values on either side of the line,
    set to True. If False, then the user should indicate the continuum by just two points on either side of the line.
    Default is False.
    :return:
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
        print "With continuum defined by the ranges:"
        print "{0:.0f} - {1:.0f}".format(cont_ranges[0].x[0], cont_ranges[0].x[-1])
        print "\nWe find:"
        print "ISP: p = {0:.3f} +/- {1:.3f} | q = {2:.3f} +/- {3:.3f} | u = {4:.3f} +/- {5:.3f}" .format(pisp, pisp_r,qisp, qisp_r,uisp,uisp_r)
        print "Continuum: p = {0:.3f} +/- {1:.3f} | q = {2:.3f} +/- {3:.3f} | u = {4:.3f} +/- {5:.3f}" .format(pcont, pcont_r,qcont, qcont_r,ucont,ucont_r)

    return pol_isp, pol_cont
    

#pol_isp, pol_cont = isp_emline('dc_11hs_ep1.pol','dc_11hs_ep1.flx' )




