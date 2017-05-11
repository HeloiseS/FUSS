"""
20 -  March - 2017 / H. F. Stevance / fstevance1@sheffield.ac.uk

stat.py is a sub-module of the FUSS package. It does statistics stuff.

Pre-requisites:
---------------
numpy, math, matplotlib.patches, scipy.odr

Functions (by sections):
------------------------
For input details check docstrings of individual functions.

# ######### Covariance Matrix ######### #
cov_el(args): cov_el calculates and returns the value of 1 element of the covariance matrix
cov_mat(args): Puts the covariance matrix together and returns it - uses cov_el().

# ######### Principal Components Analysis ########## #
pca(args): Does PCA, requires cov_mat(). Finds the ellipse of best fit to the data. Returns the axis ratio,
the angle of the dominant axis and the angle of the orthogonal axis in DEGREES.

draw_ellipse(args): Draw an ellipse. Nothing to do with statistics but if you do PCA to find the ellipse of best
fit to your data you'll want to plot an ellipse.
Returns: ellipse. To use as input in add_artist(args) [a matplotlib.axes function]

# ########## Pearson Correlation ############# #
pearson(args): Returns pearson coefficient. Uses cov_el(args)

# ########## ODR FITS ################### #
func(args): lienar 2D function
odr_fit(args): fits 2D data with a line using Orthogonal Distance Regression 

"""
from __future__ import division
import numpy as np
import math as m
from matplotlib.patches import Ellipse
from scipy.odr import ODR, Model, Data, RealData, odr, Output

# ############################# Covariance Matrix #################### #


def cov_el(j, k, q, u, q_r, u_r):
    """
    cov_el provides value of 1 element of the covariance matrix
    :param j: col number
    :param k: row number
    /!\ q, u and their errors need to have the same dimensions
    :param q: q
    :param q_r: error on q
    :param u: u
    :param u_r: error on u
    :return: element
    """
    N = len(q)
    if j == 0:
        xj=q
        dxj=q_r
    elif j == 1:
        xj=u
        dxj=u_r

    if k == 0:
        xk=q
        dxk=q_r
    elif k == 1:
        xk=u
        dxk=u_r

    xj_avg = np.average(xj, weights=1/dxj**2)
    xk_avg = np.average(xk, weights=1/dxk**2)
    num = []  # numerator
    den = []  # denominator

    for i in range(len(q)):
        num_i = (1/dxj[i])*(1/dxk[i]) * (xj[i]-xj_avg)*(xk[i]-xk_avg)
        den_i = (1/dxj[i])*(1/dxk[i])
        num.append(num_i)
        den.append(den_i)

    element = (np.sum(num) * N) / (np.sum(den) * (N-1))

    return element


def cov_mat(q, u, q_r, u_r):
    """
    Creates the covariance matrix. Requires cov_el.
    :return: The covariance matrix
    """

    coor_ls = [[0,0], [0,1], [1, 0], [1,1]]
    els = [0,1,2,3]
    for coor in coor_ls:
        j = coor[0]
        k = coor[1]
        els[coor_ls.index(coor)] = cov_el(j, k, q, u, q_r, u_r)

    cov_matrix=np.array([[els[0],els[1]],[els[2],els[3]]])
    print "Covariance Matrix:"
    print cov_matrix,'\n'
    return cov_matrix
    
# ############################ Principal Components analyis ################################### #


def pca(q, u, q_r, u_r):
    """
    Does the Principal Components Analysis. Notation for variables is consistent with Stokes parameters but could be any
    x and y coordinates and their errors instead. All arrays provided must have the same dimensions.
    :param q: Array containing q values
    :param u: Array containing u values
    :param q_r: Array containing errors on q
    :param u_r: Array containing errors on u
    :return: b/a (axis ratio), rotation angel of major axis, rotation angle of minor axis. Angles given in degrees.
    """
    matrix = cov_mat(q, u, q_r, u_r)
    evalues, evectors = np.linalg.eig(matrix)
      
    if evalues[0] > evalues[1]:
        b_a = evalues[1]/evalues[0]
        alpha_dom = 180*m.atan(evectors[1,0]/evectors[1,1])/m.pi
        alpha_ort = 180*m.atan(evectors[0,0]/evectors[0,1])/m.pi
        
    else:
        b_a = evalues[0]/evalues[1]
        alpha_dom = 180*m.atan(evectors[0,0]/evectors[0,1])/m.pi
        alpha_ort = 180*m.atan(evectors[1,0]/evectors[1,1])/m.pi

    return b_a, alpha_dom, alpha_ort    


def draw_ellipse(q,u, a, ratio, alpha_dom ):
    """
    Draws an ellipse. Duh.
    :param q: Array containing q values
    :param u: Array containing u values
    :param a: Major axis length
    :param ratio: Axis ratio (<1)
    :param alpha_dom: Rotation major axis in degrees.
    :return: ellipse. to use as input in add_artist() [a matplotlib.axes function]
    """
    centre=[np.average(q),np.average(u)]
    ellipse = Ellipse(centre, a, ratio*a, alpha_dom, facecolor = 'none', lw=2, alpha=0.8, ls='-.',zorder=1000)
    return ellipse
    
# ##################################### Pearson Correlation ###################################### #


def pearson(q, u, qr, ur):
    """
    Perform a Pearson's test on the data provided. All arrays must have the same length.
    :param q: Array containing q values
    :param u: Array containing u values
    :param q_r: Array containing errors on q
    :param u_r: Array containing errors on u
    :return: Pearson's coefficient
    """
    cov = cov_el(0,1,q, u, qr, ur)
    sq = np.sqrt(cov_el(0,0,q, u, qr, ur))
    su = np.sqrt(cov_el(1,1,q, u, qr, ur))
    return cov/(su*sq)
    
# ##################################### ODR FITS ########################################### #


def func(beta, x):
    """
    Just a linear function. Used by odr_fit but can use on its own.
    :param beta: Array [intercept, gradient]
    :param x: independent variable (generally x)
    """
    # Expression of the line that we want to fit to the data
    y = beta[0] + beta[1] * x
    return y


def odr_fit(x, xr, y, yr):
    """
    Performs an ODR fit on the data provided.
    :param x:
    :param xr:
    :param y:
    :param yr:
    :return:
    """

    data = RealData(x, y, xr, yr)
    model = Model(func)
    odr = ODR(data, model, [0, 0])

    odr.set_job(fit_type=0)  # fit_type = 0 => explicit ODR.
    output = odr.run()

    print "Line = a*x + b"
    print "a = " + str(output.beta[1]) + " +/- " + str(output.sd_beta[1])
    print "b = " + str(output.beta[0]) + " +/- " + str(output.sd_beta[0]) + "\n"

    grad = output.beta[1]
    grad_r = output.sd_beta[1]
    intercept = output.beta[0]
    intercept_r = output.sd_beta[0]
    
    return grad, grad_r, intercept, intercept_r




