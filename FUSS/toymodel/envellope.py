"""
TOYMODEL - fstevance1@sheffield.ac.uk - 22/11/2017

This module defines the tools required to create a toy model that simulates the polarisation resulting form oblate
photospheres.
"""

from __future__ import division
import numpy as np
from numpy.random import uniform
from itertools import compress
import math as m
from astropy.modeling import models
from time import time
import progressbar
from FUSS.datred import pol_deg
import cPickle as pickle

import warnings

warnings.simplefilter("error", RuntimeWarning)


class Envellope(object):
    """
    Creates a supernova envelope with a photsphore and (optional) an absorption region

    Parameters
    ----------
        a : float
            semi_major axis. Must be less than one.
        b : float
            semi_minor axis. Must be less than one.
        number_of_points : int
            number of photons generated in photosphere
        ptpl : float
            probability of a photon at the limb to have a polarisation tangential to the photosphere

    Attributes
    ----------
        a : float
            semi_major axis. Must be less than one.
        b : float
            semi_minor axis. Must be less than one.
        N : int
            number of photons generated in photosphere
        maxproba : float
            probability of a photon at the limb to have a polarisation tangential to the photosphere
        photonsx : list
            List containing the x coordinates of the photons
        photonsy : list
            List containing the y coordinates of the photons
        photonsr : list
            List containing the projected radii at the position of each photon (r^2 = x^2 + y^2)
        photonsR : list
            List containint the radius of the photosphere at the position of each photon
            (R = x^2 + b^2(a^2-x^2)/a^2)
        photonsi : list
            Intensity of each photon
        photonspa : list
            Angle of polarisation of each photon
        photonsq : list
            q values of each photon in percent
        photonsu : list
            u values of each photon in percent
        photonsproba : list
            probability of being tangentially polarised at the position of each photon
            (doesn't tell you whether the polarisation in tangentially polarised)
        mask : list of boolean
            Mask created when defining absorption regions used to mask the photons that have been absorbed
        absorb_r : float
            Radius of the circular absorption region
        centerx : float
            X coordinate of the center of the absorption region
        centery : float
            Y coordinate of the center of the absorption region
        absorb_a : float
            semi_major axis of the absorbing ellipse
        absorb_b : float
            semi_minor axis of the absorbing ellipse
        absorb_angle : int
            angle of rotation of the absorbing ellipse (not sure it's implemented yet)
        xbarmin : float
            x coordinate at the beginning of the absorption band
        xbarmax : float
            x coordinate at the end of the absorption band
        ybarmin :
            y coordinate at the beginning of the absorption band
        ybarmax :
            y coordinate at the end of the absorption band

    Methods
    -------
        _calc_pol(): calculates the pol ofthe photons. Used internally
        absorbing_circle(): creates a circular absorption region
        absorbing_ellipes(): creates an elliptical absorption region
        absorbing_band(): creates an absorbing band accross the photosphere

    """
    def __init__(self,semi_major,semi_minor,number_of_points, ptpl):
        self.a = semi_major
        self.b = semi_minor
        self.N = number_of_points
        self.maxproba = ptpl
        self.photonsx = []
        self.photonsy = [] 
        self.photonsr = []
        self.photonsR = []
        self.photonsi = []   
        self.photonspa = []
        self.photonsq = []
        self.photonsu = []
        self.photonsproba = []
        self.mask = []
        # For Absorption regions
        self.absorb_r = None
        self.centerx = None
        self.centery = None
        self.absorb_a = None
        self.absorb_b = None
        self.absorb_angle = None
        self.xbarmin = None
        self.xbarmax = None
        self.ybarmin = None
        self.ybarmax = None
        #Average pol
        self.avg_q = None
        self.avg_u = None

        i=0
        while i < self.N+1:
            xcoord = float(uniform(-self.a, self.a, size=1))
            ycoord = float(uniform(-self.b, self.b, size=1))
            condition = (xcoord**2)/(self.a**2) + (ycoord**2)/(self.b**2)
            
            if condition <=1:
                self.photonsx.append(xcoord)
                self.photonsy.append(ycoord)
                i += 1                        ## we got coordinates
        
        k = 0.5    
        for coord in zip(self.photonsx, self.photonsy):
            x = coord[0]
            y = coord[1]
            r = np.sqrt(x*x + y*y)
            self.photonsr.append(r)
            
            if self.a >= self.b:
                R = np.sqrt(x**2 + (self.b**2)*(self.a**2-x**2)/self.a**2)     
            else: 
                R = np.sqrt(y**2 + (self.a**2)*(self.b**2-y**2)/self.b**2)
                
            self.photonsR.append(R)  
            
            try:
                self.photonsi.append(1-k*(1-np.sqrt(1-((r**2)/(R**2)))))
            except RuntimeWarning:
                print "Bad value in the sqrt, intensity set to 0. You've been warned"
                print self.photonsR, self.photonsr
                self.photonsi.append(0)

        self._calc_pol()

    def _calc_pol(self):
        for x, y, r, R in zip(self.photonsx, self.photonsy, self.photonsr, self.photonsR):
            pol_cond = uniform(0,100,size=1)
            proba = ((r/R)**2)*self.maxproba
            self.photonsproba.append(proba)
            #proba=100
            if pol_cond <= proba:
                grad = - (self.b**2 * x)/(self.a**2 * y) #polarisation tangential to the photosphere.
                angle = m.atan(grad)*180/m.pi
                               
                if angle < 0.0:
                    angle += 180

            else:
                angle=float(uniform(0,180,size=1)) #random orientation

                
            self.photonspa.append(angle)
            #### Q U ####     
            q = m.cos(2*(angle-90)*m.pi/180)*100 #0 P.A is North direction so -90  
            u = m.sin(2*(angle-90)*m.pi/180)*100
            self.photonsq.append(q)  
            self.photonsu.append(u)               
            
    def absorbing_circle(self, radius, centerx, centery):
        #self.photons_left = []
        self.mask=[]
        
        for x, y in zip(self.photonsx, self.photonsy):
            cond = np.sqrt((centerx-x)**2 + (centery-y)**2)
            if cond < radius:
                self.mask.append(False)
            else:
                self.mask.append(True)

        self.photonsx=list(compress(self.photonsx,self.mask))
        self.photonsy=list(compress(self.photonsy,self.mask))
        self.photonsi=list(compress(self.photonsi,self.mask))
        self.photonsq=list(compress(self.photonsq,self.mask))
        self.photonsu=list(compress(self.photonsu,self.mask))

    def absorbing_ellipse(self, absorb_a, absorb_b, (center), angle):
        self.absorb_a = absorb_a
        self.absorb_b = absorb_b
        self.centerx = center[0]
        self.centery = center[1]
        self.absorb_angle = angle
        self.mask=[]
        
        for coord in zip(self.photonsx, self.photonsy):
            cond = (center[0]-coord[0])**2/self.absorb_a**2 + (center[1]-coord[1])**2/self.absorb_b**2
            if cond < 1:
                self.mask.append(False)
            else:
                self.mask.append(True)

        self.photonsx=list(compress(self.photonsx,self.mask))
        self.photonsy=list(compress(self.photonsy,self.mask))
        self.photonsi=list(compress(self.photonsi,self.mask))
        self.photonsq=list(compress(self.photonsq,self.mask))
        self.photonsu=list(compress(self.photonsu,self.mask))
    
    def absorbing_band(self,rangex=None,rangey=None):
        self.mask=[]
        if rangey is None:
            self.xbarmin=rangex[0]
            self.xbarmax=rangex[1]
            for coord in self.photonsx:
                if coord<self.xbarmax and coord>self.xbarmin:
                    self.mask.append(False)
                else:
                    self.mask.append(True)
        elif rangex is None:
            self.ybarmin=rangey[0]
            self.ybarmax=rangey[1]       
            for coord in self.photonsy:
                if coord<self.ybarmax and coord>self.ybarmin:
                    self.mask.append(False)
                else:
                    self.mask.append(True)

        self.photonsx=list(compress(self.photonsx,self.mask))
        self.photonsy=list(compress(self.photonsy,self.mask))
        self.photonsi=list(compress(self.photonsi,self.mask))
        self.photonsq=list(compress(self.photonsq,self.mask))
        self.photonsu=list(compress(self.photonsu,self.mask))

    def avg_pol(self):
        self.avg_q = np.average(self.photonsq,
                    weights= self.photonsi)
        self.avg_u =  np.average(self.photonsu,
                    weights = self.photonsi)


def pickle_phot(order_mag_photon, ptpl, semi_axes, num_iteration = 0): # Tested
    start = time()
    print "\nRUNNING... \n"
    i=0
    if order_mag_photon > 10:
        return "This order of magnitude is very big. If you are sure you want to do this, change the code."

    num_photon = 10**order_mag_photon
    semi_major, semi_minor = semi_axes
    axis_ratio = semi_minor / semi_major

    print "\n PTPL ", ptpl, "/ semi-major axis =",semi_major,"/ semi-minor axis =", semi_minor,"\n"
    q_ls=[]
    u_ls=[]
    theta_ls=[]
    bar = progressbar.ProgressBar(max_value=num_iteration+1)
    while i < num_iteration+1:
        bar.update(i)
        photosphere = Envellope(semi_major,semi_minor,num_photon,ptpl)
        print "Made the photosphere, now need to pickle it"

        q = np.average(photosphere.photonsq,
                        weights= photosphere.photonsi)

        u = np.average(photosphere.photonsu,
                        weights = photosphere.photonsi)
        print "Q = ",q,"/ U = ", u
        print "Theta = ", 0.5*m.atan2(u,q)*180/m.pi
        q_ls.append(q)
        u_ls.append(u)
        theta_ls.append(0.5*m.atan2(u,q)*180/m.pi)

        output_file = "photosphere"+str(order_mag_photon)+"_ax"+str(int(100*axis_ratio))+"_ptpl"+str(int(ptpl*10))+"_iter"+str(i+1)+".pkl"
        with open(output_file, 'wb') as output:
            pickle.dump(photosphere, output, pickle.HIGHEST_PROTOCOL)

        del photosphere

        i+=1
    bar.update(i)

    return "\nTotal run time: ", (time() - start)/60, " min / ", time()-start, "seconds\n"
