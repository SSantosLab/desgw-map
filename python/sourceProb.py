import numpy as np
import os
import scipy.stats
import mags
import hp2np

license="""
   Copyright (C) 2014 James Annis

   This program is free software; you can redistribute it and/or modify it
   under the terms of version 3 of the GNU General Public License as
   published by the Free Software Foundation.

   More to the points- this code is science code: buggy, barely working,
   with little or no documentation. Science code in the the alpine fast 
   & light style. (Note the rate at which people who stand on the
   summit of K2 to successfully make it down.)

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
"""

class map(object):
    """
    """
    def __init__(self, observation, lumModel="gaussian") :
        data_dir = os.environ["DESGW_DATA_DIR"] 

        self.limitingMag = observation.maglim
        self.limits      = observation.limits
        self.lumModel    = lumModel
        self.delMag = 0.1
        self.minAbsMag = -10.0
        self.maxAbsMag = -13.0
        self.modelAbsoluteMagnitude = -11.1
        min = 1.0   # Mpc
        max = 100.0 # Mpc
        self.minDistance = min
        self.maxDistance = max
        # we'll do this mpc by mpc
        self.searchDistance = np.arange(min, max+1, 1.0)

        # get the P_recognition
        self.pra = observation.pra
        self.pdec = observation.pdec
        self.precog = observation.precog

        # keep the answer
        self.probMap = ""

        # for plotting contours, keep the intermediate data around
        self.xi, self.yi, self.zi = ["","",""]
        self.lastCtype = ""
        
    # one can choose to plot the ligo contours
    # or the ligo*obs_prob contours
    #   type="ligo"   type="ls"
    #           houranlge chooses HA instead of RA projected into mcbryde coords
    def plotLigoContours(self, obs, type="ligo", whiteLine="False", hourangle=False) :
        import matplotlib
        import matplotlib.pyplot as plt
        con_levels=10
        if self.zi == "" or self.lastCtype != type:
            print "\t calculating contours for type = ",type
            if hourangle == False :
                xmin = obs.x.min(); xmax = obs.x.max()
                ymin = obs.y.min(); ymax = obs.y.max()
            else :
                xmin = obs.hx.min(); xmax = obs.hx.max()
                ymin = obs.hy.min(); ymax = obs.y.max()
            xi=np.linspace(xmin, xmax, 500)
            yi=np.linspace(ymin, ymax, 500)
            if type == "ligo" :
                probMap = obs.map
            if type == "ls" :
                probMap = obs.map*self.probMap
            if hourangle == False :
                x = obs.x; y = obs.y
            else :
                x = obs.hx; y = obs.hy
            zi = matplotlib.mlab.griddata(x,y,probMap,xi,yi)
            self.xi, self.yi, self.zi = xi,yi,zi
            self.lastCtype = type
        if not whiteLine :
            plt.contour(self.xi,self.yi,self.zi,con_levels,linewidths=3,colors="k")
        plt.contour(self.xi,self.yi,self.zi,con_levels,linewidths=0.66,colors="w")
        

    def calculateProb(self) :
        # bookkeeping
        self.zi= ""
        # realwork
        searchDistance = self.searchDistance
        limitingMag = self.limitingMag 
        telescopeLimits = self.limits
        lumModel = self.lumModel
# radial
            # for now, lets skip the source distance probability
        #PsourceDistance = self.sourceDistance(searchDistance, 
        #    model="gaussian", distanceMean=70, distanceSigma=10.)
        PvolumeWeight = self.stellarProbabilty(searchDistance, model="volume") 
        # for testing
        self.PvolumeWeight = PvolumeWeight
            # for now, lets skip the source distance probability
        #self.PsourceDistance = PsourceDistance
# radial feeds into radial per map pixel
            # for now, lets skip the source distance probability
        #Pdistance = PsourceDistance * PvolumeWeight
        Pdistance = PvolumeWeight
        # ok. The sims don't want a volume weight- the simulation is for
        # an object -at- a particular distance.
        Pdistance = np.ones(PvolumeWeight.size)
        Pdistance = Pdistance/Pdistance.sum()
        prob = self.probWeSeeIt(searchDistance, limitingMag, Pdistance, lumModel)
        # for testing
        self.Pdistance = Pdistance
        self.prob = np.copy(prob)
# probability of recognition propto star density
        prob = prob* self.precog

# finally, eliminate every where the telescope cannot point
        prob = prob* telescopeLimits

        self.probMap = prob
        print "\t probMap: made "
        return 
        
    # 
    #  probWeSeeIt  = P_m * P_dist * P_d
    #   P_m is the probability at a given abs mag (or brighter) given the model
    #   P_dist the probability the source is at a given distance, comes in as cdf from near to far
    #   P_d is the probability of detection given ap mag (1 or 0)
    #
    # if the model is "distance", then the return is a distance map, usefully.
    # 
    def probWeSeeIt(self, searchDistance, limitingMag, Pdistance, model="gaussian") :
        import sys

        # under some models, this resets minAbsMag and maxAbsMag
        self.sourceMagnitudeProbabilityInit(model) 

        delMag    = self.delMag 
        minAbsMag = self.minAbsMag 
        maxAbsMag = self.maxAbsMag 

        # limitingMag is the healpix map
        print "\t source abs mag: ",
        sumProb = np.zeros(limitingMag.size)
        for absm in np.arange(maxAbsMag, minAbsMag, delMag) :
            print absm, " ",
            sys.stdout.flush()
            probOfAbsMag = self.sourceMagnitudeProbability(model, absm) 
            # loop over distance checking app mag brighter than mag limit
            for i in range(0,searchDistance.size) :
                apparentMag = self.apparentMagnitudeFromAbs(searchDistance[i], absm)
                # find the pixels in the map such that the apparent mag
                # is brighter than the limiting magnitude
                ix = apparentMag <= limitingMag 
                # the aim here is to write each abs mag probility, in order
                # so that the brightest absMag brighter than the apparentMag 
                #   sets  the probability
                if model == "distance" :
                    # calculate the distance limit
                    sumProb[ix] = searchDistance[i]
                else :
                    # calculate the probabilites
                    #print "\nprobOfAbsMag: ", probOfAbsMag
                    #print "Pdistance[i]: ", Pdistance[i]
                    sumProb[ix] = probOfAbsMag * Pdistance[i]
        print " "
        return sumProb

    def sourceMagnitudeProbabilityInit(self, model) :
        if model == "uniform"  : print "\t source model: uniform"
        elif model == "gaussian" : print "\t source model: gaussian"
        elif model == "gaussianAM" : print "\t source model: gaussianAM"
        elif model == "delta" : print "\t source model: delta function"
        elif model == "known" : print "\t source model: known model"
        elif model == "distance" : print "\t source model: delta function \
            calculating distance instead of probability"
        else : raise Exception("no such source model {}".format(model))
        if model == "delta" or model == "distance" :
            self.minAbsMag = -11.0
            self.maxAbsMag = -11.1
        if model == "known" :
            self.minAbsMag = self.modelAbsoluteMagnitude+self.delMag
            self.maxAbsMag = self.modelAbsoluteMagnitude

        return ""
    def sourceMagnitudeProbability(self, model, absoluteMagnitude) :
        # what is the probability the source has this absolute magnitude or brighter
        if model == "uniform"  : 
            probOfAbsMag = self.uniformMagnitude(absoluteMagnitude)
        elif model == "gaussian" : 
            probOfAbsMag = self.gaussianMagnitude(absoluteMagnitude)
        elif model == "gaussianAM" : 
            probOfAbsMag = self.gaussianMagnitudeAM(absoluteMagnitude)
        elif model == "delta" : 
            probOfAbsMag = self.delta(absoluteMagnitude)
        elif model == "distance" : 
            probOfAbsMag = self.delta(absoluteMagnitude)
        elif model == "known" : 
            probOfAbsMag = self.delta(absoluteMagnitude)
        return probOfAbsMag


    # aiming at the idea that if the abs mag is this or less luminous,
    # then yea; cdf goes from 0 at high lum to 1 at low lum
    # 
    # We'll use a guassian about the analytic kilonova model
    #   what is the probability the object has this mag or brighter
    #   So, caclulate the cumulative distribution function,
    #   relying on brighter = more negative number, and cdf goes from 
    #   small to large indepenent variable
    # 
    def gaussianMagnitudeAM(self, absoluteMagnitude) :
        mean = -11.1
        sigma = 0.1
        pdf= scipy.stats.norm(loc=mean,scale=sigma)
        probability = pdf.cdf(absoluteMagnitude)
        return probability
    def delta(self, absoluteMagnitude) :
        # taking it from the object means that the absolute magnitude
        # can be changed and everything recomputed after the creation 
        # of the object
        mean = self.modelAbsoluteMagnitude
        sigma = 0.00001
        # the role of epsilon is get the cdf past the 0.5 inflexion point
        epsilon = sigma * 10
        #print ""
        #print mean, sigma
        #print ""
        pdf= scipy.stats.norm(loc=mean,scale=sigma)
        probability = pdf.cdf(absoluteMagnitude+epsilon);
        return probability
    # 
    # We'll use a guassian fit to the i-band mags in the Barnes and Kasen table 1. 
    #   what is the probability the object has this mag or brighter
    #   So, caclulate the cumulative distribution function.
    # 
    def gaussianMagnitude(self, absoluteMagnitude) :
        mean = -14.05
        sigma = 1.34
        pdf= scipy.stats.norm(loc=mean,scale=sigma)
        probability = pdf.cdf(absoluteMagnitude)
        return probability

    # 
    # We won't assume complete ignorance about source magnitudes-
    # just uniform probability between -12 and -15, in 0.1 mag bins 
    #   So, caclulate the cumulative distribution function.
    # 
    def uniformMagnitude(self, absoluteMagnitude) :
        delMag    = self.delMag 
        minAbsMag = self.minAbsMag 
        maxAbsMag = self.maxAbsMag 
        raise Exception("error - not yet converted to cdf, gigo")
        if (absoluteMagnitude <= minAbsMag) & (absoluteMagnitude >= maxAbsMag) :
            nbins  = (maxAbsMag-minAbsMag)/delMag
            probability = 1.0/nbins
        else :
            probability = 0.0
        return probability

    # 
    # prior probability of source distance
    #
    def sourceDistance  (self, searchDistance, model="uniform", 
        distanceMean=50., distanceSigma=10.) :
        if model == "uniform" :
            probability = self.uniformDistance(searchDistance)
        elif model == "gaussian" :
            probability = self.gaussianDistance(searchDistance, distanceMean, distanceSigma)
        else :
            raise Exception ("{} is not a model that is implemented".format(model))
        return probability

    # 
    # Maximum ignorance of the source distance
    # uniform prior between 1.0 and 100 Mpc
    #
    def uniformDistance(self, searchDistance) :
        min =self.minDistance 
        max = self.maxDistance
        ix = np.nonzero( (searchDistance < min) | (searchDistance > max) )
        probability = searchDistance*0.0 + 1.0
        probability[ix] = 0.0
        probalility = probability/probability.sum()
        return probability

    # 
    # gaussian probability on a distance
    #
    def gaussianDistance(self, searchDistance, mean, sigma) :
        min =self.minDistance 
        max = self.maxDistance
        ix = np.nonzero( (searchDistance < min) | (searchDistance > max) )
        pdf= scipy.stats.norm(loc=mean,scale=sigma)
        probability = pdf.pdf(searchDistance)
        probalility = probability/probability.sum()
        return probability

    def stellarProbabilty(self, searchDistance, model="volume") :
        if model == "volume" :
            probability = self.volumeWeight(searchDistance)
        else :
            raise Exception ("{} is not a model that is implemented".format(model))
        probability = probability/probability.sum()
        return probability
    #
    # model based on known galaxies
    #
    def galaxyCatWeight(self,) :
        pass
        return ""

    # 
    # Assume the probability for event is
    # proportional to volume at that distance
    #   this is, for 1 Mpc slices, (theta*d)^2
    #   divided by (1/3) theta maxDistance^3
    #
    def volumeWeight (self, searchDistance) :
        probability = searchDistance**2/(3.*self.maxDistance**3)
        # now turn this into a cumulative function, from small to large distance
        cdf = np.cumsum(probability)/probability.sum()
        return cdf

    def apparentMagnitudeFromAbs(self, searchDistance, absMag) :
        distanceModulus = 5*np.log10(searchDistance*1e6/10.)
        apparentMag  = absMag + distanceModulus
        return apparentMag

