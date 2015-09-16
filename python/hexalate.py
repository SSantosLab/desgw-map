license="""
   Copyright (C) 2015 James Annis

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

import numpy as np
import decam2hp

#
# Example usage
#
#    obs.resetTime(mjd+time)
#    obs.limitMag("i", exposure=180)
#
#    sm=sourceProb.map(obs, lumModel="known");
#    models_at_t = modelsAtTimeT (models, time)
#    abs_mag = models_at_t[0]
#    sm.modelAbsoluteMagnitude = abs_mag
#    sm.searchDistance = np.array([distance,])
#    sm.calculateProb()
#
#    hexVals, rank = cutAndHexalate (obs, sm, raHexen, decHexen)

#
# Given an obs object and a sm object, along with a map of hex centers,
# this routine returns the sum of the probability inside of DECam camera outlines
# for all hexes in the map, along with a rank of the hexes sorted from large to small.
#
# domain knowlege- the blanco cuts keep everything inside an hour angle range
def hexalateNHexes (obs, sm, nHexes, allskyDesHexes) :
    raHexen, decHexen, hexVals, rank = cutAndHexalate(obs, sm, allskyDesHexes)
    raHexen = raHexen[0:nHexes]
    decHexen = decHexen[0:nHexes]
    hexVals = hexVals[0:nHexes]
    rank = rank[0:nHexes]
    return raHexen, decHexen, hexVals, rank

# be aware that while most of my ra,dec are in degrees,
# those in obs and sm are in radians
def cutAndHexalate (obs, sm, allskyDesHexes="../data/all-sky-hexCenters.txt") :
    verbose = False
    obsHourAngle = obs.ha*360./(2*np.pi)
    obsRa        = obs.ra*360./(2*np.pi)
    obsDec       = obs.dec*360./(2*np.pi)
    # based on blanco horizen limits
    ix = (abs(obsHourAngle) <= 83. ) & (obsDec < 43.)

    raHexen, decHexen = getHexCenters (allskyDesHexes)
    ix2 = decHexen < 43.

    probabilities = obs.map*sm.probMap
    if verbose  :
        print "\t cutAndHexalate probabilities sum",probabilities.sum()

    hexVals = np.zeros(raHexen.size)
    hexVals[ix2] = decam2hp.hexalateMap(obsRa[ix],obsDec[ix], probabilities[ix],
        raHexen[ix2], decHexen[ix2])
    if verbose  :
        print "hexVals max", hexVals.max()
    rank=np.argsort(hexVals); 
    rank = rank[::-1];# sort from large to small by flipping natural argsort order
    return raHexen, decHexen, hexVals, rank

def getHexCenters (allskyDesHexes = "../data/all-sky-hexCenters.txt") :
    ra,dec = np.genfromtxt(allskyDesHexes, unpack=True)
    return ra,dec

