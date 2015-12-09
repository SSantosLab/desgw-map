import numpy as np
from scipy.spatial import cKDTree

license="""
   Copyright (C) 2014 James Annis

   This program is free software; you can redistribute it and/or modify it
   under the terms of version 3 of the GNU General Public License as
   published by the Free Software Foundation.

   More to the points- this code is science code: buggy, barely working,
   with little or no documentation. Science code in the the alpine fast 
   & light style.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
"""

def blancoLimits (ha, dec) :
    import matplotlib.path
    transmission = np.zeros(ha.size)

    # get the limits
    tel_ha,tel_dec = blancoHorizonLimits()
    blanco = matplotlib.path.Path(zip(tel_ha,tel_dec))

    # cut away hopless areas in ha, dec
    near_ix = (abs(ha) <= 6) & ( dec < 45)
    # if there is nothing near, no reason to check further
    if np.all(~near_ix) :
        return transmission

    # find the answer
    inside_ix = blanco.contains_points( zip(ha,dec) )

    # the issue here is that near_ix is large, though few are True
    # thus inside_ix is small, and even less are True.
    # One wants to adjust near_ix using
    # the truth values of inside_ix.
    near_ix[near_ix] = inside_ix
    # simple, but not obvious

    transmission[near_ix] = 1.0
    return transmission


# read from http://www.ctio.noao.edu/noao/content/Horizon-Limits
def blancoHorizonLimits () :
    ha = [-78., -70., -60., -50., -40., -30., -20., -10.,   0.,  10.,  20.,  30.,  40.,  50.,  60.,
  70.,  78.,  78.,  78.,  78.,  78.,  78.,  78.,  78.,  70.,  65.,  60.,  55.,  50.,  46.,
  40.,  32.,  30.,  25.,  20.,  15.,  10.,   5.,   0.,  -5., -10., -15., -20., -25., -30.,
 -32., -40., -46., -50., -55., -60., -65., -70., -78., -78., -78., -78., -78., -78., -78.,
 -78.]

    dec = [
 -89., -89., -89., -89., -89., -89., -89., -89., -89., -89., -89., -89., -89., -89., -89.,
 -89., -89., -80., -70., -60., -50., -40., -30., -28., -10.,  -1.,   5.,  11.,  17.,  20.,
  25.,  30.,  31.,  33.,  34.,  35.,  36.,  37.,  37.,  37.,  36.,  35.,  34.,  33.,  31.,
  30.,  25.,  20.,  17.,  11.,   5.,  -1., -10., -28., -30., -40., -50., -60., -70., -80.,
 -89.]

    ha = np.array(ha)*2*np.pi/360.
    dec = np.array(dec)*2*np.pi/360.

    return ha,dec

#
# obsolete codes
#

#        tel_dir           = data_dir
#        tel_file          = "blanco-limits-ha-dec.fits"
#        ha,dec,ttlim      = telescope.loadLimits(tel_dir, tel_file)
def loadLimits(dir = "/home/s1/annis/daedalean/desgw-map/data/",
        file = "blanco-limits-ha-dec.fits") :
    print "\t loading telescope limits ",dir+file
    import pyfits
    f=pyfits.open(dir+file)
    t = f[1].data;
    ha = t.field("ha")
    dec = t.field("dec")
    transmission = t.field("transmission")
    ix = np.nonzero(transmission < 0.1)
    transmission[ix] = 0.0

    return ha, dec, transmission

def findLimits (ha, dec, tree, trans) :
    coords = zip(ha,dec)
    transmission = []
    for i in range(0, ha.size) :
        testHA  = coords[i][0]
        testDec = coords[i][1]
        if (testDec > 40 or abs(testHA) > 80 or testDec < -88.5) :
            # clearly out of the ok region
            transmission.append( 10.0**-100.)
        elif testDec <= -30  :
            # else inside the square box region if at dec < -30
            transmission.append( 1.0 )
        elif ( (testDec <= 10.) and abs(testHA) <= 60. ) :
            # else work way up the wedding cake of rectangles  1
            transmission.append( 1.0 )
        elif ( (testDec <= 20.) and abs(testHA) <= 40. ) :
            # else work way up the wedding cake of rectangles  2
            transmission.append( 1.0 )
        elif ( (testDec <= 30.) and abs(testHA) <= 20. ) :
            # else work way up the wedding cake of rectangles  3
            transmission.append( 1.0 )
        else :
            # guess not, must do the calculation
            distAndCoord = tree.query(coords[i])
            tcoord = distAndCoord[1]
            transmission.append(trans[tcoord])
    transmission = np.array(transmission)
    return transmission

def tree (ha, dec) :
    coords = zip(ha,dec)
    tree = cKDTree(coords)
    return tree

