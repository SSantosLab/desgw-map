import numpy as np
import matplotlib.path

# the intent is that these two main routines,
#   hexesOnMap & hexalateMap
# are widely useful ways to deal with DECam and DES
# data with healpix maps
# nside=32 ~ camera   3.36 sq-degree
# nside=256 ~ ccd     0.052 sq-degrees
# nside=512 resolves ccds

# keeps count of times a hex overlies a map pixel
# camera outline for the hexes given
def hexesOnMap(ra,dec, raHexen, decHexen) :
    print "    checking {} hexes : ".format(raHexen.size),
    count = np.zeros(ra.size)
    for i in range(0,raHexen.size) :
        ix = radecInHex( raHexen[i], decHexen[i], ra, dec)
        count[ix] += 1
    print " "
    return count

# if one wants to know the sum of the ligo probability in
# a set of observed hexes, hexalateMap is the routine to use.
#
# return the sum of the map vals inside the hexMap hexes
#
def hexalateMap(ra,dec, vals, raHexen, decHexen, verbose=1, useCircle=False, radius=1.0) :
    if verbose : print "\t hexalateMap \t nhex = {},".format(raHexen.size),
    if verbose: print " npix = {}".format(ra.size)
    hexVal = np.zeros(raHexen.size)
    for i in range(0,raHexen.size) :
        if useCircle :
            ix = radecInCircle( raHexen[i], decHexen[i], ra, dec, radius=radius)
        else :
            ix = radecInHex( raHexen[i], decHexen[i], ra, dec)
        try  :
            hexVal[i] = vals[ix].sum()
        except Exception: 
            print "why are there exceptions in hexalateMap?"
            hexVal[i] = 0
    return hexVal

#
# Use a circle of area pi sq-degrees as a reasonable
# approximation to the camera outline.
#
def radecInCircle (raCenter, decCenter, ra, dec, radius=1.0) :
    raDist  = (ra  - raCenter )
    decDist = (dec - decCenter)
    raDist  = raDist * np.cos(dec * 2 * np.pi/360.)
    distance = np.sqrt(raDist**2 + decDist**2)

    ix = distance <= radius

    return ix

# return an index that is true/exists if ra,dec is inside
# the camera outline for the hex 
def radecInHex ( raCenter, decCenter, ra, dec) :
    # cut away hopless areas in ra, dec
    cosDec = np.cos(decCenter*2*np.pi/360.)
    near_ix = (abs(ra-raCenter)*cosDec <= 1.1) & (dec-decCenter <= 1.1)
    ## just doing the dec sped it up by x2
    #near_ix = (dec-decCenter < 1.5)
    # if there is nothing near, no reason to check further
    if np.all(~near_ix) :
        return near_ix

    # find the answer
    hex = hexPath(raCenter, decCenter)
    inside_ix = radecInHexPath( hex, ra[near_ix], dec[near_ix])
    
    # the issue here is that near_ix is of size ~1,000,000
    # although only ~3000 are True, so inside_ix is of size ~3,000
    # and about ~200 are True. One wants to adjust near_ix using
    # the truth values of inside_ix.
    near_ix[near_ix] = inside_ix
    # simple, but not obvious

    return near_ix
    
def radecInHexPath ( hex_path, ra, dec) :
    ix = hex_path.contains_points( zip(ra,dec) )
    return ix

def hexPath (raCenter, decCenter) :
    ra,dec = cameraOutline(raCenter, decCenter)
    hex = matplotlib.path.Path(zip(ra,dec))
    return hex

def cameraOutline (raCenter, decCenter) :
    ra  = [-0.47256, -0.47256, 0, 0.47256, 0.47256, 0.63008, 0.63008,
        0.7876, 0.7876, 0.94512, 0.94512, 0.94512, 1.10264, 1.10264,
        1.10264, 0.94512, 0.94512, 0.94512, 0.7876, 0.7876, 0.63008,
        0.63008, 0.47256, 0.47256, 0, -0.47256, -0.47256, -0.63008, -0.63008,
        -0.7876, -0.7876, -0.94512, -0.94512, -0.94512, -1.10264, -1.10264,
        -1.10264, -0.94512, -0.94512, -0.94512, -0.7876, -0.7876, -0.63008,
        -0.63008, -0.47256]
    dec = [ 0.824266666667, 0.98912, 0.98912, 0.98912, 
        0.824266666667, 0.824266666667, 0.659413333333,
        0.659413333333, 0.49456, 0.49456, 0.329706666667, 
        0.164853333333, 0.164853333333, 0.0,
        -0.164853333333, -0.164853333333, -0.329706666667, 
        -0.49456, -0.49456, -0.659413333333, -0.659413333333,
        -0.824266666667, -0.824266666667, -0.98912, -0.98912, 
        -0.98912, -0.824266666667, -0.824266666667, -0.659413333333,
        -0.659413333333, -0.49456, -0.49456, -0.329706666667, 
        -0.164853333333, -0.164853333333, 0.0,
        0.164853333333, 0.164853333333, 0.329706666667, 0.49456, 
        0.49456, 0.659413333333, 0.659413333333,
        0.824266666667, 0.824266666667]
    ra = np.array(ra)
    dec=np.array(dec)

    dec = decCenter + dec
    ra = raCenter + ra/np.cos(dec*2*np.pi/360.)
    return ra,dec

