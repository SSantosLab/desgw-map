# countour_sum:
#   find the number of pixels that contain some total
#   sum of vals (if val is normed to 1,
#   this is a percentage) 
# area
#   one can convert this to area
#
import healpy as hp
import numpy as np
def area(ra,dec,vals, threshold, nside) :
    npix = contour_sum(ra,dec,vals,threshold)
    area_per_pixel = hp.nside2pixarea(nside)*((360./2/np.pi)**2)
    area = npix*area_per_pixel
    return area

def contour_sum(ra, dec, vals, threshold) :
    ix = np.argsort(vals)[::-1] ;# descending order
    ra_sort = ra[ix]
    dec_sort = dec[ix]
    vals_sort = vals[ix]

    cumulative = 0
    i = 0
    if vals.sum() < threshold :
        print "contour_sim: vals.sum() < threshold, npix set to max. sum= ", vals.sum()
        npixels = ra.size
    else :
        while cumulative < threshold :
            cumulative +=  vals_sort[i]
            i += 1
        npixels =ra_sort[:i].size 
    return npixels

