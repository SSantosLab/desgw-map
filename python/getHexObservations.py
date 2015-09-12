import numpy as np
import jsonMaker
#
#       the routine mapsAtTimeT.oneDayOfTotalProbability
#           breaks the night into slots of 32 minutes duration
#       the routine mapsAtTimeT.probabilityMapSaver
#           makes maps for each of these slots
#       the routine getHexObservations.observing
#           once told how many slots can be used
#           returns "hoursObserving" containing ra,dec,prob,time
#               of the hexes to observe each slot
#       the routine getHexObservations.observingStats
#           places the ra,dec,prob,times per slot onto single  lists
#           tells sum(prob) on the list
#       the routine getHexObservations.observingPlot
#           needs to be told simNumber, "slot", data_dir, nslots
#           and then will plot the sky limitMag+ligo contour
#           and all hexes (in the night) on the map for that slot
#

# ==== prep the observation by calculating necessary maps
#       distance is a problem- the question is whether we want a horizon
#       distance or the estimated distance from LIGO?
#       >>> distance = 75. ;# Mpc, as an estimated horizon distance
def prepare(skymap, mjd, trigger_id, outfolder,
        distance=75, exposure_length=180.) :
    import mapsAtTimeT
    import mags
    import modelRead
    import healpy as hp
    import hp2np
    # === prep the maps
    ligo = hp.read_map(skymap)
    ra,dec,ligo = hp2np.map2np(ligo,256, fluxConservation=True)
    obs = mags.observed(ra,dec,ligo, mjd)
    obs.limitMag("i",exposure=exposure_length)
    # ==== get the neutron star explosion models
    models = modelRead.getModels()
    # ==== calculate maps during a full night of observing
    probs,times = mapsAtTimeT.oneDayOfTotalProbability(obs,mjd,distance,models)
    mapsAtTimeT.probabilityMapSaver (obs, trigger_id, mjd, \
        distance, models, times, probs, outfolder)
    return probs, times


# ==== figure out what to observe
def now(n_slots, mapDirectory="jack/", simNumber=13681, mapZero=0) :

    hoursObserving=observing(simNumber,n_slots,mapDirectory, mapZero=mapZero)
    ra,dec,prob,mjd,slotNumbers = observingStats(hoursObserving)
    ra,dec,prob,mjd,slotNumbers = observingRecord(
        hoursObserving, simNumber, mapDirectory)
    maxProb_slot = slotNumbers[np.argsort(prob)][0]
    
    name = eventName(mapDirectory, str(simNumber)) + ".json"
    jsonMaker.writeJson(ra,dec, jsonFilename=name)

    return maxProb_slot, slotNumbers, mjd

#   makeObservingPlots(n_slots, slotNumbers, mjd, simNumber, best_slot, mapDirectory) 

def eventName(data_dir, event) :
    name=data_dir+str(event)
    return name

# ok, I designed observing for the case
def contemplateTheDivisionsOfTime(probs, times, hoursAvailable=6) :
    n_slots = findNSlots(hoursAvailable)
    n_maps = times.size
    if n_maps == n_slots : 
        mapZero = 0
    elif n_maps < n_slots : 
        mapZero = 0
        n_slots = n_maps
    elif n_maps > n_slots :
        mapZero = findStartMap ( probs, times, n_slots )
    else :
        raise Exception ("no possible way to get here")
    return n_slots, mapZero

def findNSlots(hoursAvailable, slotDuration=32.) :
    nslots = int(round(round(hoursAvailable*60./32.)))   ;# 32 minutes/slot
    return nslots


#   where the number of slots in a night
#   was equal to the number of maps made.
# This is unrealist for two reasons:
#   the night could be longer than the time desgw wishes to allocate
#   the object could be up for less than the time desgw could allocate
# so:
    # possibilities
    # n_maps = n_slots     design case
    # n_maps < n_slots     set n_slots = n_maps
    # n_maps > n_slots     pick highest contiguous set of n_slot maps
# the first two are easy. 
# the third runs into the naming convention of the maps
def findStartMap ( probs, times, n_slots ) :
    n_maps = times.size
    mapNums = range(0,n_maps)
    # sum the probability in each allowed nslot range
    n_mapsToTProb = np.array([])
    n_mapStart = np.array([])
    for map in range(0,n_maps-n_slots) :
        ix = (mapNums >= map) & (mapNums < map+n_slots)
        n_mapsToTProb = np.append(n_mapsToTProb, probs[ix].sum() )
        n_mapStart = np.append(n_mapStart, map)
    minToMax = np.argsort(n_mapsToTProb)
    bestStart = n_mapStart[minToMax[-1]]
    bestStart = int(bestStart)
    return bestStart


# Load all n -hexVals files,
#   pick the highest probability one, 
#   put it into one of the n slots, unless that slots is maxed out
#   remove that hex from all the hexVals lists
#   do it again, untill all n time slots are full.
#       maxHexesPerSlot=4 comes from 32 minute duration slots
#       and 8 minutes/hex (izzi 2 min/image)
def observing(sim, nslots, data_dir, 
        maxHexesPerSlot = 4, mapZero = 0, verbose=0) :
    # prep the observing lists
    observingSlots = np.arange(0,nslots)
    slotsObserving = dict()
    slotsObserving["nslots"] = nslots
    for i in observingSlots :
        slotsObserving[i] = 0
        slotsObserving[i,"ra"]   = np.array([])
        slotsObserving[i,"dec"]  = np.array([])
        slotsObserving[i,"prob"] = np.array([])
        slotsObserving[i,"mjd"] = np.array([])

    # read in the hexelated probability data
    hexData = dict()
    for i in observingSlots :
        map_i = i + mapZero
        raHexen, decHexen, hexVal, rank, mjd = loadHexalatedProbabilities(
            sim, map_i, data_dir)
        print map_i, "map size= ",raHexen.size, 
        ix = np.nonzero(hexVal < 1e-6)
        raHexen, decHexen, hexVal, mjd  = \
            np.delete(raHexen, ix), \
            np.delete(decHexen, ix), \
            np.delete(hexVal, ix), \
            np.delete(mjd, ix) 
        print "; n hexes greater than 1e-6 probability= ",raHexen.size
        hexData[i] = raHexen, decHexen, hexVal, mjd

    # start the search for all max probabilities
    # we'll assume the list is less than 40,000 long, the n-sq-degrees/sky
    for n in range(0,40000) :
        # search for a single max probabilities
        maxRa, maxDec, maxProb, maxMjd, maxSlot = findMaxProbOfAllHexes(hexData, observingSlots, n, verbose) 
        maxData = maxRa,maxDec,maxProb,maxMjd,maxSlot

        # we've found the maximum probability on the lists, 
        # so add it to the obs lists # unless not possible. 
        # If the latter, delete it from that slot
        slot = maxData[4]
        if slotsObserving[slot] < maxHexesPerSlot : 
            # it is possible to make the observation, 
            # put it onto the observing lists
            slotsObserving = addObsToSlot (slotsObserving, maxData, slot)
            if verbose >= 1: print n, "slot of max:",slot
        else :
            # but if this slot of observing is full, it is not possible 
            # to make the observation,
            # so move on AFTER deleting it from the list
            hexData = deleteHexFromSlot (hexData, slot, maxProb) 
        #if verbose >= 2: 
        #   if n > 7: raise Exception("jack")
    
        # perform the necessary bookkeeping, 
        # eliminating this hex from future observing
        hexData = deleteHexFromAllSlots (hexData, observingSlots, maxRa, maxDec, verbose, n) 

        # do some summary statistics
        sumHexes = 0
        sumObs = 0
        for i in range(0,nslots) :
            sumHexes += hexData[i][0].size
            sumObs += slotsObserving[i]

        if verbose >= 2: 
            print "sumHexes =", sumHexes, 
            print "   slots left=", len(observingSlots),
            print "   slots=",observingSlots,
            print "   n_obs=",
            for i in observingSlots:
                print slotsObserving[i],
            print "   sum prob= ",
            for i in observingSlots:
                print " {:8.6f}".format( slotsObserving[i,"prob"].sum()) ,
            print ""

        # eliminate full observing slots
        observingSlots = eliminateFullObservingSlots(\
            hexData, slotsObserving, observingSlots, maxHexesPerSlot, verbose) 

        # check and see if we are done
        # two conditions: observing is full, or candidates empty
        if (len(observingSlots)==0) | (sumHexes == 0) :
            print "======================================== "
            if verbose >= 1: 
                print "n slots =", len(observingSlots)," == 0?"
                print "sumHexes = ", sumHexes, "==? 0"
            print "number of hexes observed = ", sumObs
            print "======================================== "
            return slotsObserving 

        # otherwise, go back and do it again
    
    # we've done everything on the lists, we can observe it all,
    # return this great success that will never be reached.
    return slotsObserving
    
#
# examine the statistics of the observing lists
#
def observingStats( slotsObserving, outfile="") :
    nslots = slotsObserving["nslots"]
    for i in range(0,nslots) :
        print i, "n obs= {}".format(slotsObserving[i,"ra"].size), 
        print "  sum prob= {:7.4f} %".format( 100*slotsObserving[i,"prob"].sum())
    ra,dec,prob,mjd,slots = slotsObservingToNpArrays(slotsObserving) 

    print "tot prob= {:.1f}%".format(100.*prob.sum())
    return ra,dec,prob,mjd,slots

def observingRecord(slotsObserving, simNumber, data_dir) :
    name = eventName(data_dir, str(simNumber)) + "-ra-dec-prob-mjd-slot.txt"
    ra,dec,prob,mjd,slots = slotsObservingToNpArrays(slotsObserving) 
    data = np.array([ra, dec, prob, mjd, slots]).T
    np.savetxt(name, data, "%.6f %.5f %.6f %.4f %d")
    return ra,dec,prob,mjd,slots

#     ra,dec,prob,mjd,slots = readObservingRecord(simNumber, data_dir)
def readObservingRecord(simNumber, data_dir) :
    name = eventName(data_dir, str(simNumber)) + "-ra-dec-prob-mjd-slot.txt"
    ra,dec,prob,mjd,slots = np.genfromtxt(name,unpack=True)
    return ra,dec,prob,mjd,slots

def slotsObservingToNpArrays(slotsObserving) :
    nslots = slotsObserving["nslots"]
    ra = np.array([])
    dec = np.array([])
    prob = np.array([])
    mjd = np.array([])
    slots = np.array([])
    for i in range(0,nslots) :
        ra = np.append(ra, slotsObserving[i,"ra"])
        dec = np.append(dec, slotsObserving[i,"dec"])
        prob = np.append(prob, slotsObserving[i,"prob"])
        mjd = np.append(mjd, slotsObserving[i,"mjd"])
        slots = np.append(slots, i*np.ones( (slotsObserving[i,"mjd"]).size ))
    return ra,dec,prob,mjd,slots
    

#===================================================================
#
# Read in all of the hexalated probability files
#
def loadHexalatedProbabilities(sim, slot, data_dir) :
    nameStem = eventName(data_dir, str(sim)) + "-{}".format(str(slot)) 
    name = nameStem + "-hexVals.txt"
    raHexen, decHexen, hexVal, rank, mjd = np.genfromtxt(name, unpack=True, delimiter=",")
    return raHexen, decHexen, hexVal, rank, mjd


#
# search for the single highest probability hex over all of the possible hexes
# in the hexData slots 
#
def findMaxProbOfAllHexes(hexData, observingSlots, n="", verbose = 0) :
    maxProb = -1
    for i in observingSlots :
        data = hexData[i]
        hexRa  = data[0]
        hexDec = data[1]
        hexVal = data[2]
        hexMjd = data[3]
        if hexVal.size == 0: continue
        if verbose >= 2: 
            if i == 2: print n,"====",i, "hexSize =",hexRa.size
        # now check for max prob
        newProb = hexVal.max()
        if verbose >= 4: print n,i, maxProb, ">?", newProb, "     n=",hexVal.size
        if newProb > maxProb :
            if verbose >= 1: print n,"==== new max", i, "       ",newProb , ">", maxProb
            ix = hexVal == newProb
            maxRa  = hexRa[ix]
            maxDec = hexDec[ix]
            maxVal = hexVal[ix]
            maxMjd = hexMjd[ix]
            maxData = maxRa,maxDec,maxVal,maxMjd,i
            maxProb = newProb
    if maxProb == -1 : raise Exception("no max probability found")
    return maxRa, maxDec, maxVal, maxMjd, i

# we've found a hex,slot that can be observed so add it the the observing lists
def addObsToSlot (slotsObserving, maxData, slot) :
    maxRa  = maxData[0]
    maxDec = maxData[1]
    maxVal = maxData[2]
    maxMjd = maxData[3]
    slotsObserving[slot,"ra"]   =  np.append(slotsObserving[slot,"ra"], maxRa)
    slotsObserving[slot,"dec"]  =  np.append(slotsObserving[slot,"dec"], maxDec)
    slotsObserving[slot,"prob"] =  np.append(slotsObserving[slot,"prob"], maxVal)
    slotsObserving[slot,"mjd"]   =  np.append(slotsObserving[slot,"mjd"], maxMjd)
    slotsObserving[slot] += 1
    return slotsObserving
# there can be no more observing in this slot, so this hex,slotj
# is impossible, delete it from the list.
def deleteHexFromSlot (hexData, slot, maxProb) :
    hexRa, hexDec, hexVal, hexMjd = hexData[slot] 
    ix = np.nonzero(hexVal == maxProb) 
    hexRa  = np.delete(hexRa, ix)
    hexDec = np.delete(hexDec, ix)
    hexVal = np.delete(hexVal, ix)
    hexMjd = np.delete(hexVal, ix)
    hexData[slot] = hexRa, hexDec, hexVal, hexMjd
    return hexData
# the hex,slot has made it onto an observing list, so remove the hex
# from all hex lists ( rm hex,*)
def deleteHexFromAllSlots (hexData, observingSlots, maxRa, maxDec, verbose=0, n="") :
    for i in observingSlots:
        data = hexData[i]
        hexRa  = data[0]
        hexDec = data[1]
        hexVal = data[2]
        hexMjd = data[3]
        ix = np.nonzero((hexRa == maxRa) & (hexDec == maxDec))
        if verbose >=4 : print ix, hexRa, maxRa
        if verbose >=2 : 
            ixs = np.shape(ix)[1]
            print n,"bookkeeping",i,"  nHex=",hexRa.size, "   ix.size", ixs, 
        hexRa  = np.delete(hexRa, ix)
        hexDec = np.delete(hexDec, ix)
        hexVal = np.delete(hexVal, ix)
        hexMjd = np.delete(hexMjd, ix)
        hexData[i] = hexRa, hexDec, hexVal, hexMjd

        if verbose >= 2: 
            print "\t after delete nHex=",hexRa.size,
            if ixs > 0 : print "index = ",ix[0][0]
            else: print ""
    return hexData

# it is useful to remove full observing slots from further processing,
# though they survive to be returned in the final lists
def eliminateFullObservingSlots(hexData, slotsObserving, observingSlots, maxHexesPerSlot, verbose) :
    full = []
    for i in range(0,len(observingSlots)):
        # either we've observed as much as we can, or there is nothing to see
        slot = observingSlots[i]
        if (slotsObserving[slot] >= maxHexesPerSlot) | (hexData[slot][0].size ==0): 
            full.append(i)
    if verbose >= 2: 
        print "hiding full observing slot ", 
        for f in full: print observingSlots[f],
        print ""
    observingSlots = np.delete(observingSlots, full)
    return observingSlots

# ==================================
# plotting 
# ==================================

def makeObservingPlots(nslots, slotNumbers, mjd, simNumber, best_slot,
        mapDirectory) :
    import matplotlib
    matplotlib.use("Agg"); # matplotlib.use("TkAgg") ??
    import matplotlib.pyplot as plt
    figure = plt.figure(1,figsize=(5.5*1.618,5.5))


    for i in range(0,nslots) :
        obsTime = ""
        ix = slotNumbers == i
        if np.any(ix) : 
            obsTime = mjd[ix[0]]
            print "making observingPlot-{}.png".format(i)
            plt.clf();
            observingPlot(figure,simNumber,i,mapDirectory,nslots,
                extraTitle=obsTime)
            name = str(simNumber)+"-observingPlot-{}.png".format(i)
            plt.savefig(mapDirectory+name)

# modify mcbryde to have alpha=center of plot
#   "slot" is roughly hour during the night at which to make plot
def observingPlot(figure, simNumber, slot, data_dir, nslots, extraTitle="") :
    import os
    from equalArea import mcbryde
    import matplotlib.pyplot as plt
    import healpy as hp
    import jsonMaker

    # get the planned observations
    ra,dec,prob,mjd,slotNumbers = readObservingRecord(simNumber, data_dir)
    
    # get the maps for a reasonable slot
    name = eventName(data_dir, str(simNumber)) + "-"+str(slot)
    print "\t reading ",name+"-ra.hp  etc"
    raMap     =hp.read_map(name+"-ra.hp");
    decMap    =hp.read_map(name+"-dec.hp");
    #haMap     =hp.read_map(name+"-ha.hp");
    ligoMap   =hp.read_map(name+"-map.hp");
    maglimMap =hp.read_map(name+"-maglim.hp");
    probMap  =hp.read_map(name+"-probMap.hp");
    #total = ligoMap*probMap
    raMap=raMap/(2*np.pi/360.)
    decMap=decMap/(2*np.pi/360.)

    ix = probMap > 0
    medianRA = np.median(raMap[ix])
    decMin = -90.; decMax = 40.
    raMin = medianRA -90.
    raMax = medianRA +90.
    alpha= -1*medianRA

    box=5.
    decMin = dec.min()-box
    decMax = dec.max()+box
    decMid = (decMax-decMin)/2.
    box=box*3
    boxRa = box/np.cos(decMid*2*np.pi/360.)
    raMin = ra.min()-boxRa
    raMax = ra.max()+boxRa
    raMid = (raMax-raMin)/2.
    alpha= -1*raMid

    v1 = np.array([raMin, raMax, raMax, raMin, raMin])
    v2 = np.array([decMin, decMin, decMax, decMax, decMin])
    x,y = mcbryde.mcbryde(v1, v2, alpha=alpha)
    xmin = x.min(); xmax = x.max()
    ymin = y.min(); ymax= y.max()

    xMap,yMap = mcbryde.mcbryde(raMap, decMap, alpha=alpha)
    x,y = mcbryde.mcbryde(ra, dec, alpha=alpha)
    # show me the whole plot  (extremely slow)
    #xmin = xMap.min(); xmax = xMap.max()
    #ymin = yMap.min(); ymax= yMap.max()


    ix=np.nonzero((xMap > xmin) & (xMap  <= xmax) & (yMap > ymin) & (yMap <= ymax) )
    #ix=np.nonzero((xMap > xmin) & (xMap  <= xmax) & (yMap > ymin) & (yMap <= ymax) &(maglimMap>10))


    cmap = "cubehelix_r"
    cmap = "YlGnBu"
    gridsize = 110
    gridsize = 66
    plt.clf();
    plt.hexbin(xMap[ix],yMap[ix],maglimMap[ix],
        vmin=19., vmax=23.1, gridsize=gridsize,cmap=cmap,mincnt=1); 
    plt.colorbar()
        #vmin=20., vmax=23., gridsize=100,cmap=cmap); plt.colorbar()

    plotLigoContours(xMap[ix],yMap[ix], probMap[ix]*ligoMap[ix], whiteLine=True) 

    ax = figure.add_subplot(1,1,1)
    ax=plotDecamHexen(ax, ra,dec,alpha, color="r", lw=1) 
    offsets = jsonMaker.tileOffsets()
    delRa = offsets[8][0]
    delDec = offsets[8][1]
    tdec = dec+delDec
    tra = ra + delRa/np.cos(tdec*2*np.pi/360.)
    #ax=plotDecamHexen(ax, tra,tdec,alpha, color="orange", lw=1) 
    delRa = offsets[9][0]
    delDec = offsets[9][1]
    tdec = dec+delDec
    tra = ra + delRa/np.cos(tdec*2*np.pi/360.)
    #ax=plotDecamHexen(ax, tra,tdec,alpha, color="orange", lw=1) 

    title = "i-band limiting magnitude"
    if extraTitle != "" :
        extraTitle = " mjd {:.2f}: ".format(extraTitle)
        title = extraTitle+title
    plt.title(title)
    plt.xlabel("total prob. countours at prob. max/[1.1, 3, 10, 30]")
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.axes().set_aspect('equal'); 
    #plt.axes().set_frame_on(False); 
    plt.axes().set_xticks([]); 
    plt.axes().set_yticks([])
    plt.show()

def plotDecamHexen(ax, ra,dec,alpha, color="r", lw=1) :
    import decam2hp
    import matplotlib.patches 
    import matplotlib.path 
    from equalArea import mcbryde
    nHex = ra.size
    for i in range(0,nHex) :
        hexRa,hexDec = decam2hp.cameraOutline(ra[i], dec[i])
        hexX,hexY = mcbryde.mcbryde(hexRa,hexDec, alpha=alpha)
        hex_path = matplotlib.path.Path(zip(hexX,hexY))
        hex_patch = matplotlib.patches.PathPatch(hex_path, edgecolor=color, lw=lw, fill=False)
        #hex_patch = matplotlib.patches.PathPatch(hex_path, edgecolor="w", lw=1.5, fill=False)
        ax.add_patch(hex_patch)
        #x,y=mcbryde.mcbryde(tra[i],tdec[i], alpha=alpha)
        #plt.text(x,y,"{}".format(i), ha="center", va="center", color="w")
    return ax

def plotLigoContours(x,y, vals, whiteLine=False) :
    import matplotlib
    import matplotlib.pyplot as plt
    from scipy.interpolate import griddata
    con_levels=5
    max=vals.max()
    levels=[max/30.,max/10.,max/3.,max/1.1]
    print "\t\t contours at max/1.1,3,10,30"

    xmin = x.min(); xmax = x.max()
    ymin = y.min(); ymax = y.max()
    
    coord = np.array(zip(x,y))
    xi=np.linspace(xmin, xmax, 500)
    yi=np.linspace(ymin, ymax, 500)
    xi,yi=np.meshgrid(xi,yi)
    zi = griddata(coord,vals,(xi,yi),method="cubic")
    if whiteLine==True :
        plt.contour(xi,yi,zi,con_levels,linewidths=0.66,colors="w", levels=levels)
    elif whiteLine == "red" :
        plt.contour(xi,yi,zi,con_levels,linewidths=0.66,colors="r", levels=levels)
    else :
        plt.contour(xi,yi,zi,con_levels,linewidths=3,colors="k", levels=levels)
