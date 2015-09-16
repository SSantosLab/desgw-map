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

#========================================================================
#
# main routines: these five are called in Dillon's recycler.py
#
#========================================================================

# ==== prep the observation by calculating necessary maps
#       distance is a problem- the question is whether we want a horizon
#       distance or the estimated distance from LIGO?
#       >>> distance = 75. ;# Mpc, as an estimated horizon distance
def prepare(skymap, mjd, trigger_id, data_dir,
        distance=75, exposure_length=180.) :
    import mapsAtTimeT
    import mags
    import modelRead
    import healpy as hp
    import hp2np
    # === prep the maps
    ligo = hp.read_map(skymap)
    ra,dec,ligo = hp2np.map2np(ligo,256, fluxConservation=True)
    obs = mags.observed(ra,dec,ligo, mjd, verbose=False)
    obs.limitMag("i",exposure=exposure_length)
    # ==== get the neutron star explosion models
    models = modelRead.getModels()
    # ==== calculate maps during a full night of observing
    probs,times = mapsAtTimeT.oneDayOfTotalProbability(obs,mjd,distance,models)
    mapsAtTimeT.probabilityMapSaver (obs, trigger_id, mjd, \
        distance, models, times, probs,data_dir)
    return probs, times

# if the 1% cut isn't in place in mapsAtTimeT.oneDayOfTotalProbability
# then one expects circa 40-45 maps as there is about 2/hour
#   with the 1% cut, then one expects far fewer. Perhaps zero.
def contemplateTheDivisionsOfTime(probs, times, hoursAvailable=6) :
    # if the number of slots is zero, nothing to observe or plot
    if np.size(times) == 0 : return 0,0
    if probs.sum() < 1e-9 : return 0,0
    verbose = 0
    n_slots = findNSlots(hoursAvailable)
    n_maps = times.size
    if verbose: print n_slots, n_maps
    if n_maps == n_slots : 
        mapZero = 0
    elif n_maps < n_slots : 
        mapZero = 0
        n_slots = n_maps
    elif n_maps > n_slots :
        mapZero = findStartMap ( probs, times, n_slots )
    else :
        raise Exception ("no possible way to get here")
    print "======= >>>>>>>> ========== ",
    print "n_maps = {}, n_slots = {}, mapZero = {}, prob_max = {:.6}".format(
        n_maps, n_slots, mapZero, probs.max())
    return n_slots, mapZero

# ==== figure out what to observe
def now(n_slots, mapDirectory="jack/", simNumber=13681, mapZero=0) :
    # if the number of slots is zero, nothing to observe or plot
    if n_slots == 0: 
        return 0,0,0
    # compute the observing schedule
    hoursObserving=observing(simNumber,n_slots,mapDirectory, mapZero=mapZero)
    # print stats to screen
    ra,dec,prob,mjd,slotNumbers,islots = observingStats(hoursObserving)
    # save results to the record
    observingRecord(hoursObserving, simNumber, mapDirectory)
    # write jsons and get slot number  of maximum probability
    maxProb_slot = turnObservingRecordIntoJSONs(
        ra,dec,prob,mjd,slotNumbers, simNumber, mapDirectory) 

    return maxProb_slot

#
# there are possibilities. Show them.
#
def makeObservingPlots(nslots, simNumber, best_slot, mapDirectory) :
    print "================ >>>>>>>>>>>>>>>>>>>>> =================== "
    print "makeObservingPlots(",nslots, simNumber, best_slot,mapDirectory," )"
    print "================ >>>>>>>>>>>>>>>>>>>>> =================== "
    import matplotlib
    matplotlib.use("Agg"); # matplotlib.use("TkAgg") ??
    import matplotlib.pyplot as plt
    figure = plt.figure(1,figsize=(8.5*1.618,8.5))

    # if the number of slots is zero, nothing to observe or plot
    if nslots == 0 : return 0

    # first, make the probability versus something plot
    ra,dec,prob,slotMjd,slotNumbers = readObservingRecord(
        simNumber, mapDirectory)

    probabilityPlot(figure, prob, slotNumbers, simNumber, mapDirectory) 

    # now make the hex observation plots
    counter = 1   ;# already made one
    for i in np.unique(slotNumbers) :
        i = np.int(i)
        obsTime = ""
        ix = slotNumbers == i
        if np.any(ix) : 
            ix = np.nonzero(slotNumbers == i)
            obsTime = slotMjd[ix[0]].mean()
            print "making observingPlot-{}.png".format(i)
            observingPlot(figure,simNumber,i,mapDirectory,nslots,
                extraTitle=obsTime)
            name = str(simNumber)+"-observingPlot-{}.png".format(i)
            plt.savefig(mapDirectory+name)
            counter += 1

    counter+= equalAreaPlot(figure,best_slot,simNumber,mapDirectory)

    # return the number of plots made
    return counter
#
# its a disaster, compute something
#
def nothingToObserveShowSomething(skymap, mjd, exposure_length) :
    import healpy as hp
    import hp2np
    import mags
    import sourceProb
    ligo = hp.read_map(skymap)
    ra,dec,ligo = hp2np.map2np(ligo,256, fluxConservation=True)
    obs = mags.observed(ra,dec,ligo,round(mjd,0), verbose=False)
    obs.limitMag("i",exposure=exposure_length)
    maglim = obs.maglim
    sm=sourceProb.map(obs, lumModel="known")
    sm.calculateProb()
    probMap = sm.probMap
    return ra, dec, ligo, maglim, probMap
#
# no, no, no, we actually can see something: lets see the best plots
#
#   raMap, decMap, ligoMap, maglimMap, probMap, haMap, xMap,yMap, hxMap,hyMap = readMaps(
#   ra, dec, ligo, maglim, prob, ha, x,y hx,hy = readMaps(
def readMaps(data_dir, simNumber, slot) :
    import healpy as hp
    # get the maps for a reasonable slot
    name = eventName(data_dir, str(simNumber)) + "-"+str(slot)
    print "\t reading ",name+"-ra.hp  & etc"
    raMap     =hp.read_map(name+"-ra.hp", verbose=False);
    decMap    =hp.read_map(name+"-dec.hp", verbose=False);
    haMap     =hp.read_map(name+"-ha.hp", verbose=False);
    xMap      =hp.read_map(name+"-x.hp", verbose=False);
    yMap      =hp.read_map(name+"-y.hp", verbose=False);
    hxMap     =hp.read_map(name+"-hx.hp", verbose=False);
    hyMap     =hp.read_map(name+"-hy.hp", verbose=False);
    ligoMap   =hp.read_map(name+"-map.hp", verbose=False);
    maglimMap =hp.read_map(name+"-maglim.hp", verbose=False);
    probMap   =hp.read_map(name+"-probMap.hp", verbose=False);
    haMap=haMap/(2*np.pi/360.)
    raMap=raMap/(2*np.pi/360.)
    decMap=decMap/(2*np.pi/360.)
    return raMap, decMap, ligoMap, maglimMap, probMap, \
        haMap, xMap, yMap, hxMap, hyMap

#========================================================================
# 
# support routines
# 
#========================================================================

def eventName(data_dir, event) :
    name=data_dir+str(event)
    return name

def findNSlots(hoursAvailable, slotDuration=32.) :
    verbose = 0
    if verbose:
        print hoursAvailable
        print hoursAvailable*60./32., round(hoursAvailable*60./32.)
        print int(round(hoursAvailable*60./32.))
    nslots = int(round(hoursAvailable*60./32.))   ;# 32 minutes/slot
    return nslots

# ok, I designed observing for the case
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
    mapNums = np.arange(0,n_maps)
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
    #print "================ >>>>>>>>>>>>>>>>>>>>> =================== "
    #print "================ >>>>>>>>>>>>>>>>>>>>> =================== "
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
        slotsObserving[i,"slotNum"] = np.array([]) ;#non-zero prob slots
        slotsObserving[i,"islot"] = np.array([]) ;# range(0,nslots)

    # read in the hexelated probability data
    hexData = dict()
    for i in observingSlots :
        map_i = i + mapZero
        raHexen, decHexen, hexVal, rank, mjd, slotNum = \
                loadHexalatedProbabilities( sim, map_i, data_dir)
        islot = i*np.ones(raHexen.size)
        print map_i, "map size= {};".format(raHexen.size), 

        impossible = 1e-12
        impossible = 1e-7
        ix = np.nonzero(hexVal < impossible)
        raHexen, decHexen, hexVal, mjd, slotNum, islot  = \
            np.delete(raHexen, ix), \
            np.delete(decHexen, ix), \
            np.delete(hexVal, ix), \
            np.delete(mjd, ix) , \
            np.delete(slotNum, ix), \
            np.delete(islot, ix)
        print " n hexes >{} probability=".format(str(impossible)),
        print "{:4d};".format(raHexen.size),
        print "  sum prob= {:7.4f} %".format( 100*hexVal.sum())
        hexData[i] = raHexen, decHexen, hexVal, mjd, slotNum, islot
        #print np.sort(hexVal), hexVal.sum(), 100.*hexVal.sum(),"%"

    # start the search for all max probabilities
    # we'll assume the list is less than 40,000 long, the n-sq-degrees/sky
    for n in range(0,40000) :
        # search for a single max probabilities
        maxRa, maxDec, maxProb, maxMjd, maxSlotNum, maxIslot  = \
            findMaxProbOfAllHexes(hexData, observingSlots, n, verbose) 
        maxData = maxRa,maxDec,maxProb,maxMjd,maxSlotNum, maxIslot

        # we've found the maximum probability on the lists, 
        # so add it to the obs lists # unless not possible. 
        # If the latter, delete it from that slot
        slot = maxIslot
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
        hexData = deleteHexFromAllSlots (
            hexData, observingSlots, maxRa, maxDec, verbose, n) 

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
        print i, "slotnum={} n obs= {}".format(
            slotsObserving[i,"slotNum"],slotsObserving[i,"ra"].size), 
        print "  sum prob= {:7.4f} %".format( 100*slotsObserving[i,"prob"].sum())
    ra,dec,prob,mjd,slotNum,islot = slotsObservingToNpArrays(slotsObserving) 

    print "observable prob_tot = {:.1f}%".format(100.*prob.sum())
    return ra,dec,prob,mjd,slotNum,islot

def observingRecord(slotsObserving, simNumber, data_dir) :
    name = eventName(data_dir, str(simNumber)) + "-ra-dec-prob-mjd-slot.txt"
    ra,dec,prob,mjd,slotNum,islot = slotsObservingToNpArrays(slotsObserving) 
    data = np.array([ra, dec, prob, mjd, slotNum]).T
    np.savetxt(name, data, "%.6f %.5f %.6f %.4f %d")
    return ra,dec,prob,mjd,slotNum

#     ra,dec,prob,mjd,slotNum,islot = readObservingRecord(simNumber, data_dir)
def readObservingRecord(simNumber, data_dir) :
    import os
    name = eventName(data_dir, str(simNumber)) + "-ra-dec-prob-mjd-slot.txt"
    if not os.path.exists(name) :
        ra,dec,prob,mjd,slotNum = \
            np.array(0),np.array(0),np.array(0), \
            np.array(0)
    else :
        ra,dec,prob,mjd,slotNum = np.genfromtxt(name,unpack=True)
    return ra,dec,prob,mjd,slotNum

def slotsObservingToNpArrays(slotsObserving) :
    nslots = slotsObserving["nslots"]
    ra = np.array([])
    dec = np.array([])
    prob = np.array([])
    mjd = np.array([])
    slotNum = np.array([])
    islot = np.array([])
    for i in range(0,nslots) :
        ra = np.append(ra, slotsObserving[i,"ra"])
        dec = np.append(dec, slotsObserving[i,"dec"])
        prob = np.append(prob, slotsObserving[i,"prob"])
        mjd = np.append(mjd, slotsObserving[i,"mjd"])
        slotNum = np.append(slotNum, slotsObserving[i,"slotNum"])
        islot = np.append(islot, slotsObserving[i,"islot"])
    return ra,dec,prob,mjd,slotNum,islot
    

#===================================================================
#
# Read in all of the hexalated probability files
#
def loadHexalatedProbabilities(sim, slot, data_dir) :
    nameStem = eventName(data_dir, str(sim)) + "-{}".format(str(slot)) 
    name = nameStem + "-hexVals.txt"
    raHexen, decHexen, hexVal, rank, mjd = np.genfromtxt(name, unpack=True, delimiter=",")
    slots = np.ones(raHexen.size)*slot
    return raHexen, decHexen, hexVal, rank, mjd, slots


#
# search for the single highest probability hex over all of the possible hexes
# in the hexData slots 
#
def findMaxProbOfAllHexes(hexData, observingSlots, n="", verbose = 0) :
    maxProb = -1
    for i in observingSlots :
        data = hexData[i]
        hexRa     = data[0]
        hexDec    = data[1]
        hexVal    = data[2]
        hexMjd    = data[3]
        hexMyslot = data[4]
        if hexVal.size == 0: continue
        if verbose >= 2: 
            if i == 2: print n,"====",i, "hexSize =",hexRa.size
        # now check for max prob
        newProb = hexVal.max()
        if verbose >= 4: print n,i, maxProb, ">?", newProb, "     n=",hexVal.size
        if newProb > maxProb :
            if verbose >= 1: print n,"==== new max", i, "       ",newProb , ">", maxProb
            ix = hexVal == newProb
            maxRa     = hexRa[ix]
            maxDec    = hexDec[ix]
            maxVal    = hexVal[ix]
            maxMjd    = hexMjd[ix]
            maxProb   = newProb
            maxSlot   = hexMyslot[ix]
            islot = i
    if maxProb == -1 : 
        raise Exception("no max probability found")
    return maxRa, maxDec, maxVal, maxMjd, maxSlot, islot

# we've found a hex,slot that can be observed so add it the the observing lists
def addObsToSlot (slotsObserving, maxData, slot) :
    maxRa  = maxData[0]
    maxDec = maxData[1]
    maxVal = maxData[2]
    maxMjd = maxData[3]
    maxSlotNum = maxData[4]
    maxIslot = maxData[5]
    slotsObserving[slot,"ra"]   =  np.append(slotsObserving[slot,"ra"], maxRa)
    slotsObserving[slot,"dec"]  =  np.append(slotsObserving[slot,"dec"], maxDec)
    slotsObserving[slot,"prob"] =  np.append(slotsObserving[slot,"prob"], maxVal)
    slotsObserving[slot,"mjd"]   =  np.append(slotsObserving[slot,"mjd"], maxMjd)
    slotsObserving[slot,"slotNum"]   =  np.append(slotsObserving[slot,"slotNum"], maxSlotNum)
    slotsObserving[slot,"islot"]   =  np.append(slotsObserving[slot,"islot"], maxIslot)
    slotsObserving[slot] += 1
    return slotsObserving
# there can be no more observing in this slot, so this hex,slotj
# is impossible, delete it from the list.
def deleteHexFromSlot (hexData, slot, maxProb) :
    hexRa, hexDec, hexVal, hexMjd, hexSlotNum, hexIslot = hexData[slot] 
    ix = np.nonzero(hexVal == maxProb) 
    hexRa  = np.delete(hexRa, ix)
    hexDec = np.delete(hexDec, ix)
    hexVal = np.delete(hexVal, ix)
    hexMjd = np.delete(hexMjd, ix)
    hexSlotNum = np.delete(hexSlotNum, ix)
    hexIslot   = np.delete(hexIslot, ix)
    hexData[slot] = hexRa, hexDec, hexVal, hexMjd, hexSlotNum, hexIslot
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
        hexSlotNum = data[4]
        hexIslot   = data[5]
        ix = np.nonzero((hexRa == maxRa) & (hexDec == maxDec))
        if verbose >=4 : print ix, hexRa, maxRa
        if verbose >=2 : 
            ixs = np.shape(ix)[1]
            print n,"bookkeeping",i,"  nHex=",hexRa.size, "   ix.size", ixs, 
        hexRa  = np.delete(hexRa, ix)
        hexDec = np.delete(hexDec, ix)
        hexVal = np.delete(hexVal, ix)
        hexMjd = np.delete(hexMjd, ix)
        hexSlotNum = np.delete(hexSlotNum, ix)
        hexIslot   = np.delete(hexIslot, ix)
        hexData[i] = hexRa, hexDec, hexVal, hexMjd, hexSlotNum, hexIslot

        if verbose >= 2: 
            print "\t after delete nHex=",hexRa.size,
            if ixs > 0 : print "index = ",ix[0][0]
            else: print ""
    return hexData

# it is useful to remove full observing slots from further processing,
# though they survive to be returned in the final lists
def eliminateFullObservingSlots(
        hexData, slotsObserving, observingSlots, maxHexesPerSlot, verbose) :
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


def turnObservingRecordIntoJSONs(
        ra,dec,prob,mjd,slotNumbers, simNumber, mapDirectory) :
    # write big json file
    name = jsonName(simNumber, mapDirectory)
    jsonMaker.writeJson(ra,dec, jsonFilename=name)

    # write slot json files
    for slot in np.unique(slotNumbers) :
        ix = slotNumbers == slot
        slotMJD = mjd[ix][0]  ;# just get first mjd in this slot
        name = jsonUTCName(slot, slotMJD, simNumber, mapDirectory)
        jsonMaker.writeJson(ra[ix],dec[ix], jsonFilename=name)
        
    # find slot with the maximum probability
    maxProb = -1; maxProb_slot = -1
    for slot in np.unique(slotNumbers) :
        ix = slotNumbers == slot
        sumProb = prob[ix].sum()
        if sumProb > maxProb :
            maxProb = sumProb
            maxProb_slot = slot
    maxProb_slot = np.int(maxProb_slot)

    return maxProb_slot

def jsonName(simNumber, mapDirectory) :
    name = eventName(mapDirectory, str(simNumber)) + ".json"
    return name
def jsonUTCName (slot, mjd, simNumber, mapDirectory) :
    from pyslalib import slalib
    date = slalib.sla_djcl(mjd)
    year = np.int(date[0])
    month= np.int(date[1])
    day = np.int(date[2])
    hour = np.int(date[3]*24.)
    minute = np.int( (date[3]*24.-hour)*60.  )
    time = "UTC-{}-{}-{}-{}:{}:00".format(year,month,day,hour,minute)
    slot = "-{}-".format(np.int(slot))
    name = eventName(mapDirectory, str(simNumber)) + slot + time + ".json"
    return name
     
# ==================================
# plotting 
# ==================================
def probabilityPlot(figure, prob, slotNumbers, simNumber, data_dir) :
    import matplotlib.pyplot as plt
    ax = figure.add_subplot(111)
    plotProb, plotSlot,plotN = np.array([]), np.array([]), np.array([])
    for i in np.unique(slotNumbers) :
        ix = slotNumbers == i
        if prob[ix].sum() > 0 :
            plotN = np.append(plotN, prob[ix].size)
            plotSlot = np.append(plotSlot,i)
            plotProb = np.append(plotProb,100.*prob[ix].sum())
    print "making probabilityPlot.png"
    plt.clf();
    plt.plot(plotSlot,plotProb,c="blue")
    plt.scatter(plotSlot,plotProb,c="red",s=50)
    plt.text(0.80,1.02,"total probability = {:5.1f}%".format(prob.sum()*100.),
        transform = ax.transAxes,   horizontalalignment='left',
        verticalalignment='center',)
    avghex = str( np.round(plotN.mean(),1) )
    plt.text(0.80,0.92,"n hexes per slot: {}".format(avghex),
        transform = ax.transAxes,   horizontalalignment='left',
        verticalalignment='center',)
    plt.ylim(0.0,plt.ylim()[1])
    plt.xlabel("slot number")
    plt.ylabel("probability per slot (%)")
    plt.title("sum(prob*ligo)")
    name = str(simNumber)+"-probabilityPlot.png"
    plt.savefig(data_dir+name)

def equalAreaPlot(figure,slot,simNumber,data_dir) :
    import matplotlib.pyplot as plt
    from equalArea import mcplot
    from equalArea import mcbryde
    import insideDesFootprint

    ra, dec, ligo, maglim, prob, ha, x,y, hx,hy = \
        readMaps(data_dir, simNumber, slot)
    # x,y are the mcbryde projection of ra, dec
    # hx,hy are the mcbryde projection of ha, dec
    ra, dec = x, y

    # des footprint
    desra, desdec = insideDesFootprint.getFootprintRaDec()
    desx, desy = mcbryde.mcbryde(desra, desdec)

    plt.axes().set_aspect('equal')

    name = data_dir+str(simNumber)+str(slot)+"-ligo-eq.png"
    print "making ",name
    plt.clf();mcplot.plot(ra,dec,ligo)
    plt.plot(desx,desy,color="w")
    plt.xlabel("RA");plt.ylabel("Dec")
    plt.savefig(name)

    name = data_dir+str(simNumber)+str(slot)+"-maglim-eq.png"
    print "making ",name
    plt.clf();mcplot.plot(ra,dec,maglim,vmin=17);
    plt.plot(desx,desy,color="w")
    plt.xlabel("RA");plt.ylabel("Dec")
    plt.savefig(name)

    name = data_dir+str(simNumber)+str(slot)+"-prob-eq.png"
    print "making ",name
    plt.clf();mcplot.plot(ra,dec,prob)
    plt.plot(desx,desy,color="w")
    plt.xlabel("RA");plt.ylabel("Dec")
    plt.savefig(name)

    name = data_dir+str(simNumber)+str(slot)+"-probXligo-eq.png"
    print "making ",name
    plt.clf();mcplot.plot(ra,dec,prob*ligo)
    plt.plot(desx,desy,color="w")
    plt.xlabel("RA");plt.ylabel("Dec")
    plt.savefig(name)
    # return the number of plots made
    return 4 

# modify mcbryde to have alpha=center of plot
#   "slot" is roughly hour during the night at which to make plot
def observingPlot(figure, simNumber, slot, data_dir, nslots, extraTitle="") :
    import os
    from equalArea import mcbryde
    import matplotlib.pyplot as plt
    import healpy as hp
    import jsonMaker
    import insideDesFootprint

    # get the planned observations
    ra,dec,prob,mjd,slotNumbers = readObservingRecord(simNumber, data_dir)
    
    # get the maps for a reasonable slot
    raMap, decMap, ligoMap, maglimMap, probMap, \
        haMap, xMap, yMap, hxMap, hyMap = readMaps(data_dir, simNumber, slot)

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
#    ix = slotNumbers == slot
#    ax=plotDecamHexen(ax, ra[ix],dec[ix],alpha, color="r", lw=1) 
    ax=plotDecamHexen(ax, ra,dec,alpha, color="r", lw=1) 
    ix =np.invert( insideDesFootprint.insideFootprint(ra, dec))
    ax=plotDecamHexen(ax, ra[ix],dec[ix],alpha, color="orange", lw=1) 

#    offsets = jsonMaker.tileOffsets()
#    delRa = offsets[8][0]
#    delDec = offsets[8][1]
#    tdec = dec+delDec
#    tra = ra + delRa/np.cos(tdec*2*np.pi/360.)
#    ax=plotDecamHexen(ax, tra,tdec,alpha, color="orange", lw=1) 
#    delRa = offsets[9][0]
#    delDec = offsets[9][1]
#    tdec = dec+delDec
#    tra = ra + delRa/np.cos(tdec*2*np.pi/360.)
#    ax=plotDecamHexen(ax, tra,tdec,alpha, color="orange", lw=1) 

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
