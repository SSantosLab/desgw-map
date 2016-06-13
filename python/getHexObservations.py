import numpy as np
import jsonMaker
import os
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
# main routines: these seven are called in Dillon's recycler.py
#   prepare
#   contemplateTheDivisionsOfTime   
#   now
#   economics
#   makeObservingPlots      
#   nothingToObserveShowSomething
#   readMaps   
#
#========================================================================

# ==== prep the observation by calculating necessary maps
#       distance is a problem- the question is whether we want a horizon
#       distance or the estimated distance from LIGO?
#       >>> distance = 75. ;# Mpc, as an estimated horizon distance
#
#    deltaTime = 0.0223  => 32 minute slots  for izzi 90sec exp, 4 hex/slot
#    deltaTime = 0.0446  => 62 minute slots  for izzi 90sec exp, 8 hex/slot
#    deltaTime = 0.0417  => 60 minute slots  for izz 90sec exp, 10 hex/slot
#    deltaTime = 0.0208  => 30 minute slots  for izz 90sec exp,  5 hex/slot
#       The first is flexible to observing conditions,
#       while the second is x2 faster as there are less maps to hexalate
#       The third shows a different tiling/filter complement
#
#   skipHexelate reuses existing hexelated maps, the biggest compute time
#   skipAll reuses  hexelated maps and probabilities/times, the 2nd largest time
#   saveHexalationMap saves all maps but the hexes, which it skips computing.
#       skipHexalate really should be named "reuseHexelationMap"
#   doOnlyMaxProbability = True selects the max prob from the list
#       and runs the map saving only on it
#
#   exposure_length is only used in the computation of the limiting mag map
#
#   if start_mjd =0, then the computations start at the burst_mjd,
#   else at start_mjd
#
def prepare(skymap, burst_mjd, trigger_id, data_dir, mapDir,
        distance=60., exposure_list = [90,], filter_list=["i",],
        overhead=30., maxHexesPerSlot=6,
        start_mjd = 0, skipHexelate=False, skipAll=False, 
        onlyHexesAlreadyDone="", 
        saveHexalationMap=True, doOnlyMaxProbability=False, resolution=256) :
    import mapsAtTimeT
    import mags
    import modelRead
    import healpy as hp
    import hp2np
    import os
    debug = 0
    if start_mjd == 0: start_mjd = burst_mjd

    exposure_list = np.array(exposure_list)
    filter_list = np.array(filter_list)
    ix = filter_list == "i"
    exposure_length = exposure_list[ix].sum()

    answers = slotCalculations( burst_mjd, exposure_list, overhead, 
        nHexes=maxHexesPerSlot) 
    hoursPerNight = answers["hoursPerNight"] ;# in minutes
    slotDuration = answers["slotDuration"] ;# in minutes
    deltaTime = slotDuration/(60.*24.) ;# in days

    probabilityTimesCache = os.path.join(data_dir,\
        "probabilityTimesCache_"+str(trigger_id)+".txt")
    if skipAll and not os.path.isfile(probabilityTimesCache) :
        print "=============>>>> forced to calculate probs as cache file nonexistent"
        skipAll = False
        skipHexelate = True

    if skipAll :
        print "=============>>>> ",
        print "prepare: using cached probabilities, times, and maps"
        print "\t reading ",probabilityTimesCache
        if os.stat(probabilityTimesCache).st_size == 0 :
            probs, times = np.array([0,]),np.array([0,])
            print "\t got nothing- we're skipping this one"
        else :
            data = np.genfromtxt(probabilityTimesCache, unpack=True)
            probs, times = data[0],data[1]
        return probs, times, slotDuration, hoursPerNight
        
    # ==== get the neutron star explosion models
    models = modelRead.getModels()

    # === prep the maps
    ligo = hp.read_map(skymap)
    ra,dec,ligo = hp2np.map2np(ligo, resolution, fluxConservation=True)
    obs = mags.observed(ra,dec,ligo, start_mjd, verbose=False)
    obs.limitMag("i",exposure=exposure_length)

    # ==== calculate maps during a full night of observing
    probs,times = mapsAtTimeT.oneDayOfTotalProbability(
        obs, burst_mjd, distance, models, 
        deltaTime=deltaTime, start_mjd= start_mjd,
        probTimeFile= probabilityTimesCache )
    if skipHexelate:
        print "=============>>>> prepare: using cached maps"
        return probs, times, slotDuration
    if debug :
        return  obs, trigger_id, mjd, distance, models, times, probs,data_dir
    if doOnlyMaxProbability :
        if len(probs) == 0 : return [0,],[0,],[0,],[0,]
        ix = np.argmax(probs)
        probs = [probs[ix],]
        times = [times[ix],]

    mapsAtTimeT.probabilityMapSaver (obs, trigger_id, burst_mjd, \
        distance, models, times, probs,mapDir, \
        onlyHexesAlreadyDone=onlyHexesAlreadyDone, 
        performHexalatationCalculation=saveHexalationMap)
    return probs, times, slotDuration, hoursPerNight

# ========== do simple calculations on how to divide the night
#
#   one of the questions is how many hours to devote to observations
#       hoursAvailable,  another is the slot duration
#
# if the 1% cut isn't in place in mapsAtTimeT.oneDayOfTotalProbability
# then one expects circa 40-45 maps as there is about 2/hour
#   with the 1% cut, then one expects far fewer. Perhaps zero.
#
def contemplateTheDivisionsOfTime(
        probs, times, slotDuration=30., 
            hoursPerNight= 10., hoursAvailable=6) :
    if hoursAvailable > hoursPerNight:
        hoursAvailable = hoursPerNight
        
    # if the number of slots is zero, nothing to observe or plot
    if np.size(times) == 0 : return 0,0
    if probs.sum() < 1e-9 : return 0,0
    verbose = 0
    n_slots = findNSlots(hoursAvailable,slotDuration=slotDuration)
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
    print "=============>>>>  contemplateTheDivisionsOfTime:"
    print "\t n_maps = {}, n_slots = {}, mapZero = {}, prob_max = {:.6}".format(
        n_maps, n_slots, mapZero, probs.max())
    return n_slots, mapZero

# ==== figure out what to observe
#
# Another fast routine
#   basic for this routine is how many hexes per slot
#
#   skipJson does just that, speeding the routine up considerably
#
#   if a doneFile is given, which should be
#   something with ra,dec in columns 1,2, like 
#       G184098-ra-dec-prob-mjd-slot.txt
#   then those ra,decs are interpreted as hex centers,
#   and hexes in the hexVals maps within 36" of the ra,decs
#   are removed- they are done, the question is what to observe next
#
def now(n_slots, mapDirectory="jack/", simNumber=13681, 
        mapZero=0, maxHexesPerSlot=5, 
        exposure_list = [90,90,90], filter_list=["i","z","z"],
        doneFile = "", skipJson = False ) :
    # if the number of slots is zero, nothing to observe or plot
    if n_slots == 0: 
        return 0
    # compute the observing schedule
    print "=============>>>>  now: observing"
    hoursObserving=observing(
        simNumber,n_slots,mapDirectory, mapZero=mapZero,
        maxHexesPerSlot = maxHexesPerSlot, doneFile=doneFile)
    # print stats to screen
    print "=============>>>>  now: observingStats"
    ra,dec,prob,mjd,slotNumbers,islots = observingStats(hoursObserving)
    # save results to the record
    observingRecord(hoursObserving, simNumber, mapDirectory)
    # write jsons and get slot number  of maximum probability
    maxProb_slot = maxProbabilitySlot(prob,slotNumbers)
    if not skipJson :
        print "=============>>>>  now: JSON"
        turnObservingRecordIntoJSONs(
            ra,dec,prob,mjd,slotNumbers, simNumber, 
            exposure_list=exposure_list, filter_list=filter_list, 
            mapDirectory=mapDirectory) 

    return maxProb_slot

# ===== The economics analysis
#
#   area_left is th enumber of hexes we have left to observe this season
#   days_left is the number of days left in the season
#   rate is the effective rate of triggers
#       p_gw is that for which the table cumul_table_pgw50.txt was  made.
#
def economics (simNumber, best_slot, mapDirectory, 
        area_left=200., days_left=60., rate=1/30.,  p_gw = 0.10) :
    import healpy as hp
    import cumul
    import des_optimization
    import os
    gw_data_dir = os.environ["DESGW_DATA_DIR"]
    ra, dec, ligo, maglim, prob, ha, x,y, hx,hy = \
        readMaps(mapDirectory, simNumber, best_slot)

    area_bar_p,area_bar = np.genfromtxt(
        gw_data_dir+"/area_bar_table.txt",unpack=True)
    avge_cumu_area,avge_cumu = np.genfromtxt(
        gw_data_dir+"/cumul_table_pgw10.txt",unpack=True)

    obsProb = ligo*prob
    nsides = hp.get_nside(obsProb)
    # max area viewable by Blanco at one time is 11734. sq-degrees
    max_area=11734.
    area, cum_prob  = cumul.area(ra,dec,obsProb, p_gw, nsides, max_area=max_area)
    area_to_cover_p_gw = area
    #print avge_cumu_area
    #print area
    ix = np.searchsorted(avge_cumu_area, area)
    if ix >= avge_cumu_area.size :
        fraction_of_sims_better_than_this_trigger = 1.0
    else :
        fraction_of_sims_better_than_this_trigger = avge_cumu[ix]

    prob, N_max = des_optimization.evaluate_average_event(
        area_left, days_left, rate, avge_cumu, avge_cumu_area, area_bar, area_bar_p)

    if fraction_of_sims_better_than_this_trigger < 1./N_max :
        area, cum_prob = cumul.area(ra,dec,obsProb, prob, nsides, max_area=max_area)
        if area>area_left:
            print "\t maxing out area: \t {:.3f} -> ".format( cum_prob),
            cum_prob = cumul.probability_covered(ra,dec,obsProb, area_left, nsides, max_area=max_area)
            print "{:.3f}".format(cum_prob)
            area=area_left
    else :
        print "\t ignore event"
        area = 0
        prob = 0

    probability_covered = cum_prob
    quality = fraction_of_sims_better_than_this_trigger 
    return probability_covered, area, area_to_cover_p_gw, quality

#
# ====== there are possibilities. Show them.
#
def makeObservingPlots(nslots, simNumber, best_slot, data_dir, mapDirectory) :
    print "================ >>>>>>>>>>>>>>>>>>>>> =================== "
    print "makeObservingPlots(",nslots, simNumber, best_slot,data_dir," )"
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

    probabilityPlot(figure, prob, slotNumbers, simNumber, data_dir) 

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
            observingPlot(figure,simNumber,i,mapDirectory, nslots,
                extraTitle=obsTime)
            name = str(simNumber)+"-observingPlot-{}.png".format(i)
            plt.savefig(os.path.join(mapDirectory,name))
            counter += 1
            counter+= equalAreaPlot(figure,i,simNumber,data_dir,mapDirectory)

    #counter+= equalAreaPlot(figure,best_slot,simNumber,data_dir,mapDirectory)

    # return the number of plots made
    return counter
#
# ===== its a disaster, compute something
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
#   ra, dec, ligo, maglim, prob, ha, x,y, hx,hy = readMaps(
def readMaps(mapDir, simNumber, slot) :
    import healpy as hp
    # get the maps for a reasonable slot
    name = os.path.join(mapDir, str(simNumber) + "-"+str(slot))
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
# for economics analysis
def time_cost_per_hex (nvisits, overhead, exposure_length) :
    tot_exptime = (np.array(overhead)+np.array(exposure_length)).sum
    time_cost_per_hex = nvisits * tot_exptime #sec
    return time_cost_per_hex
    
# for economics analysis
def area_left (area_per_hex, time_budget, time_cost_per_hex) :
    area_per_hex * (time_budget * 3600)/(time_cost_per_hex)
    return area_per_hex
#time_cost_per_hex = nvisits * nexposures * (overhead + exposure_length) #sec
#area_left =  area_per_hex * (time_budget * 3600)/(time_cost_per_hex)

# place holder for the code brought from desisurvey...
def hoursPerNight (mjd) :
    import mags
    night = mags.findNightDuration(mjd)
    night = night*24.
    return night


# These calculations used to be spread over hither and yon.
# Bring them together.
#
#    deltaTime = 0.0223  => 32 minute slots  for izzi 90sec exp, 4 hex/slot
#    deltaTime = 0.0446  => 62 minute slots  for izzi 90sec exp, 8 hex/slot
#    deltaTime = 0.0417  => 60 minute slots  for izz 90sec exp, 10 hex/slot
#    deltaTime = 0.0208  => 30 minute slots  for izz 90sec exp,  5 hex/slot
#       The first is flexible to observing conditions,
#       while the second is x2 faster as there are less maps to hexalate
#       The third shows a different tiling/filter complement

# Let us redfine this: a slot is the length of time it takes
# to do 6 hexes to completion. That is usually somewhere between 30 minutes
# and one hour, so close to the original defintion, and by force is an
# even number of hexes. Ok. Use n=6 for the forcing definition
def slotCalculations(mjd, exposure_lengths, overhead, nHexes = 6) :
    tot_exptime = (np.array(overhead)+np.array(exposure_lengths)).sum()
    slot_time = tot_exptime*nHexes
    slot_duration = slot_time/60. ;# in minutes
    hoursAvailable = hoursPerNight(mjd)
    answers = dict()
    answers["slotDuration"] = slot_duration
    answers["hoursPerNight"] = hoursAvailable
    return answers

#SHOULD NOT BE NEEDED ANYMORE - DILLON DEC 7TH
'''def eventName(data_dir, event) :
    name=os.path.join(data_dir,str(event))
    return name
'''
# find the number of slots per night
def findNSlots(hoursAvailable, slotDuration=32.) :
    verbose = 0
    if verbose:
        print hoursAvailable
        print hoursAvailable*60./slotDuration, round(hoursAvailable*60./slotDuration)
        print int(round(hoursAvailable*60./slotDuration))
    nslots = int(round(hoursAvailable*60./slotDuration))   ;# 32 minutes/slot
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
    return bestStart


# Load all n -hexVals files,
#   pick the highest probability one, 
#   put it into one of the n slots, unless that slots is maxed out
#   remove that hex from all the hexVals lists
#   do it again, untill all n time slots are full.
#       maxHexesPerSlot=4 comes from 32 minute duration slots
#       and 8 minutes/hex (izzi 2 min/image)
#
#   if we do zzi at 2 mins/image then 4 min/hex + 2 min/hex2 = 6 mins
#   call it 60 minute slots  and 10 hexes/slot
def observing(sim, nslots, data_dir, 
        maxHexesPerSlot = 4, mapZero = 0, verbose=0,
        doneFile = "jack2/G184098-ra-dec-prob-mjd-slot.txt") :
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
        if doneFile != "" :
            raHexen, decHexen, hexVal, rank, mjd, slotNum = \
                eliminateHexesAlreadyDone(doneFile, 
                raHexen, decHexen, hexVal, rank, mjd, slotNum )
        islot = i*np.ones(raHexen.size)
        print "\t", map_i, "map size= {};".format(raHexen.size), 

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
            print "\t======================================== "
            if verbose >= 1: 
                print "n slots =", len(observingSlots)," == 0?"
                print "sumHexes = ", sumHexes, "==? 0"
            print "\tnumber of hexes observed = ", sumObs
            print "\t======================================== "
            return slotsObserving 

        # otherwise, go back and do it again
    
    # we've done everything on the lists, we can observe it all,
    # return this great success that will never be reached.
    return slotsObserving
    
#
# examine the statistics of the observing lists
#
def observingStats( slotsObserving ) :
    nslots = slotsObserving["nslots"]
    for i in range(0,nslots) :
        print "\t",i, 
        #print "slotnum={} ".format( slotsObserving[i,"slotNum"]),
        print "n obs= {}".format( slotsObserving[i,"ra"].size), 
        print "  sum prob= {:7.4f} %".format( 100*slotsObserving[i,"prob"].sum())
    ra,dec,prob,mjd,slotNum,islot = slotsObservingToNpArrays(slotsObserving) 

    print "\tobservingStats:  ",
    print "observable prob_tot = {:.1f}%".format(100.*prob.sum())
    return ra,dec,prob,mjd,slotNum,islot

def observingRecord(slotsObserving, simNumber, data_dir) :
    name = os.path.join(data_dir, str(simNumber) + "-ra-dec-prob-mjd-slot.txt")
    ra,dec,prob,mjd,slotNum,islot = slotsObservingToNpArrays(slotsObserving) 
    data = np.array([ra, dec, prob, mjd, slotNum]).T
    np.savetxt(name, data, "%.6f %.5f %.6f %.4f %d")
    return ra,dec,prob,mjd,slotNum

#     ra,dec,prob,mjd,slotNum,islot = readObservingRecord(simNumber, data_dir)
def readObservingRecord(simNumber, data_dir) :
    import os
    name = os.path.join(data_dir, str(simNumber) + "-ra-dec-prob-mjd-slot.txt")
    if not os.path.exists(name) :
        ra,dec,prob,mjd,slotNum = \
            np.array(0),np.array(0),np.array(0), \
            np.array(0),np.array(0)
    else :
        ra,dec,prob,mjd,slotNum = np.genfromtxt(name,unpack=True,comments="#")
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
    nameStem = os.path.join(data_dir, str(sim) + "-{}".format(str(slot)))
    name = nameStem + "-hexVals.txt"
    raHexen, decHexen, hexVal, rank, mjd = np.genfromtxt(name, unpack=True, delimiter=",")
    slots = np.ones(raHexen.size)*slot
    return raHexen, decHexen, hexVal, rank, mjd, slots

# ra,dec of the already done file should be in cols 1,2 (1-based)
# as in G184098-ra-dec-prob-mjd-slot.txt
def eliminateHexesAlreadyDone(infile, ra, dec, hexVal, rank, mjd, slots) :
    doneRa, doneDec = np.genfromtxt(infile, unpack=True, usecols=(0,1))
    mask = np.ones(len(ra), dtype=bool)
    for i in range(0,doneRa.size) :
        cosdec = np.cos(doneDec[i]*np.pi/180)
        # within 36" of each other
        ix = ( ((doneRa[i]-ra)/cosdec)**2 + (doneDec[i]-dec)**2) < 0.01
        mask[ix] = False
    return ra[mask], dec[mask], hexVal[mask], rank[mask], mjd[mask], slots[mask]

def onlyReturnHexesDone(hexRa, hexDec, ra, dec, vals) :
    mask = np.zeros(len(ra), dtype=bool)
    for i in range(0,hexRa.size) :
        cosdec = np.cos(hexDec[i]*np.pi/180)
        # within 36" of each other
        ix = ( ((hexRa[i]-ra)/cosdec)**2 + (hexDec[i]-dec)**2) < 0.01
        mask[ix] = True
    return ra[mask], dec[mask], vals[mask]
    
    

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
        ra,dec,prob,mjd,slotNumbers, simNumber, 
        exposure_list, filter_list, mapDirectory) :
    seqtot =  ra.size
    seqzero = 0

    # write slot json files
    for slot in np.unique(slotNumbers) :
        ix = slotNumbers == slot
        slotMJD = mjd[ix][0]  ;# just get first mjd in this slot
        tmpname, name = jsonUTCName(slot, slotMJD, simNumber, mapDirectory)
        jsonMaker.writeJson(ra[ix],dec[ix], 
            simNumber, seqzero, seqtot, exposureList= exposure_list, 
            filterList= filter_list, jsonFilename=tmpname)

        desJson(tmpname, name, mapDirectory) 
        seqzero =+ ra[ix].size
        
def maxProbabilitySlot(prob,slotNumbers) :
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


# verbose can be 0, 1=info, 2=debug
def desJson(tmpname, name, data_dir, verbose = 1) :
    import os
    import logging
    from collections import defaultdict
    import gwwide
    gw_data_dir          = os.environ["DESGW_DATA_DIR"]
    des_json  = gw_data_dir + "all_wide_sispi_queue.json"
    log_levels = defaultdict(lambda:logging.WARNING)
    log_levels[1] = logging.INFO
    log_levels[2] = logging.DEBUG
    logging.basicConfig(filename=data_dir+"json-conversion.log",
        format='%(asctime)s %(message)s', level=verbose)

    logging.info("Begin")
    gwwide.file_gwwide(tmpname, des_json, name)

def jsonUTCName (slot, mjd, simNumber, mapDirectory) :
    from pyslalib import slalib
    date = slalib.sla_djcl(mjd)
    year = np.int(date[0])
    month= np.int(date[1])
    day = np.int(date[2])
    hour = np.int(date[3]*24.)
    minute = np.int( (date[3]*24.-hour)*60.  )
    time = "UTC-{}-{}-{}-{}:{}:00".format(year,month,day,hour,minute)
    tmpname, name = jsonName(slot, time, simNumber, mapDirectory)

    return tmpname, name
def jsonName (slot, utcString, simNumber, mapDirectory) :
    slot = "-{}-".format(np.int(slot))
    tmpname = os.path.join(mapDirectory, str(simNumber) + slot + utcString + "-tmp.json")
    name = os.path.join(mapDirectory, str(simNumber) + slot + utcString + ".json")
    return tmpname, name

def jsonFromRaDecFile(radecfile, nslots, slotZero, 
        hexesPerSlot, simNumber, mjdList, data_dir) :
    ra,dec = np.genfromtxt(radecfile, unpack=True,usecols=(0,1),comments="#")

    seqtot =  ra.size
    seqzero = 0

    # instead, just reorder the ra,dec before feeding to this routine
    #ix = np.argsort(ra)

    counter = 0
    slot = slotZero
    slotRa = np.array([])
    slotDec = np.array([])
    for i in range(0,ra.size) :
        slotRa = np.append(slotRa, ra[i])
        slotDec = np.append(slotDec, dec[i])
        counter += 1
        if counter == hexesPerSlot :
            tmpname, name = jsonName(slot, mjdList[slot-slotZero], 
                simNumber,data_dir)
            jsonMaker.writeJson(slotRa,slotDec, 
                simNumber, seqzero+(hexesPerSlot*(slot-slotZero)), 
                seqtot, jsonFilename=tmpname)
            desJson(tmpname, name, data_dir) 
            counter = 0
            slot += 1
            slotRa = np.array([])
            slotDec = np.array([])
    if counter > 0 :
        tmpname, name = jsonName(slot, mjdList[slot-slotZero], 
            simNumber,data_dir)
        jsonMaker.writeJson(slotRa,slotDec, 
            simNumber, seqzero+(hexesPerSlot*(slot-slotZero)), 
            seqtot, jsonFilename=tmpname)
        desJson(tmpname, name, data_dir) 
        
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
    plt.savefig(os.path.join(data_dir,name))

def equalAreaPlot(figure,slot,simNumber,data_dir,mapDir) :
    import matplotlib.pyplot as plt
    from equalArea import mcplot
    from equalArea import mcbryde
    import insideDesFootprint

    ra, dec, ligo, maglim, prob, ha, x,y, hx,hy = \
        readMaps(mapDir, simNumber, slot)
    # x,y are the mcbryde projection of ra, dec
    # hx,hy are the mcbryde projection of ha, dec
    ra, dec = x, y

    # des footprint
    desra, desdec = insideDesFootprint.getFootprintRaDec()
    desx, desy = mcbryde.mcbryde(desra, desdec)

    plt.axes().set_aspect('equal')

    name = os.path.join(data_dir,str(simNumber)+"-"+str(slot)+"-ligo-eq.png")
    print "making ",name
    plt.clf();mcplot.plot(ra,dec,ligo)
    plt.plot(desx,desy,color="w")
    plt.xlabel("RA");plt.ylabel("Dec")
    plt.savefig(name)

    name = os.path.join(data_dir,str(simNumber)+"-"+str(slot)+"-maglim-eq.png")
    print "making ",name
    plt.clf();mcplot.plot(ra,dec,maglim,vmin=17);
    plt.plot(desx,desy,color="w")
    plt.xlabel("RA");plt.ylabel("Dec")
    plt.savefig(name)

    name = os.path.join(data_dir,str(simNumber)+"-"+str(slot)+"-prob-eq.png")
    print "making ",name
    plt.clf();mcplot.plot(ra,dec,prob)
    plt.plot(desx,desy,color="w")
    plt.xlabel("RA");plt.ylabel("Dec")
    plt.savefig(name)

    name = os.path.join(data_dir,str(simNumber)+"-"+str(slot)+"-probXligo-eq.png")
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
    import plotMapAndHex

    # get the planned observations
    ra,dec,prob,mjd,slotNumbers = readObservingRecord(simNumber, data_dir)
    
    title = "i-band limiting magnitude"
    if extraTitle != "" :
        extraTitle = " mjd {:.2f}: ".format(extraTitle)
        title = extraTitle+title
    title = title + "      LIGO countours at max/[1.1, 3, 10, 30]"


    print "plotMapAndHex.mapAndHex(figure, ", simNumber, ",", slot, ",", data_dir, ",", nslots, ",ra,dec,", title,") "
    d=plotMapAndHex.mapAndHex(figure, simNumber, slot, data_dir, nslots, ra, dec, title) 
    return d


