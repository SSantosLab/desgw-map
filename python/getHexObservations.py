import numpy as np
import jsonMaker
import os
import obsSlots
#
#       the routine mapsAtTimeT.oneDayOfTotalProbability
#           breaks the night into slots of 32 minutes duration
#       the routine mapsAtTimeT.probabilityMapSaver
#           makes maps for each of these slots
#       the routine getHexObservations.observing
#           once told how many slots can be used
#           returns "hoursObserving" containing ra,dec,prob,time
#               of the hexes to observe each slot
#       the routine obsSlots.observingStats
#           places the ra,dec,id,prob,times per slot onto single  lists
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
        this_tiling = "", reject_hexes = "",
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

    answers = obsSlots.slotCalculations( burst_mjd, exposure_list, overhead, 
        hexesPerSlot=maxHexesPerSlot) 
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
        this_tiling = this_tiling, reject_hexes = reject_hexes,
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
    n_slots = obsSlots.findNSlots(hoursAvailable,slotDuration=slotDuration)
    n_maps = times.size
    if verbose: print n_slots, n_maps
    if n_maps == n_slots : 
        mapZero = 0
    elif n_maps < n_slots : 
        mapZero = 0
        n_slots = n_maps
    elif n_maps > n_slots :
        mapZero = obsSlots.findStartMap ( probs, times, n_slots )
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
#   this_tiling and reject_hexes are inverse in thier effect:
#   only hexes within 36" of the ra,decs of the centers of this_tiling
#   are to be observed, while hexes wihtin within 36" of the ra,decs
#   of the centers of reject_hexes are skipped
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
        trigger_type = "NS", skipJson = False ) :
    # if the number of slots is zero, nothing to observe or plot
    if n_slots == 0: 
        return 0
    # compute the observing schedule
    print "=============>>>>  now: observing"
    hoursObserving=obsSlots.observing(
        simNumber,n_slots,mapDirectory, mapZero=mapZero,
        maxHexesPerSlot = maxHexesPerSlot)
    # print stats to screen
    print "=============>>>>  now: observingStats"
    ra,dec,id,prob,mjd,slotNumbers,islots = obsSlots.observingStats(hoursObserving)
    # save results to the record
    obsSlots.observingRecord(hoursObserving, simNumber, mapDirectory)
    # write jsons and get slot number  of maximum probability
    maxProb_slot = obsSlots.maxProbabilitySlot(prob,slotNumbers)
    if not skipJson :
        print "=============>>>>  now: JSON"
        turnObservingRecordIntoJSONs(
            ra,dec,id,prob,mjd,slotNumbers, simNumber, 
            exposure_list=exposure_list, filter_list=filter_list, 
            trigger_type=trigger_type, mapDirectory=mapDirectory) 

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
    ra,dec,id,prob,slotMjd,slotNumbers = obsSlots.readObservingRecord(
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
def nothingToObserveShowSomething(simNumber, data_dir, mapDir) :
    import matplotlib.pyplot as plt
    figure = plt.figure(1,figsize=(8.5*1.618,8.5))
    slot = 0
    counter = equalAreaPlot(figure,slot,simNumber,data_dir,mapDir) 
    return clounter
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


#===================================================================
#
    
def turnObservingRecordIntoJSONs(
        ra,dec,id,prob,mjd,slotNumbers, simNumber, 
        exposure_list, filter_list, trigger_type, mapDirectory) :
    seqtot =  ra.size
    seqzero = 0

    # write slot json files
    for slot in np.unique(slotNumbers) :
        ix = slotNumbers == slot
        slotMJD = mjd[ix][0]  ;# just get first mjd in this slot
        tmpname, name = jsonUTCName(slot, slotMJD, simNumber, mapDirectory)
        jsonMaker.writeJson(ra[ix],dec[ix],id[ix],
            simNumber, seqzero, seqtot, exposureList= exposure_list, 
            filterList= filter_list, trigger_type=trigger_type, jsonFilename=tmpname)

        desJson(tmpname, name, mapDirectory) 
        seqzero =+ ra[ix].size
        
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
    time = utcFromMjd(mjd)
    tmpname, name = jsonName(slot, time, simNumber, mapDirectory)
    return tmpname, name

def utcFromMjd (mjd) :
    from pyslalib import slalib
    date = slalib.sla_djcl(mjd)
    year = np.int(date[0])
    month= np.int(date[1])
    day = np.int(date[2])
    hour = np.int(date[3]*24.)
    minute = np.int( (date[3]*24.-hour)*60.  )
    time = "UTC-{}-{}-{}-{}:{}:00".format(year,month,day,hour,minute)
    return time

def jsonName (slot, utcString, simNumber, mapDirectory) :
    slot = "-{}-".format(np.int(slot))
    tmpname = os.path.join(mapDirectory, str(simNumber) + slot + utcString + "-tmp.json")
    name = os.path.join(mapDirectory, str(simNumber) + slot + utcString + ".json")
    return tmpname, name

def jsonFromRaDecFile(radecfile, nslots, slotZero, 
        hexesPerSlot, simNumber, mjdList, trigger_type, data_dir) :
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
                seqtot, trigger_type, jsonFilename=tmpname)
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
            seqtot, trigger_type, jsonFilename=tmpname)
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
    ra,dec,id,prob,mjd,slotNumbers = obsSlots.readObservingRecord(simNumber, data_dir)
    
    title = "i-band limiting magnitude"
    if extraTitle != "" :
        #extraTitle = " mjd {:.2f}: ".format(extraTitle)
        extraTitle = " {} ".format(utcFromMjd(mjd))
        title = extraTitle+title
    title = title + "      LIGO countours at max/[1.1, 3, 10, 30]"


    print "plotMapAndHex.mapAndHex(figure, ", simNumber, ",", slot, ",", data_dir, ",", nslots, ",ra,dec,", title,") "
    d=plotMapAndHex.mapAndHex(figure, simNumber, slot, data_dir, nslots, ra, dec, title, slots=slotNumbers) 
    return d


