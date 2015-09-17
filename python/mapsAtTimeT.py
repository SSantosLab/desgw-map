import numpy as np
import healpy as hp

import sourceProb
import modelRead

#
# mapsAtTimeT
#
# This file contains routines to compute probability maps
#   for a given burst time and given observation time 
#   in days since the burst.
#
# probabilityMaps
#   Use the mags and sourceProb objects to calculate
#   the limiting magnitude and probability maps at a given time
#   in a given filter and exposure time
#   The models needed as input are found via:
#       models = modelRead.getModels()
#   and then models are ready to roll
#
# totalProbability
#   Return the total probability of in the map at a given time
#   in a given filter and exposure. 
#       This is a thin wrapper about probabilityMaps
#
# manyDaysOfTotalProbability 
#   Return a list of times and total probabilites
#   at deltaTime starting at startOfDays and ending at endOfDays
#   using totalProbabilityAtTimeT
#       The deltaTime parameter sets the duration of a "slot",
#       in the langauge of getHexObservations
#
# oneDayOfTotalProbability 
#   A thin wrapper about manyDaysOfTotalProbability to do one 24 hour period
#
# probabilityMapSaver 
#   Given a list of times and total probabilites,
#   use probabilityMapsAtTimeT (if the probability is high enough)
#   to re-make the maps and save them to disk
#

# Over 11 days
#   calculate sourceProb.map and given a model, sm.calculateProb()
#   sum the probability in the map, append that number to an array,
#   return
#
#   On the subject of deltaTime:
#    deltaTime = 1./24.  ;# once per hour
#    deltaTime = 0.0223  ;# 32 minutes, = 1 hex (i,z,z,i) 180s+30s = 8 minutes 
#       so 4 hexes per 32 minute slot. 
#       More convenient for survey strategy planning (as 8 doesnt go into 60)


def oneDayOfTotalProbability (obs, mjd, distance, models,
        deltaTime=0.0223, probTimeFile="probTime.txt") :

    # the work.
    totalProbs,times = manyDaysOfTotalProbability(
        obs, mjd, distance, models, startOfDays=0,endOfDays=1,
        deltaTime=deltaTime, probTimeFile=probTimeFile) 

    return totalProbs,times

def manyDaysOfTotalProbability (
        obs, mjdOfBurst, distance, models, 
        startOfDays=0, endOfDays=11, deltaTime=0.0223, 
        probTimeFile="probTime.txt") :
    times = []
    totalProbs = []

    # in the language of getHexObservations:
    #   each slot is 32 minutes, each slot can hold 4 hexes
    for time in np.arange(startOfDays,endOfDays,deltaTime) :
        if time < (1.5/24.) : continue
        print "================================== ",
        print "hours since Time Zero: {:.1f}".format(time*24.)
        totalProb = totalProbability(obs, mjdOfBurst, time, distance, models)
        times.append(time)
        totalProbs.append(totalProb)
    totalProbs =np.array(totalProbs)
    times = np.array(times)
    # informational
    print "total all-sky summed probability of detection (list1) and daysSinceBurst (list2)"
    print totalProbs,"\n",times
    print "===== times with total prob > 10**-2"
    ix = totalProbs > 10**-2; 
    print "total all-sky summed probability of detection (list1) and daysSinceBurst (list2)"
    print totalProbs[ix],"\n",times[ix]

    data = np.array([totalProbs[ix], times[ix]]).T
    np.savetxt(probTimeFile, data, "%f %f")
    return totalProbs[ix],times[ix]

def probabilityMaps(obs, mjdOfBurst, daysSinceBurst, distance, models,
        filter="i", exposure=180, ns_model="known") :
    obs.resetTime(mjdOfBurst+daysSinceBurst)
    obs.limitMag(filter, exposure=exposure)

    sm=sourceProb.map(obs, lumModel=ns_model);  
    models_at_t = modelRead.modelsAtTimeT (models, daysSinceBurst)
    abs_mag = models_at_t[0]
    sm.modelAbsoluteMagnitude = abs_mag
    sm.searchDistance = np.array([distance,])
    sm.calculateProb()
    return obs,sm

def totalProbability(obs, mjdOfBurst, daysSinceBurst, distance, models,
        filter="i", exposure=180, ns_model="known") :
    obs,sm = probabilityMaps(obs, mjdOfBurst, daysSinceBurst, distance,
        models, filter, exposure, ns_model)
    totalProb = (obs.map * sm.probMap).sum()
    return totalProb


# for each time in the times,
# calculate the probabilities, 
#   and hexalate
#   then save the maps
#       times,probabilities are the output of manyDaysOfTotalProbability
def probabilityMapSaver (obs, sim, mjd, distance, models, \
        times, probabilities, data_dir) :
    import hexalate
    import os
    gw_data_dir          = os.environ["DESGW_DATA_DIR"]
    hexFile = gw_data_dir + "all-sky-hexCenters.txt"

    counter = -1
    for time,prob  in zip(times, probabilities) :
        counter += 1
        if prob <= 0 : continue
        #print "probabilityMapSaver: counter, time= ", counter, time
        if time < 0.06: time = 0.06 ;# if less than 1.5 hours, set to 1.5 hours
        print "================== map save =====>>>>>>>>===== ",
        print "hours since Time Zero: {:.1f}".format(time*24.)
        obs,sm = probabilityMaps( obs, mjd, time, distance, models)

        raHexen, decHexen, hexVals, rank = hexalate.cutAndHexalate (
            obs, sm, hexFile)

        # obs.ra, obs.dec, obs.map  = ligo map
        # obs.hx, obs.hy   = mcbryde projection of houar angle, dec
        # obs.maglim = limiting mag
        # sm.prob = limiting mag convolve abmag convolve volume
        # sm.probMap = total prob map
        # hexRa,hexDec,hexVals
        nameStem = data_dir + str(sim) + "-{}".format(str(counter)) 

        name = nameStem + "-ra.hp"
        if os.path.exists(name): os.remove(name)
        hp.write_map(name, obs.ra)
        name = nameStem + "-dec.hp"
        if os.path.exists(name): os.remove(name)
        hp.write_map(name, obs.dec)
        name = nameStem + "-ha.hp"
        if os.path.exists(name): os.remove(name)
        hp.write_map(name, obs.ha)
        name = nameStem + "-hx.hp"
        if os.path.exists(name): os.remove(name)
        hp.write_map(name, obs.hx)
        name = nameStem + "-hy.hp"
        if os.path.exists(name): os.remove(name)
        hp.write_map(name, obs.hy)
        name = nameStem + "-x.hp"
        if os.path.exists(name): os.remove(name)
        hp.write_map(name, obs.x)
        name = nameStem + "-y.hp"
        if os.path.exists(name): os.remove(name)
        hp.write_map(name, obs.y)
        name = nameStem + "-map.hp"
        if os.path.exists(name): os.remove(name)
        hp.write_map(name, obs.map)
        name = nameStem + "-maglim.hp"
        if os.path.exists(name): os.remove(name)
        hp.write_map(name, obs.maglim)
        name = nameStem + "-maglim-global.hp"
        if os.path.exists(name): os.remove(name)
        hp.write_map(name, obs.maglimall)
        name = nameStem + "-prob.hp"
        if os.path.exists(name): os.remove(name)
        hp.write_map(name, sm.prob)
        name = nameStem + "-probMap.hp"
        if os.path.exists(name): os.remove(name)
        hp.write_map(name, sm.probMap)

        # where rank is to be understood as the indicies of the
        # ranked hexes in order; i.e., they have nothing to do with
        # raHexen, decHexen, hexVals except as a sorting key
        name = nameStem + "-hexVals.txt"
        if os.path.exists(name): os.remove(name)
        data = np.array([raHexen, decHexen, hexVals, rank, (rank*0)+(mjd+time)])
        np.savetxt(name,data.T,"%.6f, %.5f, %.4e, %d, %.4f")

    

# Get the saved maps for each day and hour.
def readMaps (data_dir, simNumber, slot) :
    name = data_dir + str(simNumber) + "-{}".format(str(slot)) 

    ra=hp.read_map(name+"-ra.hp");
    dec=hp.read_map(name+"-dec.hp");
    ha=hp.read_map(name+"-ha.hp");
    map=hp.read_map(name+"-map.hp");
    maglim=hp.read_map(name+"-maglim.hp");
    prob=hp.read_map(name+"-prob.hp");
    probMap=hp.read_map(name+"-probMap.hp");
    hx=hp.read_map(name+"-hx.hp");
    hy=hp.read_map(name+"-hy.hp");
    x=hp.read_map(name+"-x.hp");
    y=hp.read_map(name+"-y.hp");
    return ra, dec, map, maglim, prob, probMap, x,y, hx,hy

