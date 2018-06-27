import numpy as np
import healpy as hp
import hp2np
import sourceProb
import mags
import modelRead
import mapsAtTimeT

#
# trigger_type = "NS" or "BH"
#
def mainInjector (trigger_id, skymap, trigger_type,\
        outputDir, recycler_days_since_burst=0.3, 
        hexID_to_do=[], hexID_to_reject=[], 
        hours_used_by_NS=0, hours_used_by_BH=0,
        halfNight = False, firstHalf= True,
        quick=False, debug=True, camera='decam') :

    import os
    import yaml
    import getHexObservations
    skipEcon = False
    skipEcon = True

    with open("maininjector.yaml","r") as f:
        config = yaml.safe_load(f); 
    
    # debug
    debug = config["debug"]    

    # camera
    camera   = config["camera"]

    #resolution
    resolution = config["resolution"]

    # season parameters
    start_of_season   = config["start_of_season"]
    end_of_season     = config["end_of_season"]

    # static observation details
    overhead          = config["overhead"]
    area_per_hex      = config["area_per_hex"]

    # strategy
    exposure_length_ns= config["exposure_length_NS"]
    filter_list_ns    = config["exposure_filter_NS"]
    maxHexesPerSlot_ns= config["maxHexesPerSlot_NS"]
    exposure_length_bh= config["exposure_length_BH"]
    filter_list_bh    = config["exposure_filter_BH"]
    maxHexesPerSlot_bh= config["maxHexesPerSlot_BH"]

    # economics analysis for NS and for BH
    hoursAvailable_ns = config["time_budget_for_NS"]
    hoursAvailable_bh = config["time_budget_for_BH"]
    lostToWeather_ns  = config["hours_lost_to_weather_for_NS"]
    lostToWeather_bh  = config["hours_lost_to_weather_for_BH"]
    rate_bh           = config["rate_of_bh_in_O2"];# events/year
    rate_ns           = config["rate_of_ns_in_O2"];# events/year


    # configure strategy for the event type
    if trigger_type == "NS" :
        hoursAvailable       = hoursAvailable_ns - lostToWeather_ns - hours_used_by_NS
        rate                 = rate_ns
        exposure_length      = exposure_length_ns
        filter_list          = filter_list_ns
        maxHexesPerSlot      = maxHexesPerSlot_ns
        nepochs           = config["nepochs_NS"]
        epochs = np.array([])
        for i in range (1,nepochs+1) :
            epochs        = np.append(epochs, config["epoch{}_NS".format(i)])
        end_date          = config["enddate_NS"]
        distance          = -999
    elif trigger_type == "BH" :
        hoursAvailable       = hoursAvailable_bh - lostToWeather_bh - hours_used_by_BH
        rate                 = rate_bh
        exposure_length      = exposure_length_bh
        filter_list          = filter_list_bh 
        maxHexesPerSlot      = maxHexesPerSlot_bh
        nepochs           = config["nepochs_BH"]
        epochs = np.array([])
        for i in range (1,nepochs+1) :
            epochs        = np.append(epochs, config["epoch{}_BH".format(i)])
        end_date          = config["enddate_BH"]
        distance          = 1.0
    else :
        raise Exception(
            "trigger_type={}  ! Can only compute BH or NS".format(trigger_type))
    exposure_length   = np.array(exposure_length)
    
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
    
    if camera == "hsc":
        overhead = 20.
        area_per_hex = 1.5

    # make the maps
    if quick :
        skipAll = True
    else :
        skipAll = False
    #resolution = 256 ;# default, resolution element on order of ccd area size
    #resolution = 128 ;# roughly 4 ccds
    #resolution = 64 ;# very fast, debuging, roughly 1/4 of the camera size
    #if debug :
       # return getHexObservations.prepare(
       # skymap, trigger_id, outputDir, outputDir, distance=distance,
       # exposure_list=exposure_length, filter_list=filter_list,
       # overhead=overhead, maxHexesPerSlot=maxHexesPerSlot,
       # start_days_since_burst = recycler_days_since_burst, 
       # skipHexelate=False, skipAll = skipAll,
       # this_tiling = hexID_to_do, reject_hexes = hexID_to_reject, 
       # resolution=resolution, trigger_type=trigger_type, 
       # halfNight = halfNight, firstHalf= firstHalf, debug=debug,camera=camera)
    probs,times,slotDuration,hoursPerNight = getHexObservations.prepare(
        skymap, trigger_id, outputDir, outputDir, distance=distance,
        exposure_list=exposure_length, filter_list=filter_list,
        overhead=overhead, maxHexesPerSlot=maxHexesPerSlot,
        start_days_since_burst = recycler_days_since_burst, 
        skipHexelate=False, skipAll = skipAll,
        this_tiling = hexID_to_do, reject_hexes = hexID_to_reject,
        resolution=resolution, trigger_type=trigger_type,
        halfNight = halfNight, firstHalf = firstHalf,debug=debug, camera=camera)

#    print skymap, trigger_id, outputDir, outputDir, "distance=",distance, \
#        "exposure_list=",exposure_length, "filter_list=",filter_list, \
#        "overhead=",overhead, "maxHexesPerSlot=",maxHexesPerSlot, \
#        "start_days_since_burst = ",recycler_days_since_burst,  \
#        "skipHexelate=",False, "skipAll =", skipAll, \
#        "this_tiling =", hexID_to_do, "reject_hexes =", hexID_to_reject,  \
#        "resolution=",resolution, "trigger_type=",trigger_type, \
#        "halfNight =", halfNight, "firstHalf=", firstHalf
#
#    print probs,times,slotDuration,hoursPerNight  
                
    # figure out how to divide the night
    n_slots, first_slot = getHexObservations.contemplateTheDivisionsOfTime(
        probs, times, hoursPerNight=hoursPerNight,
        hoursAvailable=hoursAvailable)

    if quick :
        skipJson = True; 
    else :
        skipJson = False; 
    # compute the best observations
    best_slot = getHexObservations.now( 
        n_slots, mapDirectory=outputDir, simNumber=trigger_id, 
        maxHexesPerSlot=maxHexesPerSlot, mapZero=first_slot, 
        exposure_list=exposure_length, filter_list=filter_list, 
        trigger_type = trigger_type, skipJson =skipJson)


    if n_slots > 0 :
#   area_left is the number of hexes we have left to observe this season
#   T_left is the number of days left in the season
#   rate is the effective rate of triggers
#
        time_cost_per_hex = nepochs * np.sum(overhead + exposure_length) #sec 
        area_left =  area_per_hex * (hoursAvailable * 3600)/(time_cost_per_hex)
        time_left = end_of_season - start_of_season

        if skipEcon :
            print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            print "!!!!!!!!!!!!!     SKIPPING ECON ANALYSIS!  !!!!!!!!!!"
            print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            econ_prob, econ_area, quality = 0,0, 1.0
            need_area = 11734.0 
        # we may want to evaluate each event independently
        if not skipEcon: 
            # do Hsun-yu Chen's 
            print "======================================>>>>>>>>>>>>>>>>>>"
            print " economics "
            print "getHexObservations.economics (", trigger_id, ",",\
                best_slot, ", mapDirectory= \"",outputDir, "\" ,",\
                "area_left=",area_left, ", days_left=",time_left, ",rate=",rate,") "
            print "======================================>>>>>>>>>>>>>>>>>>"
            econ_prob, econ_area, need_area, quality = \
                getHexObservations.economics (trigger_id, 
                    best_slot, mapDirectory=outputDir, 
                    area_left=area_left, days_left=time_left, rate=rate) 
    
            if econ_area > 0.0 :
                hoursOnTarget = (econ_area/area_per_hex ) * (time_cost_per_hex/3600.)
    
                # figure out how to divide the night, 
                # given the new advice on how much time to spend
    
                n_slots, first_slot = getHexObservations.contemplateTheDivisionsOfTime(
                    probs, times, hoursPerNight=hoursPerNight,
                    hoursAvailable=hoursOnTarget)
    
                best_slot = getHexObservations.now( 
                    n_slots, mapDirectory=outputDir, simNumber=trigger_id, 
                    maxHexesPerSlot=maxHexesPerSlot, mapZero=first_slot, 
                    exposure_list=exposure_length, filter_list=filter_list, 
                    trigger_type = trigger_type, skipJson =skipJson)
    else :
        econ_prob, econ_area, quality = 0,0, 1.0
        need_area = 11734.0 

    if not quick :
        if n_slots > 0 :
            # make observation plots
            print "We're going to do {} slots with best slot {}".format(n_slots, best_slot)
            n_plots = getHexObservations.makeObservingPlots(n_slots, trigger_id, best_slot, outputDir, outputDir, camera, allSky = False)
            string = "$(ls -v {}-observingPlot*)  {}_animate.gif".format(trigger_id, trigger_id)
            os.system("convert  -delay 40 -loop 0  " + string)
        else :
            n_plots = getHexObservations.nothingToObserveShowSomething(trigger_id, outputDir, outputDir)

    # Lets find out how well we did in covering Ligo probability
    sum_ligo_prob = \
        getHexObservations.how_well_did_we_do(skymap, trigger_id, outputDir, camera, resolution)
    return best_slot, n_slots, first_slot, econ_prob, econ_area, need_area, quality



#  so what hexes did we do, given the json files?
def gethexIDfromJson (dir = "/data/des30.a/data/annis/des-gw/Jan4-2017-event/GW170104_night1_json/",
        verbose=True) :
    import json; import glob
    import os
    import hexalate
    cwd = os.getcwd(); print cwd
    os.chdir(dir)
    files = glob.glob("*.json")
    nfiles=len(files);
    ra = np.array([])
    dec = np.array([])
    slot = np.array([])
    for n in range(nfiles) :
        print files[n]
        fd=open(files[n],"r"); data=json.load(fd);fd.close()
        for i in range(len(data)): 
            ra = np.append(ra, data[i]["RA"])
            dec = np.append(dec, data[i]["dec"])
            slot = np.append(slot, files[n][9:11])
            if verbose: print data[i]["RA"],data[i]["dec"],files[n][9:11]

    os.chdir(cwd)
    id = hexalate.getHexId(ra,dec)
    return ra,dec,id, slot
def gethexIDfromDB() :
    import hexalate
    file="/data/des30.a/data/annis/des-gw/Jan4-2017-event/GW170104_exp_table_night1.txt"
    file="/data/des30.a/data/annis/des-gw/Jan4-2017-event/GW170104_exp_table_night1_2.txt" 
    ra,dec = np.genfromtxt(file, unpack=True,usecols=(1,2)); 
    ii=ra>180; ra[ii] = ra[ii]-360
    id = hexalate.getHexId(ra,dec)
    id = np.unique(id)
    return id
#
# shall we run the hexes from a ra-dec-id-prob-mjd-slot.tx file?
def hexID (file = "../Jan4-2017-event/GW170104-ra-dec-id-prob-mjd-slot.txt") :
    ra,dec = np.genfromtxt(file, unpack=True,usecols=(0,1) )
    id = np.genfromtxt(file, unpack=True,usecols=(2), dtype = "str" )
    # for GW170104; doing this cut does the obs right, not doing makes the gif cover better area
    ix = ra < 100; ra = ra[ix]; dec = dec[ix]; id = id[ix]
    return ra,dec,id
