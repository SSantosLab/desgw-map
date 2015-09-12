import numpy as np
import healpy as hp
import hp2np
import sourceProb
import mags
import modelRead
import mapsAtTimeT


# sims, mjds, distances, models = allMaps.veni(); 
# allMaps.vidi(sims, mjds, distances, models)

# sims, mjds, distances, models = allMaps.veni(False); 
# allMaps.vidi(sims, mjds, distances, models, False)

# ra, dec, map, maglim, maglimg, prob, probMap, hx,hy = allMaps.readem( 789064, 0)
# sim,mjd,distance,snr, hlv, cra, cdec, p0,p1,p2,p3,p4,p5,p6,p7,p8,p9 = np.genfromtxt("all-maxtimes-2015.txt", unpack=True);
#
#=================================================
#
# RUN the LIGO Simulations!!!
#
#=================================================

def selectSingleSim (
        simNumber, data_dir="/data/des30.a/data/annis/des-gw/ligo/sims/") :
    sims, mjds, distances, models = veni()
    ix = sims==simNumber
    sims, mjds, distances = sims[ix], mjds[ix], distances[ix] 
    print "found ",sims[0]
    simfile="bayestar-{:d}.fits.gz".format(sims[0])
    ligoMapFile = data_dir+simfile
    return sims, mjds, distances, ligoMapFile

#
#==== Big routine # 1: collect metadata
#
# Get the list of Ligo Maps
# Get the NS models
#    Definitively a prep into memory routine
#
def veni( do2015=True) :
    dir = "/data/des30.a/data/annis/des-gw/ligo/"
    #simsFile = "sims.list"
    #sims = np.genfromtxt(dir+simsFile, unpack=True)

    if do2015:
        simsFile = "2015_inj.txt"
    else :
        simsFile = "2016_inj.txt"
    sims, mjds, distances = np.genfromtxt(dir+simsFile, unpack=True, skiprows=40, usecols=(0,2,8))
    sims = sims.astype("int")

    models = modelRead.getModels()
    return sims, mjds, distances, models

#==== Big routine # 2: find the probabilities over 10 days
#
#  For each sim build the  mags.observed object
#  Then run mapsAtTimeT.manyDaysOfTotalProbability
#  Then run maximumProbabilityPerDay
#  Save the resulting stuff
#  Then run mapsAtTimeT.probabilityMapSaver
#       the hex values and rank are written to a file
#
def vidi(sims, mjds, distances, models, do2015=True) :
    import os.path
    if do2015:
        data_dir = "/data/des30.a/data/annis/des-gw/ligo/sims/"
        odata_dir = "/data/des30.a/data/annis/des-gw/ligo/sims-2015-out/"
    else :
        data_dir = "/data/des30.a/data/annis/des-gw/ligo/sims-2016/"
        odata_dir = "/data/des30.a/data/annis/des-gw/ligo/sims-2016-out/"
    file = "bayestar-{:d}.fits.gz"
    for sim, mjd, distance in zip(sims,mjds,distances) :
        simfile = file.format(sim)
        print "sim, distance: ", sim, distance
        name = maxtimeFilename(sim, odata_dir)
        if os.path.exists(name) : continue

        simfile = data_dir+simfile
        ligo = hp.read_map(simfile)
        ra,dec,ligo = hp2np.map2np(ligo, resolution=32)
        obs = mags.observed(ra,dec,ligo, mjd); 

        totalProbs, times  = mapsAtTimeT.manyDaysOfTotalProbability(
            obs, mjd, distance, models)
        maxtimes, maxprobs = maximumProbabilityPerDay(totalProbs, times)
        print "maxtimes= ",maxtimes
        print "maxprobs= ",maxprobs
        saven2 (maxtimes, maxprobs, sim, odata_dir)
        # make the small files
        mapsAtTimeT.probabilityMapSaver (obs, sim, mjd, distance, models, maxtimes, maxprobs, odata_dir)

#==== Big routine # 3:  higher resolution run against output of vidi
#
# This does most of what vidi does, except: 
#       the map is higher resolution,
#       the maxtimes are read in from a file (as opposed to calculated by maximumProbabilityPerDay)
#       and only mapsAtTimeT.probabilityMapSaver is run, 
#           not manyDaysOfTotalProbability and maximumProbabilityPerDay
#   Presumably the big maps are expensive and can only be run for certain times
#       the hex values and rank are written to a file
#
# make the big files
def vici(sims, mjds, distances, models, do2015=True) :
    import os.path
    import shutil
    if do2015:
        data_dir = "/data/des30.a/data/annis/des-gw/ligo/sims/"
        odata_dir = "/data/des30.a/data/annis/des-gw/ligo/sims-2015-out/"
        vdata_dir = "/data/des30.a/data/annis/des-gw/ligo/sims-2015-out-v/"
    else :
        data_dir = "/data/des30.a/data/annis/des-gw/ligo/sims-2016/"
        odata_dir = "/data/des30.a/data/annis/des-gw/ligo/sims-2016-out/"
        vdata_dir = "/data/des30.a/data/annis/des-gw/ligo/sims-2016-out-v/"
    file = "bayestar-{:d}.fits.gz"
    for sim, mjd, distance in zip(sims,mjds,distances) :
        simfile = file.format(sim)
        print "sim, distance: ", sim, distance

        new_maxtimeFile = maxtimeFilename(sim, vdata_dir)
        if os.path.exists(new_maxtimeFile) : continue
        old_maxtimeFile= maxtimeFilename(sim, odata_dir)
        shutil.copy(old_maxtimeFile, new_maxtimeFile)

        simfile = data_dir+simfile
        ligo = hp.read_map(simfile)
        ra,dec,ligo = hp2np.map2np(ligo, resolution=256)
        obs = mags.observed(ra,dec,ligo, mjd); 

        maxtimes, maxprobs = np.genfromtxt(new_maxtimeFile, unpack=True)
        print "maxtimes= ",maxtimes
        print "maxprobs= ",maxprobs
        mapsAtTimeT.probabilityMapSaver(obs, sim, mjd, distance, models, maxtimes, maxprobs, vdata_dir)

#==== helper routine around the higher resolution routine of Big routine # 3, collecting meta data?
#
# a Jack routine. Not useful therefore?
# same input as vici,
#   but just read the sims (not mjd, distance)
#   and write a file with ra,dec of the maximim likes
def vjack(sims, do2015=True) :
    import os.path
    import shutil
    if do2015:
        data_dir = "/data/des30.a/data/annis/des-gw/ligo/sims/"
        odata_dir = "/data/des30.a/data/annis/des-gw/ligo/sims-2015-out/"
        vdata_dir = "/data/des30.a/data/annis/des-gw/ligo/sims-2015-out-v/"
        savefile = "check-{}.txt".format("2015")
    else :
        data_dir = "/data/des30.a/data/annis/des-gw/ligo/sims-2016/"
        odata_dir = "/data/des30.a/data/annis/des-gw/ligo/sims-2016-out/"
        vdata_dir = "/data/des30.a/data/annis/des-gw/ligo/sims-2016-out-v/"
        savefile = "check-{}.txt".format("2016")
    file = "bayestar-{:d}.fits.gz"
    raList, decList = [],[]
    for sim in sims :
        simfile = file.format(int(sim))
        print "sim, distance: ", sim
        simfile = data_dir+simfile
        ligo = hp.read_map(simfile)
        ra,dec,ligo = hp2np.hp2np(simfile)
        max = ligo.max()
        ix = ligo == max
        ra = ra[ix]
        dec = dec[ix]
        if ra.size > 0 : ra = ra[0]
        if dec.size > 0 : dec = dec[0]
        raList.append(ra)
        decList.append(dec)
    np.savetxt(savefile, np.array([sims, raList, decList]).T, "%d %.6f %.5f")


#  sim,mjd,distance,snr, p0,p1,p2,p3,p4,p5,p6,p7,p8,p9 = np.genfromtxt("all-maxtimes-2015.txt", unpack=True);

#==== Big routine # 4:  collate information into a single file
#
# There is no such thing as vedi. It is veni, vidi, vici.
# This is concerned with the LIGO ancillary data
# which it is going to connect to....
#       the total probability for each day.
# and thus to a file. I don't know how the day maps are made....
#
# make the big files, the all-maxtimes files
def vedi2(sims, mjds, distances, do2015=True, doV=False) :
    import os.path
    import shutil
    savefile="all-maxtimes-{}.txt"
    snrFile = "/data/des30.a/data/annis/des-gw/ligo/{}_coinc.txt"
    checkFile = "/data/des30.a/data/annis/des-gw/ligo/check-{}.txt"
    if do2015:
        data_dir = "/data/des30.a/data/annis/des-gw/ligo/sims-2015-out/"
        if doV :
            data_dir = "/data/des30.a/data/annis/des-gw/ligo/sims-2015-out-v/"
        savefile = savefile.format("2015")
        snrFile = snrFile.format("2015")
        snrSim, snrNet = np.genfromtxt(snrFile, unpack=True,skiprows=35, usecols=(0,3))
        checkFile = checkFile.format("2015")
        network = np.zeros(snrSim.size)
    else :
        data_dir = "/data/des30.a/data/annis/des-gw/ligo/sims-2016-out/"
        if doV :
            data_dir = "/data/des30.a/data/annis/des-gw/ligo/sims-2016-out-v/"
        savefile = savefile.format("2016")
        snrFile = snrFile.format("2016")
        snrSim, snrNet = np.genfromtxt(snrFile, unpack=True,skiprows=36, usecols=(0,3))
        network = np.genfromtxt(snrFile, unpack=True,skiprows=36, usecols=(2),dtype=("str"))
        checkFile = checkFile.format("2016")
    check_sim, check_ra, check_dec = np.genfromtxt(checkFile, unpack=True)
    csim,rapidarea = np.genfromtxt("../c2016.txt",unpack=True)
    totalProbs = dict()
    for day in range(0,10) :
        totalProbs[day]  = []
    snr = []
    hlv = []
    cra, cdec = [],[]
    carea =[]
    for sim, mjd, distance in zip(sims,mjds,distances) :
        print "sim, distance: ", sim, distance

        new_maxtimeFile = maxtimeFilename(sim, data_dir)
        maxtimes, maxprobs = np.genfromtxt(new_maxtimeFile, unpack=True)
        snr_ix = np.nonzero(snrSim == sim)
        if network[snr_ix] == "HL" : 
            hlv.append(0)
        elif network[snr_ix] == "HLV" : 
            hlv.append(1)
        elif network[snr_ix] == "HV" : 
            hlv.append(2)
        elif network[snr_ix] == "LV" : 
            hlv.append(3)
        else :
            hlv.append(4)
        snr.append( snrNet[snr_ix] )
        check_ix = np.nonzero(check_sim == sim)
        cra.append(check_ra[check_ix])
        cdec.append(check_dec[check_ix])
        check_ix = np.nonzero(csim == sim)
        carea.append(rapidarea[check_ix])

        for day in range(0,10) :
            simfile = data_dir+str(sim)+"-"+str(day)+"-map.hp"
            if not os.path.exists(simfile) : 
                totalProb = 0.0
            else :
                ligo = hp.read_map(simfile)
                simfile = data_dir+str(sim)+"-"+str(day)+"-probMap.hp"
                decam = hp.read_map(simfile)
                ligo = ligo/ligo.sum()
                totalProb = (decam*ligo).sum()
            totalProbs[day].append(totalProb)
    snr= np.array(snr)
    hlv= np.array(hlv).astype(int)
    cra = np.array(cra)
    cdec = np.array(cdec)
    carea = np.array(carea)
    for day in range(0,10) :
        totalProbs[day] = np.array(totalProbs[day])
    data = np.array([sims, mjds, distances, snr, hlv, cra, cdec, carea, totalProbs[0], \
        totalProbs[1], totalProbs[2], totalProbs[3], totalProbs[4], \
        totalProbs[5], totalProbs[6], totalProbs[7], totalProbs[8], \
        totalProbs[9]])
    np.savetxt(savefile, data.T, "%d %.5f %.0f %.1f %d %.6f %.5f %.3f %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e")


# Get the saved maps for each day.
def readem (simNumber, day) :
    n=str(simNumber)+"-"+str(day)

    ra=hp.read_map(n+"-ra.hp"); 
    dec=hp.read_map(n+"-dec.hp"); 
    ha=hp.read_map(n+"-ha.hp"); 
    map=hp.read_map(n+"-map.hp");
    maglim=hp.read_map(n+"-maglim.hp");
    maglimg=hp.read_map(n+"-maglim-global.hp");
    prob=hp.read_map(n+"-prob.hp");
    probMap=hp.read_map(n+"-probMap.hp"); 
    hx=hp.read_map(n+"-hx.hp"); 
    hy=hp.read_map(n+"-hy.hp");
    return ra, dec, map, maglim, maglimg, prob, probMap, hx,hy

# Get the saved maps for each day and hour.
def readem2 (simNumber, day, hourno) :
    n=str(simNumber)+"-"+str(day)+"-"+str(hourno)

    ra=hp.read_map(n+"-ra.hp"); 
    dec=hp.read_map(n+"-dec.hp"); 
    ha=hp.read_map(n+"-ha.hp"); 
    map=hp.read_map(n+"-map.hp");
    maglim=hp.read_map(n+"-maglim.hp");
    prob=hp.read_map(n+"-prob.hp");
    probMap=hp.read_map(n+"-probMap.hp"); 
    hx=hp.read_map(n+"-hx.hp"); 
    hy=hp.read_map(n+"-hy.hp");
    return ra, dec, map, maglim, prob, probMap, hx,hy

# A clear cut save routine
def saven2 (maxtimes, maxprobs, simNumber, data_dir) :
    name = maxtimeFilename ( simNumber, data_dir)
    np.savetxt(name, np.array([maxtimes,maxprobs]).T, "%.5f %.5e")
    print "\t writing ", name

# A file finder routine
def maxtimeFilename ( simNumber, data_dir) :
    nameStem = data_dir + str(simNumber) 
    name = nameStem + "-maxtimes.txt"
    return name
    
#  for ten days,
#   find the maximum total probability for each day
#   return
def maximumProbabilityPerDay (totalProbs,times) :
    maxtimes = []
    maxprobs = []
    # for ten days
    for day in range(0,10) :
        # for times measured in days
        ix = (times >= day) & (times < day+(1.))
        max = totalProbs[ix].max()
        ix = np.nonzero((times >= day) & (times < day+(1.)) & (totalProbs == max))
        if totalProbs[ix].sum() == 0 :
            ix = ix[0][0]
        if totalProbs[ix].size > 1 :
            ix = ix[0][0]
        if type(ix) == np.int64 :
            t = times[ix]; tp = totalProbs[ix]
        else :
            t = times[ix][0]; tp = totalProbs[ix][0]
        maxtimes.append(t)
        maxprobs.append(tp)
    maxtimes =np.array(maxtimes)
    maxprobs = np.array(maxprobs)
    return maxtimes, maxprobs

#===================================================================
#

# ==================================
# plotting 
# ==================================

def histProbs () :
    import os
    os.chdir("/data/des30.a/data/annis/des-gw/ligo/sims-2015-out")
    sim,mjd,distance,p0,p1,p2,p3,p4,p5,p6,p7,p8,p9 = np.genfromtxt("all-maxtimes-2015.txt", unpack=True);
    ix = p0 > .01; p=p0[ix]; n=p0.size-p.size; ix=p0>.33; n2=p0[ix].size; plt.clf();a=plt.hist(p,bins=100,color="k"); plt.text(0.4,23,"day 0: n with p < 0.01: {:3d}      n with p > 0.33: {:3d}".format(n,n2))
    ix = p1 > .01; p=p1[ix]; n=p1.size-p.size;ix=p1>.33; n2=p1[ix].size; a=plt.hist(p,bins=100,color="r",histtype="step"); plt.text(0.4,21,"day 1: n with p < 0.01: {:3d}      n with p > 0.33: {:3d}".format(n,n2),color="r")
    ix = p4 > .01; p=p4[ix]; n=p4.size-p.size;ix=p4>.33; n2=p4[ix].size; a=plt.hist(p,bins=100,color="g",histtype="step"); plt.text(0.4,19,"day 4: n with p < 0.01: {:3d}      n with p > 0.33: {:3d}".format(n,n2),color="g")
    ix = p9 > .01; p=p9[ix]; n=p9.size-p.size;ix=p9>.33; n2=p9[ix].size; a=plt.hist(p,bins=100,color="b"); plt.text(0.4,17,"day 9: n with p < 0.01: {:3d}      n with p > 0.33: {:3d}".format(n,n2),color="b")
    plt.ylim(0,28);plt.xlabel("probability of detection");plt.ylabel("N");plt.title("2015")
    plt.savefig("hist-2015.pdf")

    os.chdir("/data/des30.a/data/annis/des-gw/ligo/sims-2016-out")
    sim,mjd,distance,p0,p1,p2,p3,p4,p5,p6,p7,p8,p9 = np.genfromtxt("all-maxtimes-2016.txt", unpack=True);
    ix = p0 > .01; p=p0[ix]; n=p0.size-p.size; ix=p0>.33; n2=p0[ix].size; plt.clf();a=plt.hist(p,bins=100,color="k"); plt.text(0.4,23,"day 0: n with p < 0.01: {:3d}      n with p > 0.33: {:3d}".format(n,n2))
    ix = p1 > .01; p=p1[ix]; n=p1.size-p.size;ix=p1>.33; n2=p1[ix].size; a=plt.hist(p,bins=100,color="r",histtype="step"); plt.text(0.4,21,"day 1: n with p < 0.01: {:3d}      n with p > 0.33: {:3d}".format(n,n2),color="r")
    ix = p4 > .01; p=p4[ix]; n=p4.size-p.size;ix=p4>.33; n2=p4[ix].size; a=plt.hist(p,bins=100,color="g",histtype="step"); plt.text(0.4,19,"day 4: n with p < 0.01: {:3d}      n with p > 0.33: {:3d}".format(n,n2),color="g")
    ix = p9 > .01; p=p9[ix]; n=p9.size-p.size;ix=p9>.33; n2=p9[ix].size; a=plt.hist(p,bins=100,color="b"); plt.text(0.4,17,"day 9: n with p < 0.01: {:3d}      n with p > 0.33: {:3d}".format(n,n2),color="b")
    plt.ylim(0,28);plt.xlabel("probability of detection");plt.ylabel("N");plt.title("2016")
    plt.savefig("hist-2016.pdf")
