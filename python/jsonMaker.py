import numpy as np

# Support:
#
# How to use JSONs
#   import json
#   from pprint import pprint
#   json_data=open('json_data')
#   data = json.load(json_data)
#   pprint(data)
#   json_data.close()
# 
# What is format of JSON for CTIO?
#   http://dns.ctio.noao.edu/noao/content/Preparing-Observing-Scripts
#   {"count": 3,
#       "expType": "object",
#       "object": "Sgr",
#       "filter": "g",
#       "RA": 281.0,
#       "dec": -29.866,
#       "expTime": 60
#       },

#
# This code wants a list ra,decs and will write a JSON file
# to cause the Blanco to observe them.
# ToDo: what is a hex_id? -get from master json file. 
# ToDo: put trigger_id and night into seq id
#
#        exposureList = [90,90,90,90], 
#        filterList = ["i","z","z","i"],
#        tilingList = [9,9,10,10], 
#        jsonFilename="des-gw.json") :
def writeJson(ra,dec, seqid="none", seqnum=0, seqtot=0,
        exposureList = [90,90,90], 
        filterList = ["i","z","z"],
        tilingList = [9,9,10], 
        jsonFilename="des-gw.json") :
    offsets = tileOffsets()
    fd = open(jsonFilename,"w")
    fd.write("[\n")

    size = ra.size
    nexp = np.size(exposureList)
    for i in range(0,size) :
        for j in range(0,nexp) :
            seqnum +=1
#not clobbered
#count, seqid, seqnum, seqtot,note,comment
            tiling = tilingList[j]
            filter = filterList[j]
            exp = exposureList[j]
            offsets[tiling]
            #print tiling, offsets[tiling]
            delRa = offsets[tiling][0]
            delDec = offsets[tiling][1]
            tra = ra[i]
            tdec = dec[i]
            tdec = tdec+delDec
            tra = tra + delRa/np.cos(tdec*2*np.pi/360.)
            comment = "DESGW: LIGO event followup with DES/DECam: hex ranked {}, filter {}, tile {}".format(i,filter, tiling)
            intra = np.int(np.round(np.int(ra[i]*10.)/10.)*10.)
            intdec = np.int(np.round(np.int(dec[i]*10.)/10.)*10.)
            signDec = "+"; 
            if tdec < 0: signDec = "-"
            object = "DES wide hex {:3d}{}{:3d}".format(intra, signDec, intdec, tiling)

            fd.write("{")
            fd.write(" \"expType\" : \"object\",\n")
            fd.write("  \"object\" : \"{}\",\n".format(object))
            fd.write("  \"seqid\" : \"{}\",\n".format(seqid))
            fd.write("  \"seqnum\" : \"{:d}\",\n".format(int(seqnum)))
            fd.write("  \"seqtot\" : \"{:d}\",\n".format(int(seqtot)))
            fd.write("  \"expTime\" : {:d},\n".format(int(exp)))
            fd.write("  \"wait\" : \"False\",\n")
            fd.write("  \"count\" : \"1\",\n")
            fd.write("  \"note\" : \"Added to queue from desgw json file, not obstac\",\n")
            fd.write("  \"filter\" : \"{}\",\n".format(filter))
            fd.write("  \"program\" : \"des gw\",\n")
            fd.write("  \"RA\" : {:.6f},\n".format(tra))
            fd.write("  \"dec\" : {:.5f},\n".format(tdec))
            fd.write("  \"comment\" : \"{}\"\n".format(comment)) 
            # note lack of comma for end
            fd.write("}")
            if (i == size-1) and ( j == nexp-1) :
                pass
            else :
                fd.write(",")
            fd.write("\n")
    fd.write("]\n")
    fd.close()

# production offsets are in DES docdb 7269
# production offsets are one based, not zero based:
def tileOffsets() :
    offsets = dict()
# blessed offsets 0
    offsets[0] = 0.000, 0.000
# we'll use jta's 8,9 (as if we were doing DES year 5)
# production offsets 9
    offsets[9] = 0.76668, 0.4227
# production offsets 10
    offsets[10] = -0.0479175, 0.388884
# production offsets 10
    offsets[11] = -0.5257, 0.7222
# production offsets 17,18
    offsets[17] = -1.1388, 0.0166
    offsets[18] =  0.0484, -0.6725
    return offsets
