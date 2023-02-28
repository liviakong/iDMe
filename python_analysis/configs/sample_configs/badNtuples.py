import uproot
import numpy as np
import json
import sys
from XRootD import client
from tqdm import tqdm
from multiprocessing import Process, Value, Array, Manager
import multiprocessing as mp
from ROOT import TTree, TFile

def checkFile(ntuple):
    #try:
    good = True
    with TFile.Open(ntuple,"read") as f:
        t = f.Get("ntuples/outT")
        if not t.GetTree():
            print("bad tree!")
            good = False
        if good:
            for i in range(0,t.GetEntries()):
                error = t.GetEntry(i)
                if error < 0:
                    print("bad event!")
                    good = False
                    break
    return good
    #except:
    #    print("exception!")
    #    return False

inputJson = sys.argv[1]
with open(inputJson) as f:
    samples = json.load(f)

for samp in samples:
    print(f"Running on {samp['name']}")
    loc = samp["location"]
    nFiles = -1
    blacklist = [] # track files that can't be opened
    sum_wgt = []
    if '.root' in loc:
        if not checkFile(loc):
            blacklist.append(loc)
    else:
        xrdClient = client.FileSystem("root://cmseos.fnal.gov")
        if type(loc) != list:
            status, flist = xrdClient.dirlist(loc)
            fullList = ["root://cmsxrootd.fnal.gov/"+loc+"/"+item.name for item in flist if '.root' in item.name]
        else:
            fullList = []
            for l in loc:
                status, flist = xrdClient.dirlist(l)
                fullList.extend(["root://cmsxrootd.fnal.gov/"+l+"/"+item.name for item in flist if '.root' in item.name])
        for f in tqdm(fullList):
            if not checkFile(f):
                blacklist.append(f)

    print('Blacklisted {0} files in {1}'.format(len(blacklist),samp['name']))
    samp['blacklist'] = blacklist

    #checkpointing
    with open(inputJson,'w') as f:
        json.dump(samples,f,indent=4)
