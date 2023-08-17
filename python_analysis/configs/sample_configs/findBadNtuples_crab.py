import uproot
import numpy as np
import json
import sys
from XRootD import client
from tqdm import tqdm
from multiprocessing import Process, Value, Array, Manager
import multiprocessing as mp

inputJson = sys.argv[1]
with open(inputJson) as f:
    samples = json.load(f)

for samp in samples:
    print(f"Running on {samp['name']}")
    loc = samp["location"]
    xrdClient = client.FileSystem("root://cmseos.fnal.gov")
    if type(loc) != list:
        status, flist = xrdClient.dirlist(loc)
        fullList = ["root://cmsxrootd.fnal.gov/"+loc+"/"+item.name for item in flist if '.root' in item.name]
    else:
        fullList = []
        for l in loc:
            status, flist = xrdClient.dirlist(l)
            fullList.extend(["root://cmsxrootd.fnal.gov/"+l+"/"+item.name for item in flist if '.root' in item.name])

    crab_job = [k for k in fullList[0].split("/") if "crab_iDMe_" in k][0]
    print("Crab job : "+crab_job)
    with open("/uscms/home/sbrightt/nobackup/iDM/iDMe_analysis/CMSSW_10_6_26/src/iDMe/AODSkimmer/crab/submissions_miniNtuplizer/"+crab_job+"/failedJobs.txt") as f:
        failedJobs = f.read().splitlines()
    failedJobs = [j+".root" for j in failedJobs]

    blacklist = []
    for rootFile in tqdm(fullList):
        if rootFile.split("_")[-1] in failedJobs:
            blacklist.append(rootFile.split("/")[-1])

    print('Blacklisted {0} files in {1}'.format(len(blacklist),samp['name']))
    if 'blacklist' in samp.keys():
        samp['blacklist'].extend(blacklist)
        samp['blacklist'] = list(set(samp['blacklist']))
    else:
        samp['blacklist'] = blacklist

    #checkpointing
    with open(inputJson,'w') as f:
        json.dump(samples,f,indent=4)
