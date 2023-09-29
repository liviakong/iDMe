from XRootD import client
import json
import sys
import subprocess
import numpy as np
import re
import datetime as dt
import os

mode = str(sys.argv[1])
year = str(sys.argv[2])
alpha = str(sys.argv[3])

if mode != "sig" and mode != "bkg":
    print("Invalid mode: use sig or bkg")
    exit
if len(sys.argv) < 3:
    print("Invalid input arguments!")
    print("Usage: python makeSignalConfigs.py mode[sig/bkg] year alpha")
    exit

xrdClient = client.FileSystem("root://cmseos.fnal.gov")

if mode == "sig":
    prefix = "/store/group/lpcmetx/iDMe//Samples/Ntuples/signal_v2/"

    status, points = xrdClient.dirlist(prefix+year)
    points = [item.name for item in points]
    output = []
    for p in points:
        mchi = float(p.split("_")[0].split("-")[1].replace("p","."))
        dmchi = float(p.split("_")[1].split("-")[1].replace("p","."))
        status, lifetimes = xrdClient.dirlist(prefix+year+"/"+p)
        lifetimes = [l.name for l in lifetimes]
        
        for l in lifetimes:
            ct = int(l.split("-")[1])
            info = {}
            info["location"] = prefix+year+"/"+p+"/"+l+"/"  
            info["Mchi"] = mchi
            info["dMchi"] = dmchi
            info["ctau"] = ct
            info["name"] = "sig_Mchi-{0}_dMchi-{1}_ct-{2}".format(info["Mchi"],info["dMchi"],info["ctau"])
            info["sum_wgt"] = 0.0
            info["type"] = "signal"
            info["year"] = int(year)
            info["alphaD"] = alpha
            info["xsec"] = 0.0
            rootFiles = [rf.name for rf in xrdClient.dirlist(info["location"])[1] if '.root' in rf.name]
            info["nFiles"] = len(rootFiles)
            output.append(info)

    out_json = "signal_{0}_{1}.json".format(year,alpha)
    with open(out_json,"w") as outfile:
        json.dump(output,outfile,indent=4)
else:
    prefix = "/store/group/lpcmetx/iDMe//Samples/Ntuples/background/"
    status, bkgs = xrdClient.dirlist(prefix+year)
    bkgs = [bkg.name for bkg in bkgs]
    for bkg in bkgs:
        base_dir = prefix+year+"/"+bkg
        rootFiles = subprocess.run(['eos','root://cmseos.fnal.gov/','find',base_dir,'-name','*.root'],stdout=subprocess.PIPE).stdout.decode('utf-8').splitlines()
        rootFiles = [r for r in rootFiles if '.root' in r]
        fileDirs = ["/".join(f.split("/")[:-1])+"/" for f in rootFiles]
        fileDirs = list(set(fileDirs)) # list of unique file directories
        
        names = fileDirs.copy()
        names = [n.split("/") for n in names]
        reference = names[0].copy()
        # removing common file path prefix from all paths
        for r in reference:
            stop = False
            for n in names:
                if r not in n:
                    stop = True
                    break
            if stop:
                break
            names = [n[1:] for n in names]
        # only need first part of path
        names = [n[0] for n in names] # head of path, i.e. qcd ht bin, for all unique file-holding dirs

        most_recent = ""
        now = dt.datetime.now()
        min_dt = dt.timedelta.max.total_seconds()
        for fd in fileDirs:
            d = re.findall("202\d_\d\d_\d\d-",fd)[0]
            delta = (now - dt.datetime.strptime(d,"%Y_%m_%d-")).total_seconds()
            if delta < min_dt:
                min_dt = delta
                most_recent = d
        print(most_recent)
        names = [names[i] for i in range(len(names)) if most_recent in fileDirs[i]]
        fileDirs = [f for f in fileDirs if most_recent in f]
        name_loc_map = {n:[] for n in set(names)}
        for i in range(len(names)):
            name_loc_map[names[i]].append(fileDirs[i])
        output = []
        for n in name_loc_map.keys():
            info = {}
            info["name"] = n
            info["location"] = name_loc_map[n]
            info["sum_wgt"] = 0.0
            info["type"] = "bkg"
            info["year"] = int(year)
            info["xsec"] = 0.0
            nFiles=0
            for loc in name_loc_map[n]:
                nFiles += len([rf.name for rf in xrdClient.dirlist(loc)[1] if '.root' in rf.name])
            info["nFiles"] = nFiles
            output.append(info)

        out_json = "bkg_{}_{}.json".format(year,bkg)
        with open(out_json,"w") as outfile:
            json.dump(output,outfile,indent=4)
