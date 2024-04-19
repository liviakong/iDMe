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
prefix = str(sys.argv[4])
name = str(sys.argv[5])

if mode != "sig" and mode != "bkg" and mode != "data":
    print("Invalid mode: use sig/bkg/data")
    exit
if len(sys.argv) < 3:
    print("Invalid input arguments!")
    print("Usage: python makeSignalConfigs.py mode[sig/bkg/data] year alpha")
    exit

xrdClient = client.FileSystem("root://cmseos.fnal.gov")

if mode == "sig":
    status, points = xrdClient.dirlist(f"{prefix}/{year}/")
    points = [item.name for item in points]
    output = []
    for p in points:
        mchi = float(p.split("_")[0].split("-")[1].replace("p","."))
        dmchi = float(p.split("_")[1].split("-")[1].replace("p","."))
        if 'mZD' in p:
            mzd = p.split("_")[2]
        else:
            mzd = ""
        status, lifetimes = xrdClient.dirlist(f"{prefix}/{year}/{p}")
        lifetimes = [l.name for l in lifetimes]
        
        for l in lifetimes:
            ct = int(l.split("-")[1])
            info = {}
            info["location"] = f"{prefix}/{year}/{p}/{l}/"
            info["Mchi"] = mchi
            info["dMchi"] = dmchi
            info["ctau"] = ct
            if mzd != "":
                info["name"] = "sig_Mchi-{0}_dMchi-{1}_ct-{2}_{3}".format(info["Mchi"],info["dMchi"],info["ctau"],mzd)
            else:
                info["name"] = "sig_Mchi-{0}_dMchi-{1}_ct-{2}".format(info["Mchi"],info["dMchi"],info["ctau"])
            info["sum_wgt"] = 0.0
            info["type"] = "signal"
            info["year"] = int(year)
            info["alphaD"] = alpha
            info["xsec"] = 0.0
            rootFiles = [rf.name for rf in xrdClient.dirlist(info["location"])[1] if '.root' in rf.name]
            info["nFiles"] = len(rootFiles)
            output.append(info)

    out_json = "signal_{0}_{1}_{2}.json".format(name,year,alpha)
    with open(out_json,"w") as outfile:
        json.dump(output,outfile,indent=4)
elif mode == "bkg":
    status, bkgs = xrdClient.dirlist(f"{prefix}/{year}/")
    bkgs = [bkg.name for bkg in bkgs]
    output = []
    for bkg in bkgs:
        base_dir = f"{prefix}/{year}/{bkg}"
        subsamples = [d.name for d in xrdClient.dirlist(base_dir)[1]]
        for subsample in subsamples:
            target_dir = f"{base_dir}/{subsample}/"
            rootFiles = subprocess.run(['eos','root://cmseos.fnal.gov/','find','-name','*.root','-f',target_dir],stdout=subprocess.PIPE).stdout.decode('utf-8').splitlines()
            rootFiles = [r for r in rootFiles if '.root' in r]
            fileDirs = ["/".join(f.split("/")[:-1])+"/" for f in rootFiles]
            fileDirs = list(set(fileDirs)) # list of unique file directories
            
            info = {}
            info["name"] = f"{bkg}_{subsample}"
            info["location"] = fileDirs[0] if len(fileDirs) == 1 else fileDirs
            info["sum_wgt"] = 0.0
            info["type"] = "bkg"
            info["year"] = int(year)
            info["xsec"] = 0.0
            nFiles=0
            for fdir in fileDirs:
                nFiles += len([rf.name for rf in xrdClient.dirlist(fdir)[1] if '.root' in rf.name])
            info["nFiles"] = nFiles
            output.append(info)

    out_json = "bkg_{0}_{1}.json".format(year,name)
    with open(out_json,"w") as outfile:
        json.dump(output,outfile,indent=4)
elif mode == "data":
    status, samples = xrdClient.dirlist(f"{prefix}/{year}/")
    samples = [samp.name for samp in samples]
    output = []
    for samp in samples:
        base_dir = f"{prefix}/{year}/{samp}"
        subsamples = [d.name for d in xrdClient.dirlist(base_dir)[1]]
        for subsample in subsamples:
            target_dir = f"{base_dir}/{subsample}/"
            rootFiles = subprocess.run(['eos','root://cmseos.fnal.gov/','find','-name','*.root','-f',target_dir],stdout=subprocess.PIPE).stdout.decode('utf-8').splitlines()
            rootFiles = [r for r in rootFiles if '.root' in r]
            fileDirs = ["/".join(f.split("/")[:-1])+"/" for f in rootFiles]
            fileDirs = list(set(fileDirs)) # list of unique file directories
            
            info = {}
            info["name"] = f"{samp}_{subsample}"
            info["location"] = fileDirs[0] if len(fileDirs) == 1 else fileDirs
            info["sum_wgt"] = 0.0
            info["type"] = "data"
            info["year"] = int(year)
            info["xsec"] = 0.0
            nFiles=0
            for fdir in fileDirs:
                nFiles += len([rf.name for rf in xrdClient.dirlist(fdir)[1] if '.root' in rf.name])
            info["nFiles"] = nFiles
            output.append(info)

    out_json = "data_{0}_{1}.json".format(year,name)
    with open(out_json,"w") as outfile:
        json.dump(output,outfile,indent=4)
