from XRootD import client
import json
import sys
import subprocess
import numpy as np

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
    prefix = "/store/group/lpcmetx/iDMe//Samples/Ntuples/signal/"

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
        common = []
        for r in reference:
            com = True
            for n in names:
                if r not in n:
                    com = False
            if com:
                common.append(r)
        for c in common:
            for n in names:
                n.remove(c)
        names = [n[0] for n in names]

        output = []
        for i in range(len(names)):
            info = {}
            info["name"] = names[i]
            info["location"] = fileDirs[i]
            info["sum_wgt"] = 0.0
            info["type"] = "bkg"
            info["year"] = int(year)
            info["xsec"] = 0.0
            rootFiles = [rf.name for rf in xrdClient.dirlist(info["location"])[1] if '.root' in rf.name]
            info["nFiles"] = len(rootFiles)
            output.append(info)

        out_json = "bkg_{}_{}.json".format(year,bkg)
        with open(out_json,"w") as outfile:
            json.dump(output,outfile,indent=4)


