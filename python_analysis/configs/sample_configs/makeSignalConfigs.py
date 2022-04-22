from XRootD import client
import json
import sys

year = str(sys.argv[1])
alpha = str(sys.argv[2])
if len(sys.argv) < 3:
    print("Not enough input arguments!")
    print("Usage: python makeSignalConfigs.py year alpha")

prefix = "/store/group/lpcmetx/iDMe//Samples/Ntuples/signal/"
xrdClient = client.FileSystem("root://cmseos.fnal.gov")
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
        info["sum_wgt"] = 0.0
        info["type"] = "signal"
        info["year"] = int(year)
        info["alphaD"] = alpha
        info["xsec"] = 0.0
        output.append(info)

out_json = "signal_{0}_{1}.json".format(year,alpha)
with open(out_json,"w") as outfile:
    json.dump(output,outfile,indent=4)