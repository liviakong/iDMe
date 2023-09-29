from XRootD import client
import sys
import json
inFile = sys.argv[1]
with open(inFile,"r") as fin:
    cfg = json.load(fin)
xrdClient = client.FileSystem("root://cmseos.fnal.gov")

for i,sample in enumerate(cfg):
    loc = sample['location']
    if type(loc) != list:
        status, flist = xrdClient.dirlist(loc)
        fullList = ["root://xcache/"+loc+"/"+item.name for item in flist if (('.root' in item.name) and (item.name not in sample['blacklist']))]
    else:
        fullList = []
        for l in loc:
            status, flist = xrdClient.dirlist(l)
            fullList.extend(["root://xcache/"+l+"/"+item.name for item in flist if (('.root' in item.name) and (item.name not in sample['blacklist']))])
    cfg[i]['fileset'] = fullList

with open(inFile,"w") as fout:
    json.dump(cfg,fout,indent=4)
