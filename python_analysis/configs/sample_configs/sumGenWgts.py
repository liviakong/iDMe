import uproot
import numpy as np
import json
import sys
from XRootD import client

inputJson = sys.argv[1]

with open(inputJson) as f:
    samples = json.load(f)
for samp in samples:
    loc = samp["location"]
    if '.root' in loc:
        tree = uproot.open(loc)['ntuples_gbm/genT']
        if tree.num_entries == 0:
            sum_wgt = 0
        else:
            sum_wgt = np.sum(tree['genWgt'].array())
    else:
        sum_wgt = 0
        xrdClient = client.FileSystem("root://cmseos.fnal.gov")
        status, flist = xrdClient.dirlist(loc)
        fullList = ["root://cmsxrootd.fnal.gov/"+loc+"/"+item.name for item in flist if '.root' in item.name]
        for f in fullList:
            tree = uproot.open(f)['ntuples_gbm/genT']
            if tree.num_entries == 0:
                continue
            else:
                sum_wgt += np.sum(tree['genWgt'].array())
    samp['sum_wgt'] = float(sum_wgt)

with open(inputJson,'w') as f:
    json.dump(samples,f,indent=4)